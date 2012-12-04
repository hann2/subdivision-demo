
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "half_edge.h"

#define x_mask (0x4)
#define y_mask (0x2)
#define z_mask (0x1)
#define unmask_face(bit_vector)   ((unsigned long)(bit_vector & -1))
#define unmask_corner(bit_vector)   ((unsigned long)((bit_vector >> 32) & -1))
#define get_corner(ifs, face, corner) (ifs->faces[face].corners[corner])

/*
    Sets up a recursive call to add_vertex.  This builds a half edge data structure from an indexed face set
*/
half_edge_structure_t * build_half_edge(indexed_face_set_t * ifs) {
    //our hash table will use face index concantated with corner as a key
    cfuhash_table_t * hash_table = cfuhash_new_with_initial_size(ifs->num_faces*3);

    //recursively add vertices
    half_edge_structure_t * hes = (half_edge_structure_t *)malloc(sizeof(half_edge_structure_t));
    hes->handle_edge = add_vertex(ifs, hash_table, 0, 0);

    hes->num_edges = cfuhash_num_entries(hash_table);
    assert((hes->num_edges & 1) == 0);
    hes->num_edges /= 2;

    cfuhash_destroy(hash_table);

    return hes;
}

/*
    Recursive function that builds a half_edge struct for a face and one of its corners.
    The corner represents the source vertex of the half edge
    Uses hash_table to keep track of what edges have already been processed.
*/
half_edge_t * add_vertex(indexed_face_set_t * ifs, cfuhash_table_t * hash_table, int face, int corner) {

    int degree = ifs->faces[face].degree;

    // build its half edge
    half_edge_t * h = (half_edge_t *)malloc(sizeof(half_edge_t));

    //get end_vert by going counter clockwise around face
    h->end_vert = get_corner(ifs, face, (corner+1)%degree);
    printf("End vert: %d\n", h->end_vert);
    h->left_face = face;

    // add to dictionary, key is 64 bits
    unsigned long long key  = generate_key(face, corner);
    cfuhash_put_data(hash_table, &key, 8, &h, 4, NULL);

    key  = generate_key(face, (corner+1)%degree);
    // if next is not in dictionary, it needs to be processed
    //if it is, copy it from hash table
    if (!cfuhash_get_data(hash_table, &key, 8, (void**)&(h->next), NULL)) {
        h->next = add_vertex(ifs, hash_table, face, (corner+1)%degree);
    } else {
        printf("End vert: %d\n", h->next->end_vert);
    }

    //look for face that shares these two vertices for opp
    key  = find_face(ifs, face, corner);
    // if opp is not in dictionary, it needs to be processed
    //if it is, copy it from hash table
    if (!cfuhash_get_data(hash_table, &key, 8, (void**)(&(h->opp)), NULL)) {
        h->opp = add_vertex(ifs, hash_table, unmask_face(key), unmask_corner(key));
    } else {
        printf("End vert: %d\n", h->next->end_vert);
    }

    return h;
}

/*
    generates a 64 bit key from two integers by concatenating them
*/
unsigned long long generate_key(int a, int b) {
    printf("%d", b);
    return (((unsigned long long)a) << 32) | ((unsigned long long)b);
}

/*
    takes face and which corner in that face (0, 1 or 2), and finds another face which shares this corner and the next corner around face counterclockwise
*/
unsigned long long find_face(indexed_face_set_t * ifs, int face, int s_corner) {
    //source and dest index for original face, will be flipped for the other
    int source_vert = get_corner(ifs, face, s_corner);
    int dest_vert   = get_corner(ifs, face, (s_corner+1)%(ifs->faces[face].degree));

    for (int adj_face = 0; adj_face < ifs->num_faces; adj_face++) {
        int degree = ifs->faces[adj_face].degree;
        for (int corner = 0; corner < degree; corner++) {

            if ((get_corner(ifs, adj_face, corner) == dest_vert) && (get_corner(ifs, adj_face, (corner+1)%degree) == source_vert)) {
                return generate_key(adj_face, (corner+1)%degree);
            } else if ((get_corner(ifs, adj_face, corner) == source_vert) && (get_corner(ifs, adj_face, (corner+1)%degree) == dest_vert) && (face != adj_face)) {
                //should never reach here
                printf("Face %d is not counterclockwise!\n", face);
                assert(0);
            }
        }
    }

    //should never reach here
    printf("%d %d\n", face, s_corner);
    printf("Cannot find a face with vertex %d followed by %d.  Fix your mesh.  Goodbye.\n\n", dest_vert, source_vert);
    assert(0);

    return 0;
}


half_edge_structure_t * catmull_clark_subdivide(indexed_face_set_t * ifs, int iterations, int flags) {

    half_edge_structure_t * hes = build_half_edge(ifs);

    if (iterations < 1 || iterations > 8) {
        return NULL;
    }

    //create array for [num_vertices + num_faces + num_edges] vertices
    //copy in old vertices
    double ** new_verts = (double **)realloc(ifs->vertices, (ifs->num_vertices + ifs->num_faces + hes->num_edges)*sizeof(double *));
    //create array for [num_edges] new faces
    face_t * new_faces = (face_t *)calloc(hes->num_edges, sizeof(face_t));

    int * midpoints = calculate_midpoints(ifs, new_verts);

    //calculate edge_points and put them into a hash table
    //key for an edge is smaller vertex index concantated with larger vertex index
    cfuhash_table_t * edge_points = calculate_edge_points(ifs, hes, new_verts, midpoints);

    //traverses vertices in hes and change their values in new_verts
    change_old_vertices(hes, midpoints, new_verts, flags);

    make_new_faces(ifs, edge_points, midpoints);

    //FIX THIS
    //free old faces NONTRIVIAL
    //free hash table of edge points
    cfuhash_destroy(edge_points);
    //free midpoints array
    free(midpoints);
    //FIX THIS
    //free all half edges

    //set new vertices
    ifs->vertices = new_verts;
    //set new faces
    ifs->num_faces = hes->num_edges;
    ifs->faces = new_faces;


    if (iterations > 1) {
        return catmull_clark_subdivide(ifs, iterations-1, flags);
    }

    return hes;
}

/*
    iterate through faces to calculate midpoints
    midpoint positions go into new vertices, index goes into midpoints array
*/
int * calculate_midpoints(indexed_face_set_t * ifs, double ** new_verts) {
    //create array for [num_faces] midpoints
    int * midpoints = (int *)calloc(ifs->num_faces, sizeof(int));
    for (int i = 0; i < ifs->num_faces; i++) {
        for (int dimension = 0; dimension < 3; dimension++) {
            double avg = 0;
            for (int corner = 0; corner < ifs->faces[i].degree; corner++) {
                avg += new_verts[get_corner(ifs, i, corner)][dimension];
            }
            avg /= ifs->faces[i].degree;
            new_verts[ifs->num_vertices] = (double *)calloc(3, sizeof(double));
            new_verts[ifs->num_vertices][dimension] = avg;
        }
        midpoints[i] = ifs->num_vertices;
        ifs->num_vertices++;
    }
    return midpoints;
}

/*
    iterate over old faces
    make new faces and add them to new faces array
*/
void make_new_faces(indexed_face_set_t * ifs, cfuhash_table_t * edge_points, int * midpoints) {
    for (int i = 0; i < ifs->num_faces; i++) {
        for (int corner = 0; corner < ifs->faces[i].degree; corner++) {
            face_t new_face;
            int old_degree = ifs->faces[i].degree;
            //catmull clark subdivision makes all quads
            new_face.degree = 4;
            new_face.corners = (int *)calloc(4, sizeof(int));
            new_face.corners[0] = get_corner(ifs, i, corner);
            new_face.corners[1] = get_edge(edge_points, get_corner(ifs, i, corner), get_corner(ifs, i, (corner + 1)%old_degree));
            new_face.corners[2] = midpoints[i];
            new_face.corners[3] = get_edge(edge_points, get_corner(ifs, i, corner), get_corner(ifs, i, (corner + old_degree - 1)%old_degree));
        }
    }
}


unsigned long long generate_edge_key(int point1, int point2) {
    assert(point1 != point2);
    if (point1 > point2) {
        int temp = point1;
        point1 = point2;
        point2 = temp;
    }
    return generate_key(point1, point2);
}

/*
    finds edge in hash_table between these two points, -1 if it finds nothing
*/
int get_edge(cfuhash_table_t * edge_points, int point1, int point2) {
    unsigned long long key = generate_edge_key(point1, point2);
    int * edge;
    if (cfuhash_get_data(edge_points, &key, 8, (void**)(&edge), NULL)){
        return *edge;
    }
    return -1;
}

cfuhash_table_t * calculate_edge_points(indexed_face_set_t * ifs, half_edge_structure_t * hes, double ** vertices, int * midpoints) {

    //create hash_table for [num_edges] edgepoints
    //key will be the lower position followed by the higher position, this will be chosen by the compare_positions function
    //value will be index in new_verts where data for the edgepoint can be found
    cfuhash_table_t * edge_points = cfuhash_new_with_initial_size(hes->num_edges);
    cfuhash_table_t * hash_table = cfuhash_new_with_initial_size(ifs->num_faces*3);

    edge_points_helper(hash_table, edge_points, ifs, vertices, midpoints, hes->handle_edge);

    return edge_points;
}


int match_corner(indexed_face_set_t * ifs, int face, int vertex) {
    for (int corner = 0; corner < ifs->faces[face].degree; corner++) {
        if (vertex == get_corner(ifs, face, corner)) {
            return corner;
        }
    }
    printf("Face %d doesn't have vertex %x!\n", face, vertex);
    assert(0);
}


void edge_points_helper(cfuhash_table_t * hash_table, cfuhash_table_t * edge_points, indexed_face_set_t * ifs, double ** vertices, int * midpoints, half_edge_t * he) {
    //printf("End vert: %d\n", he->end_vert);
    unsigned long long key  = generate_key(he->left_face, match_corner(ifs, he->left_face, he->end_vert));

    if (cfuhash_exists_data(hash_table, &key, 4)) {
        return;
    }

    //dont care about value, just using hash_table as a set
    int val = 0;
    cfuhash_put_data(hash_table, &key, 8, &val, 4, NULL);

    edge_points_helper(hash_table, edge_points, ifs, vertices, midpoints, he->next);
    edge_points_helper(hash_table, edge_points, ifs, vertices, midpoints, he->opp);

    //calculating stuffon the way up the tree of activations

    vertices[ifs->num_vertices] = (double *)calloc(3, sizeof(double));
    //edge point is the average of two touching face points
    for (int dimension = 0; dimension < 3; dimension++) {
       vertices[ifs->num_vertices][dimension] = (vertices[midpoints[he->left_face]][dimension] + vertices[midpoints[he->opp->left_face]][dimension])/2.;
    }

    //put index in vertices into hash table
    key = generate_edge_key(he->end_vert, he->opp->end_vert);
    cfuhash_put_data(edge_points, &key, 8, &(ifs->num_vertices), 4, NULL);

    ifs->num_vertices++;
}

//FIX THIS
void change_old_vertices(half_edge_structure_t * hes, int * midpoints, double ** vertices, int flags) {
    return;
}


//outline for doing a traversal of a half_edge datastructure
//don't call this
void traverse_graph(cfuhash_table_t * hash_table, indexed_face_set_t * ifs, half_edge_t * he){
    unsigned long long key  = generate_key(he->left_face, match_corner(ifs, he->left_face, he->end_vert));

    if (cfuhash_exists_data(hash_table, &key, 4)) {
        return;
    }

    //do stuff on the way down the tree of activations

    int val =0;
    cfuhash_put_data(hash_table, &key, 8, &val, 4, NULL);

    traverse_graph(hash_table, ifs, he->next);
    traverse_graph(hash_table, ifs, he->opp);

    //do stuff on the way up the tree of activations
}


