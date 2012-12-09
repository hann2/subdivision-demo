
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "half_edge.h"

#define unmask_face(bit_vector)   ((unsigned long)((bit_vector >> 32) & -1))
#define unmask_corner(bit_vector)   ((unsigned long)(bit_vector & -1))
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
    h->left_face = face;

    // add to dictionary, key is 64 bits
    unsigned long long key  = generate_key(face, corner);
    //printf("adress put into hash_table: %x\n", (unsigned int)(unsigned long)h);
    assert(cfuhash_put_data(hash_table, &key, sizeof(unsigned long long), h, sizeof(half_edge_t *), NULL));

    key  = generate_key(face, (corner+1)%degree);
    // if next is not in dictionary, it needs to be processed
    //if it is, copy it from hash table
    if (!cfuhash_get_data(hash_table, &key, sizeof(unsigned long long), (void**)&(h->next), NULL)) {
        h->next = add_vertex(ifs, hash_table, face, (corner+1)%degree);
    }
    //look for face that shares these two vertices for opp
    key  = find_face(ifs, face, corner);
    // if opp is not in dictionary, it needs to be processed
    //if it is, copy it from hash table
    if (!cfuhash_get_data(hash_table, &key, sizeof(unsigned long long), (void**)&(h->opp), NULL)) {
        h->opp = add_vertex(ifs, hash_table, unmask_face(key), unmask_corner(key));
    }

    return h;
}

/*
    generates a 64 bit key from two integers by concatenating them
*/
unsigned long long generate_key(int a, int b) {
    return (((unsigned long long)a) << 32) | ((unsigned long long)b);
}

/*
    takes face and which corner in that face (0, 1 or 2), and finds another face which shares this corner and the next corner around face counterclockwise
*/
unsigned long long find_face(indexed_face_set_t * ifs, int face, int s_corner) {
    //source and dest index for original face, source and dest will be flipped for the adj_face
    int source_vert = get_corner(ifs, face, s_corner);
    int dest_vert   = get_corner(ifs, face, (s_corner+1)%(ifs->faces[face].degree));

    for (int adj_face = 0; adj_face < ifs->num_faces; adj_face++) {
        int degree = ifs->faces[adj_face].degree;
        for (int corner = 0; corner < degree; corner++) {

            if ((get_corner(ifs, adj_face, corner) == dest_vert) && (get_corner(ifs, adj_face, (corner+1)%degree) == source_vert)) {
                return generate_key(adj_face, corner);
            } else if ((get_corner(ifs, adj_face, corner) == source_vert) && (get_corner(ifs, adj_face, (corner+1)%degree) == dest_vert) && (face != adj_face)) {
                //should never reach here
                printf("Face %d is not counterclockwise! Fix your mesh!\n", face);
                assert(0);
            }
        }
    }

    //should never reach here
    printf("Cannot find a face with vertex %d followed by %d.  Fix your mesh.  Goodbye.\n\n", dest_vert, source_vert);
    assert(0);

    return 0;
}


half_edge_structure_t * catmull_clark_subdivide(indexed_face_set_t * ifs, int iterations, int flags) {

    half_edge_structure_t * hes = build_half_edge(ifs);
    printf("Starting Catmull Clark Subdivision with %d faces, %d edges, %d vertices.\n", ifs->num_faces, hes->num_edges, ifs->num_vertices);

    if (iterations < 1 || iterations > 8) {
        return NULL;
    }

    //create array for [num_vertices + num_faces + num_edges] vertices
    //copy in old vertices
    ifs->vertices = (double **)realloc(ifs->vertices, (ifs->num_vertices + ifs->num_faces + hes->num_edges)*sizeof(double *));

    int * midpoints = calculate_midpoints(ifs);

    //calculate edge_points and put them into a hash table
    //key for an edge is smaller vertex index concantated with larger vertex index
    cfuhash_table_t * edge_points = calculate_edge_points(ifs, hes, midpoints);

    //traverses vertices in hes and change their values in ifs->vertices
    change_old_vertices(ifs, hes, midpoints, flags);

    face_t * new_faces = make_new_faces(ifs, edge_points, midpoints);

    //free hash table of edge points
    cfuhash_destroy(edge_points);
    //free midpoints array
    free(midpoints);
    //free all half edges
    free_half_edge(ifs, hes);
    //free old faces
    for (int i = 0; i < ifs->num_faces; i++) {
        free(ifs->faces[i].corners);
    }
    free(ifs->faces);

    //set new faces
    ifs->num_faces = hes->num_edges*2;
    ifs->faces = new_faces;
    printf("Completed Catmull Clark Subdivision with %d faces, %d edges, %d vertices.\n", ifs->num_faces, ifs->num_faces*4, ifs->num_vertices);


    if (iterations > 1) {
        return catmull_clark_subdivide(ifs, iterations-1, flags);
    }

    return hes;
}

/*
    iterate through faces to calculate midpoints
    midpoint positions go into new vertices, index goes into midpoints array
*/
int * calculate_midpoints(indexed_face_set_t * ifs) {
    //create array for [num_faces] midpoints
    int * midpoints = (int *)calloc(ifs->num_faces, sizeof(int));
    for (int i = 0; i < ifs->num_faces; i++) {
        ifs->vertices[ifs->num_vertices] = (double *)calloc(3, sizeof(double));
        for (int dimension = 0; dimension < 3; dimension++) {
            double avg = 0;
            for (int corner = 0; corner < ifs->faces[i].degree; corner++) {
                avg += ifs->vertices[get_corner(ifs, i, corner)][dimension];
            }
            avg /= ifs->faces[i].degree;
            ifs->vertices[ifs->num_vertices][dimension] = avg;
        
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
face_t * make_new_faces(indexed_face_set_t * ifs, cfuhash_table_t * edge_points, int * midpoints) {
    //create array for [num_edges] new faces
    face_t * new_faces = (face_t *)calloc(cfuhash_num_entries(edge_points)*2, sizeof(face_t));

    for (int old_face = 0, new_face = 0; old_face < ifs->num_faces; old_face++) {
        for (int corner = 0; corner < ifs->faces[old_face].degree; corner++, new_face++) {
            face_t new_face_struct;
            int old_degree = ifs->faces[old_face].degree;
            //catmull clark subdivision makes all quads
            new_face_struct.degree = 4;
            new_face_struct.corners = (int *)calloc(4, sizeof(int));
            new_face_struct.corners[0] = get_corner(ifs, old_face, corner);
            new_face_struct.corners[1] = get_edge(edge_points, get_corner(ifs, old_face, corner), get_corner(ifs, old_face, (corner + 1)%old_degree));
            new_face_struct.corners[2] = midpoints[old_face];
            new_face_struct.corners[3] = get_edge(edge_points, get_corner(ifs, old_face, corner), get_corner(ifs, old_face, (corner + old_degree - 1)%old_degree));
            //printf("Quad %d: %d, %d, %d, %d.\n", new_face, new_face_struct.corners[0], new_face_struct.corners[1], new_face_struct.corners[2], new_face_struct.corners[3]);
            for (int i =0; i<4; i++) {
                //printf("\tVert: %f, %f, %f.\n", ifs->vertices[new_face_struct.corners[i]][0], ifs->vertices[new_face_struct.corners[i]][1], ifs->vertices[new_face_struct.corners[i]][2]);
            }
            new_faces[new_face] = new_face_struct;
        }
    }
    return new_faces;
}


unsigned long long generate_edge_key(int point1, int point2) {
    assert(point1 != point2);
    if (point1 > point2) {
        int temp = point1;
        point1 = point2;
        point2 = temp;
        assert(generate_edge_key(point1, point2) == generate_key(point1, point2));
    }
    assert(point1 < point2);
    return generate_key(point1, point2);
}

/*
    finds edge in hash_table between these two points, kills program if it finds nothing
*/
int get_edge(cfuhash_table_t * edge_points, int point1, int point2) {
    unsigned long long key = generate_edge_key(point1, point2);
    int edge = 0;
    if (cfuhash_get_data(edge_points, &key, 8, (void**)&edge, NULL)){
        return edge;
    }
    assert(0);
    return -1;
}



void free_half_edge(indexed_face_set_t * ifs, half_edge_structure_t * hes) {

    cfuhash_table_t * hash_table = cfuhash_new_with_initial_size(ifs->num_faces*4);

    free_half_edge_helper(ifs, hash_table, hes->handle_edge);

    cfuhash_destroy(hash_table);
}

void free_half_edge_helper(indexed_face_set_t * ifs, cfuhash_table_t * hash_table, half_edge_t * he) {
    //printf("face: %d: (%d, %d, %d, %d)\n", he->left_face, ifs->faces[he->left_face].corners[0], ifs->faces[he->left_face].corners[1], ifs->faces[he->left_face].corners[2], ifs->faces[he->left_face].corners[3]);
    unsigned long long key  = generate_key(he->left_face, match_corner(ifs, he->left_face, he->end_vert));

    if (cfuhash_exists_data(hash_table, &key, 8)) {
        return;
    }

    //dont care about value, just using hash_table as a set
    cfuhash_put_data(hash_table, &key, 8, (void **)(0), 4, NULL);

    free_half_edge_helper(ifs, hash_table, he->next);
    free_half_edge_helper(ifs, hash_table, he->opp);

    free(he);

}

cfuhash_table_t * calculate_edge_points(indexed_face_set_t * ifs, half_edge_structure_t * hes, int * midpoints) {

    //create hash_table for [num_edges] edgepoints
    //key will be the lower position followed by the higher position, this will be chosen by the compare_positions function
    //value will be index in new_verts where data for the edgepoint can be found
    cfuhash_table_t * edge_points = cfuhash_new_with_initial_size(hes->num_edges);
    cfuhash_table_t * hash_table = cfuhash_new_with_initial_size(ifs->num_faces*3);

    edge_points_helper(hash_table, edge_points, ifs, midpoints, hes->handle_edge);

    cfuhash_destroy(hash_table);

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


void edge_points_helper(cfuhash_table_t * hash_table, cfuhash_table_t * edge_points, indexed_face_set_t * ifs, int * midpoints, half_edge_t * he) {
    unsigned long long key  = generate_key(he->left_face, match_corner(ifs, he->left_face, he->end_vert));

    if (cfuhash_exists_data(hash_table, &key, 8)) {
        return;
    }

    //dont care about value, just using hash_table as a set
    cfuhash_put_data(hash_table, &key, 8, (void **)(0), 4, NULL);

    edge_points_helper(hash_table, edge_points, ifs, midpoints, he->next);
    edge_points_helper(hash_table, edge_points, ifs, midpoints, he->opp);


    //calculating stuff on the way up the tree of activations

    key = generate_edge_key(he->end_vert, he->opp->end_vert);
    //only add vertex if we haven't already seen it
    if (!cfuhash_exists_data(edge_points, &key, 8)) {
        ifs->vertices[ifs->num_vertices] = (double *)calloc(3, sizeof(double));
        //edge point is the average of two touching face points
        for (int dimension = 0; dimension < 3; dimension++) {
            //FIX THIS, also add average of edge end points
            ifs->vertices[ifs->num_vertices][dimension] = 
                (ifs->vertices[midpoints[he->left_face]][dimension] +
                ifs->vertices[midpoints[he->opp->left_face]][dimension] +
                ifs->vertices[he->end_vert][dimension] +
                ifs->vertices[he->opp->end_vert][dimension])/2.;
        }

        //put index in vertices into hash table
        cfuhash_put_data(edge_points, &key, 8, (void *)(long)ifs->num_vertices, 4, NULL);

        ifs->num_vertices++;
    }
}


void change_old_vertices(indexed_face_set_t * ifs, half_edge_structure_t * hes, int * midpoints, int flags) {
    cfuhash_table_t * hash_table = cfuhash_new_with_initial_size(hes->num_edges);
    cfuhash_table_t * verts_touched = cfuhash_new_with_initial_size(ifs->num_faces*3);

    change_verts_helper(ifs, hash_table, verts_touched, hes->handle_edge, midpoints, flags);

    cfuhash_destroy(hash_table);
    cfuhash_destroy(verts_touched);
}

void change_verts_helper(indexed_face_set_t * ifs, cfuhash_table_t * hash_table, cfuhash_table_t * verts_touched, half_edge_t * he, int * midpoints, int flags) {
    unsigned long long key  = generate_key(he->left_face, match_corner(ifs, he->left_face, he->end_vert));

    if (cfuhash_exists_data(hash_table, &key, 8)) {
        return;
    }

    //dont care about value, just using hash_table as a set
    cfuhash_put_data(hash_table, &key, 8, (void **)(0), 4, NULL);

    change_verts_helper(ifs, hash_table, verts_touched, he->next, midpoints, flags);
    change_verts_helper(ifs, hash_table, verts_touched, he->opp, midpoints, flags);

    key = he->end_vert;

    if (!cfuhash_exists_data(verts_touched, &key, 4)) {
        //see wikipedia of catmull clark subdivision for calculating new vertices
        double face_avg [] = {0, 0, 0};
        double edge_avg [] = {0, 0, 0};

        int num_neighboring_faces = 0;
        half_edge_t * cur_edge = he;
        
        while (num_neighboring_faces == 0 || cur_edge != he) {
            for (int dimension = 0; dimension < 3; dimension++) {
                face_avg[dimension] += ifs->vertices[midpoints[he->left_face]][dimension];
                edge_avg[dimension] += ifs->vertices[he->end_vert][dimension];
                printf("face_avg %lf  edge1_avg %lf  ", face_avg[dimension], edge_avg[dimension]);
                edge_avg[dimension] += ifs->vertices[he->next->end_vert][dimension];
                printf("edge2_avg %lf\n", edge_avg[dimension]);
            }
            printf("\n\n");
            cur_edge = cur_edge->next->opp;
            num_neighboring_faces++;
        }

        //do the division for the averages
        for (int dimension = 0; dimension < 3; dimension++) {
            face_avg[dimension] /= num_neighboring_faces;
            edge_avg[dimension] /= (num_neighboring_faces*2);
        }

        printf("\nNew values for vertex %d.\n", he->end_vert);
        for (int dimension = 0; dimension < 3; dimension++) {
            double temp = ifs->vertices[he->end_vert][dimension];
            ifs->vertices[he->end_vert][dimension] = (face_avg[dimension] + 2*edge_avg[dimension] + (num_neighboring_faces - 3)*ifs->vertices[he->end_vert][dimension])/num_neighboring_faces;
            printf("%lf + 2*%lf + (%d-3)*%lf)/%d = %lf\n", face_avg[dimension], edge_avg[dimension], num_neighboring_faces, temp, num_neighboring_faces, ifs->vertices[he->end_vert][dimension]);
        }

        //add vertex to verts_touched
        //dont care about value, just using verts_touched as a set
        cfuhash_put_data(verts_touched, &key, 4, (void **)(0), 4, NULL);
    }
}


//outline for doing a traversal of a half_edge datastructure
//don't call this
void traverse_graph(cfuhash_table_t * hash_table, indexed_face_set_t * ifs, half_edge_t * he){
    unsigned long long key  = generate_key(he->left_face, match_corner(ifs, he->left_face, he->end_vert));

    if (cfuhash_exists_data(hash_table, &key, 8)) {
        return;
    }

    //do stuff on the way down the tree of activations

    int val =0;
    cfuhash_put_data(hash_table, &key, 8, &val, 4, NULL);

    traverse_graph(hash_table, ifs, he->next);
    traverse_graph(hash_table, ifs, he->opp);

    //do stuff on the way up the tree of activations
}


