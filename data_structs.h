

#ifndef _DATA_STRUCTS_H
#define _DATA_STRUCTS_H

typedef struct _face_t {
    //number of corners
    int degree;
    int * corners;
} face_t;

typedef struct _indexed_face_set_t {
    int num_vertices;
    int num_faces;
    double ** vertices;
    face_t * faces;
} indexed_face_set_t;

typedef struct _half_edge_t {
    _half_edge_t * opp;
    int end_vert;
    int left_face;
    _half_edge_t * next;
} half_edge_t;

typedef struct _half_edge_structure_t {
    int num_edges;
    half_edge_t * handle_edge;
} half_edge_structure_t;


#endif /* _DATA_STRUCTS_H */
