
#include "cfu_hash.h"
#include "data_structs.h"

#ifndef _HALF_EDGE_H
#define _HALF_EDGE_H

//build_half_edge and its helper functions
half_edge_structure_t * build_half_edge(indexed_face_set_t * ifs);
half_edge_t * add_vertex(indexed_face_set_t * ifs, cfuhash_table_t * hash_table, int face, int corner);
unsigned long long find_face(indexed_face_set_t * ifs, int face, int corner);

//catmull clark subdivision and its helpers
half_edge_structure_t * catmull_clark_subdivide(indexed_face_set_t * ifs, int iterations, int flags);
cfuhash_table_t * calculate_edge_points(indexed_face_set_t * ifs, half_edge_structure_t * hes, double **vertices, int * midpoints);
void edge_points_helper(cfuhash_table_t * hash_table, cfuhash_table_t * edge_points, indexed_face_set_t * ifs, double ** vertices, int * midpoints, half_edge_t * he);
void change_old_vertices(half_edge_structure_t * hes, int * midpoints, double ** vertices, int flags);
int get_edge(cfuhash_table_t * edge_points, int point1, int point2);
int * calculate_midpoints(indexed_face_set_t * ifs, double ** new_verts);
void make_new_faces(indexed_face_set_t * ifs, cfuhash_table_t * edge_points, int * midpoints);

//generating keys for hash tables
unsigned long long generate_key(int face, int corner);
unsigned long long generate_edge_key(int point1, int point2);

#endif /* _HALF_EDGE_H */

