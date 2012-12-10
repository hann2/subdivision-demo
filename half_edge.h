
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
cfuhash_table_t * calculate_edge_points(indexed_face_set_t * ifs, half_edge_structure_t * hes, int * midpoints);
void edge_points_helper(cfuhash_table_t * hash_table, cfuhash_table_t * edge_points, indexed_face_set_t * ifs, int * midpoints, half_edge_t * he);
void change_old_vertices(indexed_face_set_t * ifs, half_edge_structure_t * hes, int * midpoints, int flags);
void change_verts_helper(indexed_face_set_t * ifs, double ** new_verts, cfuhash_table_t * hash_table, cfuhash_table_t * verts_touched, half_edge_t * he, int * midpoints, int flags);
int get_edge(cfuhash_table_t * edge_points, int point1, int point2);
int * calculate_midpoints(indexed_face_set_t * ifs);
face_t * make_new_faces(indexed_face_set_t * ifs, cfuhash_table_t * edge_points, int * midpoints);
void free_half_edge_helper(indexed_face_set_t * ifs, cfuhash_table_t * hash_table, half_edge_t * he);
void free_half_edge(indexed_face_set_t * ifs, half_edge_structure_t * hes);
int match_corner(indexed_face_set_t * ifs, int face, int vertex);

//generating keys for hash tables
unsigned long long generate_key(int face, int corner);
unsigned long long generate_edge_key(int point1, int point2);

#endif /* _HALF_EDGE_H */

