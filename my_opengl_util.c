
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "my_opengl_util.h"
#include "my_string_lib.h"


void read_obj(indexed_face_set_t * ifs, const char * file_name) {

    int num_vertices = 0;
    int num_faces = 0;

    //current line being read
    char line[100];

    //open the file for reading
    FILE * obj = fopen (file_name, "rt");

    //go once through the file just counting faces and vertices
    while (fgets(line, 100, obj) != NULL) {
        if (str_starts_with(line, 'v')) {
            num_vertices++;
        } else if (str_starts_with(line, 'f')) {
            num_faces++;
        }
    }

    rewind(obj);

    ifs->vertices = (double **)calloc(num_vertices+1, sizeof(double *));
    ifs->faces = (face_t *)calloc(num_faces, sizeof(face_t));

    //will be incremented when parse_line is called
    ifs->num_vertices = 1;
    ifs->num_faces = 0;

    while (fgets(line, 100, obj) != NULL) {
        parse_line(ifs, line);
    }

    assert(ifs->num_vertices-1 == num_vertices);
    assert(ifs->num_faces == num_faces);

    //close .obj
    fclose(obj);
}


//parses a single line in an .obj file and adds it to the ifs data structure
void parse_line(indexed_face_set_t * ifs, char * line) {
    //index into a vertex or face
    //  ie:  v   23.231 42.124 51.2313
    //             0      1      2
    int index = 0;

    //current character in line
    char * cur_char = line;

    //buffer for reading floats/ints out of line
    char buffer[100];
    //index in buffer
    int buf_ind = 0;

    if (str_starts_with(line, 'v')) {
        //remove "v " from front of line
        cur_char += 2;
        ifs->vertices[ifs->num_vertices] = (double *)calloc(3, sizeof(double));

        while (cur_char == line || *(cur_char-1) != '\0') {
            //if you've reached the end of a number
            if (*cur_char == ' ' || *cur_char == '\0') {
                buffer[buf_ind] = '\0';

                ifs->vertices[ifs->num_vertices][index] = str_to_double(buffer);
                index++;
                //set to -1 b/c it will be incremented at end of the loop
                //we need it to be 0 for next iteration
                buf_ind = -1;
            } else {
                buffer[buf_ind] = *cur_char;
            }
            buf_ind++;
            cur_char++;
        }
        ifs->num_vertices++;
    } else if (str_starts_with(line, 'f')) {
        //remove "f  " from front of line
        cur_char += 3;
        //malloc corners and set degree
        int degree = count_chars(cur_char, ' ') + 1;
        ifs->faces[ifs->num_faces].degree = degree;
        ifs->faces[ifs->num_faces].corners = (int *)calloc(degree, sizeof(int));

        while (cur_char == line || *(cur_char-1) != '\0') {
            //if you've reached the end of a number
            if (*cur_char == ' ' || *cur_char == '\0') {
                buffer[buf_ind] = '\0';
                int vert = str_to_int(buffer);
                assert(vert <= ifs->num_vertices);
                ifs->faces[ifs->num_faces].corners[index] = vert;
                index++;
                //set to -1 b/c it will be incremented at end of the loop
                //we need it to be 0 for next iteration
                buf_ind = -1;
            } else {
                buffer[buf_ind] = *cur_char;
            }
            buf_ind++;
            cur_char++;
        }
        assert(index == degree);
        ifs->num_faces++;
    }

}






