
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include "jello.h"
#include "half_edge.h"
#include "my_opengl_util.h"

#define NUM_OBJ_FILES 3
#define access_ifs(face, corner, dimension) (ifs[cur_model]->vertices[ifs[cur_model]->faces[face].corners[corner]][dimension])


const char * obj_files [NUM_OBJ_FILES] = {"i_3d.obj", "donut.obj", "cube.obj"};
indexed_face_set_t * ifs [NUM_OBJ_FILES];
half_edge_structure_t * hes [NUM_OBJ_FILES];
int cur_model = 2;
int subdivides = 1;


int main(int argc, char * argv[]) {
    load_models();

    glutInit(&argc, (char**)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutCreateWindow("JELLO");

    init();
    
    glutKeyboardFunc(keyboard);
    glutDisplayFunc(display);
    glutTimerFunc(100, update, 60);
    
    glutMainLoop();

    return 0;
}

void free_models() {

}

void load_models() {
    for (int i = 2; i < NUM_OBJ_FILES; i++) {
        printf("\nBuilding %s data structures.\n", obj_files[i]);
        ifs[i] = (indexed_face_set_t * )malloc(sizeof(indexed_face_set_t));
        read_obj(ifs[i], obj_files[i]);
        catmull_clark_subdivide(ifs[i], subdivides, 0);
    }
}

void init(void) {
    glClearColor(0.,0.,0.,1.);
    glMatrixMode(GL_PROJECTION);
    gluPerspective(40.,1.,1.,40.);
    glMatrixMode(GL_MODELVIEW);
    gluLookAt(10, 5, 10, 0, 0, 0, 0, 1, 0);
}

void display(void) {
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glColor3f(1.,0.,1.);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    for (int face = 0; face < ifs[cur_model]->num_faces; face++) {
        glBegin(GL_POLYGON);
        for (int corner = 0; corner < ifs[cur_model]->faces[face].degree; corner++) {
            glVertex3f(access_ifs(face, corner, 0),
                       access_ifs(face, corner, 1),
                       access_ifs(face, corner, 2));
        }
        glEnd();
    }

    glutSwapBuffers();
}
            

void keyboard(unsigned char key, int x, int y) {
    switch(key) {
        case 'w':
            glRotatef(5.0, 1.0, 0.0, 0.0);
            break;
        case 'a':
            glRotatef(5.0, 0.0, 1.0, 0.0);
            break;
        case 's':
            glRotatef(-5.0, 1.0, 0.0, 0.0);
            break;
        case 'd':
            glRotatef(-5.0, 0.0, 1.0, 0.0);
            break;
        case '[':
            subdivides <= 0 ? 0 : subdivides--;
            free_models();
            load_models();
            break;
        case ']':
            subdivides >= 5 ? 0 : subdivides++;
            free_models();
            load_models();
            break;
    }
    if (key > '0' && key < '1' + NUM_OBJ_FILES) {
        cur_model = (int)(key - '1');
    }
}

void update(int v) {
    glutPostRedisplay();
    glutTimerFunc(1000/60, update, v);
}



