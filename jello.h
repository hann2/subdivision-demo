



#ifndef _JELLO_H
#define _JELLO_H



int main(int argc, char* argv[]);

void init(void);

void display(void);
            
void keyboard(unsigned char key, int x, int y);

void update(int v);

void free_models();

void load_models();


#endif /* _JELLO_H */

