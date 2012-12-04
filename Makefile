
COMPILER = g++
COMPILER_OPTS = -c -g -O0 -Wall -Werror
LINKER = g++
LINKER_OPTS = -lpng -lGL -lGLU -lglut

all : jello

cfu_hash.o : cfu_hash.c cfu_hash.h
	$(COMPILER) $(COMPILER_OPTS) cfu_hash.c

my_string_lib.o : my_string_lib.c my_string_lib.h
	$(COMPILER) $(COMPILER_OPTS) my_string_lib.c

half_edge.o : half_edge.c half_edge.h my_string_lib.h
	$(COMPILER) $(COMPILER_OPTS) half_edge.c

my_opengl_util.o : my_opengl_util.c my_opengl_util.h my_string_lib.h
	$(COMPILER) $(COMPILER_OPTS) my_opengl_util.c

jello.o : jello.c jello.h
	$(COMPILER) $(COMPILER_OPTS) jello.c

jello : jello.o half_edge.o my_string_lib.o my_opengl_util.o cfu_hash.o
	$(LINKER) jello.o half_edge.o my_string_lib.o my_opengl_util.o cfu_hash.o $(LINKER_OPTS) -o jello

clean :
	-rm -f *.o jello