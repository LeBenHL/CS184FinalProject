OBJS = main.o three_d_vector.o lodepng.o particle.o marching_cube.o particle_grid.o
CC = g++
INCLUDE = -I ./ -g -pg 
FLAGS = -O3 -std=c++0x

ifeq ($(shell sw_vers 2>/dev/null | grep Mac | awk '{ print $$2}'),Mac)
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -I./include/ -I/usr/X11/include -DOSX
	LDFLAGS = -framework GLUT -framework OpenGL \
    	-L"/System/Library/Frameworks/OpenGL.framework/Libraries" \
    	-lGL -lGLU -lm -lstdc++
else
	CFLAGS = -g -DGL_GLEXT_PROTOTYPES -Iglut-3.7.6-bin
	LDFLAGS = -lglut -lGLU
endif
	
RM = /bin/rm -f 
all: main 
main: $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -o main $(OBJS) $(LDFLAGS) 
main.o: main.cpp three_d_vector.h lodepng.h particle.h
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c main.cpp -o main.o
three_d_vector.o: three_d_vector.cpp three_d_vector.h
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c three_d_vector.cpp -o three_d_vector.o
lodepng.o: lodepng.h lodepng.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c lodepng.cpp -o lodepng.o
particle.o: particle.h particle.cpp three_d_vector.h
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c particle.cpp -o particle.o
marching_cube.o: marching_cube.h marching_cube.cpp three_d_vector.h particle.h particle_grid.h
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c marching_cube.cpp -o marching_cube.o
particle_grid.o: particle_grid.h particle_grid.cpp three_d_vector.h particle.h
	$(CC) $(CFLAGS) $(INCLUDE) $(FLAGS) -c particle_grid.cpp -o particle_grid.o
clean: 
	$(RM) *.o as1
 


