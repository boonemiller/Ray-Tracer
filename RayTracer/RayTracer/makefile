CC=g++
CFLAGS= -std=c++14 -O3
DEPS = Ray.hpp bvh.hpp stb_image_write.h ./glm/glm/glm.hpp ./glm/glm/gtx/io.hpp
OBJ = Ray.o bvh.o main.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ray: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

