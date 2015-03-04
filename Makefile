CC=g++

CFLAGS = -std=c++11 -g -O3

DEPS = template-rt.cpp
OBJ = template-rt.o 
SLIBS = 
LDFLAGS = 
LDLIBS = -lstdc++ -lGL -lglut -lGLEW

all: ray

%.o: %.c* $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< 

ray: $(OBJ) $(SLIBS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ) ray *~ *.ppm [#]*[#] .\#*
