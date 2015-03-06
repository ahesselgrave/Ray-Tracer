CC=g++

CFLAGS = -std=c++11 -g -O3

DEPS = template-rt.cpp
OBJ = template-rt.o 
SLIBS = 
LDFLAGS = 
LDLIBS = -lstdc++ -lGL -lglut -lGLEW

all: template-rt

%.o: %.c* $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< 

template-rt: $(OBJ) $(SLIBS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ) template-rt *~ *.ppm [#]*[#] .\#*
