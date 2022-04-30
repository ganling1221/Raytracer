LIBRARIES = -framework OpenGL -framework GLUT -lm

COMPILER = g++

COMPILERFLAGS = -O3 $(INCLUDE)

PROGRAM = ray_tracer
SOURCE = ray_tracer.cpp
OBJECT = ray_tracer.o

.cpp.o: 
	$(COMPILER) -c $(COMPILERFLAGS) $<

all: $(PROGRAM)

$(PROGRAM): $(OBJECT)
	$(COMPILER) $(COMPILERFLAGS) -o $(PROGRAM) $(OBJECT) $(LIBRARIES)

clean:
	-rm -rf core *.o *~ "#"*"#" $(PROGRAM)
