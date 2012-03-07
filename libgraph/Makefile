AR=ar
CXX=g++
#CXXFLAGS=-O3 -msse4 -m32
CXXFLAGS=-g3 -m64
GRAPHOBJECTS=graph.o cubicgraph.o layout.o hamiltonian.o graph.o fullerenegraph.o graph_fortran.o

all: ../libgraph.a

%.o: %.cc %.hh

../libgraph.a: $(GRAPHOBJECTS)
	$(AR) rcs $@ $(GRAPHOBJECTS)

clean:
	rm -f *~ *.o \#*\#
