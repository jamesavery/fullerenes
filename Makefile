#
# Makefile for FULLERENE program
#
CXX=g++
F90=gfortran
AR=ar

CXXFLAGS= -g3 -O3 -m64 -Wall -Wno-sign-compare
FFLAGS= -g3 -O3 -m64 -Wall
# if your machine has enough memory, your gfortran is sufficiently new, and you need more then 5000 atoms
# you might need to change the gfortran compiler options:
#FFLAGS= -O3 -m64 -mcmodel=medium 

OBJECTS=main.o coord.o diag.o hamilton.o isomer.o opt.o ring.o sphere.o util.o datain.o force.o geometry.o dddihedral.o hueckel.o pentindex.o schlegel.o spiral.o volume.o
GRAPHOBJECTS= graph.o cubicgraph.o layout.o hamiltonian.o graph.o planargraph.o polyhedron.o fullerenegraph.o graph_fortran.o mgmres.o

FOBJECTS=$(patsubst %.o, build/%.o, $(OBJECTS))
COBJECTS=$(patsubst %.o, build/%.o, $(GRAPHOBJECTS))
TESTINP=$(wildcard input/*.inp)
TESTOUT=$(patsubst input/%.inp, output/%.out, $(TESTINP))
#
#
fullerene: build/config.o $(FOBJECTS) build/libgraph.a
	$(F90) $(FFLAGS) $(OPTIONS) $^ $(LIBRARIES) -o $@ -lstdc++ -lgomp

#
# ############    Definition of the subroutines    ###############
#
#-----------------------------------------------------

build/config.o: source/config.f
	$(F90) $(FFLAGS) $(OPTIONS) -c $< -o $@

build/%.o: source/%.f build/config.o
	$(F90) $(FFLAGS) $(OPTIONS) -c $< -o $@

build/%.o: libgraph/%.cc
	$(CXX) $(CXXFLAGS) $(OPTIONS) -c $< -o $@

build/%.o: contrib/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTIONS) -c $< -o $@
#-----------------------------------------------------
.PHONY: build/libgraph.a
build/libgraph.a: $(COBJECTS)
	$(AR) rcs $@ $(COBJECTS)

#-----------------------------------------------------
test-%: tests/%.cc build/libgraph.a
	$(CXX) -I${PWD} $(CXXFLAGS) -o $@ $^ 
#-----------------------------------------------------

output/%.out: input/%.inp
	./fullerene < $< > $@

tests: fullerene $(TESTOUT)

tags:
	ctags -e --c-kinds=pxd -R

clean:
	find . \( -name  "*~" -or  -name "#*#" -or -name "*.o" \) -exec rm {} \;

distclean: clean
	rm -f fullerene build/libgraph.a qmga.dat config.mod

#-----------------------------------------------------
