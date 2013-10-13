#
# Makefile for FULLERENE program
#
CXX=g++
F90=gfortran
AR=ar

CXXFLAGS= -O3 -m64 -fPIC -Wall -Wno-sign-compare -Wno-char-subscripts
#CXXFLAGS= -g3 -O3 -m64 -fPIC -Wall -Wno-sign-compare -std=c++0x
FFLAGS= -g3 -O3 -m64 -Wall 
LIBRARIES=-lstdc++ -lgomp
# if your machine has enough memory, your gfortran is sufficiently new, and you need more then 5000 atoms
# you might need to change the gfortran compiler options:
#FFLAGS= -O3 -mcmodel=medium 

# Uncomment the following lines if you want to use Intel C++ and MKL
#CXX=icpc
#CXXFLAGS+=-mkl -DHAS_MKL -DHAS_LAPACK
#LIBRARIES+=-mkl
# Uncomment the following lines if you want to use system BLAS and LAPACK
CXXFLAGS+=-DHAS_LAPACK 
LIBRARIES+=-lblas -llapack 


OBJECTS=main.o coord.o diag.o hamilton.o isomer.o opt.o ring.o sphere.o util.o datain.o force.o geometry.o dddihedral.o hueckel.o pentindex.o schlegel.o spiral.o volume.o
GRAPHOBJECTS= graph.o cubicgraph.o layout.o hamiltonian.o graph.o planargraph.o polyhedron.o fullerenegraph.o graph_fortran.o mgmres.o geometryc.o unfold.o fold.o buckygen-wrapper.o triangulation.o

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

build/%.o: contrib/%.cc
	$(CXX) $(CXXFLAGS) $(OPTIONS) -c $< -o $@
#-----------------------------------------------------
.PHONY: build/libgraph.a
build/libgraph.a: $(COBJECTS)
	$(AR) rcs $@ $(COBJECTS)

#-----------------------------------------------------
test-%: tests/%.cc build/libgraph.a
	$(CXX) -I${PWD} $(CXXFLAGS) -o $@ $^ $(LIBRARIES)
#-----------------------------------------------------
app-%: apps/%.cc build/libgraph.a
	$(CXX) -I${PWD} $(CXXFLAGS) -o $@ $^ $(LIBRARIES)
#-----------------------------------------------------

app-leapfrog: apps/leapfrog.cc build/libgraph.a build/opt-standalone.o 
	$(CXX) -I${PWD} $(CXXFLAGS) -o $@ $^ -lgfortran

output/%.out: input/%.inp
	./fullerene < $< > $@

tests: fullerene $(TESTOUT)

tags:
	ctags -e --c-kinds=pxd -R

clean:
	find . \( -name  "*~" -or  -name "#*#" -or -name "*.o" \) -exec rm {} \;

distclean: clean
	rm -f fullerene build/libgraph.a qmga.dat config.mod test-*

#-----------------------------------------------------
