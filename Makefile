#
# Makefile for FULLERENE program
#
VERSION_NUMBER=4.4

CXX=g++
F90=gfortran
AR=ar

DIRECTORIES=-DFULLERENE_ROOT=\"${PWD}\" -DFULLERENE_DATABASE=\"${PWD}/database\"
WARNINGS=-Wall -Wno-sign-compare -Wno-unused-but-set-variable -Wno-char-subscripts

CXXFLAGS= -g -O3 -std=c++11 -m64 -fPIC $(WARNINGS) $(DIRECTORIES) -DVERSION_NUMBER=$(VERSION_NUMBER) 
FFLAGS= -g -O3 -m64 -Wall -cpp -D'VERSION_NUMBER="$(VERSION_NUMBER)"'
LIBRARIES=-lstdc++ -lgomp -lgfortran
# if your machine has enough memory, your gfortran is sufficiently new, and you need more then 5000 atoms
# you might need to change the gfortran compiler options:
#FFLAGS= -O3 -mcmodel=medium 

# Uncomment the following lines if you want to use Intel C++ and MKL
#CXX=icpc
#CXXFLAGS+=-mkl -DHAS_MKL -DHAS_LAPACK
#LIBRARIES+=-mkl
# Uncomment the following lines if you want to use system BLAS and LAPACK
#CXXFLAGS+=-DHAS_LAPACK 
#LIBRARIES+=-L${BLASDIR} -L${LAPACK_DIR} -mkl
#uncomment that following lines to use gsl (gnu scientific library)
#CXXFLAGS+=-DHAS_GSL
#LIBRARIES+=-lgsl


OBJECTS=main.o coord.o hamilton.o isomer.o opt.o ring.o sphere.o util.o datain.o geometry.o hueckel.o pentindex.o schlegel.o spiral.o volume.o
GRAPHOBJECTS= graph.o cubicgraph.o layout.o hamiltonian.o graph.o planargraph.o polyhedron.o polyhedron-optimize.o fullerenegraph.o graph_fortran.o mgmres.o geometryc.o unfold.o fold.o buckygen-wrapper.o triangulation.o symmetry.o isomerdb.o spherical-harmonic.o y3table.o
GRAPHFOBJECTS=geometry.o force.o diag.o dddihedral.o config.o opt-standalone.o

FOBJECTS=$(patsubst %.o, build/%.o, $(OBJECTS))
COBJECTS=$(patsubst %.o, build/%.o, $(GRAPHOBJECTS))
FLIBOBJECTS=$(patsubst %.o, build/%.o, $(GRAPHFOBJECTS))

TESTINP=$(wildcard input/*.inp)
TESTOUT=$(patsubst input/%.inp, output/%.out, $(TESTINP))
#
#
fullerene: $(FOBJECTS) build/libgraph.a
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
.PHONY: build/libgraph.a build/libgraph.so
build/libgraph.a: $(COBJECTS) $(FLIBOBJECTS)
	$(AR) rcs $@ $(COBJECTS) $(FLIBOBJECTS)

build/libgraph.so: $(COBJECTS) $(FLIBOBJECTS)
	c++ -shared -o $@ $(COBJECTS) $(FLIBOBJECTS)
#-----------------------------------------------------
test-%: tests/%.cc build/libgraph.a
	$(CXX) -I${PWD} $(CXXFLAGS) -o $@ $^ $(LIBRARIES)
#-----------------------------------------------------
app-%: apps/%.cc build/libgraph.a
	$(CXX) -I${PWD} $(CXXFLAGS) -o $@ $^ $(LIBRARIES)
#-----------------------------------------------------

output/%.out: input/%.inp
	./fullerene < $< > $@

tests: fullerene $(TESTOUT)

tags:
	ctags -e --c-kinds=pxd -R

clean:
	find . \( -name  "*~" -or  -name "#*#" -or -name "*.o" \) -exec rm {} \;

distclean: clean
	rm -f fullerene build/libgraph.a qmga.dat config.mod test-* app-*

#-----------------------------------------------------
