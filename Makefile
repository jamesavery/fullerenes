#
# Makefile for FULLERENE program
#
VERSION_NUMBER=4.5

#CXX=clang++
CXX=g++
F90=gfortran
AR=ar

DIRECTORIES=-DFULLERENE_ROOT=\"${PWD}\" -DFULLERENE_DATABASE=\"${PWD}/database\"
WARNINGS=-Wall -Wno-sign-compare -Wno-unused-but-set-variable -Wno-char-subscripts

#CXXFLAGS=-I/opt/local/include -L/opt/local/lib -g -std=c++11 -fPIC -mcmodel=large $(WARNINGS) $(DIRECTORIES) -DVERSION_NUMBER=$(VERSION_NUMBER) -fno-exceptions
CXXFLAGS=-I/opt/local/include -L/opt/local/lib -O2 -std=c++11 -fPIC -mcmodel=large $(WARNINGS) $(DIRECTORIES) -DVERSION_NUMBER=$(VERSION_NUMBER) -fno-exceptions
#FFLAGS=-cpp -L/opt/local/lib -Wall -fPIC -g  -mcmodel=large -D'VERSION_NUMBER="$(VERSION_NUMBER)"'
FFLAGS=-cpp -L/opt/local/lib -Wall -fPIC -O2  -mcmodel=large -D'VERSION_NUMBER="$(VERSION_NUMBER)"'
LIBRARIES= -lgfortran -lstdc++ -lm
#-lstdc++ 
# if your machine has enough memory, your gfortran is sufficiently new, and you need more then 5000 atoms
# you might need to change the gfortran compiler options:

# Uncomment the following lines if you want to use Intel C++ and MKL
#CXX=icpc
#CXXFLAGS+=-mkl -DHAS_MKL -DHAS_LAPACK
#LIBRARIES+=-mkl
#Uncomment the following lines if you want to use system BLAS and LAPACK
CXXFLAGS+=-DHAS_LAPACK 
LIBRARIES+=-llapack -lblas
#uncomment that following lines to use gsl (gnu scientific library)
CXXFLAGS+=-DHAS_GSL
LIBRARIES+=-lgsl


OBJECTS=main.o coord.o hamilton.o isomer.o opt.o ring.o sphere.o util.o datain.o geometry.o hueckel.o pentindex.o schlegel.o spiral.o volume.o
GRAPHOBJECTS= graph.o cubicgraph.o layout.o planargraph.o polyhedron.o polyhedron-optimize.o fullerenegraph.o graph_fortran.o mgmres.o geometryc.o unfold.o fold.o buckygen-wrapper.o triangulation.o symmetry.o isomerdb.o spherical-harmonic.o y3table.o layout-optimize.o delaunay.o spiralc.o planargraph-io.o polyhedron-io.o
GRAPHFOBJECTS=geometry.o force.o diag.o dddihedral.o config.o opt-standalone.o

FOBJECTS=$(patsubst %.o, build/%.o, $(OBJECTS))
COBJECTS=$(patsubst %.o, build/%.o, $(GRAPHOBJECTS))
FLIBOBJECTS=$(patsubst %.o, build/%.o, $(GRAPHFOBJECTS))

TESTINP=$(wildcard input/*.inp)
TESTOUT=$(patsubst input/%.inp, output/%.out, $(TESTINP))
#
#
fullerene: $(FOBJECTS) build/libgraph.a
	$(F90) $(FFLAGS) $(OPTIONS) $^ $(LIBRARIES) -o $@ -lstdc++
#-lstdc++ 

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
play-%: playground/%.cc build/libgraph.a
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
	rm -rf fullerene build/libgraph.a build/libgraph.so qmga.dat config.mod test-* app-* play-*

#-----------------------------------------------------
