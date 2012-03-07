#
# Makefile for FULLERENE program
#
DIR=./source
#FCOMP= g77 -w
#FCOMP= gfortran -w -O3 -m32
FCOMP=gfortran -g3 -O3 -m64
OBJECTS= main.o coord.o diag.o hamilton.o isomer.o opt.o ring.o sphere.o util.o datain.o force.o hueckel.o pentindex.o schlegel.o spiral.o volume.o
TESTINP=$(wildcard *.inp)
TESTOUT=$(patsubst %.inp, %.out, $(TESTINP))
#
#
fullerene: $(OBJECTS) libgraph.a
	$(FCOMP) $(OPTIONS) $^ $(LIBRARIES) -o $@ -lstdc++

#
# ############    Definition of the subroutines    ###############
#
#-----------------------------------------------------
%.o: $(DIR)/%.f
	$(FCOMP) $(OPTIONS) -c $<
#-----------------------------------------------------
.PHONY: libfullerenegraph.a
libgraph.a: 
	cd libgraph && $(MAKE) 
#-----------------------------------------------------
test-%: tests/%.cc libgraph.a
	$(CXX) -I${PWD} $(CXXFLAGS) -o $@ $^ 
#-----------------------------------------------------


output/%.out: input/%.inp
	./fullerene < $< > $@

tests: fullerene $(TESTOUT)


clean:
	rm -f *~ \#*\# *.o *.a
	cd libgraph && make clean

#-----------------------------------------------------
