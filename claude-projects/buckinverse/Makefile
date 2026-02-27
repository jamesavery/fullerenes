CXX = g++
CXXFLAGS = -std=c++20 -O2 -g -Wall -Wextra

# Project root
ROOT = ../..
INCLUDE = -I$(ROOT)/include -I.
LIBDIR = $(ROOT)/build2/src/c++
LIBS = -L$(LIBDIR) -lfullerenes -Wl,-rpath,$(LIBDIR)

# Fortran runtime (needed by libfullerenes)
FLIBS = -lgfortran

all: test_phase1 test_phase2

buckinverse.o: buckinverse.cc buckinverse.hh
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

test_phase1: test_phase1.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test_phase2: test_phase2.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test: test_phase1 test_phase2
	./test_phase1
	./test_phase2

clean:
	rm -f *.o test_phase1 test_phase2

.PHONY: all test clean
