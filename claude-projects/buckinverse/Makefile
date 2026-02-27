CXX = g++
CXXFLAGS = -std=c++20 -O2 -g -Wall -Wextra

# Project root
ROOT = ../..
INCLUDE = -I$(ROOT)/include -I.
LIBDIR = $(ROOT)/build2/src/c++
LIBS = -L$(LIBDIR) -lfullerenes -Wl,-rpath,$(LIBDIR)

# Fortran runtime (needed by libfullerenes)
FLIBS = -lgfortran

all: test_phase1 test_phase2 test_phase3 test_scale test_fast test_roundtrip

buckinverse.o: buckinverse.cc buckinverse.hh
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

test_phase1: test_phase1.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test_phase2: test_phase2.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test_phase3: test_phase3.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test_scale: test_scale.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test_fast: test_fast.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test_roundtrip: test_roundtrip.cc buckinverse.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) $(FLIBS)

test: test_phase1 test_phase2 test_phase3
	./test_phase1
	./test_phase2
	./test_phase3

clean:
	rm -f *.o test_phase1 test_phase2 test_phase3 test_scale test_fast test_roundtrip

.PHONY: all test clean
