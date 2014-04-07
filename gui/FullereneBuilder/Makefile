SWIG=swig204 -w365,366,367,368,124,362,503 
# -debug-template
# -debug-tmsearch
# -debug-top 1,2,3,4

# Python include directory is OS typedefdependent
PYTHONINCLUDE= /usr/include/python2.7

#CXX=clang -x c++ -O2
CXX=c++
#CXX=icpc

DIRECTORIES=-DFULLERENE_ROOT=${PWD}/../.. -DFULLERENE_DATABASE=\"${PWD}/../../database\"
WARNINGS=-Wall -Wno-sign-compare -Wno-unused-but-set-variable -Wno-char-subscripts

CXXFLAGS=-g -O3 -std=c++11 -m64 -march=native -fPIC $(WARNINGS) $(DIRECTORIES) -DVERSION_NUMBER=$(VERSION_NUMBER) -I../../
CXXLINK=g++ $(CXXFLAGS) -Wl,-rpath=${PWD}/../../build/


pymodule:  python-bindings_wrap.cxx FullereneSelect.o fullerenegui.py
	${CXXLINK} -o _fullwrap.so -shared python-bindings_wrap.cxx FullereneSelect.o ${PWD}/../../build/libgraph.so -lgfortran -lgsl -lblas -llapack -I${PYTHONINCLUDE}


clean:
	rm -f *~ \#*\# *.o *.oct _*.so test *_wrap.cc *_wrap.cxx octave-core space.py *.pyc


# Intermediate targets
python-bindings_wrap.cxx: python-bindings.i FullereneSelect.hh ../../build/libgraph.a
	${SWIG} -python -c++ python-bindings.i

fullerenegui.py: fullerenegui.ui
	pyuic4 -w $^ > $@
