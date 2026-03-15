add_library(fullerenes SHARED
  graph.cc cubicgraph.cc planargraph.cc fullerenegraph.cc triangulation.cc polyhedron.cc
  spiralc.cc planargraph-io.cc polyhedron-io.cc debug.cc fold.cc  unfold.cc geometryc.cc
  delaunay.cc  delaunay_old.cc  graph_fortran.cc isomerdb.cc layout.cc layout2d.cc argparser.cc
  polyhedron-optimize.cc symmetry.cc eisenstein.cc deltahedron.cc buckinverse.cc ../contrib/mgmres.cc ${MORE_SOURCES}
)

#Ensure that CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is defined
if(NOT CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
  #Terminate
  message(FATAL_ERROR "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is not defined")
endif()

find_library(GOMP_LIB gomp HINTS /usr/lib/gcc/x86_64-linux-gnu/11)
if(GOMP_LIB)
  target_link_libraries(fullerenes ${GOMP_LIB} fortran_opt buckygen)
else()
  target_link_libraries(fullerenes fortran_opt buckygen)
endif()



