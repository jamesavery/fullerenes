add_library(fullerenes SHARED
  graph.cc cubicgraph.cc planargraph.cc fullerenegraph.cc triangulation.cc polyhedron.cc
  spiralc.cc planargraph-io.cc polyhedron-io.cc debug.cc fold.cc  unfold.cc geometryc.cc
  delaunay.cc  graph_fortran.cc isomerdb.cc layout.cc layout-optimise.cc
  polyhedron-optimise.cc symmetry.cc eisenstein.cc ../contrib/mgmres.cc fast_dualize.cc
)

if(NOT CMAKE_CUDA_COMPILER)
  add_library(fullerenes SHARED isomer_batch.cc)
endif()

target_link_libraries(fullerenes fortran_opt buckygen gfortran)


