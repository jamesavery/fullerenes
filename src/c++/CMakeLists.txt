add_library(fullerenes SHARED
  graph.cc cubicgraph.cc planargraph.cc fullerenegraph.cc triangulation.cc polyhedron.cc spiralc.cc planargraph-io.cc polyhedron-io.cc
  debug.cc fold.cc  unfold.cc geometryc.cc delaunay.cc  graph_fortran.cc isomerdb.cc layout.cc layout-optimize.cc
  polyhedron-optimize.cc symmetry.cc eisenstein.cc ../contrib/mgmres.cc)

target_link_libraries(fullerenes fortran_opt buckygen gfortran)


