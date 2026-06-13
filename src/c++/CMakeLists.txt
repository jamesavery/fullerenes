add_library(fullerenes SHARED
  graph.cc cubicgraph.cc planargraph.cc fullerenegraph.cc triangulation.cc polyhedron.cc
  spiralc.cc planargraph-io.cc polyhedron-io.cc debug.cc fold.cc  unfold.cc geometryc.cc
  delaunay.cc  delaunay_embed3d.cc  delaunay_old.cc  delaunay_alexandrov12.cc  graph_fortran.cc isomerdb.cc layout.cc layout2d.cc argparser.cc
  polyhedron-optimize.cc symmetry.cc eisenstein.cc eisenstein_raster.cc barycentric.cc
  delaunay_strip.cc delaunay_unfold.cc delaunay_metric.cc eisenstein_paint.cc
  deltahedron.cc deltahedron-io.cc buckinverse.cc ../contrib/mgmres.cc ${MORE_SOURCES}
)

#Ensure that CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is defined
if(NOT CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
  #Terminate
  message(FATAL_ERROR "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is not defined")
endif()

find_package(OpenMP)
target_link_libraries(fullerenes fortran_opt buckygen ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
if(OpenMP_CXX_FOUND)
  target_link_libraries(fullerenes OpenMP::OpenMP_CXX)
endif()



