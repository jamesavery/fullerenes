add_library(fullerenes SHARED
  graph.cc cubicgraph.cc planargraph.cc fullerenegraph.cc triangulation.cc polyhedron.cc
  spiralc.cc planargraph-io.cc polyhedron-io.cc debug.cc fold.cc  unfold.cc  seam_step.cc geometryc.cc
  delaunay.cc  delaunay_embed3d.cc  delaunay_alexandrov.cc  dense_linalg.cc  graph_fortran.cc isomerdb.cc layout.cc layout2d.cc argparser.cc
  polyhedron-optimize.cc symmetry.cc eisenstein.cc eisenstein_raster.cc eisenstein_atlas.cc barycentric.cc
  delaunay_strip.cc delaunay_unfold.cc delaunay_metric.cc eisenstein_paint.cc eisenstein_paint_geometry.cc wu_forcefield.cc
  deltahedron.cc deltahedron-io.cc deltahedron-cones.cc curvature-flow.cc buckinverse.cc ../contrib/mgmres.cc ${MORE_SOURCES}
)

#Ensure that CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is defined
if(NOT CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
  #Terminate
  message(FATAL_ERROR "CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is not defined")
endif()

find_package(OpenMP)
# BLAS/LAPACK deliberately NOT linked: the last LAPACK call in the C++ core
# (planargraph's dgetrf_) migrated to the in-house BLAS-free dense_linalg,
# and the Fortran solvers never used them.  See dense_linalg.hh for why the
# library must not silently regain a BLAS dependency.
#
# fortran_opt deliberately NOT linked either: the last C++ -> Fortran call
# (optimized_geometry -> sa_optff_) was replaced by the native Wu force
# field (wu_forcefield.cc); the legacy Fortran is now reachable only from
# the legacy fullerene_program, and the cross-validation parity tests
# (claude-projects/unfortran/tests) link the Fortran archives directly.
target_link_libraries(fullerenes buckygen)
if(OpenMP_CXX_FOUND)
  target_link_libraries(fullerenes OpenMP::OpenMP_CXX)
endif()



