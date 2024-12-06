if(ENABLE_CUDA)
  set(MORE_SOURCES fast_dualize.cc ) # For some reason, fast_dualize.cc includes GPU stuff. TODO: remove.
else()
  set(MORE_SOURCES isomer_batch.cc device_io.cc)  # Don't build CPU-only IsomerBatch if it needs to interact with CUDA.
endif()


add_library(fullerenes SHARED
  graph.cc cubicgraph.cc planargraph.cc fullerenegraph.cc triangulation.cc polyhedron.cc
  spiralc.cc planargraph-io.cc polyhedron-io.cc debug.cc fold.cc  unfold.cc geometryc.cc
  delaunay.cc  graph_fortran.cc isomerdb.cc layout.cc layout-optimize.cc argparser.cc
  polyhedron-optimize.cc symmetry.cc eisenstein.cc ../contrib/mgmres.cc ${MORE_SOURCES}
)

#Set up gfortran library (don't rely on the default)
find_library(GFORTRAN_LIB gfortran  PATHS /usr/lib/gcc/x86_64-linux-gnu/11)

if (NOT GFORTRAN_LIB)
  message(FATAL_ERROR "libgfortran not found")
endif()

target_link_libraries(fullerenes fortran_opt buckygen ${GFORTRAN_LIB})


