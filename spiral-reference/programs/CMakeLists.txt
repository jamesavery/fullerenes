set(PROGRAMS gs-to-graph gs-to-polyhedron graph-to-gs gs-examples)

foreach(program ${PROGRAMS})
  set(program_source ${program}.cc)
  add_executable(${program} ${program_source})
  target_link_libraries(${program} libgraph ${libs})
endforeach(program ${PROGRAMS})

