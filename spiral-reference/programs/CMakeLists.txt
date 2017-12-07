set(PROGRAMS gs-to-graph gs-to-polyhedron)

foreach(program ${PROGRAMS})
  set(program_source ${program}.cc)
  add_executable(${program} ${program_source})
  target_link_libraries(${program} libgraph)
endforeach(program ${PROGRAMS})


