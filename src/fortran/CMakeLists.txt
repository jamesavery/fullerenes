enable_language(Fortran)

#add_library(fortran_spiral   STATIC spiral-standalone.f)
#add_library(fortran_hamilton STATIC hamilton-sa.f)


add_library(fortran_config config.f)
add_library(fortran_opt    opt-standalone.f geometry.f force.f dddihedral.f diag.f)
add_library(fortran_lib
  coord.f
  datain.f    
  hamilton.f
  hueckel.f
  isomer.f
  opt.f
  pentindex.f
  ring.f
  schlegel.f
  sphere.f
  spiral.f
  util.f
  volume.f
  )

add_executable(fullerene_program main.f)
target_link_libraries(fortran_opt fortran_config)
target_link_libraries(fortran_lib fortran_config)
target_link_libraries(fullerene_program fortran_config fortran_opt fortran_lib fullerenes)

set(CMAKE_Fortran_LINK_EXECUTABLE "")
set_property(TARGET fortran_config PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET fortran_lib    PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET fortran_opt    PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET fullerene_program PROPERTY LINKER_LANGUAGE Fortran)


install(TARGETS fullerene_program DESTINATION bin)
