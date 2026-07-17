enable_language(Fortran)

#add_library(fortran_spiral   STATIC spiral-standalone.f)
#add_library(fortran_hamilton STATIC hamilton-sa.f)
# find_package(OpenMP REQUIRED)

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
target_link_libraries(fortran_opt PRIVATE fortran_config)
target_link_libraries(fortran_lib PRIVATE fortran_config)
# Static-archive order matters now that libfullerenes no longer carries the
# fortran_opt objects: fortran_lib (spiral.f, opt.f) pulls tred2l/tqlil and
# the geometry/force kernels from fortran_opt, so fortran_opt must follow it.
target_link_libraries(fullerene_program fortran_lib fortran_opt fortran_config fullerenes ${OpenMP_CXX_LIBRARIES})

# When using Apple Clang for C++ and Homebrew gfortran for Fortran, the Fortran
# linker needs the macOS SDK library path to find libc++ (propagated from the
# C++ shared library).
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
  execute_process(COMMAND xcrun --show-sdk-path OUTPUT_VARIABLE _sdk_path OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_sdk_path)
    target_link_options(fullerene_program PRIVATE "-L${_sdk_path}/usr/lib" "-F${_sdk_path}/System/Library/Frameworks")
  endif()
endif()

set_property(TARGET fortran_config PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET fortran_lib    PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET fortran_opt    PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET fullerene_program PROPERTY LINKER_LANGUAGE Fortran)


install(TARGETS fullerene_program DESTINATION bin)
