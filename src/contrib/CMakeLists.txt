# The enumerator's vertex ceiling (dual vertices, so the largest fullerene
# is C(2*MAXN-4)): buckygen's own default is 152 (C300).  Its space is
# O(MAXN^2) static; the planar_code transfer to the wrapper caps it at 255.
set(BUCKYGEN_MAXN 152 CACHE STRING "buckygen's maximum dual vertex count (largest fullerene C(2n-4)); at most 255")
add_library(buckygen buckygen-wrapper.cc buckygen-1.0/buckygen.c)
set_source_files_properties(buckygen-1.0/buckygen.c PROPERTIES COMPILE_DEFINITIONS "MAIN=buckygen_main")
# On the whole target: buckygen-wrapper.cc includes buckygen.c and sizes its
# transfer buffers by MAXN, so both compilations must agree.
target_compile_definitions(buckygen PRIVATE MAXN=${BUCKYGEN_MAXN})
# Self-contained target for parent projects linking buckygen directly
# (buckygen-wrapper.hh lives under include/fullerenes/).
target_include_directories(buckygen PUBLIC ${Fullerenes_SOURCE_DIR}/include)
