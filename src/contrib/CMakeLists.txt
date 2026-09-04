add_library(buckygen buckygen-wrapper.cc buckygen-1.0/buckygen.c)
set_source_files_properties(buckygen-1.0/buckygen.c PROPERTIES COMPILE_DEFINITIONS "MAIN=buckygen_main")
# Self-contained target for parent projects linking buckygen directly
# (buckygen-wrapper.hh lives under include/fullerenes/).
target_include_directories(buckygen PUBLIC ${Fullerenes_SOURCE_DIR}/include)
