add_library(buckygen buckygen-wrapper.cc buckygen-1.0/buckygen.c)
set_source_files_properties(buckygen-1.0/buckygen.c PROPERTIES COMPILE_DEFINITIONS "MAIN=buckygen_main")
