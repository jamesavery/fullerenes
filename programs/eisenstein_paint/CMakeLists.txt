# programs/eisenstein_paint -- CLI tools driving the Eisenstein-paint
# pipeline: full-isomer-space sweep validator + TikZ visualisation
# drivers for the per-cell + per-arc + iDT-unfolding scenes.

# Binaries are prefixed eisenstein_paint_* to avoid collision with
# the rest of programs/.

add_executable(eisenstein_paint_validate         validate.cc)
add_executable(eisenstein_paint_tikz_cell        tikz_cell.cc)
add_executable(eisenstein_paint_tikz_paint       tikz_paint.cc)
add_executable(eisenstein_paint_tikz_strip       tikz_strip.cc)
add_executable(eisenstein_paint_tikz_iDT_unfold  tikz_iDT_unfold.cc)
add_executable(eisenstein_paint_tikz_conflict    tikz_conflict.cc)

# validate.cc sweeps isomer space under "#pragma omp parallel for", so it needs
# the OpenMP flag on its OWN command line: libfullerenes links OpenMP PRIVATEly
# (see src/c++/CMakeLists.txt) and so does not hand it down.  Without this the
# sweep silently degrades to one thread.
if(OpenMP_CXX_FOUND)
  target_link_libraries(eisenstein_paint_validate OpenMP::OpenMP_CXX)
endif()
