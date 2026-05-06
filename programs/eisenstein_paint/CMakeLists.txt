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
