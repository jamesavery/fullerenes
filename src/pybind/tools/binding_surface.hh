#pragma once
// binding_surface.hh — the single translation unit the binding generator dumps
// with `clang -Xclang -ast-dump=json`. It includes every library header whose
// classes appear in tools/spec.toml, so the AST contains the full method surface
// of the view hierarchy the bindings harvest from.

#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/isomerdb.hh"
