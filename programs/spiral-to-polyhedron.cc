// spiral-to-polyhedron: 3D embedding of a fullerene from its spiral.
//
// Supports two initial-geometry sources (Tutte/spherical zero-order vs. the
// Deltahedron extension-path optimizer) and two output representations
// (cubic fullerene polyhedron vs. dual triangulation).
//
// Dual <-> cubic conventions:
//   dual  -> cubic  : cubic vertices at triangle centroids, scaled by sqrt(3)
//                     so cubic edges match dual edges in the flat limit.
//   cubic -> dual   : dual vertices at face centers,        scaled by 1/sqrt(3).
// The product of the two scalings is 1, so repeated dual<->cubic conversions
// preserve overall scale.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/buckinverse.hh"

using namespace std;

enum OutputMode { OUT_CUBIC, OUT_DUAL };
enum InitGeom   { INIT_ZERO, INIT_DELTAHEDRON };

static void print_help(const char* prog)
{
  fprintf(stderr,
    "Usage: %s [options] <spiral> <output-file>\n"
    "\n"
    "Generate a 3D embedding of a fullerene from its spiral nomenclature.\n"
    "Output file format is determined by the filename extension\n"
    "(mol2, xyz, cc1, turbomole, gaussian, obj, ascii).\n"
    "\n"
    "Output representation:\n"
    "  --cubic           Cubic fullerene polyhedron (carbon cage) [default]\n"
    "  --dual            Dual triangulation (one vertex per 5/6-ring)\n"
    "\n"
    "Initial geometry:\n"
    "  --deltahedron     Extension-path Deltahedron optimizer on the dual\n"
    "                    (equilateral-triangle embedding; Alexandrov's theorem) [default]\n"
    "  --zero-order      Tutte layout + spherical projection\n"
    "\n"
    "Refinement:\n"
    "  --optimize        Run the cubic force-field optimizer after placement [default]\n"
    "  --no-optimize     Emit the raw initial geometry without force-field refinement\n"
    "                    (force-field refinement is only applied to cubic output)\n"
    "\n"
    "Other:\n"
    "  -h, --help        Show this help and exit\n"
    "\n"
    "Spiral syntax: <prefix>-[<spiral indices>]-<suffix>   (suffix: fullerene / cage / fulleroid)\n"
    "\n"
    "Examples:\n"
    "  %s 'C60-[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene' C60.mol2\n"
    "  %s --dual 'C60-[1,7,9,11,13,15,18,20,22,24,26,32]-fullerene' C60-dual.mol2\n"
    "  %s --zero-order --no-optimize 'C240-[1,2,3,4,5,6,121,122,237,238,239,240]-fullerene' C240-raw.mol2\n",
    prog, prog, prog, prog);
}

int main(int ac, char **av)
{
  OutputMode out_mode = OUT_CUBIC;
  InitGeom   init     = INIT_DELTAHEDRON;
  bool       do_opt   = true;
  vector<string> pos;

  for (int i = 1; i < ac; i++) {
    string a = av[i];
    if      (a == "-h" || a == "--help")       { print_help(av[0]); return 0; }
    else if (a == "--dual")                    out_mode = OUT_DUAL;
    else if (a == "--cubic" || a == "--fullerene") out_mode = OUT_CUBIC;
    else if (a == "--deltahedron")             init = INIT_DELTAHEDRON;
    else if (a == "--zero-order")              init = INIT_ZERO;
    else if (a == "--optimize")                do_opt = true;
    else if (a == "--no-optimize")             do_opt = false;
    else if (!a.empty() && a[0] == '-') {
      fprintf(stderr, "Unknown option: %s\n\n", a.c_str());
      print_help(av[0]);
      return 2;
    }
    else pos.push_back(a);
  }

  if (pos.size() != 2) {
    fprintf(stderr, "Expected <spiral> and <output-file>.\n\n");
    print_help(av[0]);
    return 1;
  }

  const string spiral_name = pos[0];
  const string out_path    = pos[1];

  const double sqrt3     = sqrt(3.0);
  const double inv_sqrt3 = 1.0 / sqrt3;

  spiral_nomenclature fsn(spiral_name);

  if (init == INIT_DELTAHEDRON) {
    // Reduce the dual triangulation to a seed extension path, then run the
    // per-step optimized pipeline to obtain an equilateral-triangle embedding.
    Triangulation T(fsn);
    buckinverse::ReducibleDual rd(T);
    buckinverse::ExtensionPath ep = rd.reduceToExtensionPath(20);
    if (ep.seed == buckinverse::SeedType::NotASeed) {
      fprintf(stderr, "Failed to reduce dual to a seed; cannot run the Deltahedron optimizer.\n");
      return 3;
    }

    Deltahedron D = Deltahedron::fromExtensionPathOptimized(ep);
    if (D.final_opt_result != OptResult::CONVERGED) {
      fprintf(stderr, "Warning: Deltahedron optimizer stopped with status %s "
                      "(angle_relerr=%.3e, n_concave=%d); emitting result anyway.\n",
              opt_result_name(D.final_opt_result),
              D.max_angle_relerr(), D.count_concave());
    }

    if (out_mode == OUT_DUAL) {
      // The deltahedron IS the dual polyhedron; write it directly.
      Polyhedron P(static_cast<const PlanarGraph&>(D),
                   vector<coord3d>(D.points.begin(), D.points.end()),
                   /*face_max=*/3);
      Polyhedron::to_file(P, out_path);
    } else {
      // Cubic from dual: place cubic vertices at triangle centroids, scaled by sqrt(3).
      // D.dual_graph() numbers cubic vertices in the same order as D.triangles().
      PlanarGraph cubic = D.dual_graph();
      auto tris = D.triangles();
      vector<coord3d> pts(cubic.N);
      const double s = sqrt3 / 3.0;  // centroid * sqrt(3)
      for (int u = 0; u < cubic.N; u++) {
        const auto& t = tris[u];
        pts[u] = (D.points[t[0]] + D.points[t[1]] + D.points[t[2]]) * s;
      }
      Polyhedron P(cubic, pts);
      if (do_opt) P.optimize();
      Polyhedron::to_file(P, out_path);
    }
  } else {
    // Zero-order initial geometry: Tutte layout + spherical projection on the cubic.
    FullereneGraph g(fsn);
    Polyhedron P(g, g.zero_order_geometry());
    if (do_opt) P.optimize();

    if (out_mode == OUT_CUBIC) {
      Polyhedron::to_file(P, out_path);
    } else {
      // Dual from cubic: face centers, scaled by 1/sqrt(3).
      Polyhedron D = P.dual();  // places dual vertices at face centers
      for (node_t u = 0; u < D.N; u++) D.points[u] *= inv_sqrt3;
      D.face_max = 3;
      Polyhedron::to_file(D, out_path);
    }
  }

  return 0;
}
