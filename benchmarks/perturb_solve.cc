// Symmetry-breaking perturbation experiment for the 14 known-failure
// fullerenes whose default-initial-r PALC converges to a degenerate
// (drum-cap or near-flat) fixed point.  The hypothesis (test of which
// is the purpose of this binary) is that those metrics admit a unique
// non-degenerate Alexandrov polytope (Alexandrov's theorem applies
// since fullerene metrics are non-degenerate), but uniform-r PALC
// preserves the metric's symmetry and lands in the symmetric flat
// basin.  A small per-vertex perturbation should break the symmetry
// and let PALC find the non-degenerate solution.
//
// Usage:
//   perturb_solve <spiral_name> [--seeds N] [--epsilons e1,e2,...]
//
// For each (seed, ε) pair, run:
//   r_init = 2·R_max · (1 + ε · u_v)   with u_v uniform in [-1, +1]^12
// then call solve() and report the outcome.  At least one combination
// converging to OK with vol_norm > 0.01 confirms the hypothesis.

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <spiral_name> "
                    "[--seeds N] [--epsilons e1,e2,...]\n", argv[0]);
    return 1;
  }
  string name = argv[1];
  int n_seeds = 5;
  vector<double> epsilons = {0.0, 0.005, 0.01, 0.05, 0.1};
  for (int i = 2; i < argc; i++) {
    string a = argv[i];
    if (a == "--seeds" && i+1 < argc) n_seeds = atoi(argv[++i]);
    else if (a == "--epsilons" && i+1 < argc) {
      epsilons.clear();
      string s = argv[++i];
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(',', pos);
        epsilons.push_back(stod(s.substr(pos, comma - pos)));
        pos = (comma == string::npos) ? s.size() : comma + 1;
      }
    }
  }

  spiral_nomenclature sn(name);
  Triangulation T(sn);
  auto D = DelaunayTriangulation::compute(T);

  // Default-initial radius scale (2·R_max).
  double R_max = 0;
  for (int h = 0; h < D.nh; h += 2)
    if (D.alive(h)) R_max = max(R_max, D.he_length[h]);
  double R0 = 2.0 * R_max;

  printf("# %s\n", name.c_str());
  printf("# nv=%d  R_max=%.4f  R0=%.4f\n", D.nv, R_max, R0);
  printf("# eps        seed   status                       vol_norm        kappa\n");

  // Three-mode perturbation:
  //   (a) random uniform [-ε, +ε] per-vertex
  //   (b) "split" — first 6 vertices ×(1+ε), other 6 ×(1−ε)
  //   (c) "single-vertex" — one vertex ×(1+ε), others unperturbed.
  //       Tries each vertex 0..11.
  printf("# eps     mode        seed   status                       vol_norm        kappa\n");

  auto run = [&](double eps, const char* mode, int seed,
                  const vector<double>& r0) {
    AlexandrovSolver solver;
    solver.D = D;
    solver.r_init_override = r0;
    solver.solve();
    printf("  %-7.4f %-11s %-6d %-28s %-15.6e %-12.3e\n",
           eps, mode, seed,
           AlexandrovSolver::status_str(solver.stats_status),
           solver.stats_volume_norm, solver.stats_final_kappa);
  };

  // Mode (d): scale-only sweep — vary R0 itself, no per-vertex
  // variation.  Tests whether different overall scale opens a
  // different basin even when symmetry is preserved.
  for (double scale : {0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0}) {
    vector<double> r0(D.nv, scale * R_max);
    run(scale, "scale", 0, r0);
  }

  for (double eps : epsilons) {
    if (eps == 0.0) {
      vector<double> r0(D.nv, R0);
      run(eps, "default", 0, r0);
      continue;
    }
    // Mode (a): random uniform.
    for (int seed = 0; seed < n_seeds; seed++) {
      vector<double> r0(D.nv, R0);
      mt19937 rng((uint32_t)(seed * 1000 + (int)(eps * 10000)));
      uniform_real_distribution<double> u(-1.0, 1.0);
      for (int v = 0; v < D.nv; v++) r0[v] = R0 * (1.0 + eps * u(rng));
      run(eps, "random", seed, r0);
    }
    // Mode (b): half-half split (use seed 0).
    {
      vector<double> r0(D.nv);
      for (int v = 0; v < D.nv; v++)
        r0[v] = R0 * (v < D.nv/2 ? (1.0 + eps) : (1.0 - eps));
      run(eps, "split", 0, r0);
    }
    // Mode (c): single-vertex bump on each vertex.
    for (int v0 = 0; v0 < D.nv; v0++) {
      vector<double> r0(D.nv, R0);
      r0[v0] = R0 * (1.0 + eps);
      run(eps, "single-v", v0, r0);
    }
  }
  return 0;
}
