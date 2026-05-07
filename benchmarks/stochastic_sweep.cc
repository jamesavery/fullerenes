// In-PALC stochastic-perturbation experiment (Direction 1).  For each
// (eps, seed) pair, run AlexandrovSolver with
// stochastic_perturbation_eps = eps applied between PALC steps.
//
// Reads spiral names from stdin (one per line), writes one row per
// (isomer, eps, seed) triple.
//
// Usage:
//   stochastic_sweep [--eps e1,e2,...] [--seeds N]
//
// Defaults: eps ∈ {0, 0.001, 0.005, 0.01, 0.05}, seeds = 5.

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
  vector<double> eps_list = {0.0, 0.001, 0.005, 0.01, 0.05};
  int n_seeds = 5;
  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    if (a == "--eps" && i+1 < argc) {
      eps_list.clear();
      string s = argv[++i];
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(',', pos);
        eps_list.push_back(stod(s.substr(pos, comma - pos)));
        pos = (comma == string::npos) ? s.size() : comma + 1;
      }
    } else if (a == "--seeds" && i+1 < argc) {
      n_seeds = atoi(argv[++i]);
    }
  }

  printf("# eps        seed   status                       vol_norm        kappa     flips  spiral\n");
  string name;
  while (getline(cin, name)) {
    if (name.empty() || name[0] == '#') continue;
    Triangulation T;
    DelaunayTriangulation D;
    try {
      spiral_nomenclature sn(name);
      T = Triangulation(sn);
      D = DelaunayTriangulation::compute(T);
    } catch (...) {
      printf("  -          -      %-28s -               -         -      %s\n",
             "PARSE_ERROR", name.c_str());
      continue;
    }

    for (double eps : eps_list) {
      int seeds_to_run = (eps == 0.0) ? 1 : n_seeds;
      for (int seed = 1; seed <= seeds_to_run; seed++) {
        AlexandrovSolver solver;
        solver.D = D;
        solver.stochastic_perturbation_eps = eps;
        solver.stochastic_seed = (uint32_t)seed;
        solver.solve();
        printf("  %-10.4f %-6d %-28s %-15.6e %-9.2e %-6d %s\n",
               eps, seed,
               AlexandrovSolver::status_str(solver.stats_status),
               solver.stats_volume_norm, solver.stats_final_kappa,
               solver.stats_flips, name.c_str());
      }
    }
    printf("\n");
  }
  return 0;
}
