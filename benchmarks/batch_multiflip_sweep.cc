// Synchronized batch multi-flip experiment (Direction 2).  Sweep
// palc_batch_multiflip_threshold; report outcome per (isomer, threshold).
//
// Usage: stdin = spiral names, optional --thresholds t1,t2,... (defaults
// 0, 0.001, 0.005, 0.01, 0.05, 0.1).

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  vector<double> thresholds = {0.0, 0.001, 0.005, 0.01, 0.05, 0.1};
  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    if (a == "--thresholds" && i+1 < argc) {
      thresholds.clear();
      string s = argv[++i];
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(',', pos);
        thresholds.push_back(stod(s.substr(pos, comma - pos)));
        pos = (comma == string::npos) ? s.size() : comma + 1;
      }
    }
  }

  printf("# threshold   status                       vol_norm        kappa     flips  spiral\n");
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
      printf("  -          %-28s -               -         -      %s\n",
             "PARSE_ERROR", name.c_str());
      continue;
    }
    for (double th : thresholds) {
      AlexandrovSolver solver;
      solver.D = D;
      solver.palc_batch_multiflip_threshold = th;
      solver.solve();
      printf("  %-10.4f %-28s %-15.6e %-9.2e %-6d %s\n",
             th,
             AlexandrovSolver::status_str(solver.stats_status),
             solver.stats_volume_norm, solver.stats_final_kappa,
             solver.stats_flips, name.c_str());
    }
    printf("\n");
  }
  return 0;
}
