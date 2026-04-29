// Coefficient sweep for the strict-interior margin schedule (Task #28).
//
// For each input spiral name (one per line on stdin), and each margin
// coefficient c, run solve() with palc_interior_margin = c.  The
// in-loop margin schedule is c·t — vanishing at the homotopy limit.
// Tabulate (status, vol_norm, kappa) per (isomer, c).
//
// Usage:
//   echo SPIRAL_NAME | margin_sweep [--coeffs c1,c2,...]
//
// Defaults: c ∈ {0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1}.

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
  vector<double> coeffs = {0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1};
  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    if (a == "--coeffs" && i+1 < argc) {
      coeffs.clear();
      string s = argv[++i];
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(',', pos);
        coeffs.push_back(stod(s.substr(pos, comma - pos)));
        pos = (comma == string::npos) ? s.size() : comma + 1;
      }
    }
  }

  printf("# coeff      status                      vol_norm        kappa     flips  steps  spiral\n");
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
      printf("  -          %-27s %-15s %-9s %s\n",
             "PARSE_ERROR", "-", "-", name.c_str());
      continue;
    }

    for (double c : coeffs) {
      AlexandrovSolver solver;
      solver.D = D;
      solver.palc_interior_margin = c;
      solver.solve();
      printf("  %-10.4f %-27s %-15.6e %-9.2e %-6d %-6d %s\n",
             c,
             AlexandrovSolver::status_str(solver.stats_status),
             solver.stats_volume_norm, solver.stats_final_kappa,
             solver.stats_flips, solver.stats_steps,
             name.c_str());
    }
    printf("\n");
  }
  return 0;
}
