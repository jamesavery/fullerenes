// Two-pass deflated-homotopy experiment (Direction 4 form (a)).
//
// Pass 1: run solve() normally; if the validator detects a near-degenerate
// fixed point (FAIL_NOT_SIMPLE or FAIL_VOLUME_DEGENERATE), record
// r* = solver.r as the empirical drum-cap point.
// Pass 2: re-run solve() with deflation toward r*, sweeping `strength`.
//
// For each (isomer, strength) pair report the post-pass-2 status,
// vol_norm, and kappa.
//
// Usage: stdin = spiral names, optional --strengths s1,s2,...
// (defaults: 0.001, 0.01, 0.1, 1.0, 10.0).

#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
  vector<double> strengths = {0.001, 0.01, 0.1, 1.0, 10.0};
  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    if (a == "--strengths" && i+1 < argc) {
      strengths.clear();
      string s = argv[++i];
      size_t pos = 0;
      while (pos < s.size()) {
        size_t comma = s.find(',', pos);
        strengths.push_back(stod(s.substr(pos, comma - pos)));
        pos = (comma == string::npos) ? s.size() : comma + 1;
      }
    }
  }

  printf("# pass strength   status                       vol_norm        kappa     spiral\n");
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
      printf("  -    -          %-28s -               -         %s\n",
             "PARSE_ERROR", name.c_str());
      continue;
    }

    // Pass 1: normal solve.
    AlexandrovSolver pass1;
    pass1.D = D;
    pass1.solve();
    printf("  1    %-10s %-28s %-15.6e %-9.2e %s\n", "-",
           AlexandrovSolver::status_str(pass1.stats_status),
           pass1.stats_volume_norm, pass1.stats_final_kappa, name.c_str());

    // If pass 1 succeeded, no deflation needed.
    if (pass1.stats_status == AlexandrovSolver::ValidationStatus::OK) {
      continue;
    }
    // Capture r* from pass 1.
    vector<double> r_star = pass1.r;

    // Pass 2: deflation sweep.
    for (double strength : strengths) {
      AlexandrovSolver pass2;
      pass2.D = D;
      pass2.r_deflate_target = r_star;
      pass2.deflate_strength = strength;
      pass2.solve();
      printf("  2    %-10.4g %-28s %-15.6e %-9.2e %s\n", strength,
             AlexandrovSolver::status_str(pass2.stats_status),
             pass2.stats_volume_norm, pass2.stats_final_kappa, name.c_str());
    }
    printf("\n");
  }
  return 0;
}
