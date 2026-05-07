// Full-population AlexandrovSolver sweep over a range of carbon counts.
// Enumerates ALL fullerene isomers in [N_min, N_max] via buckygen
// (general, not IPR-only), runs AlexandrovSolver on each, and writes
// per-failure TSV rows + a per-N OK/fail summary.
//
// OpenMP-parallel within each N (each thread holds a private buckygen
// queue chunk).  Line-buffered output for kill-9 safety.
//
// Usage:
//   scan_full_range <N_min> <N_max> --out failures.tsv --summary summary.txt
//                   [--ipr]
//
// Output (failures.tsv) — one row per non-OK isomer:
//   N  status  vol_norm  kappa  flips  steps  buckygen_idx  spiral_name
// Output (summary.txt) — per-N row:
//   N  total  ok  fail_kappa  fail_not_simple  fail_volume  fail_other
//                                                            fail_build

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <atomic>
#include <chrono>
#include <mutex>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using ValidationStatus = AlexandrovSolver::ValidationStatus;

struct PerNCount {
  int total = 0, ok = 0;
  int fail_kappa = 0, fail_not_simple = 0, fail_volume = 0;
  int fail_other = 0, fail_build = 0;
};

static const char* status_short(ValidationStatus s) {
  switch (s) {
    case ValidationStatus::OK:                       return "OK";
    case ValidationStatus::FAIL_KAPPA_NOT_CONVERGED: return "FAIL_KAPPA";
    case ValidationStatus::FAIL_NOT_SIMPLE:          return "FAIL_NOT_SIMPLE";
    case ValidationStatus::FAIL_RECONSTRUCT:         return "FAIL_RECONSTRUCT";
    case ValidationStatus::FAIL_VOLUME_DEGENERATE:   return "FAIL_VOLUME_DEGENERATE";
    case ValidationStatus::FAIL_SELF_INTERSECTING:   return "FAIL_SELF_INTERSECTING";
    case ValidationStatus::FAIL_NOT_CONVEX:          return "FAIL_NOT_CONVEX";
  }
  return "UNKNOWN";
}

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr,
      "Usage: %s <N_min> <N_max> --out failures.tsv --summary summary.txt [--ipr]\n",
      argv[0]);
    return 1;
  }
  int N_min = atoi(argv[1]);
  int N_max = atoi(argv[2]);
  string out_path, summary_path;
  bool ipr = false;
  for (int i = 3; i < argc; i++) {
    string a = argv[i];
    if      (a == "--out"     && i+1 < argc) out_path     = argv[++i];
    else if (a == "--summary" && i+1 < argc) summary_path = argv[++i];
    else if (a == "--ipr")                   ipr          = true;
  }
  if (out_path.empty() || summary_path.empty()) {
    fprintf(stderr, "--out and --summary required\n");
    return 1;
  }

  FILE* fp_fail = fopen(out_path.c_str(), "w");
  FILE* fp_sum  = fopen(summary_path.c_str(), "w");
  if (!fp_fail || !fp_sum) {
    fprintf(stderr, "cannot open output file(s)\n");
    return 1;
  }
  setvbuf(fp_fail, nullptr, _IOLBF, 0);
  setvbuf(fp_sum,  nullptr, _IOLBF, 0);
  fprintf(fp_fail, "#N\tstatus\tvol_norm\tkappa\tflips\tsteps\tbuckygen_idx\tspiral_name\n");
  fprintf(fp_sum,  "#N\ttotal\tok\tfail_kappa\tfail_not_simple\tfail_volume\tfail_other\tfail_build\n");

  std::mutex out_mtx;

  for (int N = N_min; N <= N_max; N += 2) {
    if (N < 20 || N == 22) continue;        // no fullerenes
    auto t0 = chrono::steady_clock::now();
    PerNCount cnt;
    std::atomic<int> idx{0};

    // Enumerate buckygen serially into a vector of dual-Triangulations,
    // then process them in parallel.  Buckygen is single-threaded; the
    // solver dominates per-isomer cost so this gives near-linear speedup
    // on the parallel pass.
    auto Q = BuckyGen::start(N, ipr);
    vector<Triangulation> all;
    Graph G;
    while (BuckyGen::next_fullerene(Q, G)) {
      all.emplace_back(Triangulation(G));
    }
    BuckyGen::stop(Q);
    cnt.total = (int)all.size();

    #pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < (int)all.size(); i++) {
      // Compute canonical name and run solver.
      AlexandrovSolver solver;
      string spiral_name;
      try {
        solver.D = DelaunayTriangulation::compute(all[i]);
        // Canonical spiral name (cheap relative to solve).
        spiral_nomenclature sn(all[i],
                                spiral_nomenclature::FULLERENE,
                                spiral_nomenclature::TRIANGULATION,
                                /*rarest_special_start=*/true);
        spiral_name = "C" + to_string(N) + "-" + sn.to_string();
      } catch (...) {
        std::lock_guard<std::mutex> lk(out_mtx);
        cnt.fail_build++;
        fprintf(fp_fail, "%d\tFAIL_BUILD\t-\t-\t-\t-\t%d\t-\n", N, i);
        continue;
      }

      solver.solve();

      ValidationStatus s = solver.stats_status;
      {
        std::lock_guard<std::mutex> lk(out_mtx);
        switch (s) {
          case ValidationStatus::OK:                       cnt.ok++; break;
          case ValidationStatus::FAIL_KAPPA_NOT_CONVERGED: cnt.fail_kappa++; break;
          case ValidationStatus::FAIL_NOT_SIMPLE:          cnt.fail_not_simple++; break;
          case ValidationStatus::FAIL_VOLUME_DEGENERATE:   cnt.fail_volume++; break;
          default:                                          cnt.fail_other++; break;
        }
        if (s != ValidationStatus::OK) {
          fprintf(fp_fail, "%d\t%s\t%.6e\t%.6e\t%d\t%d\t%d\t%s\n",
                  N, status_short(s),
                  solver.stats_volume_norm,
                  solver.stats_final_kappa,
                  solver.stats_flips,
                  solver.stats_steps,
                  i, spiral_name.c_str());
        }
      }
    }

    fprintf(fp_sum, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            N, cnt.total, cnt.ok,
            cnt.fail_kappa, cnt.fail_not_simple, cnt.fail_volume,
            cnt.fail_other, cnt.fail_build);
    auto dt = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    fprintf(stderr, "N=%-3d  total=%-7d  ok=%-7d  fail=%d (kappa=%d simpl=%d vol=%d other=%d build=%d)  %.1fs\n",
            N, cnt.total, cnt.ok,
            cnt.total - cnt.ok,
            cnt.fail_kappa, cnt.fail_not_simple, cnt.fail_volume,
            cnt.fail_other, cnt.fail_build, dt);
  }

  fclose(fp_fail);
  fclose(fp_sum);
  return 0;
}
