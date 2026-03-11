// Benchmark: full pipeline with method + angle-tolerance variations
// Tests method x step_angle_tol configurations per isomer, run in parallel.
// Primary convergence metric: max per-angle relative error + convexity.
//
// Usage: bench_quality_pipeline N [options]
//   --methods LBFGS,ST           Comma-separated methods (default: LBFGS,ST)
//   --step-tols 1e-1,1e-2,1e-3  Step angle tolerances (comma-separated)
//   --final-tol 1e-3             Final angle tolerance
//   --stride 1                   Buckygen stride (default: 1 = all isomers)
//   --max-enum 0                 Max buckygen index (0 = all)
//   --csv PATH                   Output CSV path (required)
//   --work-factor 20             Budget = factor * Nv^3 (default: 20)
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/buckinverse.hh"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <thread>
#include <mutex>
#include <sstream>

using namespace std;
using namespace buckinverse;

struct Config {
  OptMethod method;
  OptMethod final_method;
  double step_tol;
  const char* method_name;
  const char* final_method_name;
};

struct Quality {
  double ang_max_re, ang_mean_re, edge_relerr, ang_std;
  double h_min;
  int n_concave, iters;
  double gmax_L, ms;
  bool conv;
  int n_energy, n_grad, n_hv;
};

Quality measure_pipeline(const ExtensionPath& ep, const Config& cfg,
                         double final_tol, long long work_budget) {
  auto t0 = chrono::high_resolution_clock::now();

  // Use unreachable gmax*L tolerance so angle_tol dominates convergence
  double grad_step_tol = 1e-15;
  double grad_final_tol = 1e-15;

  Deltahedron D = Deltahedron::fromExtensionPathOptimized(
      ep, nullptr, nullptr,
      cfg.method, grad_step_tol, grad_final_tol, work_budget,
      cfg.step_tol, final_tol, cfg.final_method);
  double ms = chrono::duration<double,milli>(chrono::high_resolution_clock::now()-t0).count();

  int Nv = D.N;

  // Edge relative error (mean)
  auto edges = D.undirected_edges();
  double L = 0;
  for(auto& e : edges) L += coord3d::dist(D.points[e.first], D.points[e.second]);
  L /= edges.size();
  double edge_err_sum = 0;
  for(auto& e : edges)
    edge_err_sum += fabs(coord3d::dist(D.points[e.first], D.points[e.second]) - L);
  double edge_relerr = edge_err_sum / (edges.size() * L);

  // Angle statistics
  double target = M_PI / 3.0;
  double ang_err_sum = 0, ang_ref_sum = 0;
  vector<double> angles;
  for(int v = 0; v < Nv; v++){
    int deg = (int)D.degree(v);
    for(int j = 0; j < deg; j++){
      int a = D.nbrs(v)[j];
      int b = D.nbrs(v)[(j+1)%deg];
      coord3d ea = D.points[a] - D.points[v];
      coord3d eb = D.points[b] - D.points[v];
      double cos_th = ea.dot(eb) / (ea.norm() * eb.norm());
      cos_th = max(-1.0, min(1.0, cos_th));
      double th = acos(cos_th);
      angles.push_back(th * 180.0 / M_PI);
      ang_err_sum += fabs(th - target);
      ang_ref_sum += target;
    }
  }
  double ang_mean_re = ang_err_sum / ang_ref_sum;
  double ang_mean_deg = accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
  double ang_var = 0;
  for(double a : angles) ang_var += (a - ang_mean_deg) * (a - ang_mean_deg);
  double ang_std = sqrt(ang_var / angles.size());

  // h_min (from final reflect state)
  double h_min = 1e30;
  int n_concave = 0;
  for(int v = 0; v < Nv; v++){
    int deg = (int)D.degree(v);
    coord3d centroid(0,0,0);
    for(int nb : D.nbrs(v)) centroid = centroid + D.points[nb];
    centroid = centroid * (1.0 / deg);
    coord3d fan_normal(0,0,0);
    for(int j = 0; j < deg; j++){
      coord3d e1 = D.points[D.nbrs(v)[j]] - D.points[v];
      coord3d e2 = D.points[D.nbrs(v)[(j+1)%deg]] - D.points[v];
      fan_normal = fan_normal + e1.cross(e2);
    }
    double nn = fan_normal.norm();
    if(nn < 1e-15) continue;
    coord3d nhat = fan_normal / nn;
    double h = (D.points[v] - centroid).dot(nhat);
    h_min = min(h_min, h);
    if(h < 0) n_concave++;
  }

  // Convergence: angle-based
  bool conv = (D.final_angle_relerr < final_tol && D.final_n_concave == 0);

  return {D.final_angle_relerr, ang_mean_re, edge_relerr, ang_std,
          h_min, n_concave, D.iterations_used, D.final_gmax_L, ms, conv,
          D.n_energy_evals, D.n_grad_evals, D.n_hv_evals};
}

// Parse comma-separated doubles
vector<double> parse_doubles(const char* s) {
  vector<double> out;
  istringstream iss(s);
  string tok;
  while(getline(iss, tok, ',')) out.push_back(atof(tok.c_str()));
  return out;
}

// Parse comma-separated method names
vector<pair<OptMethod,const char*>> parse_methods(const char* s) {
  vector<pair<OptMethod,const char*>> out;
  istringstream iss(s);
  string tok;
  while(getline(iss, tok, ',')) {
    if(tok == "CG")    out.push_back({OptMethod::CG, "CG"});
    else if(tok == "LBFGS") out.push_back({OptMethod::LBFGS, "LBFGS"});
    else if(tok == "ST")    out.push_back({OptMethod::STEIHAUG, "ST"});
    else { fprintf(stderr, "Unknown method: %s\n", tok.c_str()); exit(1); }
  }
  return out;
}

struct Accum {
  vector<double> ang_max_re, ang_mean_re, edge_re, ang_std_v, h_min_v, gmax_v, work_v, ms_v;
  vector<int> conc_v, iter_v;
  int n_conv = 0, n_total = 0;
  double total_ms = 0;

  void add(const Quality& q, int Nv) {
    ang_max_re.push_back(q.ang_max_re);
    ang_mean_re.push_back(q.ang_mean_re);
    edge_re.push_back(q.edge_relerr);
    ang_std_v.push_back(q.ang_std);
    h_min_v.push_back(q.h_min);
    gmax_v.push_back(q.gmax_L);
    conc_v.push_back(q.n_concave);
    iter_v.push_back(q.iters);
    ms_v.push_back(q.ms);
    long long w = (long long)q.n_energy + 2LL * q.n_grad + 7LL * q.n_hv;
    work_v.push_back((double)w);
    if(q.conv) n_conv++;
    n_total++;
    total_ms += q.ms;
  }
};

template<typename T>
T pct(vector<T> v, double p) {
  sort(v.begin(), v.end());
  return v[min((int)(v.size()*p), (int)v.size()-1)];
}

void print_accum(const char* method, double step_tol, const Accum& a) {
  printf("%-6s step_tol=%.0e: %d/%d converged (%.1f%%)\n", method, step_tol,
         a.n_conv, a.n_total, 100.0*a.n_conv/max(1, a.n_total));
  printf("  ang_max_re: med=%.2e  p95=%.2e  max=%.2e\n",
         pct(a.ang_max_re, 0.5), pct(a.ang_max_re, 0.95), pct(a.ang_max_re, 1.0));
  printf("  ang_mean_re: med=%.2e  p95=%.2e  max=%.2e\n",
         pct(a.ang_mean_re, 0.5), pct(a.ang_mean_re, 0.95), pct(a.ang_mean_re, 1.0));
  printf("  edge_re : med=%.2e  p95=%.2e  max=%.2e\n",
         pct(a.edge_re, 0.5), pct(a.edge_re, 0.95), pct(a.edge_re, 1.0));
  printf("  ang_std : med=%.2e  p95=%.2e  max=%.2e  (degrees)\n",
         pct(a.ang_std_v, 0.5), pct(a.ang_std_v, 0.95), pct(a.ang_std_v, 1.0));
  printf("  h_min   : med=%.4f  p05=%.4f  min=%.4f\n",
         pct(a.h_min_v, 0.5), pct(a.h_min_v, 0.05), pct(a.h_min_v, 0.0));
  printf("  concave : med=%d  max=%d\n",
         pct(a.conc_v, 0.5), pct(a.conc_v, 1.0));
  printf("  iters   : med=%d  p95=%d  max=%d\n",
         pct(a.iter_v, 0.5), pct(a.iter_v, 0.95), pct(a.iter_v, 1.0));
  printf("  work    : med=%.0f  p95=%.0f  max=%.0f\n",
         pct(a.work_v, 0.5), pct(a.work_v, 0.95), pct(a.work_v, 1.0));
  printf("  gmax*L  : med=%.2e  p95=%.2e  max=%.2e\n",
         pct(a.gmax_v, 0.5), pct(a.gmax_v, 0.95), pct(a.gmax_v, 1.0));
  printf("  time    : med=%.1fms  p95=%.1fms  total=%.1fs\n",
         pct(a.ms_v, 0.5), pct(a.ms_v, 0.95), a.total_ms/1000);
}

int main(int argc, char** argv) {
  if(argc < 2) {
    fprintf(stderr, "Usage: bench_quality_pipeline N [options]\n");
    fprintf(stderr, "  --methods LBFGS,ST     Step methods (default: LBFGS,ST)\n");
    fprintf(stderr, "  --final-method CG      Final optimize method (default: same as step)\n");
    fprintf(stderr, "  --step-tols 1e-1,...   Step angle tolerances\n");
    fprintf(stderr, "  --final-tol 1e-3       Final angle tolerance\n");
    fprintf(stderr, "  --stride 1             Buckygen stride\n");
    fprintf(stderr, "  --max-enum 0           Max index (0=all)\n");
    fprintf(stderr, "  --csv PATH             Output CSV (required)\n");
    fprintf(stderr, "  --work-factor 20       Budget = factor * Nv^3\n");
    return 1;
  }

  int N = atoi(argv[1]);

  // Defaults
  const char* methods_str = "LBFGS,ST";
  const char* final_method_str = nullptr;  // nullptr = same as step method
  const char* step_tols_str = "1e-1,1e-2,1e-3";
  double final_tol = 1e-3;
  int stride = 1;
  int max_enum = 0;
  const char* csv_path = nullptr;
  int work_factor = 400;

  // Parse options (stride=0 treated as 1)
  for(int i = 2; i < argc; i++) {
    if(!strcmp(argv[i], "--methods") && i+1 < argc) methods_str = argv[++i];
    else if(!strcmp(argv[i], "--final-method") && i+1 < argc) final_method_str = argv[++i];
    else if(!strcmp(argv[i], "--step-tols") && i+1 < argc) step_tols_str = argv[++i];
    else if(!strcmp(argv[i], "--final-tol") && i+1 < argc) final_tol = atof(argv[++i]);
    else if(!strcmp(argv[i], "--stride") && i+1 < argc) stride = atoi(argv[++i]);
    else if(!strcmp(argv[i], "--max-enum") && i+1 < argc) max_enum = atoi(argv[++i]);
    else if(!strcmp(argv[i], "--csv") && i+1 < argc) csv_path = argv[++i];
    else if(!strcmp(argv[i], "--work-factor") && i+1 < argc) work_factor = atoi(argv[++i]);
    else { fprintf(stderr, "Unknown option: %s\n", argv[i]); return 1; }
  }

  if(!csv_path) { fprintf(stderr, "Error: --csv is required\n"); return 1; }
  if(stride <= 0) stride = 1;

  auto methods = parse_methods(methods_str);
  auto step_tols = parse_doubles(step_tols_str);

  // Parse final method override (if any)
  pair<OptMethod,const char*> final_method_override = {OptMethod::CG, nullptr};
  if(final_method_str) {
    auto fm = parse_methods(final_method_str);
    if(!fm.empty()) final_method_override = fm[0];
  }

  // Build configs: cross product of methods x step_tols
  vector<Config> configs;
  for(auto& [m, name] : methods)
    for(double st : step_tols) {
      OptMethod fm = final_method_override.second ? final_method_override.first : m;
      const char* fm_name = final_method_override.second ? final_method_override.second : name;
      configs.push_back({m, fm, st, name, fm_name});
    }

  int n_configs = (int)configs.size();
  int Nv = N/2 + 2;
  long long work_budget = (long long)work_factor * Nv * Nv;

  printf("Pipeline benchmark: C%d (Nv=%d), stride=%d, final_tol=%.0e, work=%d*Nv^2=%lld\n",
         N, Nv, stride, final_tol, work_factor, work_budget);
  printf("Methods: %s, step_tols: %s, %d configs\n", methods_str, step_tols_str, n_configs);
  if(max_enum > 0) printf("Max enum: %d\n", max_enum);
  printf("\n");

  FILE* csv = fopen(csv_path, "w");
  if(!csv) { fprintf(stderr, "Cannot open %s\n", csv_path); return 1; }
  fprintf(csv, "N,idx,method,step_tol,ang_max_re,ang_mean_re,edge_re,ang_std,h_min,concave,iters,work,gmax_L,ms,conv\n");

  vector<Accum> accums(n_configs);

  auto Q = BuckyGen::start(N, false);
  Graph G;
  int idx = 0, tested = 0;

  while(BuckyGen::next_fullerene(Q, G)) {
    if(idx % stride == 0) {
      ReducibleDual rd(G);
      ExtensionPath ep = rd.reduceToExtensionPath(20);

      vector<Quality> results(n_configs);
      vector<thread> threads(n_configs);

      for(int c = 0; c < n_configs; c++) {
        threads[c] = thread([&, c]() {
          results[c] = measure_pipeline(ep, configs[c], final_tol, work_budget);
        });
      }
      for(int c = 0; c < n_configs; c++) threads[c].join();

      for(int c = 0; c < n_configs; c++) {
        accums[c].add(results[c], Nv);
        const Quality& q = results[c];
        long long w = (long long)q.n_energy + 2LL * q.n_grad + 7LL * q.n_hv;
        fprintf(csv, "%d,%d,%s,%.0e,%.6e,%.6e,%.6e,%.6e,%.6f,%d,%d,%lld,%.6e,%.3f,%d\n",
                N, idx, configs[c].method_name, configs[c].step_tol,
                q.ang_max_re, q.ang_mean_re, q.edge_relerr, q.ang_std,
                q.h_min, q.n_concave,
                q.iters, w, q.gmax_L, q.ms, q.conv ? 1 : 0);
      }
      fflush(csv);

      tested++;
      if(tested % 500 == 0)
        fprintf(stderr, "  %d isomers tested (idx=%d)...\n", tested, idx);
    }
    idx++;
    if(max_enum > 0 && idx >= max_enum) break;
  }
  BuckyGen::stop(Q);
  fclose(csv);

  printf("\nResults for C%d (%d isomers, %d configs each):\n", N, tested, n_configs);
  printf("Work budget: %d*Nv^2 = %lld per optimize() call\n", work_factor, work_budget);
  printf("Final angle tol: %.0e\n\n", final_tol);

  for(int c = 0; c < n_configs; c++) {
    char label[64];
    if(strcmp(configs[c].method_name, configs[c].final_method_name) != 0)
      snprintf(label, sizeof(label), "%s+%s", configs[c].method_name, configs[c].final_method_name);
    else
      snprintf(label, sizeof(label), "%s", configs[c].method_name);
    print_accum(label, configs[c].step_tol, accums[c]);
    printf("\n");
  }

  return 0;
}
