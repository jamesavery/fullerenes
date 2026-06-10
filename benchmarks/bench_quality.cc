// Benchmark: CG vs L-BFGS vs Steihaug — triangle quality + Gaussian curvature
// Outputs summary stats to stdout, per-isomer CSV to file (if 4th arg given)
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/deltahedron.hh"
#include "fullerenes/buckinverse.hh"
#include <cstdio>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;
using namespace buckinverse;

struct Quality {
  double ang_relerr, edge_relerr, ang_std, h_min;
  double K_relerr, K_std, K_total;
  int n_concave, iters;
  double gmax_L, ms;
  bool conv;
  int n_energy, n_grad, n_hv;
};

Quality measure(Deltahedron& D, OptMethod method, const vector<coord3d>& init_pts, int budget_mult) {
  D.opt_method = method;
  D.opt_k_flat = 0;
  D.set_points(vector<coord3d>(init_pts));
  int Nv = D.N;
  // Convert iteration budget to work budget: each iter ~ N work (one gradient)
  long long max_work = (long long)budget_mult * Nv * Nv;

  auto t0 = chrono::high_resolution_clock::now();
  bool conv = (OptResult::CONVERGED == D.optimize(D.points, 0, 1e-10, {}, max_work));
  double ms = chrono::duration<double,milli>(chrono::high_resolution_clock::now()-t0).count();

  auto edges = D.undirected_edges();
  double L = 0;
  for(auto& e : edges) L += coord3d::dist(D.points[e.first], D.points[e.second]);
  L /= edges.size();

  double edge_err_sum = 0;
  for(auto& e : edges)
    edge_err_sum += fabs(coord3d::dist(D.points[e.first], D.points[e.second]) - L);
  double edge_relerr = edge_err_sum / (edges.size() * L);

  double target = M_PI / 3.0;
  double ang_err_sum = 0, ang_ref_sum = 0;
  vector<double> angles;
  vector<double> angle_sums(Nv, 0.0);

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
      angle_sums[v] += th;
    }
  }
  double ang_relerr = ang_err_sum / ang_ref_sum;
  double ang_mean = accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
  double ang_var = 0;
  for(double a : angles) ang_var += (a - ang_mean) * (a - ang_mean);
  double ang_std = sqrt(ang_var / angles.size());

  double K_err_sum = 0, K_ref_sum = 0, K_total = 0;
  vector<double> K_devs;
  for(int v = 0; v < Nv; v++){
    int deg = (int)D.degree(v);
    double K = 2.0 * M_PI - angle_sums[v];
    double K_target = 2.0 * M_PI - deg * M_PI / 3.0;
    K_total += K;
    double dev = K - K_target;
    K_devs.push_back(dev);
    K_err_sum += fabs(dev);
    K_ref_sum += fabs(K_target);
  }
  double K_relerr = (K_ref_sum > 1e-30) ? K_err_sum / K_ref_sum : K_err_sum;
  double K_dev_mean = accumulate(K_devs.begin(), K_devs.end(), 0.0) / K_devs.size();
  double K_dev_var = 0;
  for(double d : K_devs) K_dev_var += (d - K_dev_mean) * (d - K_dev_mean);
  double K_std_val = sqrt(K_dev_var / K_devs.size());

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

  return {ang_relerr, edge_relerr, ang_std, h_min,
          K_relerr, K_std_val, K_total,
          n_concave, D.iterations_used, D.final_gmax_L, ms, conv,
          D.n_energy_evals, D.n_grad_evals, D.n_hv_evals};
}

struct Accum {
  vector<double> ang_re, edge_re, ang_std_v, h_min_v, gmax_v, work_v;
  vector<double> K_re, K_std_v, K_total_v;
  vector<int> conc_v, iter_v;
  int n_conv = 0, n_total = 0;
  double total_ms = 0;

  void add(const Quality& q, int Nv) {
    ang_re.push_back(q.ang_relerr);
    edge_re.push_back(q.edge_relerr);
    ang_std_v.push_back(q.ang_std);
    h_min_v.push_back(q.h_min);
    gmax_v.push_back(q.gmax_L);
    conc_v.push_back(q.n_concave);
    iter_v.push_back(q.iters);
    K_re.push_back(q.K_relerr);
    K_std_v.push_back(q.K_std);
    K_total_v.push_back(q.K_total);
    long long w = (long long)q.n_energy + (long long)Nv * q.n_grad + (long long)Nv * q.n_hv;
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

void print_accum(const char* name, const Accum& a) {
  printf("%-8s: %d/%d converged (%.2f%%)\n", name, a.n_conv, a.n_total,
         100.0*a.n_conv/a.n_total);
  printf("  ang_re  : med=%.2e  p95=%.2e  max=%.2e\n",
         pct(a.ang_re, 0.5), pct(a.ang_re, 0.95), pct(a.ang_re, 1.0));
  printf("  edge_re : med=%.2e  p95=%.2e  max=%.2e\n",
         pct(a.edge_re, 0.5), pct(a.edge_re, 0.95), pct(a.edge_re, 1.0));
  printf("  ang_std : med=%.2e  p95=%.2e  max=%.2e  (degrees)\n",
         pct(a.ang_std_v, 0.5), pct(a.ang_std_v, 0.95), pct(a.ang_std_v, 1.0));
  printf("  K_re    : med=%.2e  p95=%.2e  max=%.2e  (Gaussian curv relerr)\n",
         pct(a.K_re, 0.5), pct(a.K_re, 0.95), pct(a.K_re, 1.0));
  printf("  K_std   : med=%.2e  p95=%.2e  max=%.2e  (K deviation std)\n",
         pct(a.K_std_v, 0.5), pct(a.K_std_v, 0.95), pct(a.K_std_v, 1.0));
  printf("  K_total : med=%.6f  min=%.6f  max=%.6f  (target=%.6f)\n",
         pct(a.K_total_v, 0.5), pct(a.K_total_v, 0.0), pct(a.K_total_v, 1.0), 4*M_PI);
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
  printf("  time    : total=%.1fs  mean=%.1fms\n", a.total_ms/1000, a.total_ms/a.n_total);
}

int main(int argc, char** argv) {
  int N = 70;
  int stride = 1;
  int budget_mult = 48;
  const char* csv_path = nullptr;
  int max_enum = 0;  // 0 = unlimited
  if(argc > 1) N = atoi(argv[1]);
  if(argc > 2) stride = atoi(argv[2]);
  if(argc > 3) budget_mult = atoi(argv[3]);
  if(argc > 4) csv_path = argv[4];
  if(argc > 5) max_enum = atoi(argv[5]);

  printf("Triangle quality benchmark: C%d, stride=%d, budget=%d*Nv", N, stride, budget_mult);
  if(max_enum > 0) printf(", first %d isomers", max_enum);
  printf("\n\n");

  FILE* csv = nullptr;
  if(csv_path) {
    csv = fopen(csv_path, "w");
    fprintf(csv, "N,idx,method,ang_re,edge_re,ang_std,K_re,K_std,h_min,concave,iters,work,gmax_L,ms,conv\n");
  }

  Accum acc_cg, acc_lb, acc_st;

  auto Q = BuckyGen::start(N, false);
  Graph G;
  int idx = 0, tested = 0;

  while(BuckyGen::next_fullerene(Q, G)) {
    if(idx % stride == 0) {
      ReducibleDual rd(G);
      ExtensionPath ep = rd.reduceToExtensionPath(20);

      Deltahedron D0 = Deltahedron::fromExtensionPathOptimized(ep);
      int Nv = D0.N;
      vector<coord3d> init(D0.points.begin(), D0.points.end());

      Deltahedron D1(D0), D2(D0), D3(D0);
      Quality q_cg = measure(D1, OptMethod::CG, init, budget_mult);
      Quality q_lb = measure(D2, OptMethod::LBFGS, init, budget_mult);
      Quality q_st = measure(D3, OptMethod::STEIHAUG, init, budget_mult);

      acc_cg.add(q_cg, Nv);
      acc_lb.add(q_lb, Nv);
      acc_st.add(q_st, Nv);

      if(csv) {
        auto emit = [&](const char* mname, const Quality& q) {
          long long w = (long long)q.n_energy + (long long)Nv * q.n_grad + (long long)Nv * q.n_hv;
          fprintf(csv, "%d,%d,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6f,%d,%d,%lld,%.6e,%.3f,%d\n",
                  N, idx, mname, q.ang_relerr, q.edge_relerr, q.ang_std,
                  q.K_relerr, q.K_std, q.h_min, q.n_concave,
                  q.iters, w, q.gmax_L, q.ms, q.conv ? 1 : 0);
        };
        emit("CG", q_cg);
        emit("LBFGS", q_lb);
        emit("ST", q_st);
        if(tested % 100 == 0) fflush(csv);
      }

      tested++;
      if(tested % 500 == 0)
        fprintf(stderr, "  %d isomers tested...\n", tested);
    }
    idx++;
    if(max_enum > 0 && idx >= max_enum) break;
  }
  BuckyGen::stop(Q);
  if(csv) fclose(csv);

  printf("Results for C%d (%d isomers):\n\n", N, tested);
  print_accum("CG", acc_cg);
  printf("\n");
  print_accum("L-BFGS", acc_lb);
  printf("\n");
  print_accum("Steihaug", acc_st);

  return 0;
}
