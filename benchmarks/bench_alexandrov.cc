// Validation benchmark for the Alexandrov embedding solver.
// Enumerates all isomers of size N via buckygen, runs the solver on each,
// and reports pass/fail counts plus quality metrics (max edge-length
// error, angle-defect curvature, face-crossings, orientation).
//
// Usage: bench_alexandrov N

#include "fullerenes/delaunay_alexandrov12.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/spiral.hh"

#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>
#include <set>
#include <sstream>
#include <algorithm>
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;

static bool segment_triangle_intersect(
    const coord3d& p, const coord3d& q,
    const coord3d& v0, const coord3d& v1, const coord3d& v2) {
  coord3d e1 = v1 - v0, e2 = v2 - v0, d = q - p;
  coord3d h = d.cross(e2);
  double a = e1.dot(h);
  if (fabs(a) < 1e-15) return false;
  double f = 1.0 / a;
  coord3d s = p - v0;
  double u = f * s.dot(h);
  if (u < -1e-8 || u > 1+1e-8) return false;
  coord3d qv = s.cross(e1);
  double v = f * d.dot(qv);
  if (v < -1e-8 || u + v > 1+1e-8) return false;
  double t = f * e2.dot(qv);
  return (t > 1e-8 && t < 1 - 1e-8);
}

static int count_crossings(const DelaunayTriangulation& D, const vector<coord3d>& x) {
  struct F { int u,v,w; };
  vector<F> faces;
  for (int f = 0; f < D.nf; f++) {
    if (D.f_he[f] < 0) continue;
    int h0 = D.f_he[f], h1 = D.he_next[h0], h2 = D.he_next[h1];
    faces.push_back({D.he_origin[h0], D.he_origin[h1], D.he_origin[h2]});
  }
  int nf = faces.size(), count = 0;
  for (int i = 0; i < nf; i++)
    for (int j = i+1; j < nf; j++) {
      set<int> vi = {faces[i].u, faces[i].v, faces[i].w};
      if (vi.count(faces[j].u) || vi.count(faces[j].v) || vi.count(faces[j].w)) continue;
      bool cross = false;
      for (int a = 0; a < 3 && !cross; a++) {
        int fa[3] = {faces[i].u, faces[i].v, faces[i].w};
        cross = segment_triangle_intersect(x[fa[a]], x[fa[(a+1)%3]],
                  x[faces[j].u], x[faces[j].v], x[faces[j].w]);
      }
      for (int a = 0; a < 3 && !cross; a++) {
        int fb[3] = {faces[j].u, faces[j].v, faces[j].w};
        cross = segment_triangle_intersect(x[fb[a]], x[fb[(a+1)%3]],
                  x[faces[i].u], x[faces[i].v], x[faces[i].w]);
      }
      if (cross) count++;
    }
  return count;
}

static void scan(int N, double scale, AlexandrovSolver::Continuation method) {
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  int total = 0;
  int n_ok = 0, n_fail = 0, n_cross = 0, n_badori = 0;
  double worst_edge = 0, worst_curv = 0;
  long long total_steps = 0, total_newton = 0, total_flips = 0;
  vector<double> times_ms;

  // Record every failing isomer to a flushed file (stdout is unreliable here:
  // buckygen's fork() duplicates the stdout buffer, which both double-prints
  // and can drop mid-run lines). One ring-spiral per line; reconstructible.
  string fails_path = "bench_alexandrov_C" + to_string(N) + "_fails.txt";
  FILE* fails = nullptr;
  bool fails_open_failed = false;
  auto record_fail = [&](const char* kind, int idx, const string& rspi) {
    if (!fails && !fails_open_failed) {
      fails = fopen(fails_path.c_str(), "w");
      if (!fails) {
        fails_open_failed = true;
        fprintf(stderr, "warning: cannot open %s; failures not recorded\n",
                fails_path.c_str());
      }
    }
    if (fails) { fprintf(fails, "%s C%d #%d %s\n", kind, N, idx, rspi.c_str());
                 fflush(fails); }
  };

  // Solve in parallel chunks: pull CHUNK triangulations serially from
  // buckygen, then solve them in parallel via OpenMP.  Per-isomer times
  // are gathered into a per-chunk array and appended after each parallel
  // region, so they are preserved for regression checks.
  const int CHUNK = 1024;
  auto t_wall_start = chrono::high_resolution_clock::now();

  while (true) {
    vector<Triangulation> batch;
    batch.reserve(CHUNK);
    Triangulation T;
    while ((int)batch.size() < CHUNK && BuckyGen::next_fullerene(Q, T))
      batch.push_back(T);
    if (batch.empty()) break;
    int nb = batch.size();
    vector<double> chunk_times(nb, 0.0);

    #if defined(_OPENMP)
    #pragma omp parallel for schedule(dynamic) \
      reduction(+:n_ok,n_fail,n_cross,n_badori,total_steps,total_newton,total_flips) \
      reduction(max:worst_edge,worst_curv)
    #endif
    for (int b = 0; b < nb; b++) {
      auto D = DelaunayTriangulation::compute(batch[b]);

      // Scale-invariance probe: multiply the metric by `scale` before solving.
      // The solver should give scale-invariant pass/fail + metrics.
      if (scale != 1.0)
        for (int h = 0; h < D.nh; h++)
          if (D.alive(h)) D.he_length[h] *= scale;

      AlexandrovSolver solver;
      solver.D = D;
      solver.continuation = method;
      auto t0 = chrono::high_resolution_clock::now();
      auto coords = solver.solve();
      auto t1 = chrono::high_resolution_clock::now();
      chunk_times[b] = chrono::duration<double, milli>(t1 - t0).count();

      total_steps  += solver.stats_steps;
      total_newton += solver.stats_newton_total;
      total_flips  += solver.stats_flips;

      int idx = total + b;

      if (coords.empty()) {
        n_fail++;
        std::ostringstream oss;
        oss << batch[b].get_general_spiral();
        #pragma omp critical
        record_fail("FAIL", idx, oss.str());
        continue;
      }

      double max_edge = 0;
      for (int h = 0; h < solver.D.nh; h += 2) {
        if (!solver.D.alive(h)) continue;
        int u = solver.D.he_origin[h], v = solver.D.dest(h);
        double target = solver.D.he_length[h];
        double actual = (coords[u] - coords[v]).norm();
        max_edge = max(max_edge, fabs(actual - target) / target);
      }
      worst_edge = max(worst_edge, max_edge);

      int cx = count_crossings(solver.D, coords);
      if (cx > 0) n_cross++;

      double vol = 0;
      for (int f = 0; f < solver.D.nf; f++) {
        if (solver.D.f_he[f] < 0) continue;
        int h0 = solver.D.f_he[f], h1 = solver.D.he_next[h0], h2 = solver.D.he_next[h1];
        vol += coords[solver.D.he_origin[h0]].dot(
          coords[solver.D.he_origin[h1]].cross(coords[solver.D.he_origin[h2]]));
      }
      if (vol < 0) n_badori++;

      double max_curv = 0;
      for (int v = 0; v < solver.D.nv; v++) {
        if (solver.D.v_out[v] < 0) continue;
        double angle_sum = 0;
        int h0 = solver.D.v_out[v], h = h0;
        do {
          int a = solver.D.dest(h);
          int hcw = solver.D.cw(h);
          int b2 = solver.D.dest(hcw);
          coord3d ea = coords[a] - coords[v], eb = coords[b2] - coords[v];
          double na = ea.norm(), nb = eb.norm();
          if (na > 1e-15 && nb > 1e-15)
            angle_sum += acos(clamp(ea.dot(eb)/(na*nb), -1.0, 1.0));
          h = hcw;
        } while (h != h0);
        max_curv = max(max_curv, fabs(2*M_PI - angle_sum - M_PI/3));
      }
      worst_curv = max(worst_curv, max_curv);

      bool pass = cx == 0 && vol > 0 && max_curv < 0.01;
      if (pass) n_ok++;
      else {
        std::ostringstream oss;
        oss << "cross=" << cx << " vol=" << vol << " curv=" << max_curv
            << " edge=" << max_edge*100 << "% rspi=" << batch[b].get_general_spiral();
        #pragma omp critical
        record_fail("GFAIL", idx, oss.str());
      }
    }

    times_ms.insert(times_ms.end(), chunk_times.begin(), chunk_times.end());
    total += nb;
    fprintf(stderr, "C%d: %d scanned\n", N, total);
    fflush(stdout);
  }
  BuckyGen::stop(Q);
  if (fails) fclose(fails);

  auto t_wall_end = chrono::high_resolution_clock::now();
  double wall_s = chrono::duration<double>(t_wall_end - t_wall_start).count();

  sort(times_ms.begin(), times_ms.end());
  int nt = times_ms.size();
  double total_ms = 0; for (double t : times_ms) total_ms += t;

  printf("SUMMARY C%d: %d isomers, %d pass, %d conv_fail, %d cross, %d badori\n",
         N, total, n_ok, n_fail, n_cross, n_badori);
  printf("  worst_edge=%.4f%% worst_curv=%.6f rad\n", worst_edge*100, worst_curv);
  #if defined(_OPENMP)
  printf("  wall: %.1fs (%d threads, %.1f isomers/s)\n",
         wall_s, omp_get_max_threads(), total / wall_s);
  #else
  printf("  wall: %.1fs (1 thread, %.1f isomers/s)\n",
         wall_s, total / wall_s);
  #endif
  printf("  per-isomer cpu: total=%.1fs  mean=%.2fms  median=%.2fms  p99=%.2fms  max=%.2fms\n",
         total_ms/1000, total_ms/total,
         nt > 0 ? times_ms[nt/2] : 0,
         nt > 0 ? times_ms[min(nt-1, (int)(nt*0.99))] : 0,
         nt > 0 ? times_ms[nt-1] : 0);
  printf("  steps: total=%lld  mean=%.1f  newton=%lld  flips=%lld\n",
         total_steps, (double)total_steps/total, total_newton, total_flips);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s N [scale] [natural|palc]\n", argv[0]);
    return 1;
  }
  double scale = (argc >= 3) ? atof(argv[2]) : 1.0;
  auto method = AlexandrovSolver::Continuation::NATURAL;
  if (argc >= 4) {
    string m = argv[3];
    if (m == "palc") method = AlexandrovSolver::Continuation::PALC;
    else if (m != "natural")
      fprintf(stderr, "warning: unknown method '%s'; using natural\n", m.c_str());
  }
  printf("method: %s, scale: %g\n",
         method == AlexandrovSolver::Continuation::PALC ? "palc" : "natural",
         scale);
  fflush(stdout);  // flush before scan()'s buckygen fork() duplicates the buffer
  scan(atoi(argv[1]), scale, method);
  return 0;
}
