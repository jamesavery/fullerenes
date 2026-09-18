// Full-isomer-space validation of the .geo format (GEO-FORMAT.md) on the cubic
// Alexandrov embeddings of all C_N isomers.
//
// write: every isomer's cubic polytope (realize_cubic: its 20..60 cone positions and
// the kappa = 0 iDT on them) is centred on its bounding box and appended, in buckygen
// order, to three archives -- 8-bit and 16-bit fixed point with a scale per record,
// and an f64 reference -- plus a sidecar of (isomer index, fingerprint of the exact
// record), 16 bytes per record. Resumable: the archives and the sidecar are cut back
// to their common prefix and the enumeration continues after it.
//
// check: reads every record of the three archives back and reports
//   - integrity: checksums, every record decoding, identical connectivity in all
//     three archives, and the f64 records matching their write-time fingerprints;
//   - coordinates: the error of u8/u16 against the f64 reference, and its ratio to
//     the per-record bound sigma/2 (which must never exceed 1);
//   - geometry: each record rebuilt as a DelaunayTriangulation with chord lengths.
//     The kappa = 0 cells are flat pieces of the polytope, so its cone-angle defects
//     are the cubic metric's p * pi/15 (12 degrees per pentagon at the cone). p is
//     read off the f64 record, whose defect error is the solver's own precision.
//
// Usage: geo_validate write|check|all N dir [--threads T] [--stop-after K]
//   --stop-after K  end the write phase once K isomers are enumerated (to test resuming)

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/eisenstein_paint_geometry.hh"
#include "fullerenes/geo-format.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/triangulation.hh"

#include <algorithm>
#include <bit>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

#define XXH_INLINE_ALL
#include "../src/contrib/xxhash/xxhash.h"

using namespace std;
using namespace eisenstein_paint;

namespace {

// ---------------------------------------------------------------------------
// Archives
// ---------------------------------------------------------------------------

enum { U8, U16, F64, N_ARCHIVES };
constexpr const char* archive_name[N_ARCHIVES] = {"u8", "u16", "f64"};

// N = 60: a fullerene has at most 60 pentagon vertices. Degrees 3..66 cover 3..59, the
// degrees of any simple triangulation on at most 60 vertices (realize_cubic refuses a
// non-simplicial iDT), so no isomer can exceed the capacities.
geo_options archive_options(int a) {
  geo_options o;
  o.type          = a == F64 ? geo_type::F64 : geo_type::FIXED;
  o.width         = a == U8 ? 8 : a == U16 ? 16 : 0;
  o.triangulation = true;
  o.record_n      = true;
  o.N             = 60;
  o.deg_min       = 3;
  o.deg_bits      = 6;
  return o;
}

struct SideEntry { uint64_t isomer, fingerprint; };

string path_of(const string& dir, int N, const string& what) {
  return dir + "/C" + to_string(N) + "." + what;
}

FILE* open_rw(const string& path) {
  const int fd = ::open(path.c_str(), O_RDWR | O_CREAT, 0666);
  FILE* f = fd >= 0 ? ::fdopen(fd, "r+b") : nullptr;
  if (!f) { perror(path.c_str()); exit(1); }
  return f;
}

uint64_t file_size(FILE* f) {
  fflush(f);
  struct stat st;
  fstat(fileno(f), &st);
  return uint64_t(st.st_size);
}

uint64_t count_of(FILE* f) { return file_size(f) == 0 ? 0 : geo::read_header(f).count; }

// Keep the first `keep` records: subtract the dropped records' hashes and rewrite count
// and checksum in one write. The dropped bytes become the ignored tail (sec. 3).
void drop_tail(FILE* f, uint64_t keep) {
  geo_header H = geo::read_header(f);
  if (H.count <= keep) return;
  const uint64_t seed = geo::header_seed(H);
  vector<uint8_t> bytes(H.record_size());
  for (uint64_t i = keep; i < H.count; i++) {
    fseeko(f, off_t(H.record_offset(i)), SEEK_SET);
    if (fread(bytes.data(), 1, bytes.size(), f) != bytes.size()) { fprintf(stderr, "drop_tail: short read\n"); exit(1); }
    H.checksum -= geo::record_hash(seed, i, bytes);
  }
  H.count = keep;
  const auto hb = geo::header_bytes(H);
  fseeko(f, 16, SEEK_SET);
  if (fwrite(hb.data() + 16, 1, 16, f) != 16 || fflush(f) != 0) { fprintf(stderr, "drop_tail: write failed\n"); exit(1); }
}

// ---------------------------------------------------------------------------
// Records
// ---------------------------------------------------------------------------

// D and x as the .geo record DelaunayTriangulation::to_geo writes, re-derived here
// (GEO-FORMAT.md sec. 7.4): rows counter-clockwise from v_out.
geo_record record_of(const DelaunayTriangulation& D, const vector<coord3d>& x) {
  geo_record r;
  r.n = D.nv;
  r.x = x;
  r.degree.assign(D.nv, 0);
  vector<int> slot(D.nh, -1), order;
  for (int v = 0; v < D.nv; v++) {
    int h = D.v_out[v];
    do { slot[h] = int(order.size()); order.push_back(h); r.degree[v]++; h = D.ccw(h); } while (h != D.v_out[v]);
  }
  r.twin.resize(order.size());
  for (size_t i = 0; i < order.size(); i++) r.twin[i] = slot[D.twin(order[i])];
  return r;
}

uint64_t fingerprint(const geo_record& r) {
  XXH3_state_t s;
  XXH3_INITSTATE(&s);
  XXH3_64bits_reset(&s);
  XXH3_64bits_update(&s, &r.n, sizeof r.n);
  XXH3_64bits_update(&s, r.x.data(), r.x.size() * sizeof(coord3d));
  XXH3_64bits_update(&s, r.degree.data(), r.degree.size() * sizeof(int));
  XXH3_64bits_update(&s, r.twin.data(), r.twin.size() * sizeof(int));
  return XXH3_64bits_digest(&s);
}

struct Realized {
  optional<DelaunayTriangulation> D;
  vector<coord3d> x;     // centred on the bounding box
  string error;
};

Realized realize(const Triangulation& T) {
  Realized out;
  try {
    CubicPolytope P = realize_cubic(T);
    const int n = int(P.cone_pos.size());
    P.D.compact_vertices();
    if (P.D.nv != n) throw runtime_error("iDT has " + to_string(P.D.nv) + " live vertices for " + to_string(n) + " cones");
    coord3d lo = P.cone_pos[0], hi = P.cone_pos[0];
    for (const auto& p : P.cone_pos)
      for (int k = 0; k < 3; k++) { lo[k] = min(lo[k], p[k]); hi[k] = max(hi[k], p[k]); }
    const coord3d mid = lo * 0.5 + hi * 0.5;
    for (const auto& p : P.cone_pos) out.x.push_back(p - mid);
    out.D = std::move(P.D);
  } catch (const exception& e) {
    ostringstream spiral;
    spiral << T.get_general_spiral();
    out.error = string(e.what()) + " | " + spiral.str();
  }
  return out;
}

// ---------------------------------------------------------------------------
// write
// ---------------------------------------------------------------------------

void write_phase(int N, const string& dir, uint64_t stop_after) {
  fflush(nullptr);
  auto Q = BuckyGen::start(N, false, false);     // forks: open the outputs afterwards

  FILE* archives[N_ARCHIVES];
  for (int a = 0; a < N_ARCHIVES; a++) archives[a] = open_rw(path_of(dir, N, string(archive_name[a]) + ".geo"));
  FILE* side = open_rw(path_of(dir, N, "sidecar"));
  FILE* failures = fopen(path_of(dir, N, "failures.txt").c_str(), "ab");

  // Resume: the common prefix of the three archives and the sidecar.
  uint64_t done = file_size(side) / sizeof(SideEntry);
  for (FILE* f : archives) done = min(done, count_of(f));
  for (FILE* f : archives) if (count_of(f) > done) drop_tail(f, done);
  if (ftruncate(fileno(side), off_t(done * sizeof(SideEntry))) != 0) { perror("sidecar"); exit(1); }
  uint64_t next_isomer = 0;
  if (done > 0) {
    SideEntry last;
    fseeko(side, off_t((done - 1) * sizeof(SideEntry)), SEEK_SET);
    if (fread(&last, sizeof last, 1, side) != 1) { fprintf(stderr, "sidecar: short read\n"); exit(1); }
    next_isomer = last.isomer + 1;
  }
  fseeko(side, 0, SEEK_END);

  const int64_t total = IsomerDB::number_isomers(N);
  fprintf(stderr, "C%d write: %lld isomers; resuming after %llu stored records (isomer %llu)\n",
          N, (long long)total, (unsigned long long)done, (unsigned long long)next_isomer);

  Triangulation T;
  uint64_t isomer = 0;
  while (isomer < next_isomer && BuckyGen::next_fullerene(Q, T)) isomer++;

  const auto t0 = chrono::steady_clock::now();
  uint64_t stored = done, failed = 0, realized_now = 0;
  int threads = 1;
#if defined(_OPENMP)
  threads = omp_get_max_threads();
#endif
  const size_t batch_size = 64 * size_t(threads);
  bool more = isomer == next_isomer;
  while (more) {
    vector<Triangulation> batch;
    while (batch.size() < batch_size && isomer + batch.size() < stop_after
           && (more = BuckyGen::next_fullerene(Q, T)))
      batch.push_back(T);
    if (batch.empty()) break;

    vector<Realized> results(batch.size());
    #pragma omp parallel for schedule(dynamic)
    for (size_t b = 0; b < batch.size(); b++) results[b] = realize(batch[b]);

    for (size_t b = 0; b < batch.size(); b++, isomer++) {
      const Realized& R = results[b];
      if (!R.D) {
        failed++;
        fprintf(failures, "C%d isomer %llu: %s\n", N, (unsigned long long)isomer, R.error.c_str());
        fflush(failures);
        continue;
      }
      for (int a = 0; a < N_ARCHIVES; a++)
        if (!DelaunayTriangulation::to_geo(*R.D, R.x, archives[a], true, archive_options(a))) {
          fprintf(stderr, "C%d isomer %llu: writing the %s archive failed\n", N, (unsigned long long)isomer,
                  archive_name[a]);
          exit(1);
        }
      const SideEntry e{isomer, fingerprint(record_of(*R.D, R.x))};
      if (fwrite(&e, sizeof e, 1, side) != 1) { perror("sidecar"); exit(1); }
      stored++;
    }
    fflush(side);
    realized_now += batch.size();
    const double s = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    fprintf(stderr, "C%d write: %llu/%lld isomers, %llu stored, %llu failed this run, %.0f/s\n", N,
            (unsigned long long)isomer, (long long)total, (unsigned long long)stored, (unsigned long long)failed,
            realized_now / s);
  }
  BuckyGen::stop(Q);

  const double seconds = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
  printf("C%d write: %llu isomers enumerated (%lld expected), %llu records stored, %llu failed this run; "
         "%.1f s for %llu isomers (%.0f/s, %d threads)\n",
         N, (unsigned long long)isomer, (long long)total, (unsigned long long)stored, (unsigned long long)failed,
         seconds, (unsigned long long)realized_now, realized_now / seconds, threads);
  if (isomer < stop_after && int64_t(isomer) != total) printf("C%d write: WARNING enumerated %llu isomers, expected %lld\n",
                                       N, (unsigned long long)isomer, (long long)total);
  for (FILE* f : archives) fclose(f);
  fclose(side);
  fclose(failures);
}

// ---------------------------------------------------------------------------
// check
// ---------------------------------------------------------------------------

struct RecordCheck {
  string error;                       // non-empty: a record failed to decode
  int n = 0, dmin = 0, dmax = 0;
  bool topology = false, lossless = false, p_valid = false;
  double err[N_ARCHIVES] = {};        // max |x - x_f64| (f64: 0)
  double ratio[N_ARCHIVES] = {};      // err / (sigma/2)
  double sigma[N_ARCHIVES] = {};
  double edge_rel[N_ARCHIVES] = {};   // max relative chord-length error
  double defect_err[N_ARCHIVES] = {}; // max |defect - p pi/15|
  string dt_error[N_ARCHIVES];        // chord-metric DelaunayTriangulation refused
};

float record_sigma(FILE* f, const geo_header& H, uint64_t i) {
  uint8_t b[4];
  fseeko(f, off_t(H.record_offset(i)), SEEK_SET);
  if (fread(b, 1, 4, f) != 4) return NAN;
  uint32_t v = uint32_t(b[0]) | uint32_t(b[1]) << 8 | uint32_t(b[2]) << 16 | uint32_t(b[3]) << 24;
  return std::bit_cast<float>(v);
}

vector<double> cone_defects(FILE* f, uint64_t i, const vector<coord3d>& x, string& error) {
  try {
    const DelaunayTriangulation D = DelaunayTriangulation::from_geo(f, i,
        [&](const DelaunayTriangulation& G, int h) { return (x[G.he_origin[h]] - x[G.dest(h)]).norm(); },
        [](int) { return 3; });
    vector<double> defect(D.nv);
    for (int v = 0; v < D.nv; v++) defect[v] = 2 * M_PI - D.vertex_angle_sum(v);
    return defect;
  } catch (const exception& e) {
    error = e.what();
    return {};
  }
}

RecordCheck check_record(FILE* files[N_ARCHIVES], const geo_header H[N_ARCHIVES], uint64_t i, uint64_t expected_fp) {
  RecordCheck c;
  geo_record r[N_ARCHIVES];
  try {
    for (int a = 0; a < N_ARCHIVES; a++) r[a] = geo::read_record(files[a], i);
  } catch (const exception& e) {
    c.error = e.what();
    return c;
  }
  const geo_record& ref = r[F64];
  c.n = ref.n;
  c.dmin = *min_element(ref.degree.begin(), ref.degree.end());
  c.dmax = *max_element(ref.degree.begin(), ref.degree.end());
  c.lossless = fingerprint(ref) == expected_fp;
  c.topology = true;
  for (int a : {U8, U16})
    c.topology = c.topology && r[a].n == ref.n && r[a].degree == ref.degree && r[a].twin == ref.twin;
  if (!c.topology) return c;

  // Coordinates against the reference, and chord lengths edge by edge.
  const vector<int> off = geo::row_offsets(ref.degree);
  vector<int> origin(off.back());
  for (int v = 0; v < ref.n; v++) fill(origin.begin() + off[v], origin.begin() + off[v + 1], v);
  for (int a : {U8, U16}) {
    c.sigma[a] = record_sigma(files[a], H[a], i);
    for (int v = 0; v < ref.n; v++)
      for (int k = 0; k < 3; k++) c.err[a] = max(c.err[a], fabs(r[a].x[v][k] - ref.x[v][k]));
    c.ratio[a] = c.err[a] / (c.sigma[a] / 2);
    for (int h = 0; h < off.back(); h++) {
      const int t = ref.twin[h];
      if (t < h) continue;
      const double L0 = (ref.x[origin[h]] - ref.x[origin[t]]).norm();
      const double L = (r[a].x[origin[h]] - r[a].x[origin[t]]).norm();
      c.edge_rel[a] = max(c.edge_rel[a], fabs(L - L0) / L0);
    }
  }

  // Cone-angle defects: multiples of pi/15, read off the reference.
  vector<double> defect[N_ARCHIVES];
  for (int a = 0; a < N_ARCHIVES; a++) defect[a] = cone_defects(files[a], i, r[a].x, c.dt_error[a]);
  if (defect[F64].empty()) return c;
  vector<int> p(ref.n);
  int p_sum = 0;
  c.p_valid = true;
  for (int v = 0; v < ref.n; v++) {
    p[v] = int(llround(defect[F64][v] * 15 / M_PI));
    p_sum += p[v];
    c.p_valid = c.p_valid && p[v] >= 1 && p[v] <= 3;
  }
  c.p_valid = c.p_valid && p_sum == 60;
  for (int a = 0; a < N_ARCHIVES; a++) {
    if (defect[a].empty()) continue;
    for (int v = 0; v < ref.n; v++) c.defect_err[a] = max(c.defect_err[a], fabs(defect[a][v] - p[v] * M_PI / 15));
  }
  return c;
}

struct Summary { double min = NAN, median = NAN, p99 = NAN, max = NAN; };

Summary summarize(vector<double> v) {
  Summary s;
  if (v.empty()) return s;
  sort(v.begin(), v.end());
  s.min = v.front();
  s.median = v[v.size() / 2];
  s.p99 = v[min(v.size() - 1, size_t(0.99 * v.size()))];
  s.max = v.back();
  return s;
}

void check_phase(int N, const string& dir) {
  const auto t0 = chrono::steady_clock::now();
  vector<string> paths;
  geo_header H[N_ARCHIVES];
  uint64_t sizes[N_ARCHIVES];
  bool verified[N_ARCHIVES];
  for (int a = 0; a < N_ARCHIVES; a++) {
    paths.push_back(path_of(dir, N, string(archive_name[a]) + ".geo"));
    FILE* f = fopen(paths[a].c_str(), "rb");
    if (!f) { perror(paths[a].c_str()); exit(1); }
    H[a] = geo::read_header(f);
    sizes[a] = file_size(f);
    verified[a] = geo::verify(f);
    fclose(f);
  }
  vector<SideEntry> side;
  {
    ifstream in(path_of(dir, N, "sidecar"), ios::binary);
    SideEntry e;
    while (in.read(reinterpret_cast<char*>(&e), sizeof e)) side.push_back(e);
  }
  const uint64_t count = H[F64].count;
  const bool counts_agree = H[U8].count == count && H[U16].count == count && side.size() == count;

  vector<RecordCheck> checks(counts_agree ? count : 0);
  #pragma omp parallel
  {
    FILE* files[N_ARCHIVES];
    for (int a = 0; a < N_ARCHIVES; a++) files[a] = fopen(paths[a].c_str(), "rb");
    #pragma omp for schedule(dynamic, 256)
    for (uint64_t i = 0; i < checks.size(); i++) checks[i] = check_record(files, H, i, side[i].fingerprint);
    for (FILE* f : files) fclose(f);
  }
  const double seconds = chrono::duration<double>(chrono::steady_clock::now() - t0).count();

  // Aggregate.
  uint64_t decode_failures = 0, topology_mismatch = 0, not_lossless = 0, p_invalid = 0, over_bound[N_ARCHIVES] = {};
  uint64_t dt_failures[N_ARCHIVES] = {};
  string first_error, first_dt_error[N_ARCHIVES];
  vector<double> n_values, err[N_ARCHIVES], ratio[N_ARCHIVES], sigma[N_ARCHIVES], edge[N_ARCHIVES], defect[N_ARCHIVES];
  int dmin = INT_MAX, dmax = 0;
  for (const RecordCheck& c : checks) {
    if (!c.error.empty()) { if (decode_failures++ == 0) first_error = c.error; continue; }
    n_values.push_back(c.n);
    dmin = min(dmin, c.dmin);
    dmax = max(dmax, c.dmax);
    topology_mismatch += !c.topology;
    not_lossless += !c.lossless;
    if (!c.topology) continue;
    p_invalid += !c.p_valid;
    for (int a = 0; a < N_ARCHIVES; a++) {
      if (!c.dt_error[a].empty()) { if (dt_failures[a]++ == 0) first_dt_error[a] = c.dt_error[a]; }
      else if (c.dt_error[F64].empty()) defect[a].push_back(c.defect_err[a]);
      if (a == F64) continue;
      err[a].push_back(c.err[a]);
      ratio[a].push_back(c.ratio[a]);
      sigma[a].push_back(c.sigma[a]);
      edge[a].push_back(c.edge_rel[a]);
      over_bound[a] += !(c.ratio[a] <= 1 + 1e-9);
    }
  }

  printf("\n=== C%d check (%.1f s) ===\n", N, seconds);
  printf("records: %llu u8, %llu u16, %llu f64, %zu sidecar entries -- %s\n", (unsigned long long)H[U8].count,
         (unsigned long long)H[U16].count, (unsigned long long)count, side.size(),
         counts_agree ? "agree" : "DISAGREE");
  const Summary ns = summarize(n_values);
  printf("cones per record: %.0f..%.0f (median %.0f); iDT degrees %d..%d (capacity %d..%llu)\n", ns.min, ns.max,
         ns.median, dmin, dmax, H[F64].opt.deg_min, (unsigned long long)H[F64].deg_max());
  printf("%-4s %8s %12s %9s %7s\n", "", "R (B)", "file (MB)", "B/isomer", "verify");
  for (int a = 0; a < N_ARCHIVES; a++) {
    const bool size_ok = sizes[a] == H[a].record_offset(H[a].count);
    printf("%-4s %8llu %12.2f %9.1f %7s%s\n", archive_name[a], (unsigned long long)H[a].record_size(),
           sizes[a] / 1e6, count ? double(sizes[a]) / count : 0.0, verified[a] ? "ok" : "FAILED",
           size_ok ? "" : "  (file size != 32 + count*R)");
  }
  printf("decode failures: %llu%s%s\n", (unsigned long long)decode_failures, first_error.empty() ? "" : " -- ",
         first_error.c_str());
  printf("connectivity identical in all three archives: %llu mismatches\n", (unsigned long long)topology_mismatch);
  printf("f64 records equal to the write-time originals (fingerprints): %llu mismatches\n",
         (unsigned long long)not_lossless);
  printf("\ncoordinates vs f64, per record (bond lengths):\n");
  printf("%-4s %10s %10s %10s %10s   %9s %9s %9s   %s\n", "", "err max", "median", "p99", "min", "sigma min",
         "median", "max", "worst err/(sigma/2)");
  for (int a : {U8, U16}) {
    const Summary e = summarize(err[a]), s = summarize(sigma[a]), r = summarize(ratio[a]);
    printf("%-4s %10.3e %10.3e %10.3e %10.3e   %9.3e %9.3e %9.3e   %.6f (%llu above 1)\n", archive_name[a], e.max,
           e.median, e.p99, e.min, s.min, s.median, s.max, r.max, (unsigned long long)over_bound[a]);
  }
  printf("\nrelative chord-length error, per record:\n");
  for (int a : {U8, U16}) {
    const Summary e = summarize(edge[a]);
    printf("%-4s max %.3e  median %.3e  p99 %.3e\n", archive_name[a], e.max, e.median, e.p99);
  }
  printf("\ncone-angle defects vs p*pi/15 from the chord-metric iDT (rad; step pi/15 = %.4f):\n", M_PI / 15);
  printf("pentagon counts read off f64 invalid (p outside 1..3 or sum != 60): %llu\n", (unsigned long long)p_invalid);
  for (int a : {F64, U16, U8}) {
    const Summary d = summarize(defect[a]);
    printf("%-4s max %.3e  median %.3e  p99 %.3e  | chord metric refused: %llu%s%s\n", archive_name[a], d.max,
           d.median, d.p99, (unsigned long long)dt_failures[a], first_dt_error[a].empty() ? "" : " -- ",
           first_dt_error[a].c_str());
  }
  fflush(stdout);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s write|check|all N dir [--threads T] [--stop-after K]\n", argv[0]);
    return 1;
  }
  const string mode = argv[1];
  const int N = atoi(argv[2]);
  const string dir = argv[3];
  uint64_t stop_after = UINT64_MAX;
  for (int a = 4; a < argc; a++) {
    if (string(argv[a]) == "--stop-after" && a + 1 < argc) {
      stop_after = strtoull(argv[++a], nullptr, 10);
    } else if (string(argv[a]) == "--threads" && a + 1 < argc) {
#if defined(_OPENMP)
      omp_set_num_threads(atoi(argv[++a]));
#else
      ++a;
      fprintf(stderr, "--threads ignored: built without OpenMP\n");
#endif
    } else {
      fprintf(stderr, "unknown option %s\n", argv[a]);
      return 1;
    }
  }
  if (mode != "write" && mode != "check" && mode != "all") { fprintf(stderr, "unknown mode %s\n", mode.c_str()); return 1; }
  filesystem::create_directories(dir);
  if (mode != "check") write_phase(N, dir, stop_after);
  if (mode != "write") check_phase(N, dir);
  return 0;
}
