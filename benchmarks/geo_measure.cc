// Measurement for the compact binary .geo format: per isomer, the raw quantities
// that decide its coordinate framing and graph capacities, for the two Alexandrov
// embeddings we want to archive -- the 12-cone dual polytope and the 20..60-cone
// cubic polytope.  One JSON line per isomer; geo_measure_report.py turns the lines
// into bytes-per-record vs. quantization-error tables.
//
// Per embedding: cone count n, live half-edges A, degree range, self-loops and
// parallel edges of the kappa=0 iDT, the T-bar arc count (dual only), and the
// extents that set a fixed-point scale -- Rc (max |x - centroid| per axis), and the
// bounding-box half-extents, descending, in the given frame (hb) and in the principal
// frame (hp).  Also the graph diameters of the dual and cubic graphs.
//
// Usage: geo_measure N out.jsonl [--chunks C --chunk-index I] [--max M]
//   --chunks / --chunk-index  buckygen's own splitting, to sample a large isomer space
//   --max M                   stop after M isomers of this chunk
// Resumes: complete lines already in out.jsonl are kept, and that many isomers skipped.

#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/delaunay_alexandrov.hh"
#include "fullerenes/eisenstein_paint_geometry.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/triangulation.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace eisenstein_paint;

static string json_escape(const string& s) {
  string r;
  for (char c : s) {
    if (c == '"' || c == '\\') { r += '\\'; r += c; }
    else if (c == '\n')        r += "\\n";
    else if ((unsigned char)c < 0x20) r += ' ';
    else r += c;
  }
  return r;
}

static int graph_diameter(const GraphView& G) {
  int diam = 0;
  vector<int> dist(G.N);
  vector<node_t> queue(G.N);
  for (node_t s = 0; s < G.N; s++) {
    fill(dist.begin(), dist.end(), -1);
    size_t head = 0, tail = 0;
    dist[s] = 0; queue[tail++] = s;
    while (head < tail) {
      node_t u = queue[head++];
      for (node_t v : G.nbrs(u))
        if (dist[v] < 0) { dist[v] = dist[u] + 1; queue[tail++] = v; }
    }
    if (tail != size_t(G.N)) throw runtime_error("graph_diameter: disconnected graph");
    diam = max(diam, dist[queue[tail - 1]]);
  }
  return diam;
}

// The JSON object body for one realized polytope: its iDT D (cones = vertices 0..n-1)
// and the cone positions x.
static string polytope_stats(const DelaunayTriangulation& D, const vector<coord3d>& x) {
  const int n = x.size();
  for (int v = 0; v < D.nv; v++)
    if ((D.v_out[v] >= 0) != (v < n))
      throw runtime_error("iDT live vertices are not exactly the cones 0.." + to_string(n - 1));

  vector<int> deg(n, 0);
  set<pair<int,int>> pairs;
  int A = 0, loops = 0, multi = 0;
  for (int h = 0; h < D.nh; h++) {
    if (!D.alive(h)) continue;
    A++;
    deg[D.he_origin[h]]++;
    if (h & 1) continue;
    int u = D.he_origin[h], v = D.dest(h);
    if (u == v) loops++;
    else if (!pairs.insert({min(u, v), max(u, v)}).second) multi++;
  }
  auto [dmin, dmax] = minmax_element(deg.begin(), deg.end());

  coord3d c(0, 0, 0);
  for (const auto& p : x) c += p;
  c /= double(n);
  double Rc = 0;
  coord3d lo = x[0], hi = x[0];
  matrix3d cov;
  for (const auto& p : x) {
    const coord3d d = p - c;
    for (int i = 0; i < 3; i++) {
      Rc = max(Rc, fabs(d[i]));
      lo[i] = min(lo[i], p[i]);
      hi[i] = max(hi[i], p[i]);
      for (int j = 0; j < 3; j++) cov(i, j) += d[i] * d[j];
    }
  }
  double hb[3];
  for (int i = 0; i < 3; i++) hb[i] = (hi[i] - lo[i]) / 2;
  sort(hb, hb + 3, greater<double>());

  const matrix3d C = cov.eigensystem().second;   // row k = k-th principal axis
  coord3d plo(INFINITY, INFINITY, INFINITY), phi(-INFINITY, -INFINITY, -INFINITY);
  for (const auto& p : x) {
    const coord3d q = C * (p - c);
    for (int i = 0; i < 3; i++) { plo[i] = min(plo[i], q[i]); phi[i] = max(phi[i], q[i]); }
  }
  double hp[3];
  for (int i = 0; i < 3; i++) hp[i] = (phi[i] - plo[i]) / 2;
  sort(hp, hp + 3, greater<double>());

  char buf[512];
  snprintf(buf, sizeof buf,
           "\"n\":%d,\"A\":%d,\"dmin\":%d,\"dmax\":%d,\"loops\":%d,\"multi\":%d,"
           "\"Rc\":%.9g,\"hb\":[%.9g,%.9g,%.9g],\"hp\":[%.9g,%.9g,%.9g]",
           n, A, *dmin, *dmax, loops, multi, Rc, hb[0], hb[1], hb[2], hp[0], hp[1], hp[2]);
  return buf;
}

static string measure_isomer(int N, long idx, const Triangulation& T) {
  ostringstream line;
  line << "{\"N\":" << N << ",\"i\":" << idx;
  bool failed = false;

  try {
    const DualPolytope P = realize_dual(sorted_dual(T));
    vector<int> labels(P.D.nv);
    for (int v = 0; v < P.D.nv; v++) labels[v] = v;
    const auto tbar = AlexandrovSolver::polytope_tesselation(P.D, P.r, labels);
    int tbar_arcs = 0;
    for (const auto& cell : tbar.cells) tbar_arcs += cell.size();
    line << ",\"dual\":{" << polytope_stats(P.D, P.cone_pos)
         << ",\"tbar_arcs\":" << tbar_arcs << ",\"diam\":" << graph_diameter(T) << "}";
  } catch (const exception& e) {
    failed = true;
    line << ",\"dual_err\":\"" << json_escape(e.what()) << "\"";
  }

  try {
    const CubicPolytope P = realize_cubic(T);
    line << ",\"cubic\":{" << polytope_stats(P.D, P.cone_pos)
         << ",\"diam\":" << graph_diameter(T.dual_graph()) << "}";
  } catch (const exception& e) {
    failed = true;
    line << ",\"cubic_err\":\"" << json_escape(e.what()) << "\"";
  }

  if (failed) {
    ostringstream rspi;
    rspi << T.get_general_spiral();
    line << ",\"rspi\":\"" << json_escape(rspi.str()) << "\"";
  }
  line << "}\n";
  return line.str();
}

// Keep the complete lines of a previous run (dropping a torn last line) and count them.
static long resume_count(const string& path) {
  if (!filesystem::exists(path)) return 0;
  ifstream in(path, ios::binary);
  const string text((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
  const size_t end = text.rfind('\n');
  const size_t keep = (end == string::npos) ? 0 : end + 1;
  if (keep != text.size()) filesystem::resize_file(path, keep);
  return count(text.begin(), text.begin() + keep, '\n');
}

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr, "Usage: %s N out.jsonl [--chunks C --chunk-index I] [--max M]\n", argv[0]);
    return 1;
  }
  const int N = atoi(argv[1]);
  const string out_path = argv[2];
  size_t chunks = 1, chunk_index = 0;
  long max_isomers = -1;
  for (int a = 3; a < argc; a++) {
    const string opt = argv[a];
    if (a + 1 >= argc) { fprintf(stderr, "%s: missing value\n", opt.c_str()); return 1; }
    if      (opt == "--chunks")      chunks      = strtoul(argv[++a], nullptr, 10);
    else if (opt == "--chunk-index") chunk_index = strtoul(argv[++a], nullptr, 10);
    else if (opt == "--max")         max_isomers = strtol(argv[++a], nullptr, 10);
    else { fprintf(stderr, "unknown option %s\n", opt.c_str()); return 1; }
  }
  if (chunks == 0 || chunk_index >= chunks) { fprintf(stderr, "need chunk-index < chunks\n"); return 1; }

  const long done = resume_count(out_path);
  if (max_isomers >= 0 && done >= max_isomers) {
    fprintf(stderr, "C%d: %ld isomers already in %s\n", N, done, out_path.c_str());
    return 0;
  }

  // Open the output only after buckygen has forked, so no stdio buffer is duplicated.
  fflush(nullptr);
  auto Q = BuckyGen::start(N, false, false, chunk_index, chunks);
  FILE* out = fopen(out_path.c_str(), "ab");
  if (!out) { perror(out_path.c_str()); BuckyGen::stop(Q); return 1; }

  const int CHUNK = 256;
  long idx = 0;
  Triangulation T;
  while (idx < done && BuckyGen::next_fullerene(Q, T)) idx++;

  bool more = (idx == done);
  while (more) {
    vector<Triangulation> batch;
    while (int(batch.size()) < CHUNK && (max_isomers < 0 || idx + long(batch.size()) < max_isomers)
           && (more = BuckyGen::next_fullerene(Q, T)))
      batch.push_back(T);
    if (batch.empty()) break;

    vector<string> lines(batch.size());
    #pragma omp parallel for schedule(dynamic)
    for (size_t b = 0; b < batch.size(); b++)
      lines[b] = measure_isomer(N, idx + b, batch[b]);

    for (const auto& l : lines) fwrite(l.data(), 1, l.size(), out);
    fflush(out);
    idx += batch.size();
    fprintf(stderr, "C%d: %ld isomers\n", N, idx);
    if (max_isomers >= 0 && idx >= max_isomers) break;
  }
  BuckyGen::stop(Q);
  fclose(out);
  return 0;
}
