// delaunay_metric.cc -- intrinsic surface metric on a metric delta-complex.
//
// Per-cone-pair geodesic distances and geodesics on the piecewise-flat surface
// represented by a DelaunayTriangulation (DCEL + edge lengths).  Priority-
// bounded BFS over triangle strips unfolded into the Eisenstein plane: exact in
// Z[w] for the simple distances, floating point only in the APSP-on-roots step.
// Mirror of TriangulationView's surface-metric methods (graphview.hh), but
// intrinsic -- the iDT IS the cone graph, so no node subset and no source mesh.
//
// Promoted from the delta-complex surface-metric project, validated on every
// fullerene dual C20-C160 (211,203,353 isomers, 0 failures).  See
// claude-projects/delta-complex/DELTA-COMPLEX-SURFACE-METRIC.md for the design.
// Cone-to-cone only: self-geodesics (loops based at one cone) are not computed
// -- the convex sector walk cannot wrap >pi around a cone to close a loop.

#include "fullerenes/delaunay.hh"
#include "fullerenes/eisenstein.hh"
#include "fullerenes/matrix.hh"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace {  // file-internal walk machinery

using simple_geodesic = DelaunayTriangulation::simple_geodesic;
using geodesic        = DelaunayTriangulation::geodesic;

// Squared length of half-edge h, as an exact integer (Loeschian on a valid iDT).
int Lsq_of(const DelaunayTriangulation& D, int h) {
  double L = D.he_length[h];
  return (int)std::llround(L * L);
}

// The three squared edge lengths of the triangle entered through h_in:
//   abs2 = |AB|^2 (h_in itself), acs2 = |AC|^2 (prev CCW edge), bcs2 = |BC|^2 (next).
struct FaceLengths { int abs2, acs2, bcs2; };
FaceLengths face_lengths_through(const DelaunayTriangulation& D, int h_in) {
  int h_next = D.he_next[h_in];
  int h_prev = D.he_next[h_next];
  return { Lsq_of(D, h_in), Lsq_of(D, h_prev), Lsq_of(D, h_next) };
}

// Parent-linked sequence of placed cone positions, for the cone-on-line check:
// a placement is invalid if a prior cone lies on the open segment origin -> C.
struct ConeChain { int parent; Eisenstein pos; };
bool chain_shadows(const std::vector<ConeChain>& chain, int idx, Eisenstein C) {
  while (idx >= 0) {
    if (is_on_open_segment(chain[idx].pos, C)) return true;
    idx = chain[idx].parent;
  }
  return false;
}

// BFS front: one strip face awaiting expansion.
struct Front {
  int        h_in;            // half-edge through which the face is entered
  Eisenstein A, B;            // positions of h_in's two corners in the unfolding
  Sector     sector;          // cumulative cone-sector from the seed cone
  int        chain_idx;       // index into the cone chain (parent linkage)
  long long  priority;        // = front_priority(A, B)
};
struct FrontGreater {
  bool operator()(const Front& x, const Front& y) const { return x.priority > y.priority; }
};
long long front_priority(Eisenstein A, Eisenstein B) {
  return std::min(A.norm2(), B.norm2());
}

// Per-seed walk result: squared distance + Eisenstein displacement to each cone.
struct WalkResult {
  enum class Code { Ok, SeedChiralityBarrier, StateCapHit };
  Code                    code = Code::Ok;
  std::vector<long long>  dist;     // squared distance per cone; LLONG_MAX if unreached
  std::vector<Eisenstein> disp;     // displacement per cone (valid iff dist < LLONG_MAX)
  long long               n_states = 0;
};

// Safety cap on processed BFS states.  A well-bounded walk never hits it.
constexpr long long kWalkStateCap = 2'000'000;

enum class StepKind { Continue, Drained, BoundedOut, CapHit, ApexChirality };

// One pop-and-expand BFS step; mutates pq, chain, R.dist/R.disp, R.n_states.
// Throws DegenerateTriangle (via place_third_eis) on a degenerate face.
StepKind step_walk(const DelaunayTriangulation& D, int chirality, long long bound_sq,
                   std::priority_queue<Front, std::vector<Front>, FrontGreater>& pq,
                   std::vector<ConeChain>& chain, WalkResult& R, long long state_cap) {
  if (pq.empty())                   return StepKind::Drained;
  if (pq.top().priority > bound_sq) return StepKind::BoundedOut;
  Front fr = pq.top(); pq.pop();
  if (++R.n_states > state_cap)     return StepKind::CapHit;

  FaceLengths fl = face_lengths_through(D, fr.h_in);
  std::optional<Eisenstein> apex =
      place_third_eis(fr.A, fr.B, fl.abs2, fl.acs2, fl.bcs2, chirality);
  if (!apex) return StepKind::ApexChirality;   // off-lattice for this chirality (soft prune)

  const Eisenstein C   = *apex;
  const long long  n_C = C.norm2();
  const int h_next = D.he_next[fr.h_in];
  const int h_prev = D.he_next[h_next];
  const int v_apex = D.he_origin[h_prev];

  if (fr.sector.contains(C) && !chain_shadows(chain, fr.chain_idx, C)
      && v_apex >= 0 && v_apex < (int)R.dist.size() && n_C < R.dist[v_apex]) {
    R.dist[v_apex] = n_C;
    R.disp[v_apex] = C;
  }

  const int new_chain_idx = (int)chain.size();
  chain.push_back({fr.chain_idx, C});

  auto push_child = [&](int twin, Eisenstein A_ch, Eisenstein B_ch) {
    if (!D.alive(twin) || D.he_face[twin] < 0) return;
    long long prio = front_priority(A_ch, B_ch);
    if (prio > bound_sq) return;
    std::optional<Sector> entry = entry_sector(A_ch, B_ch);
    if (!entry) return;
    Sector child = fr.sector.narrow_with(entry->R, entry->L);
    if (child.is_empty()) return;
    pq.push({twin, A_ch, B_ch, child, new_chain_idx, prio});
  };
  push_child(h_next ^ 1, C,    fr.B);
  push_child(h_prev ^ 1, fr.A, C);
  return StepKind::Continue;
}

// Geodesic distances from cone u_start via one seed config (h_seed, B_seed, chirality).
// Throws std::logic_error if B_seed does not realize the seed edge length, and
// DegenerateTriangle (from place_third_eis) on a degenerate face -- both deep
// invariants that cannot occur on a valid iDT.
WalkResult walk_from_seed(const DelaunayTriangulation& D, int u_start,
                          int h_seed, Eisenstein B_seed, int chirality,
                          long long bound_sq, long long state_cap = kWalkStateCap) {
  using Code = WalkResult::Code;
  WalkResult R;
  R.dist.assign(D.nv, LLONG_MAX);
  R.disp.assign(D.nv, Eisenstein(0, 0));

  if (B_seed.norm2() != Lsq_of(D, h_seed))
    throw std::logic_error(
        "DelaunayTriangulation surface metric: seed |B|^2 != Lsq(h_seed) at h="
        + std::to_string(h_seed));

  FaceLengths fl = face_lengths_through(D, h_seed);
  std::optional<Eisenstein> apex_seed =
      place_third_eis({0, 0}, B_seed, fl.abs2, fl.acs2, fl.bcs2, chirality);
  if (!apex_seed) { R.code = Code::SeedChiralityBarrier; return R; }

  R.dist[u_start] = 0;
  R.disp[u_start] = Eisenstein(0, 0);
  std::priority_queue<Front, std::vector<Front>, FrontGreater> pq;
  std::vector<ConeChain> chain;
  chain.push_back({-1, B_seed});
  pq.push({h_seed, {0, 0}, B_seed, ccw_order(B_seed, *apex_seed), /*chain_idx=*/0, /*priority=*/0});

  for (;;) {
    StepKind k = step_walk(D, chirality, bound_sq, pq, chain, R, state_cap);
    if (k == StepKind::ApexChirality) continue;                 // soft prune; keep going
    if (k == StepKind::CapHit) { R.code = Code::StateCapHit; break; }
    if (k == StepKind::Drained || k == StepKind::BoundedOut) break;
  }
  return R;
}

// The outgoing half-edges of cone u that bound a real face: u's angular sectors.
std::vector<int> fan_half_edges(const DelaunayTriangulation& D, int u) {
  std::vector<int> hs;
  const int h0 = D.v_out[u];
  if (h0 < 0) return hs;
  int h = h0;
  do { if (D.alive(h) && D.he_face[h] >= 0) hs.push_back(h); h = D.cw(h); } while (h != h0);
  return hs;
}

// Self-contained BFS bound: the iDT's edge-length-weighted graph diameter,
// squared, plus slack.  Surface distance <= graph distance, so a squared
// surface distance never exceeds this -- every cone pair is reachable.
long long metric_bound(const DelaunayTriangulation& D) {
  const int n = D.nv;
  matrix<double> W(n, n, std::numeric_limits<double>::infinity());
  for (int u = 0; u < n; u++) W(u, u) = 0.0;
  for (int h = 0; h < D.nh; h++) {
    if (!D.alive(h)) continue;
    int u = D.he_origin[h], v = D.he_origin[h ^ 1];
    if (u == v || u < 0 || v < 0 || u >= n || v >= n) continue;
    if (D.he_length[h] < W(u, v)) W(u, v) = D.he_length[h];
  }
  matrix<double> G = W.APSP(false);
  double dmax = 0.0;
  for (std::size_t i = 0; i < G.size(); i++)
    if (std::isfinite(G[i])) dmax = std::max(dmax, G[i]);
  return (long long)std::ceil(dmax * dmax) + 4;
}

struct SimpleMetric {
  matrix<long long>       square;   // square(u,v) = squared simple distance; LLONG_MAX if none
  matrix<simple_geodesic> geo;      // realizing simple geodesic (u's frame)
};

// For every ordered cone pair, the shortest simple geodesic over every seed
// config at u (each fan half-edge x each sector-0 rep of its squared length x
// both chiralities).  Throws on a deep-invariant failure -- non-Loeschian seed,
// degenerate face, or a BFS state-cap trip -- none of which occur on a valid iDT.
SimpleMetric compute_simple(const DelaunayTriangulation& D, long long bound_sq) {
  const int n = D.nv;
  SimpleMetric S{ matrix<long long>(n, n, LLONG_MAX),
                  matrix<simple_geodesic>(n, n, simple_geodesic{}) };
  for (int u = 0; u < n; u++) S.square(u, u) = 0;   // geo(u,u) stays trivial

  for (int u = 0; u < n; u++)
    for (int h : fan_half_edges(D, u)) {
      std::vector<Eisenstein> Bs = sector0_reps_of_norm(Lsq_of(D, h));
      if (Bs.empty())
        throw std::logic_error(
            "DelaunayTriangulation surface metric: non-Loeschian seed Lsq="
            + std::to_string(Lsq_of(D, h)) + " at h=" + std::to_string(h));
      for (Eisenstein B : Bs)
        for (int chi : {+1, -1}) {
          WalkResult R = walk_from_seed(D, u, h, B, chi, bound_sq);
          if (R.code == WalkResult::Code::StateCapHit)
            throw std::logic_error(
                "DelaunayTriangulation surface metric: BFS state cap hit at u="
                + std::to_string(u) + " h=" + std::to_string(h));
          for (int v = 0; v < n; v++)
            if (R.dist[v] < S.square(u, v)) {
              S.square(u, v) = R.dist[v];
              S.geo(u, v)    = simple_geodesic{ R.disp[v], h, B };   // stamp the winning seed
            }
        }
    }
  return S;
}

}  // anonymous namespace

// ============================================================================
// Public DelaunayTriangulation surface-metric methods
// ============================================================================

matrix<long long> DelaunayTriangulation::simple_square_surface_distances() const {
  return compute_simple(*this, metric_bound(*this)).square;
}

matrix<DelaunayTriangulation::simple_geodesic>
DelaunayTriangulation::simple_geodesics() const {
  return compute_simple(*this, metric_bound(*this)).geo;
}

DelaunayTriangulation::geodesic
DelaunayTriangulation::compose_simple_geodesics(const std::vector<int>& path,
                                                const matrix<simple_geodesic>& simple) {
  geodesic g;
  if (path.size() < 2) return g;
  g.segments.reserve(path.size() - 1);
  for (std::size_t k = 0; k + 1 < path.size(); k++)
    g.segments.push_back(simple(path[k], path[k + 1]));
  return g;
}

matrix<double>
DelaunayTriangulation::surface_distances(matrix<geodesic>* geodesics_out) const {
  SimpleMetric S = compute_simple(*this, metric_bound(*this));

  // sqrt of the integer simple distances (unreachable -> +inf), then APSP
  // smoothing to enforce the triangle inequality across intermediate cones.
  matrix<double> H(S.square.m, S.square.n, 0.0);
  for (std::size_t i = 0; i < S.square.size(); i++)
    H[i] = (S.square[i] < LLONG_MAX) ? std::sqrt((double)S.square[i])
                                     : std::numeric_limits<double>::infinity();

  if (geodesics_out == nullptr)
    return H.APSP(false).square_elementwise();

  const int n = S.square.m;
  APSPResult<double> apsp = H.APSP_with_paths();
  *geodesics_out = matrix<geodesic>(n, n, geodesic{});
  for (int u = 0; u < n; u++)
    for (int v = 0; v < n; v++)
      (*geodesics_out)(u, v) =
          compose_simple_geodesics(reconstruct_path(apsp.preds, u, v), S.geo);
  return apsp.dist.square_elementwise();
}

matrix<DelaunayTriangulation::geodesic>
DelaunayTriangulation::surface_geodesics() const {
  matrix<geodesic> G(0, 0, geodesic{});
  surface_distances(&G);
  return G;
}
