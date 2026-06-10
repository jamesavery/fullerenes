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
// Default scope: cone-to-cone.  Pass calculate_self_geodesics = true to also
// compute self-geodesics (closed loops based at one cone).  Self mode lifts
// the R.dist[u_start] pin, records seed-edge self-loops at seed setup, and
// bypasses the convex-sector containment test for u_start apex placements
// at recording time so the wrap-around closure flavour is also captured.

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
// `target_label`: -1 for cone-to-cone; the u_start label in self mode (the
// apex's containment test relaxes when v_apex == target_label so wrap-around
// closures are recorded; the chain-shadow test is *not* relaxed -- a closure
// whose line passes through a prior chain cone is a composed loop, picked up
// by APSP).
// Throws DegenerateTriangle (via place_third_eis) on a degenerate face.
StepKind step_walk(const DelaunayTriangulation& D, int chirality, long long bound_sq,
                   int target_label,
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

  const bool is_closure = (target_label >= 0 && v_apex == target_label && n_C > 0);
  const bool sector_ok  = is_closure || fr.sector.contains(C);
  if (sector_ok && !chain_shadows(chain, fr.chain_idx, C)
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
// `target_label`: -1 for cone-to-cone (R.dist[u_start] is pinned to 0); u_start
// for self mode (the pin is lifted, the seed-edge self-loop is recorded at seed
// setup, and step_walk relaxes the sector test for u_start apex placements).
// Throws std::logic_error if B_seed does not realize the seed edge length, and
// DegenerateTriangle (from place_third_eis) on a degenerate face -- both deep
// invariants that cannot occur on a valid iDT.
WalkResult walk_from_seed(const DelaunayTriangulation& D, int u_start,
                          int h_seed, Eisenstein B_seed, int chirality,
                          long long bound_sq, int target_label = -1,
                          long long state_cap = kWalkStateCap) {
  using Code = WalkResult::Code;
  WalkResult R;
  R.dist.assign(D.nv, LLONG_MAX);
  R.disp.assign(D.nv, Eisenstein(0, 0));

  if (B_seed.norm2() != Lsq_of(D, h_seed))
    throw std::logic_error(
        "DelaunayTriangulation surface metric: seed |B|^2 != Lsq(h_seed) at h="
        + std::to_string(h_seed));

  // Seed-time self-loop record (self mode): run BEFORE the seed-apex check so
  // a chirality-barrier early return doesn't discard a valid seed-edge closure
  // -- B_seed is a Loeschian representative of the seed-edge length regardless
  // of whether the adjacent face's third vertex lands on the lattice.
  if (target_label == u_start && D.dest(h_seed) == u_start) {
    R.dist[u_start] = B_seed.norm2();
    R.disp[u_start] = B_seed;
  } else if (target_label != u_start) {
    // Cone-to-cone API pin (the BFS would otherwise record valid self-geodesics).
    R.dist[u_start] = 0;
    R.disp[u_start] = Eisenstein(0, 0);
  }

  FaceLengths fl = face_lengths_through(D, h_seed);
  std::optional<Eisenstein> apex_seed =
      place_third_eis({0, 0}, B_seed, fl.abs2, fl.acs2, fl.bcs2, chirality);
  if (!apex_seed) { R.code = Code::SeedChiralityBarrier; return R; }
  std::priority_queue<Front, std::vector<Front>, FrontGreater> pq;
  std::vector<ConeChain> chain;
  chain.push_back({-1, B_seed});
  pq.push({h_seed, {0, 0}, B_seed, ccw_order(B_seed, *apex_seed), /*chain_idx=*/0, /*priority=*/0});

  for (;;) {
    StepKind k = step_walk(D, chirality, bound_sq, target_label, pq, chain, R, state_cap);
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
// both chiralities).  Self mode (calculate_self_geodesics): each walk runs with
// target_label = u so the diagonal is filled with the closed-geodesic squared
// length when one fits within bound_sq, or sealed with bound_sq (the "no closure
// in this bound" sentinel).  Throws on a deep-invariant failure -- non-Loeschian
// seed, degenerate face, or a BFS state-cap trip -- none of which occur on a
// valid iDT.
SimpleMetric compute_simple(const DelaunayTriangulation& D, long long bound_sq,
                            bool calculate_self_geodesics) {
  const int n = D.nv;
  SimpleMetric S{ matrix<long long>(n, n, LLONG_MAX),
                  matrix<simple_geodesic>(n, n, simple_geodesic{}) };
  if (!calculate_self_geodesics)
    for (int u = 0; u < n; u++) S.square(u, u) = 0;   // geo(u,u) stays trivial

  for (int u = 0; u < n; u++) {
    const int target = calculate_self_geodesics ? u : -1;
    for (int h : fan_half_edges(D, u)) {
      std::vector<Eisenstein> Bs = sector0_reps_of_norm(Lsq_of(D, h));
      if (Bs.empty())
        throw std::logic_error(
            "DelaunayTriangulation surface metric: non-Loeschian seed Lsq="
            + std::to_string(Lsq_of(D, h)) + " at h=" + std::to_string(h));
      for (Eisenstein B : Bs)
        for (int chi : {+1, -1}) {
          WalkResult R = walk_from_seed(D, u, h, B, chi, bound_sq, target);
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
  }

  // Seal the diagonal: cones with no closed geodesic found get a sentinel
  // strictly greater than any recordable squared length (the BFS records iff
  // priority <= bound_sq).  bound_sq + 1 avoids ambiguity against a real
  // closure of length exactly bound_sq and stays positive when bound_sq == 0.
  // Off-diagonal "unreached" entries stay at LLONG_MAX.
  if (calculate_self_geodesics)
    for (int u = 0; u < n; u++)
      if (S.square(u, u) == LLONG_MAX) S.square(u, u) = bound_sq + 1;
  return S;
}

}  // anonymous namespace

// ============================================================================
// Public DelaunayTriangulation surface-metric methods
// ============================================================================

// Self mode walks need to reach the (2*diameter) bound for the wrap-around
// closure flavour to be visible, mirroring TriangulationView's M *= 2.
static long long bound_for(const DelaunayTriangulation& D, bool calculate_self_geodesics) {
  long long b = metric_bound(D);
  return calculate_self_geodesics ? 4 * b : b;
}

matrix<long long>
DelaunayTriangulation::simple_square_surface_distances(bool calculate_self_geodesics) const {
  return compute_simple(*this, bound_for(*this, calculate_self_geodesics),
                        calculate_self_geodesics).square;
}

matrix<DelaunayTriangulation::simple_geodesic>
DelaunayTriangulation::simple_geodesics(bool calculate_self_geodesics) const {
  return compute_simple(*this, bound_for(*this, calculate_self_geodesics),
                        calculate_self_geodesics).geo;
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
DelaunayTriangulation::surface_distances(bool calculate_self_geodesics,
                                         matrix<geodesic>* geodesics_out) const {
  SimpleMetric S = compute_simple(*this, bound_for(*this, calculate_self_geodesics),
                                  calculate_self_geodesics);

  // sqrt of the integer simple distances (unreachable -> +inf), then APSP
  // smoothing to enforce the triangle inequality across intermediate cones.
  // In self mode the diagonal carries the closed-geodesic squared length (or
  // bound sentinel); the APSP relaxes it to min(closed_loop, 2*nearest_neighbor).
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
    for (int v = 0; v < n; v++) {
      // Self-mode diagonal: reconstruct_path(preds, u, u) returns [u], which
      // compose_simple_geodesics turns into an empty geodesic{}.  Emit the
      // single recorded closing simple_geodesic directly so the diagonal
      // matches the squared distance flowing through APSP.
      if (u == v && calculate_self_geodesics && S.geo(u, u).g.norm2() > 0) {
        (*geodesics_out)(u, v).segments.push_back(S.geo(u, u));
        continue;
      }
      (*geodesics_out)(u, v) =
          compose_simple_geodesics(reconstruct_path(apsp.preds, u, v), S.geo);
    }
  return apsp.dist.square_elementwise();
}

matrix<DelaunayTriangulation::geodesic>
DelaunayTriangulation::surface_geodesics(bool calculate_self_geodesics) const {
  matrix<geodesic> G(0, 0, geodesic{});
  surface_distances(calculate_self_geodesics, &G);
  return G;
}
