// DeltahedronView<double>::find_cones -- locate curvature-concentration cones on a
// triangulation surface (e.g. the 12 fullerene pentagon centres of a tomogram shell).
// Declaration + SurfaceCone/ConeFinderParams: include/fullerenes/graphview.hh.
//
// A closed genus-0 surface carries total Gaussian curvature 4*pi (Gauss-Bonnet,
// chi=2); a fullerene shell concentrates it in 12 cones of pi/3 each. This finds
// those cones as surface points, each with the angle-defect mass it integrates.
//
// Algorithm (the "M2_taubin" method):
//   0. Taubin-smooth a COPY of the points (DeltahedronView::taubin_smooth) to
//      suppress marching-cubes staircase noise. The input mesh is untouched.
//   1. Signed discrete angle defect K (DeltahedronView::angle_defects) on the copy.
//   2. Greedily pack disjoint R disks at strict 1-ring local maxima of the signed
//      disk-integral detector (cosine-tapered top-hat, maintained INCREMENTALLY as
//      disks are claimed -- f[i] = sum_j W[i,j] K[j] [j unclaimed], updated by a
//      sparse AXPY over each claimed column, collapsing C convolutions to ~one),
//      each accepted iff its hard R-disk integral of signed K clears threshold.
//   3. Mean-shift each centre to its (positive-part) K-mass centroid, then project
//      onto the ORIGINAL (unsmoothed) mesh.
//
// Distances are Euclidean (chordal), a good geodesic proxy on near-isotropic shells;
// spatial queries use a uniform cell list (cell = R). Detection and the acceptance
// gate use the SIGNED defect (positive/negative noise cancels in flat regions); the
// mean-shift centroid uses the positive part (a signed centroid is ill-defined).

#include "fullerenes/deltahedron.hh"
#include "fullerenes/delaunay.hh"
#include "fullerenes/geometry.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace {

// ---------------------------------------------------------------- cell list

// Uniform grid over a point set for Euclidean radius queries. Cell size is the
// disk radius R, so a radius-r query scans the +-ceil(r/cell) cell shell.
struct CellList {
  std::span<const coord3d> P;            // the points (caller owns; must outlive)
  coord3d  lo;
  double   cell;
  int      nx, ny, nz;
  std::vector<int> start;   // size nx*ny*nz + 1 (prefix sums)
  std::vector<int> items;   // point indices, bucketed by cell

  CellList(std::span<const coord3d> pts, double cell_size) : P(pts), cell(cell_size) {
    coord3d hi;
    lo = hi = pts.empty() ? coord3d(0,0,0) : pts[0];
    for (const auto& p : pts)
      for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], p[d]); hi[d] = std::max(hi[d], p[d]); }
    // The cell is the query radius, but clamp it up so the grid never has far more
    // cells than points (a pathologically small radius on a large mesh would
    // otherwise overflow nx*ny*nz). A coarser grid stays correct: a radius-r query
    // scans ceil(r/cell) cell shells either way, it just buckets more per cell.
    const double span = std::max({hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2], 0.0});
    cell = std::max(cell, span / std::max(1.0, std::cbrt(8.0 * (double)pts.size())));
    nx = std::max(1, (int)((hi[0]-lo[0])/cell) + 1);
    ny = std::max(1, (int)((hi[1]-lo[1])/cell) + 1);
    nz = std::max(1, (int)((hi[2]-lo[2])/cell) + 1);

    const size_t ncells = (size_t)nx*ny*nz;            // size_t: nx*ny*nz can exceed INT_MAX
    start.assign(ncells + 1, 0);
    for (const auto& p : pts) start[cell_of(p) + 1]++;
    for (size_t c = 0; c < ncells; c++) start[c+1] += start[c];
    items.resize(pts.size());
    std::vector<int> cursor(start.begin(), start.end() - 1);
    for (int i = 0; i < (int)pts.size(); i++) items[cursor[cell_of(pts[i])]++] = i;
  }

  int axis_cell(double x, int lo_d, int n) const {
    int c = (int)((x - lo[lo_d]) / cell);
    return std::min(std::max(c, 0), n - 1);
  }
  size_t cell_of(const coord3d& p) const {
    return ((size_t)axis_cell(p[0],0,nx)*ny + axis_cell(p[1],1,ny))*nz + axis_cell(p[2],2,nz);
  }
  // Longest axis span: a query of this radius covers the whole grid.
  double extent() const { return (double)std::max({nx,ny,nz}) * cell; }

  // Visit every point index whose cell lies within ceil(r/cell) of q's cell.
  template <class F>
  void query(const coord3d& q, double r, F&& f) const {
    const int nc = (int)std::ceil(r / cell);
    const int cx = axis_cell(q[0],0,nx), cy = axis_cell(q[1],1,ny), cz = axis_cell(q[2],2,nz);
    for (int ix = std::max(0,cx-nc); ix <= std::min(nx-1,cx+nc); ix++)
      for (int iy = std::max(0,cy-nc); iy <= std::min(ny-1,cy+nc); iy++)
        for (int iz = std::max(0,cz-nc); iz <= std::min(nz-1,cz+nc); iz++) {
          const size_t c = ((size_t)ix*ny + iy)*nz + iz;
          for (int k = start[c]; k < start[c+1]; k++) f(items[k]);
        }
  }

  // Visit each point j with true Euclidean distance d = |P[j] - q| < r, as f(j, d).
  template <class F>
  void query_within(const coord3d& q, double r, F&& f) const {
    query(q, r, [&](int j) { double d = (P[j]-q).norm(); if (d < r) f(j, d); });
  }
};

// ---------------------------------------------------------------- geometry

// Max distance from a face centroid to its own vertices, over all faces. This is
// the slack that makes the centroid cell-list query provably find the nearest
// face: a face whose centroid is farther than best_dist + tri_radius_max from the
// query cannot have a closer surface point, so once the search shell reaches that
// radius the current best is the global nearest.
double max_triangle_radius(std::span<const coord3d> P, const std::vector<tri_t>& tris,
                           const std::vector<coord3d>& centroids) {
  double r = 0;
  for (size_t i = 0; i < tris.size(); i++)
    for (int k = 0; k < 3; k++) r = std::max(r, (P[tris[i][k]] - centroids[i]).norm());
  return r;
}

// Closest point on the triangulation to `p`, as (triangle, barycentric, pos). The
// centroid cell-list query is grown until best_dist + tri_radius_max is covered --
// the bound that proves no unscanned face can be closer -- so the true nearest face
// is always found, not merely the nearest within the first non-empty shell. (The
// common near-surface query terminates at the first radius-R shell; the expansion
// only triggers for a far query, e.g. reprojecting a SMOOTHED-surface point onto the
// ORIGINAL mesh by a displacement that exceeds R.)
struct Projection { int tri; std::array<double,3> bary; coord3d pos; };

Projection project_to_surface(const coord3d& p, std::span<const coord3d> P,
                              const std::vector<tri_t>& tris, const CellList& face_cells,
                              double R, double tri_radius_max) {
  if (tris.empty()) return {0, {1,0,0}, P.empty() ? coord3d() : P[0]};
  Projection best{0, {1,0,0}, P[tris[0][0]]};
  double best_d2 = std::numeric_limits<double>::infinity();
  for (double r = R; ; r *= 2) {
    face_cells.query(p, r, [&](int ti) {
      ClosestPoint cp = closest_point_on_triangle(p, tris[ti], P);
      if (cp.dist2 < best_d2) { best_d2 = cp.dist2; best = {ti, cp.bary, cp.pos}; }
    });
    if (best_d2 < std::numeric_limits<double>::infinity() &&
        r >= std::sqrt(best_d2) + tri_radius_max) break;   // proven nearest
    if (r > face_cells.extent()) break;                    // whole mesh scanned
  }
  return best;
}

// Gaussian-windowed mean-shift to a K-mass centre (positive-part, unclaimed
// mass) on the work mesh. Returns the refined surface point and convergence.
struct MeanShift { Projection at; bool converged; int iters; };

// @anchor meanshift-center
// @pre    nonneg_mass: all_of(K_pos, [](double k){ return k >= 0; })
// @variant ||c_proj - c||   -- the reprojected step length, until < tol or budget
// @post   on converged:  result.converged && result.iters <= max_iters
// @post   on empty mass: !result.converged && result.iters < max_iters  (the 3*sigma
//                        window held no positive unclaimed K, so the centre cannot move)
// @post   on budget:     !result.converged && result.iters == max_iters
MeanShift meanshift_center(std::span<const coord3d> P, const std::vector<tri_t>& tris,
                           const CellList& vert_cells, const CellList& face_cells,
                           const std::vector<double>& K_pos,
                           const std::vector<uint8_t>& alive,
                           coord3d seed, double R, int max_iters, double tri_radius_max) {
  const double sigma = R/2.0, cutoff = 3.0*sigma, tol = 1e-3*sigma;
  const double inv_2s2 = 1.0/(2.0*sigma*sigma);
  Projection cur = project_to_surface(seed, P, tris, face_cells, R, tri_radius_max);
  for (int it = 0; it < max_iters; it++) {
    coord3d num; double den = 0;
    vert_cells.query_within(cur.pos, cutoff, [&](int j, double d) {
      if (!alive[j] || K_pos[j] == 0) return;
      double w = std::exp(-d*d*inv_2s2) * K_pos[j];
      den += w; num += P[j]*w;
    });
    if (den < 1e-30) return {cur, false, it};
    Projection nxt = project_to_surface(num/den, P, tris, face_cells, R, tri_radius_max);
    if ((nxt.pos - cur.pos).norm() < tol) return {nxt, true, it+1};
    cur = nxt;
  }
  return {cur, false, max_iters};
}

// ---------------------------------------------- R-disk taper-stencil operator

// The detector operator W[i,j] = taper(|i-j|) for |i-j| < R (a cosine-tapered
// top-hat), as CSR. apply computes f = W x; subtract_column removes a column.
// @anchor taper-stencil
struct Stencil {
  std::vector<size_t> row_start;
  std::vector<int>    col;
  std::vector<double> wgt;

  std::vector<double> apply(const std::vector<double>& x) const {     // f = W x
    std::vector<double> f(row_start.size()-1, 0.0);
    for (size_t i = 0; i+1 < row_start.size(); i++)
      for (size_t k = row_start[i]; k < row_start[i+1]; k++) f[i] += wgt[k]*x[col[k]];
    return f;
  }
  void subtract_column(std::vector<double>& f, int j, double s) const { // f -= s W[:,j]
    for (size_t k = row_start[j]; k < row_start[j+1]; k++) f[col[k]] -= wgt[k]*s;
  }
};

Stencil build_taper_stencil(const CellList& cells, std::span<const coord3d> P, double R) {
  const int N = (int)P.size();
  const double r_taper = 0.8*R;                           // top-hat plateau; cosine roll-off in (0.8R, R]
  Stencil W;
  W.row_start.assign(N+1, 0);
  // Two-pass CSR build: pass 1 counts each row's disk size into the prefix sums,
  // pass 2 fills the cosine-tapered weights. Exact sizing, no reallocs, no estimate.
  for (int i = 0; i < N; i++) {
    size_t cnt = 0;
    cells.query_within(P[i], R, [&](int, double) { cnt++; });
    W.row_start[i+1] = W.row_start[i] + cnt;
  }
  W.col.resize(W.row_start[N]);
  W.wgt.resize(W.row_start[N]);
  for (int i = 0; i < N; i++) {
    size_t k = W.row_start[i];
    cells.query_within(P[i], R, [&](int j, double d) {
      W.col[k] = j;
      W.wgt[k] = d <= r_taper ? 1.0
               : 0.5*(1.0 + std::cos(M_PI*(d - r_taper)/(R - r_taper)));
      k++;
    });
  }
  return W;
}

// The detector field: f = W (K . alive), built once and maintained incrementally
// as vertices are claimed (alive 1 -> 0). The detector's evolving state.
// @anchor detector-field
struct DetectorField {
  Stencil W;
  std::vector<double> K, f;
  std::vector<uint8_t> alive;                            // 1 = unclaimed

  DetectorField(const CellList& cells, std::span<const coord3d> P,
                std::vector<double> K_, double R)
    : W(build_taper_stencil(cells, P, R)), K(std::move(K_)), f(W.apply(K)),
      alive(K.size(), 1) {}

  void claim(int j) {
    if (!alive[j]) return;
    alive[j] = 0;
    W.subtract_column(f, j, K[j]);
  }
};

// The unclaimed vertex of strictly-largest field among its unclaimed 1-ring,
// with f >= threshold; -1 if none. The packer's candidate-ranking step. NOTE: the
// strictness is deliberate, but it means a perfectly symmetric mesh where two
// adjacent vertices share a bit-identical peak field (e.g. a regular icosahedron)
// has NO strict maximum there and that cone is missed -- a non-issue for real,
// noisy shells, where smoothing breaks such ties.
int best_local_max(const DetectorField& d, const DeltahedronView<double>& mesh,
                   double threshold) {
  int best = -1; double best_f = threshold;
  for (int i = 0; i < (int)d.f.size(); i++) {
    if (!d.alive[i] || d.f[i] < threshold) continue;
    bool is_max = true, has_alive_nbr = false;
    for (node_t j : mesh[i]) {
      if (!d.alive[j]) continue;
      has_alive_nbr = true;
      if (d.f[j] >= d.f[i]) { is_max = false; break; }
    }
    // A vertex whose entire 1-ring is already claimed is not a strict local
    // max (matches the reference: an empty unclaimed neighbour set disqualifies).
    if (is_max && has_alive_nbr && d.f[i] >= best_f) { best_f = d.f[i]; best = i; }
  }
  return best;
}

// ---------------------------------------------------------------- canonical order

// Deterministic output order: strongest cone first (largest integrated_K), with a
// lexicographic position tiebreak. The cone-finder is deterministic, so within one
// build each cone's integrated_K is bit-reproducible and the order is stable; the
// position tiebreak resolves exact K ties. No axis/centre/polar assumption -- the
// cones are geometric points with no underlying graph to canonicalise against.
// (Across builds/platforms two near-equal-K cones could swap before the tiebreak
// engages; consumers needing a physical-site key should match by position.)
void canonical_cone_order(std::vector<SurfaceCone>& cones) {
  std::sort(cones.begin(), cones.end(), [](const SurfaceCone& a, const SurfaceCone& b) {
    if (a.integrated_K != b.integrated_K) return a.integrated_K > b.integrated_K;
    for (int d = 0; d < 3; d++) if (a.position[d] != b.position[d]) return a.position[d] < b.position[d];
    return false;
  });
}

} // namespace

template<>
std::vector<SurfaceCone>
DeltahedronView<double>::find_cones(const ConeFinderParams& params,
                                   std::vector<coord3d>* smoothed_points) const {
  if (!(params.R_pent > 0) || !std::isfinite(params.R_pent))
    throw std::invalid_argument("find_cones: R_pent must be a finite value > 0");
  if (!(params.threshold_frac >= 0 && params.threshold_frac <= 1))
    throw std::invalid_argument("find_cones: threshold_frac must be in [0,1]");
  // Finite-geometry precondition: a NaN/inf vertex would silently corrupt the cell
  // list, the detector field, and the output sort (NaN breaks strict-weak ordering).
  for (const coord3d& p : points)
    if (!(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2])))
      throw std::invalid_argument("find_cones: mesh has a non-finite vertex coordinate");

  const double R = params.R_pent, pi_3 = M_PI/3.0;
  const double threshold = params.threshold_frac * pi_3;
  const std::vector<tri_t> tris = triangles();
  if (N == 0 || tris.empty()) return {};               // nothing to find on an empty mesh

  // stage 0: Taubin-smoothed COPY (this mesh untouched); detection runs on it.
  Deltahedron work(static_cast<const TriangulationView&>(*this), points);
  work.taubin_smooth(params.taubin_iters);
  std::vector<coord3d> smoothed(work.points.begin(), work.points.end());
  std::span<const coord3d> P_smooth = smoothed;
  std::span<const coord3d> P_orig = points;
  if (smoothed_points) *smoothed_points = smoothed;

  // stage 1: signed angle defect on the smoothed mesh (reuse tris; positive part for mean-shift).
  std::vector<double> K = work.angle_defects(tris);
  std::vector<double> K_pos(N);
  for (int i = 0; i < N; i++) K_pos[i] = std::max(K[i], 0.0);

  // spatial indices (built once): vertices, and the two meshes' face centroids.
  CellList vert_cells(P_smooth, R);
  std::vector<coord3d> work_centroids(tris.size()), orig_centroids(tris.size());
  for (size_t i = 0; i < tris.size(); i++) {
    work_centroids[i] = tris[i].centroid(P_smooth);
    orig_centroids[i] = tris[i].centroid(P_orig);
  }
  CellList work_face_cells(work_centroids, R), orig_face_cells(orig_centroids, R);
  const double work_tri_rmax = max_triangle_radius(P_smooth, tris, work_centroids);
  const double orig_tri_rmax = max_triangle_radius(P_orig,   tris, orig_centroids);

  // detector field, maintained incrementally as disks are claimed.
  DetectorField detector(vert_cells, P_smooth, std::move(K), R);

  // stage 2/3: greedy packing -- pick best local max, refine, gate, record, claim.
  std::vector<SurfaceCone> cones;
  for (int c = 0; c < params.max_centres; c++) {
    int best = best_local_max(detector, *this, threshold);
    if (best < 0) break;

    MeanShift ms = meanshift_center(P_smooth, tris, vert_cells, work_face_cells, K_pos,
                                    detector.alive, P_smooth[best], R, params.meanshift_iters, work_tri_rmax);

    // the unclaimed vertices within R of the refined centre -- one query, used for
    // both the hard top-hat acceptance gate (signed K sum) and the claim.
    std::vector<int> disk;
    vert_cells.query_within(ms.at.pos, R, [&](int j, double) { if (detector.alive[j]) disk.push_back(j); });
    double f_at_c = 0; for (int j : disk) f_at_c += detector.K[j];
    if (f_at_c < threshold) { detector.claim(best); continue; }   // anti-stall, retry

    // report on the ORIGINAL geometry, then claim the disk and the seed (the latter
    // in case mean-shift drifted the centre so far that `best` fell outside the disk,
    // which would otherwise leave it re-selectable on the next iteration).
    Projection op = project_to_surface(ms.at.pos, P_orig, tris, orig_face_cells, R, orig_tri_rmax);
    cones.push_back({op.tri, op.bary, op.pos, f_at_c, ms.converged, ms.iters});
    for (int j : disk) detector.claim(j);
    detector.claim(best);
  }

  canonical_cone_order(cones);
  return cones;
}

// The intrinsic surface metric: the iDT whose edge lengths are this mesh's Euclidean
// edge lengths. The lone step that reads 3D coordinates; the curvature flow downstream
// (curvature-flow.hh) is purely intrinsic.
template<>
DelaunayTriangulation DeltahedronView<double>::intrinsic_metric() const {
  Triangulation T(static_cast<const GraphView&>(*this));   // combinatorial part, oriented
  auto length = [&](node_t u, node_t v) { return (points[u] - points[v]).norm(); };
  return DelaunayTriangulation::from_intrinsic_metric(T, length);
}
