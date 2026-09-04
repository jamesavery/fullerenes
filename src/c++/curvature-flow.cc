// curvature-flow -- see include/fullerenes/curvature-flow.hh for the method. Purely
// intrinsic on the DCEL metric (DelaunayTriangulation); the disk partition and the
// cone insertion are DelaunayTriangulation methods, the 3D->metric seed is
// DeltahedronView::intrinsic_metric().

#include "fullerenes/curvature-flow.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ------------------------------------------------------------ curvature

// The weight of a disk member under a mode -- the only thing that differs between the
// four prescriptions. (v, d) is the member and its geodesic distance from its cone.
double weight(CurvatureMode mode, std::pair<int,double> member, int cone,
              double sigma, const std::vector<double>& K) {
  auto [v, d] = member;
  switch (mode) {
    case CurvatureMode::Shape:    return std::max(K[v], 0.0);
    case CurvatureMode::Uniform:  return 1.0;
    case CurvatureMode::Gaussian: return std::exp(-d*d / (2.0*sigma*sigma));
    case CurvatureMode::Point:    return v == cone ? 1.0 : 0.0;
  }
  return 0.0;
}

// Spread pi/3 (or, with explicit `deficits`, disk i's deficits[i]) over each (disjoint)
// disk in proportion to the mode's weight, then rescale to Sigma K* = 4pi -- a no-op up
// to rounding when the per-disk shares already sum to 4pi (the equal-pi/3 case, or a
// `deficits` input that sums to 4pi by construction), kept only to kill FP drift; skipped
// when there are no cones. The four modes ARE this one operation; only weight() changes.
// `deficits[i]` is disk i's budget: geodesic_disks returns disks 1:1 with its `sources`
// in order (disks[i].source == sources[i], its ORDER CONTRACT), and confine_curvature
// seeds those sources in `cones` order, so deficits[i] <-> cones[i].
std::vector<double> prescribe_curvature(const DelaunayTriangulation& D,
                                        const std::vector<GeodesicDisk>& disks,
                                        CurvatureMode mode, double R,
                                        std::span<const double> deficits) {
  const std::vector<double> K = mode == CurvatureMode::Shape ? D.curvature()
                                                             : std::vector<double>{};
  const double sigma = R / std::sqrt(2.0*std::log(100.0));   // 99% of the 2D mass within R

  std::vector<double> kstar(D.nv, 0.0);
  for (size_t di = 0; di < disks.size(); di++) {
    const GeodesicDisk& disk = disks[di];
    const double per_disk = deficits.empty() ? M_PI/3.0 : deficits[di];
    std::vector<double> w(disk.members.size());
    for (size_t i = 0; i < w.size(); i++) w[i] = weight(mode, disk.members[i], disk.source, sigma, K);
    const double total = std::accumulate(w.begin(), w.end(), 0.0);
    if (total <= 0.0)
      throw std::runtime_error("confine_curvature: a disk carries no curvature weight in this mode");
    for (size_t i = 0; i < w.size(); i++) kstar[disk.members[i].first] += per_disk * w[i] / total;
  }

  const double sum = std::accumulate(kstar.begin(), kstar.end(), 0.0);
  if (sum > 0.0) for (double& k : kstar) k *= 4.0*M_PI / sum;
  return kstar;
}

// ------------------------------------------------------------ cone insertion

// Floor on barycentric coordinates: keeps the cone strictly inside its face, so a
// detector-clamped on-edge/vertex centre cannot make a zero-length spoke or sliver.
constexpr double kMinBary = 1e-3;

// Spoke lengths {P->a, P->b, P->c} from a barycentric interior point P to the corners
// of a triangle with the given edge lengths: lay the triangle flat and measure
// |P - corner|. Metric-preserving since P stays inside the flat face.
std::array<double,3> cone_spoke_lengths(double L_ab, double L_bc, double L_ca,
                                        std::array<double,3> bary) {
  const coord2d A(0,0), B(L_ab, 0), C = place_apex(L_ab, L_bc, L_ca);
  const coord2d P = B*bary[1] + C*bary[2];        // bary[0]*A vanishes (A at origin)
  return {(P-A).norm(), (P-B).norm(), (P-C).norm()};
}

// The half-edge a->b of the face (a, b, c) (CCW). Throws if absent -- e.g. an earlier
// cone already consumed it (cones too close to insert independently).
int face_half_edge(const DelaunayTriangulation& D, int a, int b, int c) {
  for (int h : D.incident(a))
    if (D.dest(h) == b && D.dest(D.he_next[h]) == c) return h;
  throw std::runtime_error("confine_curvature: cone face not found (cone patches overlap?)");
}

// Insert a flat vertex at the cone's (interior-nudged) barycentric position by a 1->3
// split of its face. Returns the new vertex id.
int insert_cone(DelaunayTriangulation& D, const ConeSite& cone) {
  const int h0 = face_half_edge(D, cone.face[0], cone.face[1], cone.face[2]);   // a->b
  const int h1 = D.he_next[h0], h2 = D.he_next[h1];
  const auto spoke = cone_spoke_lengths(D.he_length[h0], D.he_length[h1], D.he_length[h2],
                                        interior_bary(cone.bary));
  return D.split_face(h0, spoke);
}

} // namespace

// Exposed (declared in curvature-flow.hh): nudge a barycentric triple strictly interior (each coord
// >= kMinBary) and renormalize. Both insert_cone (above) and warm-start placement use it so the two
// agree on where a cone actually sits.
std::array<double,3> interior_bary(std::array<double,3> b) {
  for (double& w : b) w = std::max(w, kMinBary);
  const double s = b[0] + b[1] + b[2];
  return {b[0]/s, b[1]/s, b[2]/s};
}

ConePrescription confine_curvature(DelaunayTriangulation surface,
                                   const std::vector<ConeSite>& cones,
                                   double R, CurvatureMode mode, DiskMetric disk_metric,
                                   std::span<const double> deficits) {
  if (!deficits.empty()) {
    if (deficits.size() != cones.size())
      throw std::runtime_error("confine_curvature: deficits.size() must equal cones.size()");
    double sum = 0; for (double d : deficits) sum += d;
    if (std::fabs(sum - 4.0*M_PI) > 1e-6)
      throw std::runtime_error("confine_curvature: deficits sum to " + std::to_string(sum)
                               + ", not 4*pi (Gauss-Bonnet)");
  }
  std::vector<int> cone_ids;
  for (const ConeSite& cone : cones) cone_ids.push_back(insert_cone(surface, cone));
  std::vector<GeodesicDisk> disks = surface.geodesic_disks(cone_ids, R, disk_metric);
  std::vector<double>       kstar = prescribe_curvature(surface, disks, mode, R, deficits);
  return { std::move(surface), std::move(cone_ids), std::move(kstar) };
}
