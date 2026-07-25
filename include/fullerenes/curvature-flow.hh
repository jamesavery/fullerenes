#pragma once

// curvature-flow -- the intrinsic half of the prescribed-curvature flow.
//
// Given a closed genus-0 fullerene-shell surface (as its intrinsic metric) and a set
// of detected cones (centre + geodesic radius R), confine all 4*pi of Gaussian
// curvature (Gauss-Bonnet) into the geodesic R-disks around the cones -- pi/3 per
// disk, zero outside -- and emit the prescribed curvature field K* on the cone-refined
// metric. The conformal solve that realises this metric is a separate, later step.
//
// The surface is an intrinsic DCEL (DelaunayTriangulation): half-edge topology +
// per-edge length, no 3D positions and no maximum vertex degree. A 3D mesh enters the
// pipeline through exactly ONE step -- DeltahedronView::intrinsic_metric() seeds the
// edge lengths from coordinates -- and never again: cone insertion, the geodesic
// disks, and the curvature prescription are all purely intrinsic (lengths only).
//
// The flow is three operations:
//   1. insert a flat vertex at each cone (DelaunayTriangulation::split_face),
//   2. partition the surface into geodesic R-disks -- one multi-source front from all
//      cones at once, so each vertex is claimed by its nearest cone and the disks are
//      DISJOINT (DelaunayTriangulation::geodesic_disks, a geodesic Voronoi clipped at
//      R; DiskMetric: edge-graph Dijkstra, fast; or fast-marching, more exact),
//   3. confine 4*pi of curvature into the disks.
//
// Step 3 is one operation with four faces: every mode distributes pi/3 over each disk
// in proportion to a per-vertex WEIGHT, differing only in the weight field --
//   shape:    w = max(K_current, 0)      (follow the shell's own curvature)
//   uniform:  w = 1                       (spread evenly over the disk)
//   gaussian: w = exp(-d^2 / 2 sigma^2)   (geodesic bump, sigma = R/sqrt(2 ln 100))
//   point:    w = [v is the cone]         (a single delta at the centre)
// -- then a global rescale to Sigma K* = 4*pi.

#include "fullerenes/delaunay.hh"   // DelaunayTriangulation, DiskMetric, tri_t

#include <array>
#include <span>
#include <vector>

// How pi/3 is distributed WITHIN each geodesic disk (see header comment); the
// confinement to the disks is the same across modes.
enum class CurvatureMode { Shape, Uniform, Gaussian, Point };

// A cone's combinatorial location: the face (its 3 CCW vertices) and the barycentric
// position of the cone within it. Deliberately minimal -- decouples the intrinsic flow
// from the 3D cone finder's SurfaceCone (graphview.hh).
struct ConeSite { tri_t face; std::array<double,3> bary; };

// Where confine_curvature actually inserts a cone, given its requested barycentric: the triple
// nudged strictly interior (each coord >= 1e-3) and renormalised, so a cone detected on a triangle
// edge/vertex still lands inside the face. Exposed so a warm-start 3D placement can match the
// inserted cone's intrinsic position rather than the un-nudged detector centre.
std::array<double,3> interior_bary(std::array<double,3> bary);

// The prescribed curvature on the cone-refined metric.
struct ConePrescription {
  DelaunayTriangulation surface;   // DCEL metric with the cone vertices inserted
  std::vector<int>      cones;     // the inserted cone vertex ids
  std::vector<double>   kstar;     // prescribed Gaussian curvature per vertex; Sigma = 4*pi
};

// Insert a flat vertex at each cone, partition the surface into disjoint geodesic
// R-disks (disk_metric), and confine all 4*pi of Gaussian curvature into them (mode
// sets the within-disk distribution). Purely intrinsic; `surface` is consumed into the
// returned ConePrescription (the caller moves it in).
//
// `deficits`, if given, overrides the default equal pi/3-per-disk split: deficits[i] is
// the total curvature confined into cones[i]'s disk (in place of the fixed pi/3), so
// disks may carry unequal shares -- e.g. a measured cubic metric whose cones carry a
// varying number of incident pentagons. Must be empty or the same size as `cones`, and
// (when non-empty) must sum to 4*pi (Gauss-Bonnet) -- both are checked, throwing
// otherwise. The Sigma K* = 4*pi rescale still runs afterwards, a no-op up to FP drift
// since the deficits already sum to 4*pi (as they do by construction for the
// per-pentagon-count case).
ConePrescription confine_curvature(DelaunayTriangulation surface,
                                   const std::vector<ConeSite>& cones,
                                   double R, CurvatureMode mode, DiskMetric disk_metric,
                                   std::span<const double> deficits = {});
