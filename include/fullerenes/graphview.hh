#pragma once

// View hierarchy for the graph class system.
//
// All view types are trivially copyable (TC) -- just spans + scalars.
// This enables zero-copy interop with SYCL kernels, numpy arrays, and
// batch GPU pipelines.
//
// The view hierarchy provides full polymorphism:
//   GraphView -> PlanarGraphView -> CubicGraphView -> FullereneGraphView
//   GraphView -> PlanarGraphView -> TriangulationView -> FullereneDualView
//
// Geometry views (PolyhedronView, DeltahedronView) add a span<coord3d>
// for vertex coordinates.
//
// Algorithms are compiled once on the view types. Owned<View> (owned.hh)
// adds vector storage for any view type.

#include "fullerenes/dense_graph.hh"
#include "fullerenes/geometry.hh"

// --- Forward declarations for owned types ---
template<typename View> struct Owned;

struct GraphView;
struct PlanarGraphView;
struct CubicGraphView;
struct FullereneGraphView;
struct TriangulationView;
struct FullereneDualView;
struct PolyhedronView;
struct DeltahedronView;

using Graph         = Owned<GraphView>;
using PlanarGraph   = Owned<PlanarGraphView>;
using Triangulation = Owned<TriangulationView>;
using FullereneDual = Owned<FullereneDualView>;
using Polyhedron    = Owned<PolyhedronView>;

// Thin wrappers (defined after Owned in their own headers):
struct CubicGraph;
struct FullereneGraph;
struct Deltahedron;

// ---------------------------------------------------------------------------
// GraphView: general undirected/directed graph with adjacency lists.
// ---------------------------------------------------------------------------
struct GraphView : Spanify::RSRAdjacencyView<node_t> {
    using RSRAdjacencyView::RSRAdjacencyView;
    static constexpr uint8_t default_dmax = 10;
    GraphView() = default;
};

// ---------------------------------------------------------------------------
// PlanarGraphView: planar graph with oriented embedding.
// ---------------------------------------------------------------------------
struct PlanarGraphView : GraphView {
    using GraphView::GraphView;
    static constexpr uint8_t default_dmax = 6;
};

// ---------------------------------------------------------------------------
// CubicGraphView: 3-regular planar graph.
// ---------------------------------------------------------------------------
struct CubicGraphView : PlanarGraphView {
    using PlanarGraphView::PlanarGraphView;
    static constexpr uint8_t default_dmax = 3;
};

// ---------------------------------------------------------------------------
// FullereneGraphView: fullerene graph (3-regular, 12 pentagons, rest hex).
// ---------------------------------------------------------------------------
struct FullereneGraphView : CubicGraphView {
    using CubicGraphView::CubicGraphView;
    static constexpr uint8_t default_dmax = 3;
};

// ---------------------------------------------------------------------------
// TriangulationView: planar triangulation (max degree 6).
// ---------------------------------------------------------------------------
struct TriangulationView : PlanarGraphView {
    using PlanarGraphView::PlanarGraphView;
    static constexpr uint8_t default_dmax = 6;
};

// ---------------------------------------------------------------------------
// FullereneDualView: dual of a fullerene (triangulation with 12 degree-5
// vertices, rest degree 6).
// ---------------------------------------------------------------------------
struct FullereneDualView : TriangulationView {
    using TriangulationView::TriangulationView;
    static constexpr uint8_t default_dmax = 6;
};

// ---------------------------------------------------------------------------
// PolyhedronView: planar graph with 3D vertex coordinates.
// ---------------------------------------------------------------------------
struct PolyhedronView : PlanarGraphView {
    std::span<coord3d> points;
    static constexpr uint8_t default_dmax = 10;

    PolyhedronView() = default;

    // Construct from adjacency view + coordinate span.
    PolyhedronView(const PlanarGraphView& g, std::span<coord3d> pts)
        : PlanarGraphView(g), points(pts) {}

    // Full view constructor (adjacency + coordinates).
    PolyhedronView(node_t N, int dmax,
                   std::span<node_t> neighbours, std::span<uint8_t> deg,
                   std::span<coord3d> pts,
                   std::span<uint8_t> twin = {})
        : PlanarGraphView(N, dmax, neighbours, deg, twin), points(pts) {}
};

// ---------------------------------------------------------------------------
// DeltahedronView: triangulation with 3D vertex coordinates (equilateral
// triangle embedding).
// ---------------------------------------------------------------------------
struct DeltahedronView : TriangulationView {
    std::span<coord3d> points;
    static constexpr uint8_t default_dmax = 6;

    DeltahedronView() = default;

    // Construct from adjacency view + coordinate span.
    DeltahedronView(const TriangulationView& t, std::span<coord3d> pts)
        : TriangulationView(t), points(pts) {}

    // Full view constructor (adjacency + coordinates).
    DeltahedronView(node_t N, int dmax,
                    std::span<node_t> neighbours, std::span<uint8_t> deg,
                    std::span<coord3d> pts,
                    std::span<uint8_t> twin = {})
        : TriangulationView(N, dmax, neighbours, deg, twin), points(pts) {}
};

// ---------------------------------------------------------------------------
// TC verification: all view types must be trivially copyable.
// ---------------------------------------------------------------------------
static_assert(std::is_trivially_copyable_v<GraphView>,
    "GraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<PlanarGraphView>,
    "PlanarGraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<CubicGraphView>,
    "CubicGraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<FullereneGraphView>,
    "FullereneGraphView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<TriangulationView>,
    "TriangulationView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<FullereneDualView>,
    "FullereneDualView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<PolyhedronView>,
    "PolyhedronView must be trivially copyable");
static_assert(std::is_trivially_copyable_v<DeltahedronView>,
    "DeltahedronView must be trivially copyable");

// Verify correct hierarchy relationships
static_assert(std::is_base_of_v<GraphView, PlanarGraphView>);
static_assert(std::is_base_of_v<PlanarGraphView, CubicGraphView>);
static_assert(std::is_base_of_v<CubicGraphView, FullereneGraphView>);
static_assert(std::is_base_of_v<PlanarGraphView, TriangulationView>);
static_assert(std::is_base_of_v<TriangulationView, FullereneDualView>);
static_assert(std::is_base_of_v<PlanarGraphView, PolyhedronView>);
static_assert(std::is_base_of_v<TriangulationView, DeltahedronView>);
