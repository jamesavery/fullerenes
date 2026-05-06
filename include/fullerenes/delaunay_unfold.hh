#pragma once

#include "fullerenes/eisenstein.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"

#include <array>
#include <iosfwd>
#include <map>
#include <utility>
#include <vector>

// =====================================================================
// iDT face unfolding into a single Eisenstein-lattice frame.
//
// Diagnostic / visualisation only.  Lays out every live iDT cell into
// one Eisenstein-plane drawing via spanning-tree DFS over face
// adjacency:
//   - Each arc u -> v is either "placed" (the cell to its left has a
//     global placement) or "open".
//   - Pop an open arc, place the cell to its left with the affine
//     transform that pins its two shared cones to the already-placed
//     global positions.
//   - If a cone gets placed at a different position than a previous
//     placement, count it as a "tear" -- this is intrinsic to delta-
//     complex unfoldings (Gauss-Bonnet) and cannot be avoided, only
//     re-routed by choosing a different DFS spanning tree.
//
// Independent of any specific embedding pipeline: the input is a
// vector<CellPlacement>, one per iDT face, that describes each cell's
// corner positions in its own local Eisenstein frame plus the
// (optional) lattice points in that frame.  Pipelines that produce
// such per-cell data (e.g. delaunay-fillin) build the CellPlacement
// vector and pass it in.
// =====================================================================

// One iDT cell's combinatorics + Eisenstein-lattice placement.
//
// Caller-supplied per cell:
//   * cell_id                 -- D's face id
//   * c0, c1, c2              -- T_sorted corner vertex IDs, CCW
//   * P0, P1, P2              -- corner positions in the cell's own
//                                local Eisenstein frame
//   * lattice_points          -- Eisenstein lattice points strictly
//                                inside / on the cell, each tagged
//                                with its T_sorted vertex id.  May be
//                                empty if the caller doesn't need
//                                lattice-point output in the unfolding.
//   * ok                      -- true iff the placement is meaningful;
//                                false slots are skipped by unfold_iDT
struct CellPlacement {
    int        cell_id = -1;
    int        c0 = -1, c1 = -1, c2 = -1;
    Eisenstein P0{0, 0};
    Eisenstein P1{0, 0};
    Eisenstein P2{0, 0};
    std::vector<std::pair<Eisenstein, int>> lattice_points;
    bool       ok = false;
};

// Result of unfold_iDT.
struct LatticeUnfolding {
    // All DISTINCT global positions seen for each cone across placed
    // cells.  On a delta-complex with cone angle deficits some cones
    // end up at 2+ positions (= spanning-tree cross-edges force "tears"
    // in the unfolding); the global net carries all of them.
    std::map<int, std::vector<Eisenstein>> cone_positions;

    // For each placed iDT cell: corners (c0, c1, c2) and their global
    // positions in the unfolding's spanning-tree frame, plus the
    // lattice_points translated into those global coords.
    struct UnfoldedCell {
        int cell_id;
        int c0, c1, c2;
        Eisenstein P0, P1, P2;
        int parent_cell_id = -1;     // -1 = seed
        std::vector<std::pair<Eisenstein, int>> entries;  // (global pos, vertex_id)
    };
    std::vector<UnfoldedCell> cells;
    int n_cells = 0;
    int n_tears = 0;       // #(cone, cell) pairs where cone ends up at
                           // a position different from its first sighting
};

// Unfold all live iDT cells into a single Eisenstein-lattice frame via
// a spanning-tree DFS over cell adjacency.  Every live cell with
// ok == true is placed exactly once (via the arc that first reaches
// it in the DFS), so n_cells == count of live ok cells.  Cone
// deficits manifest as "tears": a cone may end up at multiple
// positions across placed cells.
//
// `start_cell_id == -1` means pick the first ok cell.  `dfs_max_seeds`
// controls how many different seed cells to try (1 = single DFS;
// more = try additional seeds and keep the one minimising tear count).
LatticeUnfolding unfold_iDT(const DelaunayTriangulation& D,
                            const std::vector<CellPlacement>& cells,
                            int start_cell_id = -1,
                            int dfs_max_seeds = 0);

// TikZ picture of a LatticeUnfolding plus a text debug report (list of
// duplicate-position hex, missing hex, and per-cell placements).
void dump_lattice_unfolding_tikz(const LatticeUnfolding& U,
                                 const Triangulation& T_sorted,
                                 std::ostream& tikz_os,
                                 std::ostream& debug_os,
                                 double scale = 0.4);

// Local "F + its 3 neighbours" tikz visual, with each neighbour
// unfolded into F's frame across the shared edge.
//
// neighbours[k] is the cell on the OTHER side of F's k-th edge
// (k = 0, 1, 2 for edges c0c1 / c1c2 / c2c0).  `valid = false` means
// "no neighbour" (e.g. dead face) -- distinct from placement.ok = false
// which would mean placement attempt failed.
struct NeighbourCell {
    CellPlacement placement;
    bool          valid = false;
};
void dump_neighbour_unfolding_tikz(const CellPlacement& F,
                                   const std::array<NeighbourCell, 3>& neighbours,
                                   const Triangulation& T_sorted,
                                   std::ostream& os,
                                   double scale = 0.5);
