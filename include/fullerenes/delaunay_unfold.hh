#pragma once

#include "fullerenes/eisenstein.hh"
#include "fullerenes/eisenstein_paint_tables.hh"   // ParamTablesView
#include "fullerenes/triangulation.hh"
#include "fullerenes/delaunay.hh"

#include <iosfwd>
#include <map>
#include <utility>
#include <vector>

// =====================================================================
// iDT face unfolding into a single Eisenstein-lattice frame.
//
// Diagnostic / visualisation only.  Lays out live iDT cells into one
// Eisenstein-plane drawing via spanning-tree DFS over face adjacency:
//   - Each half-edge h is either "placed" (the cell to its left has a
//     global placement) or "open".
//   - Pop an open half-edge, place the cell to its left with the
//     lattice isometry that pins the shared edge -- identified by its
//     CYCLE SLOT in the cell's half-edge 3-cycle, never by corner
//     label (labels can repeat on delta-complex cells) -- to the
//     already-placed global endpoint positions of its two cones.
//   - If a cone gets placed at a different position than a previous
//     placement, count it as a "tear" -- this is intrinsic to delta-
//     complex unfoldings (Gauss-Bonnet) and cannot be avoided, only
//     re-routed by choosing a different DFS spanning tree.
//
// Per-cell input is a ParamTablesView (eisenstein_paint's flat chart
// tables): cell f contributes its CCW corner ids (corner_ids), corner
// positions in its own local Eisenstein frame (frame_points), and its
// claimed lattice points (cell_entries); dead or uncharted slots
// (!cell_live) are skipped.  Build one with eisenstein_paint::
// parametrize and pass view(); V must be the parametrization of the
// SAME complex D (unfold_iDT refuses V.nf != D.nf).
// =====================================================================

// Result of unfold_iDT.
struct LatticeUnfolding {
    // Sorted labels < n_cones are cones (copied from the chart tables;
    // THE cone predicate for readers of this result).
    int n_cones = 0;

    // All DISTINCT global positions seen for each cone across placed
    // cells.  On a delta-complex with cone angle deficits some cones
    // end up at 2+ positions (= spanning-tree cross-edges force "tears"
    // in the unfolding); the global net carries all of them.
    std::map<int, std::vector<Eisenstein>> cone_positions;

    // For each placed iDT cell: corners (c0, c1, c2) and their global
    // positions in the unfolding's spanning-tree frame, plus the
    // lattice points translated into those global coords.
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

// Unfold live iDT cells into a single Eisenstein-lattice frame via a
// spanning-tree DFS over cell adjacency.  Each cell is placed AT MOST
// once (via the half-edge that first reaches it); a cell is SKIPPED
// when no orientation-preserving gluing of its entering edge exists,
// i.e. when the local edge vector and the first-sighting global one
// are not lattice ASSOCIATES: after a tear, the label-keyed anchors
// can be ANY equal-norm vector -- a different norm (the tear crosses
// this edge), a different ideal-class orbit (no lattice alignment at
// all; the smallest norm with several classes is 49), or the
// conjugate orbit (only a REFLECTED alignment exists, and a
// development never mirrors -- a mirrored cell is never drawn).  A
// skipped cell's outgoing edges are not pushed, which is how subtrees
// strand.  (The anchors themselves always exist: the entering edge's
// endpoints are corners of the already-placed twin cell -- an
// invariant, violated = throw.)  So n_cells <= count of live charted
// cells; equality holds on tear-free unfoldings and is NOT guaranteed
// in general.  Cone deficits manifest as "tears": a cone may end up
// at multiple positions across placed cells.
//
// `start_cell_id == -1` means pick the first charted cell.
// `dfs_max_seeds` controls how many different seed cells to try
// (0 or 1 = single DFS; more = try additional seeds and keep the one
// placing the most cells, ties broken by minimal tear count).
//
// @pre V.nf == D.nf -- throws otherwise (a cell-count check; that V
// really is THE parametrization of D is the caller's contract, gated
// upstream by build_atlas's chart/cycle agreement check).
LatticeUnfolding unfold_iDT(const DelaunayTriangulation& D,
                            const eisenstein_paint::ParamTablesView& V,
                            int start_cell_id = -1,
                            int dfs_max_seeds = 0);

// TikZ picture of a LatticeUnfolding plus a text debug report (list of
// duplicate-position non-cone vertices, missing vertices, and per-cell
// placements).
void dump_lattice_unfolding_tikz(const LatticeUnfolding& U,
                                 const Triangulation& T_sorted,
                                 std::ostream& tikz_os,
                                 std::ostream& debug_os,
                                 double scale = 0.4);

// Local "F + its 3 neighbours" tikz visual, with each neighbour
// unfolded into cell f's frame across the shared edge.  The three
// neighbours and the shared-edge slot pairing are derived from D's
// half-edge structure (slot-keyed, multi-edge-sound); uncharted or
// dead neighbours are skipped.
//
// @pre V.nf == D.nf (V parametrizes D) -- throws otherwise.
void dump_neighbour_unfolding_tikz(const DelaunayTriangulation& D,
                                   const eisenstein_paint::ParamTablesView& V,
                                   int f,
                                   std::ostream& os,
                                   double scale = 0.5);
