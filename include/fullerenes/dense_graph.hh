#pragma once

#include <vector>
#include <span>
#include <array>
#include <tuple>
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <initializer_list>
#include <type_traits>
#include <utility>
#include <stdexcept>
#include <sstream>
#include <string>

#include "planar_csr.hh"
#include "fullerenes/batch/batch.hh"

namespace Spanify {

template<typename K> struct OwnedDenseGraph;  // forward declaration

// ---------------------------------------------------------------------------
// A violated precondition of embedded-edge surgery (insert_edge_at /
// remove_edge_at / the slot-addressed cores, and GraphView's vertex-pair
// wrappers).  The caller asserted something about the rotation system that
// does not hold, so the operation has no meaningful result -- the exception
// channel, not a modeled outcome (the unoriented_surface_error precedent in
// graphview.hh).  The reason is a CLOSED code, so a handler (or a test) can
// ask which contract was violated instead of parsing prose; the what()
// string additionally carries the operation and the actual row content, so
// the message reads as a diagnosis rather than a shrug.
// ---------------------------------------------------------------------------
struct graph_surgery_error : public std::logic_error {
    enum class Code {
        VertexOutOfRange,        // an endpoint is not a vertex of the graph
        SelfLoop,                // u == v (simple graphs only)
        RowFull,                 // an endpoint's degree is already dmax
        SlotOutOfRange,          // a slot outside [0, deg] / [0, deg)
        AsymmetricAdjacency,     // one direction of the edge exists, the other not
        EdgePresent,             // the edge exists where the contract needs it absent
        NotMutualPair,           // the two slots do not name each other's reverse
        NotLiveArc,              // an arcix names no live arc
        SuccessorNotNeighbour,   // a named successor is absent from the row
        TwinAbsent,              // reverse_arc on a graph with no twin table
    };
    Code code;
    long long u, v;   // the endpoints the operation was called with (-1: n/a)
    graph_surgery_error(Code code, const std::string& what,
                        long long u = -1, long long v = -1)
        : std::logic_error(what), code(code), u(u), v(v) {}
};

// ---------------------------------------------------------------------------
// The one suffix-shift word: open or close a one-slot gap in a slot-indexed
// row.  The surgery members below express BOTH their neighbour-row and
// twin-row shifts through these, and a caller that carries per-arc
// attributes across surgery applies the identical word to its own rows (the
// surgery @post names slot and pre-call degree, which are exactly these
// arguments) -- so the off-by-one arithmetic lives here once.
//
// @anchor row-open-gap
// @pre  gap_fits: int(row.size()) > old_deg
// @post entries [slot, old_deg) sit at [slot+1, old_deg+1); row[slot] is
//       stale for the caller to write
template<typename T>
void row_open_gap(std::span<T> row, int slot, int old_deg) {
    for (int i = old_deg; i > slot; --i) row[i] = row[i-1];
}

// @anchor row-close-gap
// @pre  slot_live: slot < old_deg
// @post entries [slot+1, old_deg) sit at [slot, old_deg-1); row[old_deg-1]
//       is stale for the caller to clear
template<typename T>
void row_close_gap(std::span<T> row, int slot, int old_deg) {
    for (int i = slot; i < old_deg-1; ++i) row[i] = row[i+1];
}

// ---------------------------------------------------------------------------
// RSRAdjacencyView: non-owning view over fixed-stride flat adjacency data.
//
// Each vertex has dmax slots; a degree counter tracks how many are filled.
// Supports O(dmax) insert/remove at arbitrary cyclic positions.
// Optionally carries a twin array for O(1) reverse-arc lookup.
//
// operator[] returns a span over the row's active entries.
//
// Template parameter:
//   K -- index type: int32_t (CPU), uint16_t (GPU)
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct RSRAdjacencyView {
    using node_type = K;
    using arc_t   = std::pair<K, K>;       // {source, target}
    using arcix_t = std::pair<K, uint8_t>; // {source, position_in_row}

    static constexpr uint8_t default_dmax = 3;  // cubic
    // The named "no such slot": what vacated_corner returns for an endpoint
    // whose row the removal emptied, and what find_arc yields for an absent
    // arc.  No live arc carries it (deg <= dmax <= 255 caps live slots).
    static constexpr uint8_t no_slot = uint8_t(-1);

    K N = 0;
    int dmax = 0;
    std::span<K> neighbours;     // [N * dmax] -- neighbor IDs, row-major
    std::span<uint8_t> deg;      // [N] -- current degree of each vertex
    std::span<uint8_t> twin;     // [N * dmax] -- twin position per arc (optional)
    //   twin[u * dmax + i] = j  means: neighbours[v][j] == u,
    //   where v = neighbours[u][i]. Empty when not computed.

    // --- Constructors ---

    RSRAdjacencyView() = default;

    // View constructor: wrap existing spans.
    RSRAdjacencyView(K N, int dmax, std::span<K> neighbours, std::span<uint8_t> deg,
                     std::span<uint8_t> twin = {})
        : N(N), dmax(dmax), neighbours(neighbours), deg(deg), twin(twin) {}

    // -----------------------------------------------------------------------
    // Batchability contract (see include/fullerenes/batch/batchable.hh).
    //
    // Tuple of span-typed fields in canonical order: {neighbours, deg, twin}.
    // get_element_counts gives each field's element count for ONE batch entry
    // -- absolute counts, so N-proportional and constant-size fields ride the
    // same law:
    //   neighbours : N * dmax
    //   deg        : N
    //   twin       : N * dmax  (empty until computed)
    //
    // Derived views that add span fields (e.g. PolyhedronView adds `points`)
    // override n_fields / to_tuple / get_element_counts to extend the tuple.
    // -----------------------------------------------------------------------
    static constexpr std::size_t n_fields = 3;

    auto to_tuple() {
        return std::forward_as_tuple(neighbours, deg, twin);
    }
    auto to_tuple() const {
        return std::forward_as_tuple(neighbours, deg, twin);
    }

    static constexpr std::array<std::size_t, n_fields>
    get_element_counts(int N, int dmax) {
        return { (std::size_t)N * dmax, (std::size_t)N, (std::size_t)N * dmax };
    }

    // --- operator[] returns span over active entries ---

    std::span<K> operator[](K u) {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }
    std::span<const K> operator[](K u) const {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }

    // --- Span-based accessors ---

    std::span<const K> nbrs(K u) const {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }

    std::span<K> nbrs_mut(K u) {
        return {neighbours.data() + u * dmax, (size_t)deg[u]};
    }

    int degree(K u) const { return deg[u]; }

    // --- Twin query ---
    bool has_twin() const { return !twin.empty(); }

    // --- Vertex-pair navigation (O(degree) per call) ---

    int find(K u, K v) const {
        const K* row = neighbours.data() + u * dmax;
        for (int i = 0; i < deg[u]; ++i)
            if (row[i] == v) return i;
        return -1;
    }

    // --- Arc-index navigation (O(1) per call, except find_arc which is O(degree)) ---

    arcix_t find_arc(K u, K v) const {
        return {u, uint8_t(find(u, v))};
    }

    // Flat arc id in [0, N*dmax): the dense canonical index of the directed arc at
    // position i in u's row (the arc u -> neighbours[u*dmax+i]). It indexes both
    // `neighbours` and any arc-attribute vector<T>(N*dmax) identically, so an arc-flag
    // array replaces a hash map keyed by arc_t -- denser, allocation-free, and parallel
    // arcs / self-loops stay distinct (an (u,v) key conflates them). The arc VALUE is
    // recoverable O(1): source = id/dmax (arc_of), target = neighbours[id]. size_t so
    // the index can never overflow before the (size_t-sized) backing buffer does.
    size_t  arcid(K u, int i) const { return size_t(u) * dmax + i; }
    size_t  arcid(arcix_t a)  const { return arcid(a.first, a.second); }
    arcix_t arc_of(size_t id) const { return {K(id / dmax), uint8_t(id % dmax)}; }

    K target(arcix_t a) const {
        return neighbours[a.first * dmax + a.second];
    }

    arcix_t next(arcix_t a) const {
        auto [u, i] = a;
        uint8_t d = deg[u];
        return {u, uint8_t((i + 1) % d)};
    }

    arcix_t prev(arcix_t a) const {
        auto [u, i] = a;
        uint8_t d = deg[u];
        return {u, uint8_t((i + d - 1) % d)};
    }

    // @pre twin: has_twin()  (violation throws graph_surgery_error{TwinAbsent})
    arcix_t reverse_arc(arcix_t a) const {
        if (!has_twin())
            surgery_fail("reverse_arc", graph_surgery_error::Code::TwinAbsent,
                         a.first, K(-1));
        auto [u, i] = a;
        K v = neighbours[arcid(u, i)];
        uint8_t j = twin[arcid(u, i)];
        return {v, j};
    }

    arcix_t next_on_face(arcix_t a) const {
        return next(reverse_arc(a));
    }

    // =======================================================================
    // Embedded-edge surgery
    //
    // The operations in this section edit the rotation system as a SYMMETRIC
    // pair -- both directed arcs of an edge in one call -- and are the only
    // mutators that maintain the twin table.  What "valid twin" means is
    // defined ONCE, executably, by twin_is_valid() below; every contract in
    // this section cites it by name.
    //
    // Positions are corners: insert_edge_at(a_u, a_v) lands the new edge in
    // the corner AFTER arc a_u in u's rotation and after a_v in v's.  The
    // slot-addressed cores below them are the ONE implementation (SSoT) --
    // the corner forms and GraphView's vertex-pair wrappers all resolve to
    // them.  Simple graphs only: no self-loops, no parallel arcs (the
    // hierarchy's open question on loops/parallel arcs is Delta Complex
    // territory -- see graphview.hh).  The cores enforce this on insertion
    // and on self-loops; a pre-existing parallel arc is outside their
    // contract.
    //
    // Contract violations THROW graph_surgery_error, whose closed Code names
    // the violated clause and whose message carries the offending rows;
    // nothing is silently repositioned or skipped, and every check runs
    // BEFORE the first write (a throwing call leaves the graph unchanged).
    // Host only (throwing); device code has no caller for these.
    //
    // The per-row primitives further down (push_back / insert_at / erase_at
    // / assign_row / ...) edit one row at a time, pass through asymmetric
    // intermediates, and do NOT maintain twin -- @post of each: twin entries
    // covering the edited row (and any entry pointing into it) are STALE;
    // detectable by twin_is_valid(), repaired only by recomputing the twin.
    // =======================================================================

    using Code = graph_surgery_error::Code;

    // THE twin-validity predicate -- the definition every @pre/@post twin
    // clause in this section cites.  Empty twin is vacuously valid; else
    // every live arc's entry names a live slot of a valid target row that
    // points back at it (its reverse), both ways (involution).  The RSR
    // sibling of validate(PlanarCSR) in planar_csr.hh.  Total: survives
    // arbitrary garbage (out-of-range targets are reported false, never
    // read through).  O(N*dmax) -- a validator, not a hot-path check.
    bool twin_is_valid() const {
        if (!has_twin()) return true;
        for (K u = 0; u < N; ++u)
            for (int i = 0; i < deg[u]; ++i) {
                const K v = neighbours[arcid(u, i)];
                const int j = twin[arcid(u, i)];
                if (size_t(v) >= size_t(N))     return false;
                if (j >= deg[v])                return false;
                if (neighbours[arcid(v, j)] != u) return false;
                if (twin[arcid(v, j)] != i)     return false;
            }
        return true;
    }

    // The linear slot at which an arc inserted into the corner after a
    // lands.  Deliberately NOT next(a): next wraps slot deg-1 back to 0,
    // while the corner after the last arc is the append slot deg.
    int slot_after(arcix_t a) const { return a.second + 1; }

    // New edge u--v landing in the corner after a_u (in u's rotation) and
    // after a_v (in v's).
    // @anchor rsr-insert-edge-at
    // @pre  live_u: size_t(a_u.first) < size_t(N) && int(a_u.second) < degree(a_u.first)
    // @pre  live_v: size_t(a_v.first) < size_t(N) && int(a_v.second) < degree(a_v.first)
    // @pre  (the slot core's @pre, at slots slot_after(a_u) / slot_after(a_v))
    // @post target(result.first) == a_v.first && target(result.second) == a_u.first
    // @post twin: twin_is_valid()
    // @throws graph_surgery_error{NotLiveArc} and the slot core's codes
    std::pair<arcix_t, arcix_t> insert_edge_at(arcix_t a_u, arcix_t a_v) {
        require_live_arc("insert_edge_at", a_u);
        require_live_arc("insert_edge_at", a_v);
        return insert_edge_slots(a_u.first, slot_after(a_u),
                                 a_v.first, slot_after(a_v));
    }

    // Remove the edge carried by live arc a (and its reverse, resolved via
    // twin in O(1) when present, one find otherwise).  Returns the two
    // corners it vacated -- each named by the arc now PRECEDING the gap, in
    // post-removal coordinates -- so a removal composes into the next
    // insert_edge_at without a find.  An endpoint whose degree drops to 0
    // has no corner left; its arcix comes back with slot no_slot.
    // @anchor rsr-remove-edge-at
    // @pre  live: size_t(a.first) < size_t(N) && int(a.second) < degree(a.first)
    // @pre  (the slot core's @pre)
    // @post result == {vacated_corner(u, s_u), vacated_corner(v, s_v)}
    // @post twin: twin_is_valid()
    // @throws graph_surgery_error{NotLiveArc} and the slot core's codes
    std::pair<arcix_t, arcix_t> remove_edge_at(arcix_t a) {
        require_live_arc("remove_edge_at", a);
        const arcix_t r = has_twin() ? reverse_arc(a)
                                     : find_arc(target(a), a.first);
        return remove_edge_slots(a.first, a.second, r.first, r.second);
    }

    // --- Slot-addressed cores (the one surgery implementation) -----------
    // The corner forms above and GraphView's vertex-pair wrappers resolve
    // their addressing to linear slots and call these.  Slot semantics:
    // the new arc lands AT slot s (0 <= s <= deg); existing arcs [s, deg)
    // shift up.  Arc ids in the shifted suffixes are renamed by +1/-1 --
    // mirror caller-owned arc-attribute rows with row_open_gap /
    // row_close_gap (declared before this struct), whose arguments are
    // exactly this call's (slot, pre-call degree).

    // @anchor rsr-insert-edge-slots
    // @pre  vertices: size_t(u) < size_t(N) && size_t(v) < size_t(N)
    // @pre  distinct: u != v
    // @pre  room:     degree(u) < dmax && degree(v) < dmax
    // @pre  slots:    0 <= s_u && s_u <= degree(u) && 0 <= s_v && s_v <= degree(v)
    // @pre  absent:   find(u, v) < 0 && find(v, u) < 0
    // @pre  twin:     twin_is_valid()
    // @post result == std::pair{arcix_t{u, s_u}, arcix_t{v, s_v}} && target(result.first) == v
    // @post twin:     twin_is_valid()
    // @throws graph_surgery_error{VertexOutOfRange, SelfLoop, RowFull,
    //         SlotOutOfRange, AsymmetricAdjacency, EdgePresent}
    std::pair<arcix_t, arcix_t> insert_edge_slots(K u, int s_u, K v, int s_v) {
        require_vertices("insert_edge_slots", u, v);
        if (u == v)
            surgery_fail("insert_edge_slots", Code::SelfLoop, u, v);
        if (deg[u] >= dmax || deg[v] >= dmax)
            surgery_fail("insert_edge_slots", Code::RowFull, u, v);
        if (s_u < 0 || s_u > deg[u] || s_v < 0 || s_v > deg[v])
            surgery_fail("insert_edge_slots", Code::SlotOutOfRange, u, v);
        const bool uv = find(u, v) >= 0, vu = find(v, u) >= 0;
        if (uv != vu)
            surgery_fail("insert_edge_slots", Code::AsymmetricAdjacency, u, v);
        if (uv)
            surgery_fail("insert_edge_slots", Code::EdgePresent, u, v);

        row_open_slot(u, s_u);
        row_open_slot(v, s_v);
        neighbours[arcid(u, s_u)] = v;
        neighbours[arcid(v, s_v)] = u;
        if (has_twin()) {
            twin[arcid(u, s_u)] = uint8_t(s_v);
            twin[arcid(v, s_v)] = uint8_t(s_u);
        }
        return {{u, uint8_t(s_u)}, {v, uint8_t(s_v)}};
    }

    // @anchor rsr-remove-edge-slots
    // @pre  vertices: size_t(u) < size_t(N) && size_t(v) < size_t(N)
    // @pre  distinct: u != v
    // @pre  mutual:   0 <= s_u && s_u < degree(u) && neighbours[arcid(u, s_u)] == v
    //              && 0 <= s_v && s_v < degree(v) && neighbours[arcid(v, s_v)] == u
    // @pre  twin:     twin_is_valid()
    // @post result == std::pair{vacated_corner(u, s_u), vacated_corner(v, s_v)}
    //       evaluated BEFORE the removal
    // @post twin:     twin_is_valid()
    // @throws graph_surgery_error{VertexOutOfRange, SelfLoop, NotMutualPair}
    std::pair<arcix_t, arcix_t> remove_edge_slots(K u, int s_u, K v, int s_v) {
        require_vertices("remove_edge_slots", u, v);
        if (u == v)   // one loop occupies ONE slot; closing it twice would
            surgery_fail("remove_edge_slots", Code::SelfLoop, u, v);
        if (s_u < 0 || s_u >= deg[u] || neighbours[arcid(u, s_u)] != v ||
            s_v < 0 || s_v >= deg[v] || neighbours[arcid(v, s_v)] != u)
            surgery_fail("remove_edge_slots", Code::NotMutualPair, u, v);
        const arcix_t c_u = vacated_corner(u, s_u), c_v = vacated_corner(v, s_v);
        row_close_slot(u, s_u);
        row_close_slot(v, s_v);
        return {c_u, c_v};
    }

    // --- Surgery internals ------------------------------------------------

    // The corner left behind by removing the arc at slot s of u's row, in
    // POST-removal coordinates: the arc preceding the gap.  s == 0 wraps to
    // the (shifted) last slot -- and when the removal empties the row
    // (deg == 1 forces s == 0), that same arithmetic yields d_new-1 == -1
    // == no_slot, which is exactly the "no corner left" answer.
    // @pre  live: s < degree(u)  (called before the removal's close)
    arcix_t vacated_corner(K u, int s) const {
        const int d_new = deg[u] - 1;
        return {u, uint8_t(s == 0 ? d_new - 1 : s - 1)};
    }

    // Open slot s in u's row: one row_open_gap on the neighbour row (and the
    // twin row in lockstep), then retarget each shifted arc's REVERSE entry
    // -- the twin table itself says where the other side is, so each repair
    // is one O(1) write (the step that is impossible to do cheaply without
    // twin).  Leaves slot s's content stale; the caller writes it.
    void row_open_slot(K u, int s) {
        const int d = deg[u];
        row_open_gap(row_slots(u), s, d);
        if (has_twin()) {
            row_open_gap(twin_slots(u), s, d);
            retarget_reverses(u, s + 1, d + 1);
        }
        deg[u]++;
    }

    // Close slot s in u's row: the mirror image of row_open_slot, with the
    // same O(1)-per-slot reverse retargeting.  Clears both vacated slots to
    // their named dead values.
    void row_close_slot(K u, int s) {
        const int d = deg[u];
        row_close_gap(row_slots(u), s, d);
        if (has_twin()) {
            row_close_gap(twin_slots(u), s, d);
            retarget_reverses(u, s, d - 1);
            twin_slots(u)[d-1] = no_slot;
        }
        row_slots(u)[d-1] = K(-1);
        deg[u]--;
    }

    // Point each reverse of u's arcs in slots [lo, hi) back at its (shifted)
    // slot.  @pre the twin row of u is already shifted alongside.
    void retarget_reverses(K u, int lo, int hi) {
        for (int i = lo; i < hi; ++i)
            twin[arcid(neighbours[arcid(u, i)], twin[arcid(u, i)])] = uint8_t(i);
    }

    // --- Surgery failure channel ------------------------------------------

    void require_vertices(const char* op, K u, K v) const {
        if (size_t(u) >= size_t(N) || size_t(v) >= size_t(N))
            surgery_fail(op, Code::VertexOutOfRange, u, v);
    }

    void require_live_arc(const char* op, arcix_t a) const {
        if (size_t(a.first) >= size_t(N) || int(a.second) >= deg[a.first])
            surgery_fail(op, Code::NotLiveArc, a.first, K(-1),
                         " (arc slot " + std::to_string(int(a.second)) + ")");
    }

    // The ONE builder for every surgery diagnosis: code -> sentence, plus
    // the actual rows of whichever endpoints are in range.
    [[noreturn]] void surgery_fail(const char* op, Code c, K u, K v,
                                   const std::string& detail = {}) const {
        const char* why = "";
        switch (c) {
            case Code::VertexOutOfRange:      why = "vertex out of range"; break;
            case Code::SelfLoop:              why = "self-loop (simple graphs only)"; break;
            case Code::RowFull:               why = "row full (degree == dmax)"; break;
            case Code::SlotOutOfRange:        why = "slot out of range"; break;
            case Code::AsymmetricAdjacency:   why = "asymmetric adjacency"; break;
            case Code::EdgePresent:           why = "edge already present"; break;
            case Code::NotMutualPair:         why = "arcs are not a mutual pair "
                                                    "(asymmetric adjacency, stale twin, or wrong slots)"; break;
            case Code::NotLiveArc:            why = "not a live arc"; break;
            case Code::SuccessorNotNeighbour: why = "the named successor is not a neighbour"; break;
            case Code::TwinAbsent:            why = "no twin table (reverse_arc needs one)"; break;
        }
        std::string msg = std::string(op) + "(" + std::to_string((long long)u) + "," +
                          std::to_string((long long)v) + "): " + why + detail;
        if (size_t(u) < size_t(N)) msg += "; row of " + std::to_string((long long)u) + " = " + row_string(u);
        if (size_t(v) < size_t(N)) msg += "; row of " + std::to_string((long long)v) + " = " + row_string(v);
        throw graph_surgery_error(c, msg, (long long)u, (long long)v);
    }

    // A row as text, in operator<<'s comma-separated shape, so a diagnostic
    // and a dump of the same graph look alike.
    std::string row_string(K u) const {
        std::ostringstream s;
        s << "[";
        for (int i = 0; i < deg[u]; ++i)
            s << (i ? "," : "") << (long long)neighbours[arcid(u, i)];
        s << "]";
        return s.str();
    }

    // --- Full-width row spans (all dmax slots, live prefix + padding) -----
    // nbrs/operator[] span the LIVE prefix; the surgery shifts need one
    // slot beyond it, which is what these provide.
    std::span<K>       row_slots(K u)  { return {neighbours.data() + arcid(u, 0), (size_t)dmax}; }
    std::span<uint8_t> twin_slots(K u) { return {twin.data()       + arcid(u, 0), (size_t)dmax}; }

    // --- Per-row mutation (no reallocation) ---
    // Construction-phase vocabulary: these edit ONE row, pass through
    // asymmetric intermediates, and do NOT maintain twin.
    // @post of each: twin entries covering this row, and any twin entry
    // pointing into it, are STALE (twin_is_valid() detects this); only a
    // twin recompute repairs it.  On a graph whose twin must stay valid,
    // use the surgery operations above.

    // @throws graph_surgery_error{RowFull} on degree overflow
    void push_back(K u, K v) {
        if (deg[u] >= dmax)
            surgery_fail("push_back", graph_surgery_error::Code::RowFull, u, v);
        neighbours[arcid(u, deg[u])] = v;
        deg[u]++;
    }

    // @throws graph_surgery_error{RowFull} on degree overflow
    void insert_at(K u, K v, int pos) {
        if (deg[u] >= dmax)
            surgery_fail("insert_at", graph_surgery_error::Code::RowFull, u, v);
        row_open_gap(row_slots(u), pos, deg[u]);
        neighbours[arcid(u, pos)] = v;
        deg[u]++;
    }

    void erase_at(K u, int pos) {
        row_close_gap(row_slots(u), pos, deg[u]);
        row_slots(u)[deg[u]-1] = K(-1);
        deg[u]--;
    }

    void clear_row(K u) { deg[u] = 0; }

    void resize_row(K u, int n) {
        assert(n <= dmax);
        for (int i = deg[u]; i < n; ++i)
            neighbours[u * dmax + i] = K(-1);
        deg[u] = uint8_t(n);
    }

    void assign_row(K u, std::initializer_list<K> v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        int i = 0;
        for (K x : v) neighbours[u * dmax + i++] = x;
    }

    void assign_row(K u, const std::vector<K>& v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        for (int i = 0; i < int(v.size()); ++i)
            neighbours[u * dmax + i] = v[i];
    }

    void assign_row(K u, std::span<const K> v) {
        assert(int(v.size()) <= dmax);
        deg[u] = uint8_t(v.size());
        for (int i = 0; i < int(v.size()); ++i)
            neighbours[u * dmax + i] = v[i];
    }

    // --- Backward compat with vector<vector<>> interface ---

    int size() const { return N; }

    bool operator==(const RSRAdjacencyView& other) const {
        if (N != other.N || dmax != other.dmax) return false;
        for (K u = 0; u < N; ++u) {
            auto a = (*this)[u], b = other[u];
            if (a.size() != b.size()) return false;
            if (!std::equal(a.begin(), a.end(), b.begin())) return false;
        }
        return true;
    }

    std::vector<std::vector<K>> to_vectors() const {
        std::vector<std::vector<K>> adj(N);
        for (K v = 0; v < N; ++v) {
            adj[v].resize(deg[v]);
            for (int i = 0; i < deg[v]; ++i)
                adj[v][i] = neighbours[v * dmax + i];
        }
        return adj;
    }
};

// Backward-compat alias
template<typename K = int32_t>
using DenseGraph = RSRAdjacencyView<K>;

// TC verification: RSRAdjacencyView with default K is trivially copyable.
static_assert(std::is_trivially_copyable_v<RSRAdjacencyView<int32_t>>,
    "RSRAdjacencyView<int32_t> must be trivially copyable");
static_assert(std::is_trivially_copyable_v<RSRAdjacencyView<uint16_t>>,
    "RSRAdjacencyView<uint16_t> must be trivially copyable");

// ---------------------------------------------------------------------------
// RSRPolyhedronView<T,K>: RSR adjacency extended with per-vertex 3D coordinates.
//
// Carries the same three adjacency fields (neighbours, deg, twin) as
// RSRAdjacencyView plus a points span of std::array<T,3> (xyz).  This makes
// (graph + geometry) atomically transferable through Batch<V>/BatchQueue<V>.
//
// Canonical tuple order: {neighbours, deg, twin, points}.  The first three
// fields' element counts match RSRAdjacencyView, so a BatchView<RSRPolyhedronView>
// can be sliced into a BatchView<RSRAdjacencyView> via as_adjacency_view().
// ---------------------------------------------------------------------------
template<typename T, typename K = int32_t>
struct RSRPolyhedronView : RSRAdjacencyView<K> {
    using typename RSRAdjacencyView<K>::node_type;

    static constexpr uint8_t default_dmax = 3;  // cubic

    std::span<std::array<T,3>> points;          // [N] -- per-vertex xyz

    RSRPolyhedronView() = default;

    RSRPolyhedronView(K N, int dmax,
                      std::span<K> neighbours, std::span<uint8_t> deg,
                      std::span<std::array<T,3>> pts,
                      std::span<uint8_t> twin = {})
        : RSRAdjacencyView<K>(N, dmax, neighbours, deg, twin), points(pts) {}

    // Construct from an existing adjacency view + points span.
    RSRPolyhedronView(const RSRAdjacencyView<K>& g,
                      std::span<std::array<T,3>> pts)
        : RSRAdjacencyView<K>(g), points(pts) {}

    // -- Batchability contract (extends adjacency with `points`) ------------
    static constexpr std::size_t n_fields = 4;

    auto to_tuple() {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }
    auto to_tuple() const {
        return std::forward_as_tuple(this->neighbours, this->deg, this->twin, points);
    }

    static constexpr std::array<std::size_t, n_fields>
    get_element_counts(int N, int dmax) {
        return { (std::size_t)N * dmax, (std::size_t)N, (std::size_t)N * dmax, (std::size_t)N };
    }
};

static_assert(std::is_trivially_copyable_v<RSRPolyhedronView<float,  uint16_t>>,
    "RSRPolyhedronView<float,uint16_t> must be trivially copyable");
static_assert(std::is_trivially_copyable_v<RSRPolyhedronView<double, uint16_t>>,
    "RSRPolyhedronView<double,uint16_t> must be trivially copyable");

// ---------------------------------------------------------------------------
// Graph-prefix slicing: view a BatchView<RSRPolyhedronView<T,K>> as a
// BatchView<RSRAdjacencyView<K>> over the shared adjacency fields.
// The points span is available via points_span().
// ---------------------------------------------------------------------------
template<typename T, typename K>
batch::BatchView<RSRAdjacencyView<K>>
as_adjacency_view(batch::BatchView<RSRPolyhedronView<T,K>> pv) {
    const auto& s = pv.spans();
    typename batch::detail::span_tuple_t<RSRAdjacencyView<K>> adj_spans{
        std::get<0>(s), std::get<1>(s), std::get<2>(s)
    };
    return batch::BatchView<RSRAdjacencyView<K>>(pv.N(), pv.dmax(), pv.size(), adj_spans);
}

template<typename T, typename K>
std::span<std::array<T,3>>
points_span(batch::BatchView<RSRPolyhedronView<T,K>> pv) {
    return std::get<3>(pv.spans());
}

// ---------------------------------------------------------------------------
// OwnedDenseGraph: owning version with vector storage.
// Inherits all view/mutation methods from RSRAdjacencyView.
// ---------------------------------------------------------------------------
template<typename K = int32_t>
struct OwnedDenseGraph : RSRAdjacencyView<K> {
    std::vector<K> owned_neighbours;
    std::vector<uint8_t> owned_deg;

    void repoint() {
        this->neighbours = std::span<K>(owned_neighbours.data(), owned_neighbours.size());
        this->deg = std::span<uint8_t>(owned_deg.data(), owned_deg.size());
    }

    // Bring base push_back(K,K) into scope (otherwise hidden by push_back(vector))
    using RSRAdjacencyView<K>::push_back;

    // --- Constructors ---

    OwnedDenseGraph() = default;

    // Explicit: N vertices, given stride, all empty.
    explicit OwnedDenseGraph(K N, int dmax)
        : owned_neighbours(N * dmax, K(-1)), owned_deg(N, 0) {
        this->N = N; this->dmax = dmax; repoint();
    }

    // Fill constructor: N vertices, each initialized with the given row.
    OwnedDenseGraph(K N, const std::vector<K>& initial_row)
        : owned_neighbours(N * int(initial_row.size()), K(-1)),
          owned_deg(N, uint8_t(initial_row.size())) {
        this->N = N; this->dmax = int(initial_row.size()); repoint();
        for (K v = 0; v < N; ++v)
            for (int i = 0; i < this->dmax; ++i)
                this->neighbours[v * this->dmax + i] = initial_row[i];
    }

    // Brace-init constructor: OwnedDenseGraph{{1,2,3},{4,5,6},...}.
    OwnedDenseGraph(std::initializer_list<std::initializer_list<K>> adj) {
        this->N = K(adj.size());
        this->dmax = 0;
        for (auto& row : adj) this->dmax = std::max(this->dmax, int(row.size()));
        owned_neighbours.resize(this->N * this->dmax, K(-1));
        owned_deg.resize(this->N, 0);
        repoint();
        int v = 0;
        for (auto& row : adj) {
            this->deg[v] = uint8_t(row.size());
            int i = 0;
            for (K x : row) this->neighbours[v * this->dmax + i++] = x;
            v++;
        }
    }

    // Converting constructor from vector<vector<K>>.
    OwnedDenseGraph(const std::vector<std::vector<K>>& adj) {
        this->N = K(adj.size());
        this->dmax = 0;
        for (auto& row : adj) this->dmax = std::max(this->dmax, int(row.size()));
        owned_neighbours.resize(this->N * this->dmax, K(-1));
        owned_deg.resize(this->N, 0);
        repoint();
        for (int v = 0; v < int(adj.size()); ++v) {
            this->deg[v] = uint8_t(adj[v].size());
            for (int i = 0; i < int(adj[v].size()); ++i)
                this->neighbours[v * this->dmax + i] = adj[v][i];
        }
    }

    // Converting constructor from RSRAdjacencyView (copies data).
    OwnedDenseGraph(const RSRAdjacencyView<K>& v) {
        this->N = v.N; this->dmax = v.dmax;
        if (v.N > 0 && v.neighbours.data()) {
            owned_neighbours.assign(v.neighbours.begin(), v.neighbours.end());
            owned_deg.assign(v.deg.begin(), v.deg.end());
        }
        repoint();
    }

    // --- Rule of 5 ---

    OwnedDenseGraph(const OwnedDenseGraph& o)
        : owned_neighbours(o.owned_neighbours), owned_deg(o.owned_deg) {
        this->N = o.N; this->dmax = o.dmax; repoint();
    }

    OwnedDenseGraph(OwnedDenseGraph&& o) noexcept
        : owned_neighbours(std::move(o.owned_neighbours)), owned_deg(std::move(o.owned_deg)) {
        this->N = o.N; this->dmax = o.dmax; repoint();
        o.neighbours = {}; o.deg = {}; o.N = 0;
    }

    OwnedDenseGraph& operator=(const OwnedDenseGraph& o) {
        if (this != &o) {
            this->N = o.N; this->dmax = o.dmax;
            owned_neighbours = o.owned_neighbours;
            owned_deg = o.owned_deg;
            repoint();
        }
        return *this;
    }

    OwnedDenseGraph& operator=(OwnedDenseGraph&& o) noexcept {
        this->N = o.N; this->dmax = o.dmax;
        owned_neighbours = std::move(o.owned_neighbours);
        owned_deg = std::move(o.owned_deg);
        repoint();
        o.neighbours = {}; o.deg = {}; o.N = 0;
        return *this;
    }

    // --- Resize operations (require reallocation) ---

    void resize(int new_N) {
        owned_neighbours.resize(new_N * this->dmax, K(-1));
        owned_deg.resize(new_N, 0);
        this->N = K(new_N);
        repoint();
    }

    void push_back(const std::vector<K>& row) {
        assert(int(row.size()) <= this->dmax);
        owned_neighbours.resize((this->N + 1) * this->dmax, K(-1));
        owned_deg.push_back(uint8_t(row.size()));
        repoint();
        for (int i = 0; i < int(row.size()); ++i)
            this->neighbours[this->N * this->dmax + i] = row[i];
        this->N++;
    }

    void pop_back() {
        assert(this->N > 0);
        this->N--;
        owned_neighbours.resize(this->N * this->dmax);
        owned_deg.resize(this->N);
        repoint();
    }

    // --- Stride change ---

    OwnedDenseGraph restride(int new_dmax) const {
        OwnedDenseGraph g(this->N, new_dmax);
        for (K v = 0; v < this->N; ++v) {
            assert(this->deg[v] <= new_dmax);
            g.deg[v] = this->deg[v];
            for (int i = 0; i < this->deg[v]; ++i)
                g.neighbours[v * new_dmax + i] = this->neighbours[v * this->dmax + i];
        }
        return g;
    }
};

// --- ostream support ---

template<typename K>
std::ostream& operator<<(std::ostream& s, const RSRAdjacencyView<K>& g) {
    s << '[';
    for (K v = 0; v < g.N; ++v) {
        if (v) s << ',';
        s << '[';
        auto row = g[v];
        for (int i = 0; i < int(row.size()); ++i) {
            if (i) s << ',';
            s << row[i];
        }
        s << ']';
    }
    s << ']';
    return s;
}

// ---------------------------------------------------------------------------
// Free-function restride: works on any RSRAdjacencyView. O(N*dmax).
// ---------------------------------------------------------------------------
template<typename K>
OwnedDenseGraph<K> restride(const RSRAdjacencyView<K>& src, int new_dmax) {
    OwnedDenseGraph<K> g(src.N, new_dmax);
    for (K v = 0; v < src.N; ++v) {
        assert(src.deg[v] <= new_dmax);
        g.deg[v] = src.deg[v];
        for (int i = 0; i < src.deg[v]; ++i)
            g.neighbours[v * new_dmax + i] = src.neighbours[v * src.dmax + i];
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: RSR -> CSR (freeze). O(N).
// ---------------------------------------------------------------------------
template<typename K>
PlanarCSR<Owned, K> freeze(const RSRAdjacencyView<K>& g) {
    PlanarCSR<Owned, K> csr;
    csr.N = g.N;
    csr.offsets.resize(g.N + 1);
    csr.offsets[0] = 0;
    for (K v = 0; v < g.N; ++v)
        csr.offsets[v+1] = csr.offsets[v] + K(g.deg[v]);
    csr.n_arcs = csr.offsets[g.N];

    csr.values.resize(csr.n_arcs);
    for (K v = 0; v < g.N; ++v) {
        K base = csr.offsets[v];
        for (int i = 0; i < g.deg[v]; ++i)
            csr.values[base + i] = g.neighbours[v * g.dmax + i];
    }

    csr.twin = compute_twin(csr.N, csr.offsets, csr.values);
    return csr;
}

// ---------------------------------------------------------------------------
// Conversion: CSR -> RSR (thaw). O(N).
// Caller must choose dmax large enough for the graph's maximum degree.
// ---------------------------------------------------------------------------
template<typename K>
OwnedDenseGraph<K> thaw(const PlanarCSR<Owned, K>& csr, int dmax) {
    OwnedDenseGraph<K> g(csr.N, dmax);
    for (K v = 0; v < csr.N; ++v) {
        K d = csr.offsets[v+1] - csr.offsets[v];
        assert(d <= dmax);
        g.deg[v] = uint8_t(d);
        for (K i = 0; i < d; ++i)
            g.neighbours[v * dmax + i] = csr.values[csr.offsets[v] + i];
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: vector<vector<int>> -> RSR. O(N).
// ---------------------------------------------------------------------------
template<typename K = int32_t>
OwnedDenseGraph<K> dense_from_neighbours(int N_verts,
                                     const std::vector<std::vector<int>>& nbrs,
                                     int dmax) {
    OwnedDenseGraph<K> g(K(N_verts), dmax);
    for (int v = 0; v < N_verts; ++v) {
        assert(int(nbrs[v].size()) <= dmax);
        g.deg[v] = uint8_t(nbrs[v].size());
        for (int i = 0; i < int(nbrs[v].size()); ++i)
            g.neighbours[v * dmax + i] = K(nbrs[v][i]);
    }
    return g;
}

// ---------------------------------------------------------------------------
// Conversion: RSR -> vector<vector<int>>. O(N).
// ---------------------------------------------------------------------------
template<typename K>
std::vector<std::vector<int>> to_neighbours(const RSRAdjacencyView<K>& g) {
    std::vector<std::vector<int>> adj(g.N);
    for (K v = 0; v < g.N; ++v) {
        adj[v].resize(g.deg[v]);
        for (int i = 0; i < g.deg[v]; ++i)
            adj[v][i] = int(g.neighbours[v * g.dmax + i]);
    }
    return adj;
}

} // namespace Spanify
