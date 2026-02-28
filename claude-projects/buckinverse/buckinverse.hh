#pragma once

// Buckinverse: Invert buckygen extension operations.
// Given a fullerene dual graph (Triangulation), find a reduction path to a seed.

#include "fullerenes/triangulation.hh"
#include <cassert>
#include <optional>
#include <vector>
#include <string>

namespace buckinverse {

// =====================================================================
// Direction: DRight (CW walk) or DLeft (CCW walk)
// =====================================================================
enum class Dir { DRight, DLeft };

inline Dir flipDir(Dir d) {
    return d == Dir::DRight ? Dir::DLeft : Dir::DRight;
}

// =====================================================================
// Navigation primitives for fullerene dual graphs (degree 5-6 triangulations)
// These wrap Graph::next/prev with buckygen-specific straight/turn operations.
// All assume CW-ordered neighbour lists (is_oriented == true).
// =====================================================================

// Advance k positions CW from v around u's neighbour list.
inline node_t advanceCW(const Graph& g, node_t u, node_t v, int k) {
    int pos = g.arc_ix(u, v);
    assert(pos >= 0 && "v must be a neighbour of u");
    int d = g.degree(u);
    return g.neighbours[u][((pos + k) % d + d) % d];
}

// Straight-ahead through u when entering from `from`.
// At degree-6: always 3 positions CW from `from`.
// At degree-5: 3 positions for DRight, 2 positions for DLeft.
inline node_t straightAhead(const Graph& g, Dir dir, node_t u, node_t from) {
    int d = g.degree(u);
    if (d == 6) return advanceCW(g, u, from, 3);
    // degree 5
    if (dir == Dir::DRight) return advanceCW(g, u, from, 3);
    else return advanceCW(g, u, from, 2);
}

// Turn-ahead (bent path) through u when entering from `from`.
// DRight: 2 positions CW. DLeft: deg-2 positions CW.
inline node_t turnAhead(const Graph& g, Dir dir, node_t u, node_t from) {
    if (dir == Dir::DRight) return advanceCW(g, u, from, 2);
    else return advanceCW(g, u, from, g.degree(u) - 2);
}

// Side neighbour: the parallel-path vertex.
// DRight: prevCW(a, b). DLeft: nextCW(a, b).
inline node_t sideNbr(const Graph& g, Dir dir, node_t a, node_t b) {
    return (dir == Dir::DRight) ? g.prev(a, b) : g.next(a, b);
}

// =====================================================================
// Degree-5 vertex enumeration
// =====================================================================

inline std::vector<node_t> deg5vertices(const Graph& g) {
    std::vector<node_t> result;
    for (int v = 0; v < g.N; ++v)
        if (g.degree(v) == 5) result.push_back(v);
    return result;
}

// =====================================================================
// Expansion types
// =====================================================================

struct ExpKind {
    enum Type { L_type, B_type, F_type };
    Type type;
    int i, j;  // L_i: i >= 0, j unused. B_{i,j}: i,j >= 0. F: both unused.

    int reductionLength() const {
        switch (type) {
            case L_type: return i + 1;
            case B_type: return i + j + 2;
            case F_type: return 0;
        }
        return 0;
    }

    int newVertices() const {
        switch (type) {
            case L_type: return i + 2;
            case B_type: return i + j + 3;
            case F_type: return 5;
        }
        return 0;
    }

    std::string toString() const {
        switch (type) {
            case L_type: return "L" + std::to_string(i);
            case B_type: return "B(" + std::to_string(i) + "," + std::to_string(j) + ")";
            case F_type: return "F";
        }
        return "?";
    }

    bool operator==(const ExpKind& o) const {
        return type == o.type && i == o.i && j == o.j;
    }
};

inline ExpKind Lk(int i)          { return {ExpKind::L_type, i, 0}; }
inline ExpKind Bk(int i, int j)   { return {ExpKind::B_type, i, j}; }
inline ExpKind Fk()               { return {ExpKind::F_type, 0, 0}; }

// =====================================================================
// Expansion / Reduction triples
// =====================================================================

struct Expansion {
    ExpKind kind;
    node_t u, v;   // Starting directed edge (u is degree-5 for L/B)
    Dir dir;

    std::string toString() const {
        return kind.toString() + " (" + std::to_string(u) + "->" +
               std::to_string(v) + ") " +
               (dir == Dir::DRight ? "R" : "L");
    }
};

// A Reduction is the same triple as an Expansion but used in the inverse direction.
using Reduction = Expansion;

// =====================================================================
// Path info: main path and parallel path computed from the graph
// =====================================================================

struct PathInfo {
    std::vector<node_t> path;      // Main path vertices
    std::vector<node_t> parallel;  // Parallel path vertices
    bool valid = false;
};

// =====================================================================
// Seed graphs
// =====================================================================

enum class SeedType { C20, C28, C30, NotASeed };

// Construct seed graphs as Triangulation objects
Triangulation makeSeedC20();
Triangulation makeSeedC28();
Triangulation makeSeedC30();

// Identify which seed a small triangulation is (by vertex count, since
// the three seeds have unique vertex counts: 12, 16, 17).
// Returns NotASeed if the graph doesn't match any seed size.
SeedType identifySeed(const Triangulation& t);

inline std::string seedName(SeedType s) {
    switch (s) {
        case SeedType::C20: return "C20";
        case SeedType::C28: return "C28";
        case SeedType::C30: return "C30";
        case SeedType::NotASeed: return "not a seed";
    }
    return "?";
}

// =====================================================================
// Path computation
// =====================================================================

// Compute a straight path of numEntries vertices starting from directed
// edge (u, v) in direction dir. Returns valid=false if path self-intersects.
PathInfo computeStraightPath(const Graph& g, node_t u, node_t v,
                             Dir dir, int numEntries);

// Compute a B_{0,0} bent path starting from (u, v) in direction dir.
PathInfo computeBentZeroPath(const Graph& g, node_t u, node_t v, Dir dir);

// Compute a B_{i,j} bent path (i+j > 0) starting from (u, v) in direction dir.
PathInfo computeBentPath(const Graph& g, node_t u, node_t v,
                         Dir dir, int bi, int bj);

// =====================================================================
// Reduction surgery (Phase 3)
// =====================================================================

// Apply a reduction to produce a smaller graph.
// Returns the reduced graph. The input must be a valid oriented triangulation
// and the reduction must be a valid reduction site (as returned by allReductions).
// NOTE: This uses the expansion-site formulation (Phase 2). For arbitrary graph
// inversion, use invertReduction() instead.
Graph applyReduction(const Graph& g, const Reduction& red);

// =====================================================================
// Inversion: reduce arbitrary graphs (Phase 3)
// =====================================================================

// A validated inversion site with all vertices identified.
struct InvSite {
    ExpKind kind;
    Dir dir;
    std::vector<node_t> strip;   // strip vertices to remove
    std::vector<node_t> path;    // path vertices (remain, get reconnected)
    std::vector<node_t> tp;      // true parallel vertices (remain, get reconnected)
};

// =====================================================================
// Extension path: portable representation of seed → fullerene route
// =====================================================================

// One expansion step (insert strip vertices and rewire path/tp).
// Vertex IDs are in the final (full-size) graph's numbering.
struct ExtensionStep {
    ExpKind kind;
    Dir dir;
    std::vector<node_t> strip;   // vertices to insert
    std::vector<node_t> path;    // existing vertices (reconnected)
    std::vector<node_t> tp;      // existing vertices (reconnected)
};

// Full extension path from seed to target fullerene.
// Steps are in expansion order (seed → full graph).
struct ExtensionPath {
    SeedType seed;
    int full_N;                           // vertex count of full graph

    // Seed state: full vertex data for alive vertices after full reduction.
    // Includes inactive nbr[] positions, which hold "shadow" strip vertex IDs
    // from prior splice/snip operations, needed by unsnip during expansion.
    struct SeedVertex {
        node_t id;
        node_t nbr[6];
        uint8_t active;
    };
    std::vector<SeedVertex> seed_state;

    std::vector<ExtensionStep> steps;     // in expansion order

    std::string toString() const {
        std::string s = seedName(seed) + " [" + std::to_string(full_N) + "v, "
                        + std::to_string(seed_state.size()) + " seed verts]:";
        for (auto& st : steps) s += " " + st.kind.toString();
        return s;
    }
};

// Find all valid L0 inversion sites on the graph.
// An L0 site is a pair of adjacent degree-5 vertices (the strip) with the
// correct expansion topology around them.
std::vector<InvSite> findL0InvSites(const Graph& g);

// Find all valid F (nanotube ring) inversion sites on the graph.
// An F site is a 5-cycle of degree-6 strip vertices sitting between a
// parallel ring 5-cycle and a set of 5 outer (cap) vertices.
std::vector<InvSite> findFRingInvSites(const Graph& g);

// Find all valid inversion sites (L0, L_i, B(i,j), F) up to given max reduction length.
std::vector<InvSite> allInvSites(const Graph& g, int maxRedLen = 5);

// Apply an inversion site to produce a smaller graph.
// Returns Graph() (N==0) on failure.
Graph invertReduction(const Graph& g, const InvSite& site);

// =====================================================================
// Reduction enumeration (Phase 2)
// =====================================================================

// Enumerate all valid expansion sites (for canonical algorithm).
// Default maxRedLen=5 matches the Haskell allReductions.
std::vector<Reduction> allReductions(const Graph& g, int maxRedLen = 5);

// Follow straight-ahead from u through v until a degree-5 vertex is found.
// Returns {endpoint, prevVertex, distance} or empty optional if not found
// within maxDist steps or if a cycle is detected.
struct StraightEndpoint {
    node_t endpoint;
    node_t prev;
    int distance;
};
std::optional<StraightEndpoint> followStraightToFive(
    const Graph& g, node_t u, node_t v, int maxDist);

// =====================================================================
// ReducibleDual: O(N) reduction to seed via bitmask-over-fixed-array.
// Follows the RemainingGraph pattern (triangulation.cc) extended with
// in-place edge mutation (splice) for reduction reconnection.
// =====================================================================

struct ReducibleDual {
    static constexpr int D_MAX = 6;  // Max degree in fullerene dual

    struct Vertex {
        node_t nbr[D_MAX];  // CW-ordered neighbours (fixed positions)
        uint8_t active;      // Bitmask: bit i set iff nbr[i] is present
    };

    // A recorded reduction step: the site plus saved strip vertex states.
    // Saved states are needed because prior steps may have modified strip
    // vertices via splice (e.g., path[k] in step i becomes strip[j] later).
    struct ReductionStep {
        InvSite site;
        std::vector<std::pair<node_t, Vertex>> saved;
    };

    std::vector<Vertex> V;   // Indexed by original vertex ID (never resized)
    std::set<node_t> deg5;   // The 12 degree-5 vertices (maintained)
    int n_alive;

    // Construction: O(N)
    explicit ReducibleDual(const Graph& g);

    // Construct with given capacity, all vertices dead.
    // Used by graphFromExtensionPath to build from seed + extension path.
    explicit ReducibleDual(int capacity);

    // Queries
    int  degree(node_t u) const { return __builtin_popcount(V[u].active); }
    bool alive(node_t u)  const { return V[u].active != 0; }
    int  N()               const { return n_alive; }
    int  arc_ix(node_t u, node_t v) const;
    node_t next(node_t u, node_t v) const;
    node_t prev(node_t u, node_t v) const;

    // Navigation (buckygen-style)
    node_t advanceCW(node_t u, node_t v, int k) const;
    node_t straightAhead(Dir dir, node_t u, node_t from) const;

    // Mutations: O(1) each
    void splice(node_t u, node_t old_v, node_t new_v);  // Replace neighbour
    void snip(node_t u, node_t v);                       // Remove neighbour
    void unsnip(node_t u, node_t v);                     // Reverse of snip
    void kill(node_t v);                                  // Delete vertex
    void unkill(node_t v, const Vertex& saved);          // Reverse of kill
    void insertAfter(node_t u, node_t v, node_t after);  // Insert v after 'after' in u's CW (5→6)

    // Reduction: O(1) per call
    std::optional<InvSite> findSite(int maxRedLen = 5) const;
    void reduce(const InvSite& site);

    // Expansion: O(1) per call — inverse of reduce (round-trip, uses saved state)
    void expand(const ReductionStep& step);

    // Standalone expansion: O(1) per call — creates strip vertices from CW formulas.
    // Does NOT require saved state; works from seed + extension path alone.
    void expand(const ExtensionStep& step);

    // Full reduction loop: O(N) total
    SeedType reduceToSeed(int maxRedLen = 5);
    SeedType reduceToSeed(std::vector<ReductionStep>& path, int maxRedLen = 5);

    // Reduce and return extension path (seed → full graph direction)
    ExtensionPath reduceToExtensionPath(int maxRedLen = 5);

    // Extract compacted Graph: O(N)
    Graph toGraph() const;

private:
    // Site finders (return first valid site or nullopt)
    std::optional<InvSite> findL0Site() const;
    std::optional<InvSite> findLiSite(int maxRedLen) const;
    std::optional<InvSite> findB00Site() const;
    std::optional<InvSite> findFRingSite() const;

    // Reconnection (reduction direction)
    void reconnectStraight(const std::vector<node_t>& path,
                           const std::vector<node_t>& strip,
                           const std::vector<node_t>& tp);
    void reconnectBent(const std::vector<node_t>& path,
                       const std::vector<node_t>& strip,
                       const std::vector<node_t>& tp,
                       int bentPos, int bentLen);
    void reconnectRing(const std::vector<node_t>& ring,
                       const std::vector<node_t>& strip,
                       const std::vector<node_t>& outer);

    // Reconnection (expansion direction — reverse of reduction)
    void expandStraight(const std::vector<node_t>& path,
                        const std::vector<node_t>& strip,
                        const std::vector<node_t>& tp);
    void expandBent(const std::vector<node_t>& path,
                    const std::vector<node_t>& strip,
                    const std::vector<node_t>& tp,
                    int bentPos, int bentLen);
    void expandRing(const std::vector<node_t>& ring,
                    const std::vector<node_t>& strip,
                    const std::vector<node_t>& outer);
};

// Reconstruct a fullerene dual graph from a seed + extension path.
// This is standalone: does not require the original graph or any saved state.
Graph graphFromExtensionPath(const ExtensionPath& ep);

} // namespace buckinverse
