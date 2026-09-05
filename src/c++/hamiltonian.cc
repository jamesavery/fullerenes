#include "fullerenes/graphview.hh"

// GraphView::hamiltonian_cycle_count -- the number of undirected Hamiltonian
// cycles (contract on the declaration in graphview.hh).
//
// Darko Babic's backtracking, as implemented in the legacy Fortran
// HamiltonCyc, generalised from cubic graphs to any degree.
// @ref hamilton.f:1-183
//
// Counting each cycle once.  Every Hamiltonian cycle meets the start vertex
// s in two of its edges, {s,n_i} and {s,n_j} with i < j in s's neighbour
// order.  The search takes the arc s->n_i as its first step and lets the
// cycle close only through a LATER neighbour n_j: phase i then finds exactly
// the cycles whose lower s-edge index is i, so the phases partition the
// cycles and each is counted once, traversed in the direction that leaves s
// along its lower-indexed edge.  (In phase i the edges {s,n_k}, k < i, are
// dead and excluded; that exclusion is a prune, not what makes the phases
// disjoint.)
//
// Three prunes, each the negation of a necessary condition, so none can lose
// a cycle; they are the named predicates below.  The budget (Babic's "pass"
// rule) is what makes C110 feasible; the plain distance-pruned search is
// not.  For the record, the budget is a sound RELAXATION of the Fortran's
// two "pass" invariants, not an equivalent: a passed closing neighbour of s
// still counts its reserved edge to s here, so the branch the Fortran prunes
// outright is entered and dies one level down.  In the other direction this
// has the BFS bound, which the Fortran lacks.  Neither difference can change
// the count.
//
// The search state is mutated in place -- a backtracking undo stack, which
// style.md leaves to judgment stated here: every mutation is undone either by
// the scope that made it (Exclusion) or on the line after the recursive call
// (place_and_count), so do and undo are always adjacent.  The enumerating
// variant (the remaining hamilton.f functionality, cycle listing and IUPAC
// naming) forks at the depth == g.N leaf on the same state and prunes, and
// belongs beside this.

namespace {

// One count; constructed, run once, discarded.
//
// @anchor hamiltonian-cycle-search
// @inv closing: closing_ends_left == count_if(g.nbrs(s), [&](node_t v){ return closing_end[v] && !on_path[v]; })
// @inv budget:  all_of(indices(g.N), [&](node_t w){ return in_range(usable_degree[w], 0, g.degree(w)); })
//
// The fact the budget prune rests on, not checkable without the path: for an
// unvisited w, usable_degree[w] is an UPPER bound on the number of edges at w
// that any cycle extending the current path can still use, because it is
// decremented only for edges no such cycle can use -- the edges a step leaves
// behind (their endpoint's two cycle edges are already fixed) and the earlier
// phases' s-edges.  A bound below two therefore proves the branch dead.
struct HamiltonianCycleSearch {
  const GraphView& g;
  const node_t s = 0;            // every Hamiltonian cycle passes through it
  vector<char> on_path;
  vector<int>  usable_degree;    // edges at an unvisited vertex not yet excluded
  vector<char> closing_end;      // neighbour of s the cycle may close through
  vector<int>  dist;             // BFS distance from s on the full graph
  int closing_ends_left = 0;     // closing ends not yet on the path

  HamiltonianCycleSearch(const GraphView& g)
    : g(g), on_path(g.N, 0), usable_degree(g.N), closing_end(g.N, 0), dist(g.N) {
    for (node_t v = 0; v < g.N; v++) usable_degree[v] = g.degree(v);
    g.single_source_shortest_paths(s, dist.data());
  }

  // The edges from the vertex a step leaves to the neighbours `ws` it skips
  // (all but `taken`), excluded for the lifetime of the branch: each unvisited
  // skipped neighbour loses one usable edge.  `viable` says whether each keeps
  // the two an unvisited vertex needs (one in, one out; for a closing
  // neighbour of s the reserved edge to s is the out).  A scope, so exclusion
  // and restoration cannot come apart; on_path is the same at both ends
  // because the only vertex placed inside the scope is `taken`.
  struct Exclusion {
    HamiltonianCycleSearch& S;
    span<const node_t> ws;
    node_t taken;
    bool viable = true;
    Exclusion(HamiltonianCycleSearch& S, span<const node_t> ws, node_t taken) : S(S), ws(ws), taken(taken) {
      for (node_t w : ws) if (w != taken && !S.on_path[w] && --S.usable_degree[w] < 2) viable = false;
    }
    ~Exclusion() { for (node_t w : ws) if (w != taken && !S.on_path[w]) S.usable_degree[w]++; }
    Exclusion(const Exclusion&) = delete;
    Exclusion& operator=(const Exclusion&) = delete;
  };

  // A vertex placed as the depth-th of the path closes the cycle over the
  // g.N - depth + 1 edges that remain.
  bool reaches_start_in_time(node_t v, int depth) const { return dist[v] <= g.N - depth + 1; }

  // The cycle must end at a closing end, so the last one off the path may not
  // be consumed before the last position.
  bool would_strand_closing_end(node_t v, int depth) const {
    return closing_end[v] && depth < g.N && closing_ends_left == 1;
  }

  // The number of Hamiltonian cycles completing the path whose `depth`
  // vertices end at `last`.
  int64_t count_completions(node_t last, int depth) {
    if (depth == g.N) return closing_end[last];
    auto nb = g.nbrs(last);
    int64_t n = 0;
    for (node_t next : nb) {
      if (on_path[next] || !reaches_start_in_time(next, depth + 1)) continue;
      Exclusion excluded(*this, nb, next);
      if (!excluded.viable || would_strand_closing_end(next, depth + 1)) continue;
      n += place_and_count(next, depth + 1);
    }
    return n;
  }

  // Put v on the path as its depth-th vertex, count, take it off again.
  int64_t place_and_count(node_t v, int depth) {
    on_path[v] = 1;  closing_ends_left -= closing_end[v];
    const int64_t n = count_completions(v, depth);
    closing_ends_left += closing_end[v];  on_path[v] = 0;
    return n;
  }

  // Phase i: the cycles whose first step from s is the arc s->ns[i] and which
  // close through a later neighbour of s.
  int64_t count_cycles_closing_after(span<const node_t> ns, int i) {
    for (size_t j = 0; j < ns.size(); j++) closing_end[ns[j]] = (j > size_t(i));
    closing_ends_left = int(ns.size()) - 1 - i;
    const node_t first = ns[i];
    Exclusion excluded(*this, ns.first(i), first);   // earlier phases' s-edges, dead here
    if (!excluded.viable || !reaches_start_in_time(first, 2)) return 0;
    return place_and_count(first, 2);
  }

  int64_t count_cycles() {
    auto ns = g.nbrs(s);
    on_path[s] = 1;
    int64_t n = 0;
    for (size_t i = 0; i + 1 < ns.size(); i++) n += count_cycles_closing_after(ns, int(i));
    return n;
  }
};

}  // namespace

int64_t GraphView::hamiltonian_cycle_count() const
{
  if (N < 3) return 0;   // before the search state (and its BFS) is built
  return HamiltonianCycleSearch(*this).count_cycles();
}
