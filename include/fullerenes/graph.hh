#pragma once

#include <stdio.h>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <assert.h>
#include <span>

using namespace std;
#include "geometry.hh"
#include "auxiliary.hh"
#include "matrix.hh"
#include "graphview.hh"
#include "permutation.hh"

struct Graph : GraphView {
  using base_t = Spanify::RSRAdjacencyView<node_t>;

  // Optional owned storage (empty when Graph is a view of external memory).
  // base_t spans (neighbours, deg) always reference either owned storage
  // or external memory.
  std::vector<node_t> owned_neighbours;
  std::vector<uint8_t> owned_deg;

  void repoint() {
    neighbours = std::span<node_t>(owned_neighbours.data(), owned_neighbours.size());
    deg = std::span<uint8_t>(owned_deg.data(), owned_deg.size());
  }

  bool owns_memory() const { return !owned_neighbours.empty(); }

  // Bring base push_back(K u, K v) into scope (otherwise hidden by push_back(vector))
  using base_t::push_back;

  // --- Constructors ---

  Graph() = default;

  // View: wrap existing adjacency memory (caller manages lifetime).
  Graph(node_t N_, int dmax_, std::span<node_t> neighbours_, std::span<uint8_t> deg_) {
    N = N_; dmax = dmax_; neighbours = neighbours_; deg = deg_;
  }

  // Allocating: N vertices, dmax stride, all degrees 0.
  Graph(size_t N_, int dmax_ = GRAPH_DMAX)
      : owned_neighbours(N_ * dmax_, node_t(-1)), owned_deg(N_, 0) {
    N = node_t(N_); dmax = dmax_; repoint();
  }

  // Allocating: N vertices, each initialized to initial_row.
  Graph(size_t N_, const std::vector<node_t>& initial_row)
      : owned_neighbours(N_ * initial_row.size()), owned_deg(N_, uint8_t(initial_row.size())) {
    N = node_t(N_); dmax = int(initial_row.size()); repoint();
    for (size_t v = 0; v < N_; ++v)
      std::copy(initial_row.begin(), initial_row.end(), owned_neighbours.data() + v * dmax);
  }

  // Copy from adjacency view (copies data, owns it).
  Graph(const base_t& adj)
      : owned_neighbours(adj.neighbours.begin(), adj.neighbours.end()),
        owned_deg(adj.deg.begin(), adj.deg.end()) {
    N = adj.N; dmax = adj.dmax; repoint();
  }

  // Move from OwnedDenseGraph (takes ownership, zero-copy).
  Graph(Spanify::OwnedDenseGraph<node_t>&& adj)
      : owned_neighbours(std::move(adj.owned_neighbours)),
        owned_deg(std::move(adj.owned_deg)) {
    N = adj.N; dmax = adj.dmax; repoint();
    adj.neighbours = {}; adj.deg = {}; adj.N = 0;
  }

  // Copy from OwnedDenseGraph.
  Graph(const Spanify::OwnedDenseGraph<node_t>& adj)
      : owned_neighbours(adj.owned_neighbours),
        owned_deg(adj.owned_deg) {
    N = adj.N; dmax = adj.dmax; repoint();
  }

  // Brace-init constructor: Graph{{1,2,3},{4,5,6},...}
  Graph(std::initializer_list<std::initializer_list<node_t>> adj) {
    N = node_t(adj.size());
    dmax = 0;
    for (auto& row : adj) dmax = std::max(dmax, int(row.size()));
    owned_neighbours.resize(N * dmax, node_t(-1));
    owned_deg.resize(N, 0);
    repoint();
    int v = 0;
    for (auto& row : adj) {
      deg[v] = uint8_t(row.size());
      int i = 0;
      for (node_t x : row) neighbours[v * dmax + i++] = x;
      v++;
    }
  }

  // --- Assignment from any GraphView (deep copy adjacency) ---
  Graph& operator=(const GraphView& src) {
    if (src.neighbours.data() == owned_neighbours.data()) return *this;
    N = src.N; dmax = src.dmax;
    owned_neighbours.assign(src.neighbours.begin(), src.neighbours.end());
    owned_deg.assign(src.deg.begin(), src.deg.end());
    repoint();
    return *this;
  }

  // --- Assignment from OwnedDenseGraph ---

  Graph& operator=(Spanify::OwnedDenseGraph<node_t>&& adj) {
    N = adj.N; dmax = adj.dmax;
    owned_neighbours = std::move(adj.owned_neighbours);
    owned_deg = std::move(adj.owned_deg);
    repoint();
    adj.neighbours = {}; adj.deg = {}; adj.N = 0;
    return *this;
  }

  Graph& operator=(const Spanify::OwnedDenseGraph<node_t>& adj) {
    N = adj.N; dmax = adj.dmax;
    owned_neighbours = adj.owned_neighbours;
    owned_deg = adj.owned_deg;
    repoint();
    return *this;
  }

  // --- Rule of 5 ---

  Graph(const Graph& o)
      : GraphView(o), owned_neighbours(o.owned_neighbours), owned_deg(o.owned_deg) {
    if (!owned_neighbours.empty()) repoint();
  }

  Graph(Graph&& o) noexcept
      : GraphView(o), owned_neighbours(std::move(o.owned_neighbours)),
        owned_deg(std::move(o.owned_deg)) {
    if (!owned_neighbours.empty()) repoint();
    o.neighbours = {}; o.deg = {}; o.N = 0;
  }

  Graph& operator=(const Graph& o) {
    if (this != &o) {
      N = o.N; dmax = o.dmax;
      owned_neighbours = o.owned_neighbours;
      owned_deg = o.owned_deg;
      if (!owned_neighbours.empty()) repoint();
      else { neighbours = o.neighbours; deg = o.deg; }
    }
    return *this;
  }

  Graph& operator=(Graph&& o) noexcept {
    N = o.N; dmax = o.dmax;
    owned_neighbours = std::move(o.owned_neighbours);
    owned_deg = std::move(o.owned_deg);
    if (!owned_neighbours.empty()) repoint();
    else { neighbours = o.neighbours; deg = o.deg; }
    o.neighbours = {}; o.deg = {}; o.N = 0;
    return *this;
  }

  // --- Vertex-level reallocation (requires ownership) ---

  void resize(int new_N) {
    owned_neighbours.resize(new_N * dmax, node_t(-1));
    owned_deg.resize(new_N, 0);
    N = node_t(new_N);
    repoint();
  }

  void push_back(const std::vector<node_t>& row) {
    assert(int(row.size()) <= dmax);
    owned_neighbours.resize((N + 1) * dmax, node_t(-1));
    owned_deg.push_back(uint8_t(row.size()));
    repoint();
    for (int i = 0; i < int(row.size()); ++i)
      neighbours[N * dmax + i] = row[i];
    N++;
  }

  void pop_back() {
    assert(N > 0);
    N--;
    owned_neighbours.resize(N * dmax);
    owned_deg.resize(N);
    repoint();
  }

  // Stride change: returns new owned adjacency with different dmax.
  Spanify::OwnedDenseGraph<node_t> restride(int new_dmax) const {
    return Spanify::restride(*this, new_dmax);
  }

  // --- Methods requiring owned storage ---
  void remove_isolated_vertices();
  void remove_vertices(set<int> &sv);

  // Relabel vertices according to pi (pi[u_old] = u_new): vertex u_old's
  // adjacency row moves to slot u_new, and every arc target t is
  // relabelled to pi[t].  Requires owned memory.  (argsort is inherited
  // from GraphView.)
  void apply_permutation(const Permutation& pi);

  friend ostream& operator<<(ostream& s, const Graph& g);
};
