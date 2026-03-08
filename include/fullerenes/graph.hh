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

struct Graph : Spanify::DenseGraph<node_t> {
  using base_t = Spanify::DenseGraph<node_t>;

  // Optional owned storage (empty when Graph is a view of external memory).
  // Follows the SpanVector pattern: base_t pointers (values, deg) always
  // point to either owned storage or external memory.
  std::vector<node_t> owned_values;
  std::vector<uint8_t> owned_deg;

  void repoint() {
    values = owned_values.data();
    deg = owned_deg.data();
  }

  // Bring base push_back(K u, K v) into scope (otherwise hidden by push_back(vector))
  using base_t::push_back;

  // --- Constructors ---

  Graph() = default;

  // Allocating: N vertices, dmax stride, all degrees 0.
  Graph(size_t N_, int dmax_ = GRAPH_DMAX)
      : owned_values(N_ * dmax_, node_t(-1)), owned_deg(N_, 0) {
    N = node_t(N_); dmax = dmax_; repoint();
  }

  // Copy from adjacency view (copies data, owns it).
  Graph(const base_t& adj)
      : owned_values(adj.values, adj.values + adj.N * adj.dmax),
        owned_deg(adj.deg, adj.deg + adj.N) {
    N = adj.N; dmax = adj.dmax; repoint();
  }

  // Move from OwnedDenseGraph (takes ownership, zero-copy).
  Graph(Spanify::OwnedDenseGraph<node_t>&& adj)
      : owned_values(std::move(adj.owned_values)),
        owned_deg(std::move(adj.owned_deg)) {
    N = adj.N; dmax = adj.dmax; repoint();
    adj.values = nullptr; adj.deg = nullptr; adj.N = 0;
  }

  // Copy from OwnedDenseGraph.
  Graph(const Spanify::OwnedDenseGraph<node_t>& adj)
      : owned_values(adj.owned_values),
        owned_deg(adj.owned_deg) {
    N = adj.N; dmax = adj.dmax; repoint();
  }

  // Brace-init constructor: Graph{{1,2,3},{4,5,6},...}
  Graph(std::initializer_list<std::initializer_list<node_t>> adj) {
    N = node_t(adj.size());
    dmax = 0;
    for (auto& row : adj) dmax = std::max(dmax, int(row.size()));
    owned_values.resize(N * dmax, node_t(-1));
    owned_deg.resize(N, 0);
    repoint();
    int v = 0;
    for (auto& row : adj) {
      deg[v] = uint8_t(row.size());
      int i = 0;
      for (node_t x : row) values[v * dmax + i++] = x;
      v++;
    }
  }

  // --- Assignment from OwnedDenseGraph ---

  Graph& operator=(Spanify::OwnedDenseGraph<node_t>&& adj) {
    N = adj.N; dmax = adj.dmax;
    owned_values = std::move(adj.owned_values);
    owned_deg = std::move(adj.owned_deg);
    repoint();
    adj.values = nullptr; adj.deg = nullptr; adj.N = 0;
    return *this;
  }

  Graph& operator=(const Spanify::OwnedDenseGraph<node_t>& adj) {
    N = adj.N; dmax = adj.dmax;
    owned_values = adj.owned_values;
    owned_deg = adj.owned_deg;
    repoint();
    return *this;
  }

  // --- Rule of 5 ---

  Graph(const Graph& o)
      : base_t(o), owned_values(o.owned_values), owned_deg(o.owned_deg) {
    if (!owned_values.empty()) repoint();
  }

  Graph(Graph&& o) noexcept
      : base_t(o), owned_values(std::move(o.owned_values)),
        owned_deg(std::move(o.owned_deg)) {
    if (!owned_values.empty()) repoint();
    o.values = nullptr; o.deg = nullptr; o.N = 0;
  }

  Graph& operator=(const Graph& o) {
    if (this != &o) {
      N = o.N; dmax = o.dmax;
      owned_values = o.owned_values;
      owned_deg = o.owned_deg;
      if (!owned_values.empty()) repoint();
      else { values = o.values; deg = o.deg; }
    }
    return *this;
  }

  Graph& operator=(Graph&& o) noexcept {
    N = o.N; dmax = o.dmax;
    owned_values = std::move(o.owned_values);
    owned_deg = std::move(o.owned_deg);
    if (!owned_values.empty()) repoint();
    else { values = o.values; deg = o.deg; }
    o.values = nullptr; o.deg = nullptr; o.N = 0;
    return *this;
  }

  // --- Vertex-level reallocation (requires ownership) ---

  void resize(int new_N) {
    owned_values.resize(new_N * dmax, node_t(-1));
    owned_deg.resize(new_N, 0);
    N = node_t(new_N);
    repoint();
  }

  void push_back(const std::vector<node_t>& row) {
    assert(int(row.size()) <= dmax);
    owned_values.resize((N + 1) * dmax, node_t(-1));
    owned_deg.push_back(uint8_t(row.size()));
    repoint();
    for (int i = 0; i < int(row.size()); ++i)
      values[N * dmax + i] = row[i];
    N++;
  }

  void pop_back() {
    assert(N > 0);
    N--;
    owned_values.resize(N * dmax);
    owned_deg.resize(N);
    repoint();
  }

  // Stride change: returns new owned adjacency with different dmax.
  Spanify::OwnedDenseGraph<node_t> restride(int new_dmax) const {
    return Spanify::restride(*this, new_dmax);
  }

  // --- Graph algorithms ---

  bool insert_edge(const arc_t& e, const node_t suc_uv=-1, const node_t suc_vu=-1);
  bool remove_edge(const edge_t& e);
  bool edge_exists(const edge_t& e) const;
  void remove_isolated_vertices();
  void remove_vertices(set<int> &sv);
  void flip_all_orientations();

  int  arc_ix(node_t u, node_t v) const;
  node_t next(node_t u, node_t v) const;
  node_t prev(node_t u, node_t v) const;
  node_t next_on_face(node_t u, node_t v) const;
  node_t prev_on_face(node_t u, node_t v) const;

  bool is_consistently_oriented() const;
  bool adjacency_is_symmetric() const;
  bool has_separating_triangles() const;

  bool is_connected(const set<node_t>& subgraph = set<node_t>()) const;
  void      single_source_shortest_paths(node_t source, int *distances, size_t max_depth = INT_MAX) const;
  matrix<int> all_pairs_shortest_paths(const vector<node_t>& V,
				       const unsigned int max_depth = INT_MAX) const;
  matrix<int> all_pairs_shortest_paths(const unsigned int max_depth = INT_MAX) const;
  vector<vector<node_t> > connected_components() const;

  vector<node_t> shortest_cycle(node_t s, const int max_depth) const;
  vector<node_t> shortest_cycle(const vector<node_t> &prefix, const int max_depth) const;
  vector<int> multiple_source_shortest_paths(const vector<node_t>& sources, const unsigned int max_depth=INT_MAX) const;

  int hamiltonian_count() const;
  int hamiltonian_count(node_t current_node, vector<bool>& used_edges, vector<bool>& used_nodes, vector<node_t>& path, const vector<int>& distances) const;

  coord2d centre2d(const vector<coord2d>& layout) const;
  coord3d centre3d(std::span<const coord3d> layout) const;

  int max_degree() const;

  vector<edge_t>  undirected_edges() const;
  vector<arc_t> directed_edges()   const;

  size_t count_edges() const;

  friend ostream& operator<<(ostream& s, const Graph& g);
};
