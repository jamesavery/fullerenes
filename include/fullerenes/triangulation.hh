#pragma once

#include <functional>

#include "fullerenes/matrix.hh"
#include "fullerenes/spiral.hh"
#include "fullerenes/planargraph.hh"

class Triangulation : public Owned<TriangulationView> {
  using base_t = Owned<TriangulationView>;
public:
  typedef function<bool(Triangulation)> predicate_t;

  Triangulation() = default;
  explicit Triangulation(int N) : base_t(N) {}
  Triangulation(const GraphView& g) : base_t(g) {}
  Triangulation(const neighbours_t& neighbours) : base_t(Graph(neighbours)) {}

  Triangulation(const vector<int>& spiral_string, const jumplist_t& jumps = jumplist_t(), const bool best_effort=false);
  Triangulation(const spiral_nomenclature &fsn): Triangulation(fsn.spiral.spiral_code, fsn.spiral.jumps, true){}
};


// A THIN allocation wrapper over Owned<FullereneDualView> (which stays
// ignorant of pentagons): its whole job beyond the base is to own the
// pentagon list's backing -- constant-size and address-stable for the
// object's lifetime -- and keep the view's span pointing at it across
// construction, copy and move.  All pentagon logic, initialisation
// included, is FullereneDualView's; the converting constructor merely
// invokes the view's establishing call at its boundary.
class FullereneDual : public Owned<FullereneDualView> {
  using base_t = Owned<FullereneDualView>;

  void repoint_pentagons() { this->pentagons = std::span<node_t>(owned_pentagons); }
public:
  pentagon_storage_t<FullereneDualView> owned_pentagons{};

  FullereneDual() { repoint_pentagons(); }
  // Allocate an N-vertex dual to be filled row by row; the pentagon list is
  // stale until the producer derives it (construction phase).
  explicit FullereneDual(int N) : base_t(N) { repoint_pentagons(); }
  // Deep copy + establish: the view derives its list from the copied graph
  // (throws pentagon_error when g is not a fullerene dual).  This is the one
  // conversion path, so every construction from a foreign graph establishes.
  FullereneDual(const GraphView& g) : base_t(g) {
    repoint_pentagons();
    if (this->N > 0) this->derive_pentagons();
  }

  FullereneDual(const FullereneDual& o)
    : base_t(o), owned_pentagons(o.owned_pentagons) { repoint_pentagons(); }
  FullereneDual(FullereneDual&& o) noexcept
    : base_t(std::move(o)), owned_pentagons(o.owned_pentagons) {
    repoint_pentagons();
    o.pentagons = {};
  }
  FullereneDual& operator=(const FullereneDual& o) {
    base_t::operator=(o);
    owned_pentagons = o.owned_pentagons;
    repoint_pentagons();
    return *this;
  }
  FullereneDual& operator=(FullereneDual&& o) noexcept {
    if (this != &o) {
      base_t::operator=(std::move(o));
      owned_pentagons = o.owned_pentagons;
      repoint_pentagons();
      o.pentagons = {};
    }
    return *this;
  }

  FullereneDual(const int N, const general_spiral& rspi) : FullereneDual(N,rspi.spiral_code,rspi.jumps) {}
  FullereneDual(const int N, const vector<int>& rspi, const jumplist_t& jumps = jumplist_t()) {
    vector<int> spiral(N/2+2,6);
    for(int i: rspi) spiral[i] = 5;
    *this = Triangulation(spiral,jumps);
  }

  static vector<general_spiral> isomer_search(const Triangulation::predicate_t& predicate, size_t N, size_t print_step=0,
                                              bool IPR=false, bool only_nontrivial_symmetry=false, size_t N_chunks=1, size_t chunk_index=0);
};


class CubicPair {
  Triangulation T;
  PlanarGraph   G;
  IDCounter<tri_t> triangle_id;
  vector<vector<arc_t>> CtoD, DtoC;

  int face_start(const face_t &f){
    node_t i_m = 0;
    for(int i=0, m=INT_MAX; i<int(f.size()); i++) if(f[i] < m){ i_m = i; m = f[i]; }
    return i_m;
  }

  CubicPair(const Triangulation &T) : G(T.dual_graph()), CtoD(G.N,vector<arc_t>(3)), DtoC(T.N)
  {
    for(const auto &t: T.triangles()) triangle_id.insert(t.sorted());

    for(node_t u=0;u<T.N;u++){
      auto nu = T.nbrs(u);
      DtoC[u].resize(nu.size());

      for(size_t i=0;i<nu.size();i++){
        node_t v = nu[i];
        node_t s = nu[(i+1)%nu.size()];
        node_t t = nu[(i+nu.size()-1)%nu.size()];

        tri_t t1 = {u,v,s}, t2 = {v,u,t};

        node_t U   = triangle_id(t1.sorted()), V = triangle_id(t2.sorted());
        node_t i_V = G.arc_ix(U,V), i_U = G.arc_ix(V,U);
        node_t i_v = T.arc_ix(u,v), i_u = T.arc_ix(v,u);

        CtoD[U][i_V] = {u,i_v};
        CtoD[V][i_U] = {v,i_u};
        DtoC[u][i_v] = {U,i_V};
        DtoC[v][i_u] = {V,i_U};
      }
    }
  }
};
