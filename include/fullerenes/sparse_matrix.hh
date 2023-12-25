#pragma once
#include <span>

using std::span;

// We mostly work with light-weight classes that don't themselves allocate memory.
// The following allocator classes make the allocation and object lifetime convienient.
namespace Allocators {
  template <typename node_t>
  struct sparsity {
    // We include both row_start and row_end, as we allow for fewer actual nonzeros than the capacity,
    // in order to enable dynamic updating within the bounds of a maximum degree.
    // If the sparsity pattern is compressed, then row_end[i] == row_start[i+1].
    vector<node_t> values, row_start, row_end;
    constexpr static node_t EMPTY_NODE = -1;
  
    sparsity(size_t N=0/*number of rows*/,
	     size_t nnz=0/*total number of nonzeros*/) : values(nnz,EMPTY_NODE), row_start(N), row_end(N) {}
  };

  // The special case of rectangular sparse matrices (same number of nonzeros for every row).
  template <typename node_t>
  struct rectangular_sparsity: public sparsity<node_t> {
    using parent = sparsity<node_t>;
    using parent::values;
    using parent::row_start;
    using parent::row_end;
    
    rectangular_sparsity(size_t N=0,  size_t d_max=3) : sparsity<node_t>(N*d_max) {
      for(size_t i=0;i<N;i++){ row_start[i] = d_max*i; row_end[i] = d_max*i; }
    }

    rectangular_sparsity(const vector<vector<node_t>>& xs) : rectangular_sparsity(xs.size(),size_max(xs)) { }

  private:
    size_t size_max(const vector<vector<node_t>>& xs){
      size_t mx = 0;
      for(const auto& x: xs) mx = max(mx,xs.size());
      return mx;
    }
  };
}

namespace Views {   

  template <typename node_t>
  struct sparsity {
    span<node_t> values, row_start, row_end;
    constexpr static node_t EMPTY_NODE = -1;
    
    sparsity(Allocators::sparsity<node_t> &G):
      values(G.values), row_start(G.row_start), row_end(G.row_end) {}

    Allocators::sparsity<node_t> compress() const {
      Allocators::sparsity C(N(), number_nonzeros());
      node_t offset = 0;
      for(node_t u=0;u<N();u++){
	size_t row_size = row(u).size();
	C.row_start[u]  = offset;
	C.row_end[u]    = offset + row_size;
	
	for(int i=0;i<row_size;i++) C.values[offset+i] = values[row_start[u]+i];
	offset         += row_size;
      }
      return C;
    }
    
    span<node_t>       row(node_t u)       { return {&values[row_start[u]], &values[row_end[u]]}; };
    span<const node_t> row(node_t u) const { return {&values[row_start[u]], &values[row_end[u]]}; };  

    span<node_t> operator[](node_t u)       { return row(u); }
    span<node_t> operator[](node_t u) const { return row(u); }    
    
    
    node_t N() const { return row_start.size(); }

    // Total number of nonzeros
    node_t number_nonzeros() const {
      node_t nnz = 0;
      for(node_t u: values) nnz += (u!=EMPTY_NODE);
      return nnz;
    }
    // Maximum number of nonzeros in row u
    node_t capacity(node_t u) const {
      if(u+1<N()) return row_start[u+1]-row_start[u];
      else        return values.size() -row_start[u];
    }
    
    node_t arc_index(node_t u, node_t v) const {
      auto ru = row(u);
      for(int i=0;i<ru.size();i++) if(ru[i]==v) return i;
      return EMPTY_NODE; // u-v is not an edge in this graph
    }
    
    bool arc_exists(node_t u, node_t v) const {
      return (arc_index(u,v)>=0);
    }
    
    bool is_symmetric() const {
      for(node_t u=0;u<N();u++)
	for(auto v: row(u))
	  if(!arc_exists(u,v)) return false;
      return true;
    }
  };
}


#if 0
template <typename node_t>
struct static_graph: public sparse_boolean_matrix<node_t> {
  using arc_t = array<node_t,2>;
  
  node_t next(node_t u, node_t v) const {
    const auto ru = row(u);
    int  j   = arc_index(u,v);
    if(j>=0) return ru[(j+1)%ru.size()];
    else return -1;		// u-v is not an edge in this graph
  }

  node_t prev(node_t u, node_t v) const {
    const auto ru = row(u);
    int  j   = arc_index(u,v);
    if(j>=0) return ru[(j+ru.size()-1)%ru.size()];
    else return -1;		// u-v is not an edge in this graph
  }  

  node_t next_on_face(node_t u, node_t v) const { return prev(v,u); }
  node_t prev_on_face(node_t u, node_t v) const { return next(v,u); }
  

  int degree(node_t u) const {
    return row(u).size();	// TODO: Won't work with rectangular format (empty slots will be counted).
  }

  
  int max_degree() const {
    size_t max_deg = 0;
    for(node_t u=0;u<N();u++) max_deg = max(max_degree,degree(u));
    return max_deg;
  }
  
  bool has_separating_triangles() const
  {
    for(node_t u=0;u<N();u++){
      const auto &ru(row(u));

      for(auto &t: ru){
	node_t v = prev(u,t), w = next(u,t); // edges: u--t, u--v, u--w
	// Invariant: We always have full edges, i.e. u->v exists iff v->u exists
	if(arc_exists(t,w) && arc_exists(t,v) && arc_exists(v,w)) return true;
      }
    }
    return false;
  }

  // TODO: General graph algorihms like SSSP, MSSP, APSP, connected components, etc?
  //       (even though not all make sense on surface)

};



template <typename node_t>
struct dynamic_graph {
  using arc_t = array<node_t,2>;
  // Insert edge u-v, represented by arcs u->v and v->u. Successor to u->v around u
  // and successor to v->u around v must be known in order to retain coherent
  // orientation.
  bool insert_edge(const arc_t& uv, node_t suc_uv=-1, node_t suc_vu=-1);

  // Remove edge u-v, represented by the two arcs u->v and v->u.
  bool remove_edge(const arc_t& uv);

  bool arc_exists(const arc_t& uv) const;
};
#endif
