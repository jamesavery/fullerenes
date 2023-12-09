// TODO
// 0. In SWsearch5: Generate A[i,j] = (d,p,i) where d is smallest d with swsite, and (p,i) is swsite 
//    (Possible optimization: Only include d>1 if connecting two components.)
// 1. vector<T> input
// 2. Breadth-first traversal

// Procedure:
// 0. Generate shape of input isomer
// 1. Read in matrix A
// 2. Breadth-first traversal as follows:
//  2.0 Push current structure on stack
//  2.1 For each child
//  2.1.1 Perform swtransform
//  2.1.2 Recurse

#include <libgraph/auxiliary.hh>
#include <libgraph/triangulation.hh>
#include <libgraph/isomerdb.hh>
#include <libgraph/patches.hh>
#include <stack>
#include <queue>

#include <iostream>
#include <netcdf.h>
using namespace std;


template <typename S> class SparseMatrix: public Graph {
public:
  using Graph::Graph;
  
  map<arc_t, S> values;
  
  const S& add_edge(const arc_t& uv, const S& x) {
    insert_arc(uv,-1);
    values[uv] += x;
    return values[uv];
  }

  friend ostream& operator<<(ostream& s, const SparseMatrix& A)
  {
    vector<arc_t> edges = A.directed_edges();
    vector<S> edge_values(edges.size());

    for(int i=0;i<edges.size();i++){
	auto x = A.values.find(edges[i]);
	edge_values[i] = x->second;
      }

    s << "{" << edges << ", " << edge_values  << "}";
    return s;
  }  
};

class SWshape {
public:
  size_t N, Nisomers;
  typedef struct { int d, p, v; } d_swsite_t;
  map< arc_t, d_swsite_t > Breduced;
  vector<arc_t> spanning_tree;
  
  SWshape(const string ncfile_path){
    int ncid, varid, dimid, rv;
    size_t dim;

    rv = nc_open(ncfile_path.c_str(), NC_NOWRITE, &ncid);

    rv = nc_inq_dimid(ncid,"N",&dimid);
    rv = nc_inq_dimlen(ncid,dimid,&N);

    rv = nc_inq_dimid(ncid,"Nisomers",&dimid);
    rv = nc_inq_dimlen(ncid,dimid,&Nisomers);

    rv = nc_inq_dimid(ncid,"B_arc",&dimid);
    rv = nc_inq_dimlen(ncid,dimid,&dim);
    vector<arc_t>    Bkeys(dim);    
    vector<d_swsite_t> Bvalues(dim);

    rv = nc_inq_varid(ncid,"B_arc",&varid);
    rv = nc_get_var_int(ncid,varid,(int*)&Bkeys[0]);
    rv = nc_inq_varid(ncid,"B_values",&varid);
    rv = nc_get_var_int(ncid,varid,(int*)&Bvalues[0]);

    rv = nc_inq_dimid(ncid,"tree",&dimid);
    rv = nc_inq_dimlen(ncid,dimid,&dim);
    spanning_tree = vector<arc_t>(dim);

    rv = nc_inq_varid(ncid,"tree",&varid);
    rv = nc_get_var_int(ncid,varid,(int*)&spanning_tree[0]);

    for(auto v: Bvalues)
      cerr << vector<int>((int*)&v,(int*)&v+3) << endl;

    rv = nc_close(ncid);
  }
};

int main()
{
  SWshape("SWsearch5-C40.nc");

  return 0;
}
