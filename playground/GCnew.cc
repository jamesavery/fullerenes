#include <fullerenes/geometry.hh>
#include <fullerenes/spiral.hh>
#include <fullerenes/triangulation.hh>

typedef dedge_t arc_t; 		// TODO: Rename centrally

template <typename T> T pop(vector<T> &s){
  T& x = s.back(); s.pop_back(); return x;
}

template <typename T> size_t push(vector<T>& s, const T& x){ s.push_back(x); return s.size(); }

// When rasterizing, each triangle "owns" exactly one of its edges, for which the fractional "pixels" are included in the triangle rasterization.
// Its neighbouring triangles owns the other two edges. This is to make rasterization fit together like a puzzle without gaps or overlaps.
vector<arc_t> owned_edges(Triangulation &dual, map<arc_t,int> &Tid)
{
  size_t N = dual.triangles.size(), Nf = dual.N;
  
  vector<bool>   seen(N,false);
  vector<arc_t>  owned_arc(N);
  
  vector<int> work_stack = {0};
  map<arc_t,bool> used; 	// TODO: Replace by array lookup
  
  while(!work_stack.empty()){
    int triangle_id = pop(work_stack);
    tri_t t = dual.triangles[triangle_id];
    
    seen[triangle_id] = true;

    for(int i=0;i<3;i++){
      node_t u = t[i], v = t[(i+1)%3];
      int tid_uv = Tid[arc_t{u,v}];

      if(!seen[tid_uv]) push(work_stack, tid_uv);

      if(!used[{v,u}]){
	used[{u,v}] = true;
	owned_arc[triangle_id] = {u,v};
	break;
      }
    }
  }
  
  return owned_arc;		
}

int main(int ac, char **av)
{
  if(ac != 3){
    printf("Syntax: %s <K> <L> <spiral>\n", av[0]);
    return -1;
  }

  int K = strtol(av[1],0,0), L = strtol(av[2],0,0);

  spiral_nomenclature spiral(av[3]);
  Triangulation dual(spiral_nomenclature);

  
  return 0;
}
