#include <fullerenes/geometry.hh>
#include <fullerenes/spiral.hh>
#include <fullerenes/triangulation.hh>

typedef dedge_t arc_t; 		// TODO: Rename centrally

template <typename T> T pop(vector<T> &s){
  T& x = s.back(); s.pop_back(); return x;
}

template <typename T> size_t push(vector<T>& s, const T& x){
  //  cout << "push " << x << " onto " << s << endl;
  s.push_back(x); return s.size();
}

// TODO: map -> collection of lookup arrays 
unordered_map<arc_t, int> generate_Tid_by_arc(Triangulation &dual)
{
  unordered_map<arc_t,int> Tid_by_arc;
  for(int i=0;i<dual.triangles.size();i++){
    tri_t t = dual.triangles[i];
    for(int j=0;j<3;j++) Tid_by_arc[{t[j],t[(j+1)%3]}] = i;
  }
  // Every arc should be accounted for
  assert(get_values(Tid_by_arc).size() == (dual.N-2)*6);
  
  return Tid_by_arc;
}

// When rasterizing, each triangle "owns" exactly one of its edges, for which the fractional "pixels" are included in the triangle rasterization.
// Its neighbouring triangles owns the other two edges. This is to make rasterization fit together like a puzzle without gaps or overlaps.
//
// TODO: Reorganize the two lookup maps to flat array lookups 
//       (Most of the time goes into building Tid_by_arc and used)
vector<arc_t> compute_owned_edges(Triangulation &dual)
{
  unordered_map<arc_t,int> Tid_by_arc = generate_Tid_by_arc(dual); // TODO: Generate metadata outside in O(N) time

  size_t N = dual.triangles.size(), Nf = dual.N;
  
  vector<bool>   seen(N);
  vector<arc_t>  owned_arc(N);

  // 0. Start with triangle 0
  vector<int> work_stack = {0};
  unordered_map<arc_t,bool> used; 	// TODO: Replace by array lookup

  // 1. For each remaining triangle t in the work set
  while(!work_stack.empty()){
    int triangle_id = pop(work_stack);
    tri_t t = dual.triangles[triangle_id];

    // 2. Process t by picking first unused edge as owned edge and mark t as seen
    if(!seen[triangle_id])
      for(int i=0;i<3;i++){
    	node_t u = t[i], v = t[(i+1)%3];
    	if(!used[{v,u}]){
    	  //	  cout << triangle_id << "/" << t << " owns arc " << arc_t{u,v} << endl;
    	  used[{u,v}] = true;
    	  owned_arc[triangle_id] = {u,v};
    	  break;
    	}
      }
    seen[triangle_id] = true;        

    // 3. Now the adjacent triangles are ready to be processed
    for(int i=0;i<3;i++){
      node_t u = t[i], v = t[(i+1)%3];
      int tid_vu = Tid_by_arc[arc_t{v,u}];

      if(!seen[tid_vu]) push(work_stack, tid_vu);
    }
  }

  // TODO: Check correctness. Each triangle must own one unique edge, and all edges must be spoken for (NOT POSSIBLE, AS N_E = 3/2 N!) 
  
  return owned_arc;		
}

int main(int ac, char **av)
{
  if(ac != 4){
    printf("Syntax: %s <K> <L> <spiral>\n", av[0]);
    return -1;
  }

  int K = strtol(av[1],0,0), L = strtol(av[2],0,0);

  spiral_nomenclature spiral(av[3]);
  Triangulation dual{spiral};

  vector<arc_t> owned_edges = compute_owned_edges(dual);

  cout << owned_edges << endl;
  
  return 0;
}
