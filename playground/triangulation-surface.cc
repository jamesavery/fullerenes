#include "libgraph/auxiliary.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/matrix.hh"

vector< vector<vector<tri_t>> > semisimple_geodesics() const
{
  matrix<int> H = matrix<int>(N,N,all_pairs_shortest_paths());
  int M = *max_element(H.begin(),H.end());      // M is upper bound to path length

  vector< vector<vector<node_t>> > Paths(N*N);

  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++){

      // Note: All Eisenstein numbers of the form (a,0) or (0,b) yield same lengths
      //       as graph distance, and are hence covered by initial step. So start from 1.
      //       M is upper bound for distance, so only need to do a^2+ab+b^2 strictly less than M.
      for(int a=1; a<M;a++)
        for(int b=1; a*a + a*b + b*b < M*M; b++){
          // Check: if(gcd(a,b) != 1) continue.
          vector<node_t> path = all_of_the_line(u,i,a,b);
	  node_t v = path[path.size()-1];

	  Paths[u*N+v].push_back(path);
        }
    }
  return Paths;
}

vector< vector<int> > Triangulation::semisimple_geodesic_lengths() const
{
  matrix<int> H = matrix<int>(N,N,all_pairs_shortest_paths());
  int M = *max_element(H.begin(),H.end());      // M is upper bound to path length
  vector<vector<int>> Path_lengths(N*N);

  for(int i=0;i<H.size();i++) H[i] *= H[i]; // Work with square distances, so that all distances are integers.

  for(node_t u=0;u<N;u++)
    for(int i=0;i<neighbours[u].size();i++){

      // Note: All Eisenstein numbers of the form (a,0) or (0,b) yield same lengths
      //       as graph distance, and are hence covered by initial step. So start from 1.
      //       M is upper bound for distance, so only need to do a^2+ab+b^2 strictly less than M.
      for(int a=1; a<M;a++)
        for(int b=1; a*a + a*b + b*b < M*M; b++){
          // Check: if(gcd(a,b) != 1) continue.
          const node_t v = end_of_the_line(u,i,a,b);

          // printf("min(H(%d,%d),|(%d,%d)|^2)  = min(%d,%d)\n",
          // 	 u,v,a,b,H(u,v), a*a+a*b+b*b);
          Path_lengths[u*N+v].push_back(a*a + a*b + b*b);
        }
    }
  return Path_lengths;
}



// Given start node u0 and adjacent face F_i, lay down triangles along the the straight
// line to Eisenstein number (a,b), and report what the final node is.
//
// Assumes a,b >= 1.
// TODO: Add special cases for (a,0) and (b,0) to make more general.
// TODO: Better name.
vector<node_t> Triangulation::all_of_the_line(node_t u0, int i, int a, int b) const
{
  node_t q,r,s,t;		// Current square

  // Square one
  q = u0; 			// (0,0)
  r = neighbours[u0][i];	// (1,0)
  s = next(q,r);	// (0,1)
  t = next(s,r);	// (1,1)

  vector<node_t> path{{u0,t}};
  if(a==1 && b==1) return path;
  
  vector<int> runlengths = draw_path(max(a,b), min(a,b));

  
  auto go_north = [&](){
    const node_t S(s), T(t); // From old square
    q = S; r = T; s = next(S,T); t = next(s,r);
    path.push_back(t);
  };

  auto go_east = [&](){
    const node_t R(r), T(t); // From old square
    q = R; s = T; r = next(s,q); t = next(s,r);
    path.push_back(t);
  };

  for(int i=0;i<runlengths.size();i++){
    int L = runlengths[i];

    if(a>=b){			// a is major axis
      for(int j=0;j<L-1;j++)    go_east();
      if(i+1<runlengths.size()) go_north();
    } else {			// b is major axis
      for(int j=0;j<L-1;j++)    go_north();

      if(i+1<runlengths.size()) go_east();
    }
  }

  return path;			// End node is upper right corner.
}



int main(int ac, char **av)
{
  assert(ac>13);
  int N = strtol(av[1],0,0);
  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[strtol(av[2+i],0,0)-1] = 5;

  cout << "spiral = " << spiral << ";\n";

  Triangulation   dG(spiral);
  PlanarGraph      G(dG.dual_graph());

  assert(dG.is_consistently_oriented());

  matrix<int> Hinit(dG.N,dG.N,dG.all_pairs_shortest_paths());
  cout << "Hinit = " << Hinit << ";\n\n";

  matrix<int> Hsimple(dG.convex_square_surface_distances());
  
  cout << "Hsimple = Sqrt[" << Hsimple << "];\n\n";

  matrix<double> H = dG.surface_distances();
  
  cout << "H = " << H << ";\n\n";

  vector< vector<int> > geodesic_distances = dG.semisimple_geodesic_lengths();

  cout << "GD = " << geodesic_distances << ";\n\n";

  vector< vector<vector<node_t>> > geodesics = dG.semisimple_geodesics();

  cout << "geodesics = " << geodesics << "\n\n";
  
  return 0;
}
