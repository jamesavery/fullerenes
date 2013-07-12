#include "fullerenegraph.hh"
#include "fortran.hh"

#include <fstream>
#include <vector>
#include <list>
#include <vector>
#include <utility> //required for pair

typedef list<pair<int,int> > jumplist_t;

pair<set< face_t>, set<face_t> > FullereneGraph::compute_faces56() const 
{
  set<face_t> pentagons, hexagons;
  for(set<edge_t>::const_iterator e(edge_set.begin()); e!= edge_set.end(); e++){
    const node_t s = e->first, t = e->second;

    const vector<node_t>& ns(neighbours[t]);
    for(unsigned int i=0;i<ns.size();i++)
      if(ns[i] != s) {
	const node_t u = ns[i];
	
	// Assumes fullerene. Remove 6 and return map instead of pair to allow cubic graphs with arbitrary polygons
	face_t face(shortest_cycle(s,t,u,6)); 
	if      (face.size() == 5) pentagons.insert(face);
	else if (face.size() == 6) hexagons.insert(face);
	else {
	  fprintf(stderr,"Graph is not a fullerene: Contains face ");
	  for(unsigned int i=0;i<face.size();i++)
	    fprintf(stderr,"%d ",face[i]);
	  fprintf(stderr,"of size %d\n",int(face.size()));
	}
    }
  }
  return pair< set<face_t>,set<face_t> >(pentagons,hexagons);
}

// Creates the m-point halma-fullerene from the current fullerene C_n with n(1+m)^2 vertices. (I.e. 4,9,16,25,36,... times)
FullereneGraph FullereneGraph::halma_fullerene(const int m, const bool planar_layout) const {
  if(m<0) return FullereneGraph(*this);

  PlanarGraph dual(dual_graph(6));
  vector<face_t> triangles(dual.compute_faces_flat(3,false));
  map<edge_t,vector<node_t> > edge_nodes;
  
  set<edge_t> edgeset_new;
  node_t v_new = dual.N;

  // TODO: Everything that's made from the outer face needs to be "reversed".
  vector<coord2d> new_layout;
  if(planar_layout && layout2d.size() == N){
    new_layout.resize(dual.N);
    for(int i=0;i<dual.N;i++) new_layout[i] = dual.layout2d[i];
  }
    
  // Create n new vertices for each edge
  for(set<edge_t>::const_iterator e(dual.edge_set.begin()); e!=dual.edge_set.end(); e++){
    vector<node_t>& nodes(edge_nodes[*e]);
    for(unsigned int i=0;i<m;i++) nodes.push_back(v_new++);
    
    if(planar_layout && layout2d.size() == N) 
      for(unsigned int i=0;i<m;i++){
	double lambda = (1.0+i)*(1.0/(m+1));
	const coord2d &a(dual.layout2d[e->first]), &b(dual.layout2d[e->second]);
	new_layout.push_back(a*(1.0-lambda) + b*lambda);
      }
  }

  // For every triangle in the dual, we create and connect a halma-type grid
  for(size_t i=0;i<triangles.size();i++){
    map<edge_t,node_t> grid;
    const face_t& T(triangles[i]);
    edge_t e0(T[0],T[1]),e1(T[1],T[2]),e2(T[2],T[0]);
    const vector<node_t>& ns0(edge_nodes[e0]), ns1(edge_nodes[e1]), ns2(edge_nodes[e2]);

    // Insert original vertices
    grid[edge_t(0,0)]     = T[0];
    grid[edge_t(m+1,0)]   = T[1];
    grid[edge_t(m+1,m+1)] = T[2];
    // Insert new edge vertices
    for(size_t j=0;j<m;j++){	
      grid[edge_t(0,j+1)]   = ns0[j];
      grid[edge_t(j+1,m+1)] = ns1[j];
      grid[edge_t(j+1,j+1)] = ns2[j];
    }
    // Create and insert inner vertices
    for(int j=1;j<m;j++)
      for(int k=j+1;k<=m;k++)
	grid[edge_t(j,k)] = v_new++;

    if(planar_layout && layout2d.size() == N){
      double sqrt2inv = 1.0/sqrt(2.0);
      const coord2d &a(dual.layout2d[T[0]]), &b(dual.layout2d[T[1]]), &c(dual.layout2d[T[2]]);      
      for(int j=1;j<m;j++)
	for(int k=j+1;k<=m;k++){
	  double s = (1+j)*(1.0/(m+2)), t = k*(1.0/(m+2));
	  //fprintf(stderr,"(s,t) = (%g,%g)\n",s,t);
	  new_layout.push_back(a+((b-a)*s + (c-a)*t)*sqrt2inv);
	}
    }
      
    // Connect the vertices in the grid
    for(int j=0;j<=m;j++)
      for(int k=j+1;k<=m+1;k++){
	node_t v(grid[edge_t(j,k)]), down(grid[edge_t(j+1,k)]), 
	  left(grid[edge_t(j,k-1)]);

	edgeset_new.insert(edge_t(v,down));
	edgeset_new.insert(edge_t(v,left));
	edgeset_new.insert(edge_t(left,down));
      }
  }

  //cerr << "new_layout.size() = " << new_layout.size() << endl;

  PlanarGraph new_dual(Graph(edgeset_new), new_layout);
  //cerr << "newdual.N = " << new_dual.N << endl;

  /*
  ofstream h("output/halma2.m");

  h << "graph = " << *this << endl;
  h << "dual = " << dual << endl;
  h << "newdual = " << new_dual << endl;
  */

  PlanarGraph G(new_dual.dual_graph(3));

  /*
  h << "halma = " << G << endl;

  h.close();
  */
  return FullereneGraph(G,G.layout2d);
}

unsigned int gcd(unsigned int a, unsigned int b)
{
  unsigned int r = a % b;
  if(r == 0) return b; else return gcd(b,r);
}


// Actually works for all cubic graphs -- perhaps stick it there instead
FullereneGraph FullereneGraph::coxeter_fullerene(const unsigned int i, const unsigned int j, const bool do_layout) const
{
  FullereneGraph CG;
  return CG;
}

// Creates the next leapfrog fullerene C_{3n} from the current fullerene C_n
FullereneGraph FullereneGraph::leapfrog_fullerene(bool planar_layout) const {
  PlanarGraph dualfrog(*this);

  vector<face_t> faces(dualfrog.compute_faces_flat(6,planar_layout)); 

  //  cout << "(*leapfrog*)outer_face = " << dualfrog.outer_face << ";\n";
  //  cout << "(*leapfrog*)faces      = " << faces << ";\n";

  node_t v_new = N;

  if(planar_layout)
    dualfrog.layout2d.resize(N+faces.size());

  for(size_t i=0;i<faces.size();i++){
    const face_t& f(faces[i]);
    for(size_t j=0;j<f.size();j++)
      dualfrog.edge_set.insert(edge_t(v_new,f[j]));
    
    if(planar_layout)
      dualfrog.layout2d[v_new] = f.centroid(layout2d);

    v_new++;
  }
  dualfrog.update_from_edgeset();

  // Note that dualfrog is no longer planar, but is a triangulation of the sphere.
  // The dual of dualfrog becomes planar again.
  vector<coord2d> new_layout;
  face_t new_outer_face;

  if(planar_layout){
    // The layout of dualfrog is not planar - faces must be computed without it
    vector<face_t> triangles(dualfrog.compute_faces_flat(3,false));
    new_layout.resize(triangles.size());

    for(int i=0;i<triangles.size();i++){
      const face_t& t(triangles[i]);
      new_layout[i] = t.centroid(dualfrog.layout2d)*coord2d(1,-1);
      for(int j=0;j<3;j++) if(t[j] == N){
	  //	  cout << "Triangle number " << i << " = " << t << " belongs to outer face.\n";
	  new_outer_face.push_back(i);
	  new_layout[i] *= 2.0; // TODO: Scale to be 1.1*radius
	  break;
      }
    }
  } 

  FullereneGraph frog(dualfrog.dual_graph(3,false), new_layout);

  sort_ccw_point CCW(new_layout,new_outer_face.centroid(new_layout));
  sort(new_outer_face.begin(),new_outer_face.end(),CCW);
  frog.outer_face = new_outer_face;

  return frog;
}

void wg_connect_backward(set<edge_t> &edge_set, list<pair<node_t, int> > &ov)
{
  list< pair<node_t,int> >::iterator second_last(ov.end());
  --second_last;
  --second_last;

  edge_set.insert(edge_t(ov.back().first, second_last->first));
  --ov.back().second;
  --(second_last->second);//decrement the last but one entry
}

void wg_connect_forward(set<edge_t> &edge_set, list<pair<node_t, int> > &ov)
{
  edge_set.insert(edge_t(ov.back().first, ov.front().first));
  --ov.back().second;
  --ov.front().second;
}

// // debug only (do not remove, please [lukas])
// void pdp(list<pair<int,int> > &open_valencies){
//   for(list<pair<int,int> >::iterator it(open_valencies.begin()); it!= open_valencies.end(); ++it){
//      cout << it->first << ": " << it->second << endl;
//   }
// }

bool do_windup_general(const int n_faces,  const vector<int> &spiral,  list<pair<int,int> > jumps,  set<edge_t> &edge_set){

  // open_valencies is a list with one entry per node that has been added to the spiral but is not fully saturated yet.  The entry contains the number of the node and the number of open valencies
  list<pair<int,int> > open_valencies;

  // set up first two nodes
  open_valencies.push_back(make_pair(0,spiral[0]));
  open_valencies.push_back(make_pair(1,spiral[1]));
  //connect first two faces
  wg_connect_backward(edge_set, open_valencies);

  //iterate over atoms
  //k=0, k=1 have been done already
  for (int k=2; k<n_faces-1; ++k){
//    cout << "k: " << k << endl;

    if(jumps.size() != 0 && k == jumps.front().first){
      // perform cyclic rotation on open_valencies
      for(int i = jumps.front().second; i>1; --i){ // 1 is no jump
        open_valencies.push_back(open_valencies.front());
        open_valencies.pop_front();
      }
      jumps.pop_front();
    }

    // add node to spiral
    open_valencies.push_back(make_pair(k,spiral[k]));

    // connect k to k-1
    wg_connect_backward(edge_set, open_valencies);

    // connect k to k-2, etc
    wg_connect_forward(edge_set, open_valencies);

    // do the remaining connect forwards
    while(open_valencies.front().second==0){
      open_valencies.pop_front();
      wg_connect_forward(edge_set, open_valencies);
    }
    // do the remaining connect backwards //not neat but the most simple way to emulate 'while second_last->second==0) ...'
    while(true){
      list<pair<int,int> >::iterator second_last(open_valencies.end());
      --second_last;
      --second_last;
      if(second_last->second==0){
        open_valencies.erase(second_last);
        wg_connect_backward(edge_set, open_valencies);
      } else break;
      
    }
//    pdp(open_valencies);

    if (open_valencies.back().second == 0){//the current atom is saturated (which may only happen for the last one)
      cout << "Cage closed but faces left (or otherwise invalid spiral)" << endl;
      return 1;
    }
   
  }//iterate over atoms

  // make sure we left the spiral in a sane state
  // open_valencies must be either 5 or 6 times '1' at this stage
  if(open_valencies.size() != spiral.back()){
    cout << "Cage not closed but no faces left (or otherwise invalid spiral), wrong number of faces left" << endl;
    return 1;
  }
  for(list<pair<int,int> >::iterator it = open_valencies.begin(); it!=open_valencies.end(); ++it){
    if(it->second!=1){
      cout << "Cage not closed but no faces left (or otherwise invalid spiral), more than one valency left for at least one face" << endl;
    return 1;
    }
  }

  // add last node to spiral // the name of the last node is n_faces -1 (because it's the last one)
  open_valencies.push_back(make_pair(n_faces-1,spiral[n_faces-1]));

  for(int i=0; i<spiral.back(); ++i){
    wg_connect_forward(edge_set, open_valencies);
    open_valencies.pop_front();
  }

  return 0;// success
}// do_windup_general


// both the pentagon indices and the jumps start at 0
// n is the number of vertices
FullereneGraph::FullereneGraph(const int n, const vector<int> spiral_indices, const list<pair<int,int> > jumps) : CubicGraph() {

  assert(spiral_indices.size() == 12);

  const int n_faces = n/2 + 2;

//  printf("n_faces = %d\n",n_faces);
  vector<int> potential_spiral (n_faces,6);
  for (int i=0; i<12; ++i){
//    printf("spiral_indices[%d] = %d\n",i,spiral_indices[i]);
    potential_spiral[spiral_indices[i]] = 5;
  }

  set<edge_t> edge_set;
  
  if(do_windup_general(n_faces, potential_spiral, jumps, edge_set)){
    cerr << "No general spiral found ... aborting.  This shouldn't happen unless the input is wrong." << endl;
    abort();
  }
  
  PlanarGraph dual(edge_set);
  dual.update_from_edgeset();

  *this = dual.dual_graph(3);
  fullerene_check();
}


// pentagon indices and jumps start to count at 0
// perform a general general spiral search and return 12 pentagon indices and the jump positions + their length
void FullereneGraph::get_general_spiral_from_fg(const node_t f1, const node_t f2, const node_t f3, vector<int> &pentagon_indices, jumplist_t &jumps) const {

  //this routine expects empty containers pentagon_indices and jumps.  we make sure they *are* empty
  pentagon_indices.clear();
  jumps.clear();

  PlanarGraph dual = this->dual_graph(6);
  dual.update_from_edgeset();

  // the spiral is a string of numbers 5 and 6 and is built up during the loop
  vector<int> spiral; 

  dual.get_vertex_spiral(f1, f2, f3, spiral, jumps);
  
  // extract spiral indices from spiral
  int k=0;
  for(vector<int>::iterator it=spiral.begin(); it != spiral.end(); ++it, ++k){
    if(*it==5){
      pentagon_indices.push_back(k);
    }
  }
  assert(pentagon_indices.size()==12);

}


// perform the canonical general general spiral search and return 12 pentagon indices and the jump positions + their length
void FullereneGraph::get_canonical_general_spiral_from_fg(vector<int> &pentagon_indices, jumplist_t &jumps) const {

  vector<int> pentagon_indices_tmp;
  vector<int> spiral_tmp;
  list<pair<int,int> > jumps_tmp;
  
  //100 times 0 to make sure size() is large (so it get's overwritten in the first cycle)
  vector<int> general_spiral_bak(100,0); // FIXME

  PlanarGraph dual(dual_graph(6));
  vector<face_t> faces(dual.compute_faces_flat(3));

//  cout << "generating all spirals ";

  for(int i=0; i<faces.size(); i++){
    int permutations[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    const face_t& f = faces[i];
    for(int j=0; j<6; j++){
      pentagon_indices_tmp.clear();

      int f1 = f[permutations[j][0]], f2 = f[permutations[j][1]], f3 = f[permutations[j][2]];

      dual.get_vertex_spiral(f1, f2, f3, spiral_tmp, jumps_tmp);

      // extract spiral indices from spiral
      int k=0;
      for(vector<int>::const_iterator it=spiral_tmp.begin(); it != spiral_tmp.end(); ++it){
        if(*it==5){
          pentagon_indices_tmp.push_back(k);
        }
        ++k;
      }
      assert(pentagon_indices_tmp.size()==12);

//      printf("Face %d:%d vertices defining the face(%d,%d,%d)\n",i,j,f1,f2,f3);

      //flatten and combine:
      vector<int> general_spiral;
      for(list<pair<int,int> >::const_iterator it(jumps_tmp.begin()); it!= jumps_tmp.end(); ++it){
        general_spiral.push_back(it->first);
        general_spiral.push_back(it->second);
      }
      general_spiral.insert(general_spiral.end(), pentagon_indices_tmp.begin(), pentagon_indices_tmp.end());

      // store the shortest / lexicographically smallest one
      if(general_spiral.size() < general_spiral_bak.size() || 
			(general_spiral.size() == general_spiral_bak.size() &&
			lexicographical_compare(general_spiral.begin(), general_spiral.end(), general_spiral_bak.begin(), general_spiral_bak.end()))){
		general_spiral_bak = general_spiral;
      }
    }
  }

  assert(general_spiral_bak.size() % 2 == 0);

  // get rspi
  vector<int> rspi(general_spiral_bak.end()-12, general_spiral_bak.end());
  pentagon_indices = rspi;

//  cout << "got rspi, size: " << general_spiral_bak.size() << endl;
//  cout << "got rspi: " << rspi << endl;

  //get jumps
  jumps.clear();
  while(general_spiral_bak.size() > 12){
    jumps.push_back(make_pair(general_spiral_bak.front(), *(general_spiral_bak.begin()+1)));
    general_spiral_bak.erase(general_spiral_bak.begin(), general_spiral_bak.begin()+2);
  }
//  cout << "got jumps, size: " << jumps.size() << endl;
//  cout << "got jumps: " << jumps << endl;

}


// create a matrix that holds the topological distances between all pentagons
vector<int> FullereneGraph::pentagon_distance_mtx() const{

  PlanarGraph dual(dual_graph(6));
  
  vector<int> pentagons;
  for(int i=0; i<dual.N; i++) if(dual.neighbours[i].size() == 5) pentagons.push_back(i);
  
  // I know the following line is very wasteful but it's good enough
  vector<int> all_distances = dual.all_pairs_shortest_paths(INT_MAX);
  vector<int> mtx_vec(144,0);
  
  //cout << "pentagon distances: " << endl;
  for(int i=0; i!=12; ++i){
    for(int j=0; j!=12; ++j){
      mtx_vec[12*i+j] = all_distances[dual.N * pentagons[i] + pentagons[j]];
    }
  }
 
//  cout << pentagon_distance_mtx << endl;
  
  return mtx_vec;

//mathematica output, please do not remove (lukas)
//    cout << "{";
//    for(int i=0; i!=12; ++i){
//      cout << "{";
//      for(int j=0; j!=12; ++j){
//        cout << pentagon_distances[12*i + j];
//        if(j!=11)cout << ", ";
//      }
//      cout << "}";
//      if(i!=11)cout << ",";
//      cout << endl;
//    }
//    cout << "}" << endl;
  
}


node_t FullereneGraph::C20_edges[30][2] ={{0,13},{0,14},{0,15},{1,4},{1,5},{1,12},{2,6},{2,13},{2,18},{3,7},{3,14},{3,19},{4,10},{4,18},{5,11},{5,19},{6,10},{6,15},{7,11},{7,15},{8,9},{8,13},{8,16},{9,14},{9,17},{10,11},{12,16},{12,17},{16,18},{17,19}};

double FullereneGraph::C20_layout2d[20][2] = {{1.548,0.503},{-1.134,-0.368},{0.,1.628},{0.957,-1.317},{-1.548,0.503},{-0.957,-1.317},{0.,2.466},{1.449,-1.995},{0.302,0.416},{0.489,-0.159},{-2.345,0.762},{-1.45,-1.995},{-0.489,-0.159},{0.7,0.965},{1.133,-0.369},{2.345,0.762},{-0.302,0.416},{0.,-0.514},{-0.7,0.965},{0.,-1.192}};
