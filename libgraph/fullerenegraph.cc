#include "fullerenegraph.hh"
#include "fortran.hh"

#include <fstream>
#include <vector>
#include <deque>
#include <utility> //required for pair

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
  return G;
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

  //  cerr << "leapfrog()\n";
  vector<face_t> faces(dualfrog.compute_faces_flat(6,planar_layout)); 
  //  cerr << "leapfrog::got "<< faces.size() << " faces.\n";

  node_t v_new = N;

  if(planar_layout)
    dualfrog.layout2d.resize(N+faces.size());

  for(size_t i=0;i<faces.size();i++){
    const face_t& f(faces[i]);
    //    cerr << "Face " << i << ": " << f << endl;
    for(size_t j=0;j<f.size();j++)
      dualfrog.edge_set.insert(edge_t(v_new,f[j]));
    
    if(planar_layout)
      dualfrog.layout2d[v_new] = f.centroid(layout2d);

    v_new++;
  }
  dualfrog.update_auxiliaries();

  // Note that dualfrog is no longer planar, but is a triangulation of the sphere.
  // The dual of dualfrog becomes planar again.
  vector<coord2d> new_layout;

  if(planar_layout){
    //    cerr << "leapfrog::find_outer_face and compute faces\n";
    if(outer_face.size() < 5) outer_face = find_outer_face();

    vector<face_t> triangles(dualfrog.compute_faces_flat(3,false));

    //    cerr << "leapfrog::planar layout of " << triangles.size() << " triangles\n";
    new_layout.resize(triangles.size());

    for(int i=0;i<triangles.size();i++){
      const face_t& t(triangles[i]);
      new_layout[i] = t.centroid(dualfrog.layout2d)*coord2d(1,-1);
      if(t[0] == N || t[1] == N || t[2] == N) // triangle belongs to old outer face
	new_layout[i] *= 2.0;		      // TODO: Ensure that loop encompasses remaining graph
    }
  } 
  //  cerr << "leapfrog::dual()\n";
  FullereneGraph frog(dualfrog.dual_graph(3), new_layout);

  return frog;
}

void wg_connect(const int i, const int j, set<edge_t> &edge_set, std::vector<int> &uv)
{
  edge_set.insert(edge_t(i,j));
  ++uv[i];
  ++uv[j];
}

bool do_windup_general(const int n, const std::vector<int> &pot_spiral, std::vector<int> &pos, std::vector<int> &dist, set<edge_t> &edge_set){
  //number of used valencies per vertex (0-5(6))
  std::vector<int> used_valencies(n,0);
  //list of vertices that have open valencies
  std::vector<int> open_valencies;

  //connect first two faces
  wg_connect(0, 1, edge_set, used_valencies);

  open_valencies.push_back(0);
  open_valencies.push_back(1);

  int x; //jump distance (x=1 is no jump)
  
  //iterate over atoms
  //k=0, k=1 have been done already
  for (int k=2; k<n; ++k){
    //check if jump
    x = 0;
    if(pot_spiral[open_valencies.front()] - used_valencies[open_valencies.front()] == 2 && open_valencies.size() > 6){
      while(pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 2){
        ++x;
      }
      //two error cases
      if (pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+2)] - used_valencies[*(open_valencies.begin()+x+2)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+3)] - used_valencies[*(open_valencies.begin()+x+3)] == 1)
                {std::cerr << "There is no spiral." << std::endl;
        return 1;}
      if (pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+2)] - used_valencies[*(open_valencies.begin()+x+2)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x)] == 5)
                {std::cerr << "There is no spiral." << std::endl;
        return 1;}
      //two jump cases
      if ((pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
         pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
         pot_spiral[*(open_valencies.begin()+x)] == 5) || 
        (pot_spiral[*(open_valencies.begin()+x)] - used_valencies[*(open_valencies.begin()+x)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+1)] - used_valencies[*(open_valencies.begin()+x+1)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x+2)] - used_valencies[*(open_valencies.begin()+x+2)] == 1 &&
        pot_spiral[*(open_valencies.begin()+x)] == 6))
        ++x;
        else
        x=0;
    }

    //jump positions and distances
    if(x>1){// x=1 is no jump
      pos.push_back(k);
      dist.push_back(x);
    }

    // perform cyclic rotation on open_valencies
    for(int i = 1; i<x; ++i){
      int j = open_valencies.front();
      open_valencies.erase(open_valencies.begin());
      open_valencies.push_back(j);
    }

    //connect k to k-1
    wg_connect(open_valencies.back(), k, edge_set, used_valencies);

    //connect k to k-2, etc
    while(open_valencies.size() != 0 && pot_spiral[open_valencies.back()] - used_valencies[open_valencies.back()] == 0){
//      std::cout << pot_spiral[open_valencies.back()] << used_valencies[open_valencies.back()] << std:: endl;
      open_valencies.erase(open_valencies.end()-1);
      if(open_valencies.size() != 0 && pot_spiral[k] - used_valencies[k] != 0){
        wg_connect(open_valencies.back(), k, edge_set, used_valencies);
      }
    }

//    std::cout << k << "vii" << std::endl;
    //connect k to oldest unconnected (etc)
    if(pot_spiral[k] - used_valencies[k] != 0){
      wg_connect(open_valencies.front(), k, edge_set, used_valencies);
    }
//    std::cout << open_valencies.front() << " " << pot_spiral[open_valencies.front()] << " " << used_valencies[open_valencies.front()] << std:: endl;
    while(open_valencies.size() != 0 && pot_spiral[open_valencies.front()] - used_valencies[open_valencies.front()] == 0){
      open_valencies.erase(open_valencies.begin());
      if(open_valencies.size() != 0 && pot_spiral[k] - used_valencies[k] != 0){
        wg_connect(open_valencies.front(), k, edge_set, used_valencies);
      }
    }
    
    //append k to the open valencies (after making sure it has some valencies left)
    if (pot_spiral[k] - used_valencies[k] != 0){//the current atom is not saturated (which may only happen for the last one)
      open_valencies.push_back(k);
    }else{
      if(k + 1 != n){
        std::cout << "Fail 1 (cage closed but faces left)" << std::endl;
        return 1;
      }
    }
  }//iterate over atoms
  if(open_valencies.size() != 0){
      std::cout << "Fail 2 (cage not closed but no faces left)" << std::endl;
      return 1;
  } 
  return 0;//success
}//windup_general


FullereneGraph::FullereneGraph(const int n, const std::vector<int> spiral_indices, bool IPR, bool general) : CubicGraph() {
  if(!general){
//    int m = n/2+2;
//    int s[m], d[m*m], ipr = IPR, error = 0;
//    cerr << "Spiral constructor: " << n << ", " << face_t(spiral_indices) << endl; 
//    assert(spiral_indices.size() == 12);
// 
//    // Call Peter's fortran routine for constructing dual from spiral.
//    for(int i=0;i<n/2+2;i++) s[i] = 6;
//    for(int i=0;i<12;i++) s[spiral_indices[i]-1] = 5;
//
//    windup_(&m,&ipr,&error,s,d);
//    if(error != 0){
//      fprintf(stderr,"Spiral windup failed after %d pentagons.\n",error);
//      //    delete d;
//      N = 0;
//    } else {
//      //    cerr << " Spiral windup is successful.\n";
//      PlanarGraph dual;
//      //    printf("Dual should have %d nodes\n",m);
//      for(node_t u=0;u<m;u++)
//        for(node_t v=0;v<m;v++)
//      if(d[u*m+v] == 1)
//        dual.edge_set.insert(edge_t(u,v)); 
//
//      //    delete d;
//      cerr << "dual = " << dual << endl;
//      dual.update_auxiliaries();
//      dual.layout2d = dual.tutte_layout(-1,-1,-1,3);
//
//      *this = dual.dual_graph(3);
//      //    cerr << "dual = " << dual << endl;
//      //    cerr << "G    = " << G << endl;
//
//    }
  }
  else
  {
    assert(spiral_indices.size() == 12);

    std::vector<int> potential_spiral (n,6);
    for (int i=0; i<12; ++i){
      //std::cout << spiral_indices[i] << " " ;
      potential_spiral[spiral_indices[i]-1] = 5;//because the spiral input starts at 1 but this vector starts at 0
    }

    set<edge_t> edge_set;
    std::vector<int> jump_positions;
    std::vector<int> jump_distances;
    
    if(do_windup_general(n, potential_spiral, jump_positions, jump_distances, edge_set)){
      std::cerr << "No general spiral found either ... aborting." << std::endl;
      abort();
    }
    
    std::cout << jump_positions.size() << " jump(s) required.";
    for (std::vector<int>::iterator i(jump_positions.begin()); i<jump_positions.end(); ++i) {
      std::cout << *i+1 << ", ";//because k is relative to 0
    }
    
    for (std::vector<int>::iterator i(jump_distances.begin()); i<jump_distances.end(); ++i) {
      std::cout << *i << ", ";
    }
    std::cout << std::endl;

    PlanarGraph dual(edge_set);
    dual.update_auxiliaries();

    *this = dual.dual_graph(3);
  }
}

void gpi_connect_forward(std::deque<pair<int,int> > &open_valencies){
  --open_valencies.back().second;
  --open_valencies.front().second;
}

void gpi_connect_backward(std::deque<pair<int,int> > &open_valencies){
  --open_valencies.back().second;
  --(*(open_valencies.end() -2)).second;//decrement the last but one entry
}

// debug only
void pdp(std::deque<pair<int,int> > &open_valencies){
  for(std::deque<pair<int,int> >::iterator it(open_valencies.begin()); it!= open_valencies.end(); ++it){
     std::cout << it->first << ": " << it->second << std::endl;
  }
}

void gpi_remove_node(const node_t i, PlanarGraph &remaining_dual, std::set<int> &remaining_nodes){
  remaining_nodes.erase(i);
  //iterate over the neighbours of i and delete all connecting edges. it is legal to delete non-existing objects from a set, which will happen if the neighbour vector hasn't been updated
  for(std::vector<node_t>::iterator it = remaining_dual.neighbours[i].begin(); it != remaining_dual.neighbours[i].end(); ++it){
    remaining_dual.edge_set.erase(edge_t(i,*it));
  }
}

// perform a general general spiral search and return 12 pentagon indices and the jump positions + their length
void FullereneGraph::get_pentagon_indices(const int f1, const int f2, const int f3, std::vector<int> &pentagon_indices, std::vector<pair<int,int> >&jumps) const {

  std::cout << "entering 'get_pentagon_indices'" << std::endl;

  PlanarGraph dual = this->dual_graph(6);
  dual.update_auxiliaries();

  std::cout << "dual created: " << dual << std::endl;

  // remaining_dual is the graph that consists of all nodes that haven't been added to the graph yet
  PlanarGraph remaining_dual(dual);
  remaining_dual.update_auxiliaries();
  std::cout << "remaining dual created" << std::endl;
  // all the nodes that haven't been added yet, not ordered and starting at 0
  std::set<int> remaining_nodes;

  // the spiral is a string of numbers 5 and 6 and is built up during the loop
  std::deque<int> spiral; 
  // valencies is a list of length N and contains the valencies of each node (5 or 6)
  std::vector<int> valencies(dual.N, 0);  
  // open_valencies is a list with one entry per node that has been added to the spiral but is not fully saturated yet.  The entry contains the number of the node and the number of open valencies
  std::deque<pair<int,int> > open_valencies;

  //init of the valency-list and the set of nodes in the remaining graph
  for(int i=0; i<remaining_dual.N; ++i){
    valencies[i] = remaining_dual.neighbours[i].size();
    //std::cout << i << ": " << valencies[i]<< std::endl;
    remaining_nodes.insert(i);
  }

  // add the first three (defining) nodes
  //first node
  spiral.push_back(valencies[f1]);
  gpi_remove_node(f1, remaining_dual, remaining_nodes);
  open_valencies.push_back(make_pair(f1,valencies[f1]));

  //second node
  spiral.push_back(valencies[f2]);
  gpi_remove_node(f2, remaining_dual, remaining_nodes);
  open_valencies.push_back(make_pair(f2,valencies[f2]));
  gpi_connect_backward(open_valencies);
  //pdp(open_valencies);

  //third node
  spiral.push_back(valencies[f3]);
  gpi_remove_node(f3, remaining_dual, remaining_nodes);
  open_valencies.push_back(make_pair(f3,valencies[f3]));
  gpi_connect_backward(open_valencies);
  gpi_connect_forward(open_valencies);
  //pdp(open_valencies);

  // iterate over all nodes (of the dual) but not by their respective number
  // starting at 3 because we added 3 already
  for(int i=3; i<dual.N; ++i){

    //PlanarGraph remaining_dual_bak(remaining_dual);

    // find *the* node in dual (not the remaining_dual), that is connected to open_valencies.back() und open_valencies.front()
    // we can't search in the remaining_dual because there are some edges deleted already
    std::set<int>::iterator j=remaining_nodes.begin();
    node_t u = open_valencies.back().first, w = open_valencies.front().first;
    for( ; j!=remaining_nodes.end(); ++j){
    //   std::cout << open_valencies.back().first<< ", " << open_valencies.front().first << ", " << *seti << ", " << valencies[*seti] << std::endl;
      if(dual.edge_set.find(edge_t(u,*j)) != dual.edge_set.end() &&
         dual.edge_set.find(edge_t(w,*j)) != dual.edge_set.end()) break;
    }
    assert(j!=remaining_nodes.end());

    std::cout << "adding node " << *j << std::endl;
    spiral.push_back(valencies[*j]);
    open_valencies.push_back(make_pair(*j,valencies[*j]));
    gpi_connect_backward(open_valencies);
    gpi_connect_forward(open_valencies);

   // there are three positions in open_valencies that can be 0---one shouldn't happen, the other two cases requires interaction.
   while(true){
     bool a=0;
     assert(open_valencies.back().second!=0);//can only happen if the psiral missed a jump
     if((*(open_valencies.end()-2)).second==0){
       open_valencies.erase (open_valencies.end()-2);
       gpi_connect_backward(open_valencies);
       a=1;
     }
     if(open_valencies.front().second==0){
       open_valencies.pop_front();
       gpi_connect_forward(open_valencies);
       a=1;
     }
     if(!a) break;
   }

   node_t v = *j;
   gpi_remove_node(*j, remaining_dual, remaining_nodes); //remove all edges of which *j is part from the remaining dual
   pdp(open_valencies);

   if(!remaining_dual.is_connected()){
     std::cout << "entering reversion" << std::endl;
     //revert the last operations
     remaining_nodes.insert(v);
     //perform  cyclic rotation on open_valencies
     open_valencies.push_back(open_valencies.front());
     open_valencies.pop_front();
   } else {
     cout << "\nremaining_dual = " << remaining_dual << " is connected.\n";
   }
  }
  
  int k=0;
  for(std::deque<int>::iterator it=spiral.begin(); it != spiral.end(); ++it, ++k){
    //std::cout << j << ", " << *it << std::endl;
    if(*it==5){
      pentagon_indices.push_back(k+1);
      //std::cout << j << ", ";
    }
  }
  std::cout << std::endl;
  assert(pentagon_indices.size()==12);

}

node_t FullereneGraph::C20_edges[30][2] ={{0,13},{0,14},{0,15},{1,4},{1,5},{1,12},{2,6},{2,13},{2,18},{3,7},{3,14},{3,19},{4,10},{4,18},{5,11},{5,19},{6,10},{6,15},{7,11},{7,15},{8,9},{8,13},{8,16},{9,14},{9,17},{10,11},{12,16},{12,17},{16,18},{17,19}};

double FullereneGraph::C20_layout2d[20][2] = {{1.548,0.503},{-1.134,-0.368},{0.,1.628},{0.957,-1.317},{-1.548,0.503},{-0.957,-1.317},{0.,2.466},{1.449,-1.995},{0.302,0.416},{0.489,-0.159},{-2.345,0.762},{-1.45,-1.995},{-0.489,-0.159},{0.7,0.965},{1.133,-0.369},{2.345,0.762},{-0.302,0.416},{0.,-0.514},{-0.7,0.965},{0.,-1.192}};
