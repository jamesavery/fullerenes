#include "fullerenegraph.hh"

#include <fstream>

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
FullereneGraph FullereneGraph::halma_fullerene(const unsigned int m, const bool planar_layout) const {
  PlanarGraph dual(dual_graph(6));
  vector<face_t> triangles(dual.compute_faces_flat(3,false));
  map<edge_t,vector<node_t> > edge_nodes;
  
  set<edge_t> edgeset_new;
  node_t v_new = dual.N;

  vector<coord2d> new_layout;
  if(planar_layout){
    new_layout.resize(dual.N);
    for(int i=0;i<dual.N;i++) new_layout[i] = dual.layout2d[i];
  }
    
  // Create n new vertices for each edge
  for(set<edge_t>::const_iterator e(dual.edge_set.begin()); e!=dual.edge_set.end(); e++){
    vector<node_t>& nodes(edge_nodes[*e]);
    for(unsigned int i=0;i<m;i++) nodes.push_back(v_new++);
    
    if(planar_layout) 
      for(unsigned int i=0;i<m;i++){
	double lambda = (1.0+i)*(1.0/(m+1));
	cerr << "lambda = " << lambda << endl;
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

    if(planar_layout){
      double sqrt2inv = 1.0/sqrt(2.0);
      const coord2d &a(dual.layout2d[T[0]]), &b(dual.layout2d[T[1]]), &c(dual.layout2d[T[2]]);      
      for(int j=1;j<m;j++)
	for(int k=j+1;k<=m;k++){
	  double s = (1+j)*(1.0/(m+2)), t = k*(1.0/(m+2));
	  fprintf(stderr,"(s,t) = (%g,%g)\n",s,t);
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

  cerr << "new_layout.size() = " << new_layout.size() << endl;

  PlanarGraph new_dual(Graph(edgeset_new), new_layout);
  cerr << "newdual.N = " << new_dual.N << endl;

  ofstream h("output/halma2.m");

  h << "graph = " << *this << endl;
  h << "dual = " << dual << endl;
  h << "newdual = " << new_dual << endl;

  PlanarGraph G(new_dual.dual_graph(3));
  h << "halma = " << G << endl;

  h.close();

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

  cerr << "leapfrog()\n";
  vector<face_t> faces(dualfrog.compute_faces_flat(6,planar_layout)); 
  cerr << "leapfrog::got "<< faces.size() << " faces.\n";

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
    cerr << "leapfrog::find_outer_face and compute faces\n";
    if(outer_face.size() < 5) outer_face = find_outer_face();

    vector<face_t> triangles(dualfrog.compute_faces_flat(3,false));

    cerr << "leapfrog::planar layout of " << triangles.size() << " triangles\n";
    new_layout.resize(triangles.size());

    for(int i=0;i<triangles.size();i++){
      const face_t& t(triangles[i]);
      new_layout[i] = t.centroid(dualfrog.layout2d)*coord2d(1,-1);
      if(t[0] == N || t[1] == N || t[2] == N) // triangle belongs to old outer face
	new_layout[i] *= 2.0;		      // TODO: Ensure that loop encompasses remaining graph
    }
  } 
  cerr << "leapfrog::dual()\n";
  FullereneGraph frog(dualfrog.dual_graph(3), new_layout);

  return frog;
}

node_t FullereneGraph::C20_edges[30][2] ={{0,13},{0,14},{0,15},{1,4},{1,5},{1,12},{2,6},{2,13},{2,18},{3,7},{3,14},{3,19},{4,10},{4,18},{5,11},{5,19},{6,10},{6,15},{7,11},{7,15},{8,9},{8,13},{8,16},{9,14},{9,17},{10,11},{12,16},{12,17},{16,18},{17,19}};

double FullereneGraph::C20_layout2d[20][2] = {{1.548,0.503},{-1.134,-0.368},{0.,1.628},{0.957,-1.317},{-1.548,0.503},{-0.957,-1.317},{0.,2.466},{1.449,-1.995},{0.302,0.416},{0.489,-0.159},{-2.345,0.762},{-1.45,-1.995},{-0.489,-0.159},{0.7,0.965},{1.133,-0.369},{2.345,0.762},{-0.302,0.416},{0.,-0.514},{-0.7,0.965},{0.,-1.192}};
