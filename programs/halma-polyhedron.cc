#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/isomerdb.hh"

int testN = 80;
int testRSPI[12] = {1,2,3,4,5,6,37,38,39,40,41,42};

struct Deltahedron {
  Triangulation g;
  vector<coord3d> points;

  vector<face_t> faces(){
    const vector<tri_t> &tris = g.triangles;
    vector<face_t> faces(tris.size());
    for(int i=0;i<tris.size();i++)
      faces[i] = vector<node_t>{tris[i][0],tris[i][1],tris[i][1]};
    return faces;
  }

  void smooth(double q)
  {
    vector<coord3d> new_points(g.N);
    
    for(node_t u=0;u<g.N;u++){
      coord3d x = points[u]*(1-q);
      for(node_t v: g.neighbours[u]) x += points[v]*q/g.neighbours[u].size();
      new_points[u] = x;
    }
    points = new_points;
  }
};

vector<coord3d> operator*(const vector<coord3d>& xs, const double s)
{
  vector<coord3d> ys = xs;
  for(int i=0;i<ys.size();i++) ys[i] *= s;
  return ys;
}

vector<tri_t> operator+(const vector<tri_t>& xs, int k)
{
  vector<tri_t> ys = xs;
  for(int i=0;i<ys.size();i++)
    for(int j=0;j<3;j++) ys[i][j] += k;
  return ys;
}

Deltahedron halma_triangulation(const Deltahedron &P, int m=1) {
  // Create n new vertices for each edge
  vector<edge_t>                edges = P.g.undirected_edges();
  map< edge_t, vector<node_t> > edge_nodes;
  size_t n_vert = P.g.N;
  size_t n_tris = 2*(n_vert-2);
  size_t n_newtris = (m+1)*(m+1)*n_tris;
  size_t n_newvert = n_newtris/2+2;
  node_t v_new  = n_vert;

  
  Deltahedron Phalma{Graph(n_newvert),P.points*(m+1)};

  for(int i=0;i<edges.size();i++){
    const edge_t &e = edges[i];
    vector<node_t>& nodes(edge_nodes[edges[i]]);
    for(unsigned int i=0;i<m;i++) nodes.push_back(v_new++);

    for(unsigned int i=0;i<m;i++){
      double lambda = (1.0+i)*(1.0/(m+1));
      const coord3d &a(Phalma.points[e.first]), &b(Phalma.points[e.second]);
      Phalma.points.push_back(a*(1.0-lambda) + b*lambda);
    }
  }

  // For every triangle, we create and connect a halma-type grid
  const vector<tri_t> triangles = P.g.triangles;
  for(size_t i=0;i<triangles.size();i++){
    map<edge_t,node_t> grid;
    const tri_t& T(triangles[i]);
    edge_t e0(T[0],T[1]),e1(T[1],T[2]),e2(T[2],T[0]);
    const vector<node_t>& ns0(edge_nodes[e0]), ns1(edge_nodes[e1]), ns2(edge_nodes[e2]);

    // Insert original vertices
    grid[edge_t(0,0)]     = T[0]; // 0,0??
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

    double sqrt2inv = 1.0/sqrt(2.0);
    const coord3d &a(Phalma.points[T[0]]), &b(Phalma.points[T[1]]),
                  &c(Phalma.points[T[2]]);      
    for(int j=1;j<m;j++)
      for(int k=j+1;k<=m;k++){
	double s = (1+j)*(1.0/(m+2)), t = k*(1.0/(m+2));
	//fprintf(stderr,"(s,t) = (%g,%g)\n",s,t);
	Phalma.points.push_back(a+((b-a)*s + (c-a)*t)*sqrt2inv);
      }
      
    // Connect the vertices in the grid
    for(int j=0;j<=m;j++)
      for(int k=j+1;k<=m+1;k++){
	node_t v(grid[edge_t(j,k)]), down(grid[edge_t(j+1,k)]), 
	  left(grid[edge_t(j,k-1)]);

	Phalma.g.insert_edge(arc_t(v,down));
	Phalma.g.insert_edge(arc_t(down,left));
	Phalma.g.insert_edge(arc_t(left,v));
      }
  }

  Phalma.g.update(false);
  return Phalma;
}



int main(int ac, char **av)
{
  int N;
  jumplist_t jumps;
  vector<int> RSPI(12);
  bool from_file = false;

  if(ac==2){			// If only one argument is given, it is a filename to read initial polyhedron from
    from_file = true;
    N = 0;
    
  } else if(ac==3) {            // Pentagon indices in quoted string
    N = strtol(av[1],0,0);
    int n = strlen(av[2]);
    for(int i=0;i<n;i++) if(av[2][i] == ',') av[2][i] = ' ';
    istringstream iss(av[2]);
    for(int i=0;i<12;i++){ iss >> RSPI[i]; RSPI[i]--; }
    cout << "RSPI = " << RSPI << endl;
  } else 
    if(ac>=14) {
      N = strtol(av[1],0,0);
      for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
    } 

  // int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  // bool IPR = ac>=4? strtol(av[3],0,0) : false;

  // IsomerDB DB(N,IPR);

  // FullereneGraph g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));


  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;

  //  cerr << "rspi="<<RSPI<<"\n";
  
  Polyhedron P0;
  PlanarGraph g;
  if(from_file){
    P0 = Polyhedron::from_file(av[1]);
    N = P0.N;
    g = P0;
    g.layout2d = g.tutte_layout(-1,-1,-1,8);
  } else {
    Triangulation T(spiral,jumps);
    
    g = T.dual_graph();
    g.layout2d = g.tutte_layout();
    P0 = Polyhedron(g,g.zero_order_geometry(),6);
  }

  string basename("polyhedron-"+to_string(N));
  Polyhedron::to_file(P0,"output/"+basename+"-P0.mol2");

  Polyhedron P(P0);
  P.optimise();

  Polyhedron::to_file(P,"output/"+basename+".mol2");

  {
    P.move_to_origin();
    matrix3d If(P.principal_axes());
    P.points = If*P.points;

    Polyhedron::to_file(P,"output/"+basename+"-if.mol2");


    // ofstream pov(("output/"+basename+"-if.pov").c_str());
    // pov << P.to_povray();
    // pov.close();
  }

  // ofstream output(("output/"+basename+".m").c_str());

  // facemap_t facemap(g.compute_faces(8,true));
  // output << "g = " << g << ";\n";
  // output << "coordinates0 = " << P0.points << ";\n";
  // output << "coordinates = "  << P.points << ";\n";
  // output << "pentagons = " << facemap[5] << ";\n"
  // 	  << "hexagons  = " << facemap[6] << ";\n"
  // 	  << "RSPI = " << RSPI << ";\n";

  // output << "P0 = " << P0 << ";\n";
  // output << "P = " << P << ";\n";

  Polyhedron D(P.dual());
  // D.layout2d = D.tutte_layout();
  // D.faces    = D.compute_faces_flat(3,true);
  // D.face_max = 3;
  // //   D.optimise();
  // output << "PD = " << D << ";\n";

  Deltahedron DD{D,D.points};
  DD.smooth(.8);
  Deltahedron HD(halma_triangulation(DD,1)); // Might as well roll all this into polyhedron (but assert graph is triangulation)/
  HD = halma_triangulation(HD,1); // Might as well roll all this into polyhedron (but assert graph is triangulation)/
  HD.smooth(.8);
  HD.smooth(.8);
  Deltahedron HD2(halma_triangulation(HD,1)); // Might as well roll all this into polyhedron (but assert graph is triangulation)/
  HD2.smooth(.8);
  HD2.smooth(.8);

  Polyhedron  PHD(HD.g,HD.points,3,HD.faces());  
  //  Polyhedron  PHD2(HD2.g,HD2.points,3,HD2.faces());
  Polyhedron  PH(PHD.dual());

  //  Polyhedron  PH2(PHD2.dual(3,false));
  PH.optimise();
  //  PH2.optimise();
  //  PHD.optimise();
		  
  // output << "PHD  = " << PHD << ";\n";
  // output << "PHD2 = " << PHD2 << ";\n";
  // output << "PH  = " << PH << ";\n";
  // output << "DDg = " << DD.g << ";\n"
  // 	 << "DDpoints = " << DD.points << ";\n"
  // 	 << "DDfaces  = " << (DD.g.triangles+1) << ";\n";

  // output << "HDg = " << HD.g << ";\n"
  // 	 << "HDpoints = " << HD.points << ";\n"
  // 	 << "HDfaces  = " << (HD.g.triangles+1) << ";\n";
  
  // output << "HD2g = " << HD2.g << ";\n"
  // 	 << "HD2points = " << HD2.points << ";\n"
  // 	 << "HD2faces  = " << (HD2.g.triangles+1) << ";\n";

  //  output.close();
  

  Polyhedron::to_file(D,"output/"+basename+"-dual.mol2");
  //    Polyhedron::to_file(D,"output/"+basename+"-dual.pov");

  Polyhedron::to_file(PH,"output/"+basename+"-halma.mol2");
  // {
  //   ofstream mol2(("output/"+basename+"-halma2.mol2").c_str());
  //   mol2 << PH2.to_mol2();
  //   mol2.close() ;
  // }

  Polyhedron::to_file(PHD,"output/"+basename+"-dual-halma.mol2");
  // {
  //   ofstream mol2(("output/"+basename+"-dual-halma2.mol2").c_str());
  //   mol2 << PHD2.to_mol2();
  //   mol2.close() ;
  // }
  return 0;
}
