#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"
#include "libgraph/triangulation.hh"

int testN = 80;
int testRSPI[12] = {1,2,3,4,5,6,37,38,39,40,41,42};

// Smallest planar cubic graph with no face spiral
Graph example1()
{
  int M = 6, N = 18;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0;i<3;i++){
    // Outer triangle
    neighbours[i*M+0][0] = i*M+1;
    neighbours[i*M+0][1] = i*M+2;
    neighbours[i*M+0][2] = ((i-1)*M+1+N)%N;
    
    neighbours[i*M+1][0] = i*M+0;
    neighbours[i*M+1][1] = i*M+2;
    neighbours[i*M+1][2] = ((i+1)*M+0)%N;

    neighbours[i*M+2][0] = i*M+0;
    neighbours[i*M+2][1] = i*M+1;
    neighbours[i*M+2][2] = i*M+3;

    // Inner triangle
    neighbours[i*M+3][0] = i*M+4;
    neighbours[i*M+3][1] = i*M+5;
    neighbours[i*M+3][2] = i*M+2;
    
    neighbours[i*M+4][0] = i*M+3;
    neighbours[i*M+4][1] = i*M+5;
    neighbours[i*M+4][2] = ((i+1)*M+5)%N;

    neighbours[i*M+5][0] = i*M+3;
    neighbours[i*M+5][1] = i*M+4;
    neighbours[i*M+5][2] = ((i-1)*M+4+N)%N;
  }
  return Graph(neighbours);
}


// Tutte graph
Graph example2(){
  const int M=3, N=46;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0; i!=3; ++i){
    neighbours[i][0] = 45;
    neighbours[i][1] = 1*M+i;
    neighbours[i][2] = 2*M+i;
    
    neighbours[1*M+i][0] =     i;
    neighbours[1*M+i][1] = 3*M+i;
    neighbours[1*M+i][2] = 4*M+i;
    
    neighbours[2*M+i][0] =     i;
    neighbours[2*M+i][1] = 3*M+i;
    neighbours[2*M+i][2] = 6*M+i;
    
    neighbours[3*M+i][0] = 1*M+i;
    neighbours[3*M+i][1] = 2*M+i;
    neighbours[3*M+i][2] = 5*M+i;
    
    neighbours[4*M+i][0] = 1*M+i;
    neighbours[4*M+i][1] = 7*M+i;
    neighbours[4*M+i][2] = 8*M+i;
    
    neighbours[5*M+i][0] = 3*M+i;
    neighbours[5*M+i][1] = 8*M+i;
    neighbours[5*M+i][2] = 9*M+i;
    
    neighbours[6*M+i][0] = 2*M+i;
    neighbours[6*M+i][1] = 9*M+i;
    neighbours[6*M+i][2] = 10*M+i;
    
    neighbours[7*M+i][0] = 4*M+i;
    neighbours[7*M+i][1] = 11*M+i;
    neighbours[7*M+i][2] = 12*M+i;
    
    neighbours[8*M+i][0] = 4*M+i;
    neighbours[8*M+i][1] = 5*M+i;
    neighbours[8*M+i][2] = 12*M+i;
    
    neighbours[9*M+i][0] = 5*M+i;
    neighbours[9*M+i][1] = 6*M+i;
    neighbours[9*M+i][2] = 13*M+i;
    
    neighbours[10*M+i][0] = 6*M+i;
    neighbours[10*M+i][1] = 11*M+((i+1+3)%3);
    neighbours[10*M+i][2] = 14*M+i;
    
    neighbours[11*M+i][0] = 7*M+i;
    neighbours[11*M+i][1] = 10*M+((i-1+3)%3);
    neighbours[11*M+i][2] = 14*M+i;
    
    neighbours[12*M+i][0] = 7*M+i;
    neighbours[12*M+i][1] = 8*M+i;
    neighbours[12*M+i][2] = 13*M+i;
    
    neighbours[13*M+i][0] = 9*M+i;
    neighbours[13*M+i][1] = 12*M+i;
    neighbours[13*M+i][2] = 14*M+i;
    
    neighbours[14*M+i][0] = 10*M+i;
    neighbours[14*M+i][1] = 11*M+i;
    neighbours[14*M+i][2] = 13*M+i;
  }

  neighbours[45][0] = 0;
  neighbours[45][1] = 1;
  neighbours[45][2] = 2;

  return Graph(neighbours);
}


// smallest non-spiral graph with face sizes less/equal 6
Graph example3(){
  const int M=3, N=36;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0; i!=3; ++i){
    neighbours[i][0] = (i+1)%M;
    neighbours[i][1] = (i-1+M)%M;
    neighbours[i][2] = M+i;
    
    neighbours[1*M+i][0] =     i;
    neighbours[1*M+i][1] = 2*M+i;
    neighbours[1*M+i][2] = 3*M+i;
    
    neighbours[2*M+i][0] = 1*M+i;
    neighbours[2*M+i][1] = 3*M+(i+1)%M;
    neighbours[2*M+i][2] = 4*M+i;
    
    neighbours[3*M+i][0] = 1*M+i;
    neighbours[3*M+i][1] = 2*M+(i-1+M)%M;
    neighbours[3*M+i][2] = 5*M+i;
    
    neighbours[4*M+i][0] = 2*M+i;
    neighbours[4*M+i][1] = 6*M+i;
    neighbours[4*M+i][2] = 7*M+i;
  
    neighbours[5*M+i][0] = 3*M+i;
    neighbours[5*M+i][1] = 6*M+i;
    neighbours[5*M+i][2] = 8*M+(i-1+M)%M;
    
    neighbours[6*M+i][0] = 4*M+i;
    neighbours[6*M+i][1] = 5*M+i;
    neighbours[6*M+i][2] = 10*M+i;
    
    neighbours[7*M+i][0] = 4*M+i;
    neighbours[7*M+i][1] = 8*M+i;
    neighbours[7*M+i][2] = 9*M+i;
  
    neighbours[8*M+i][0] = 5*M+(i+1)%M;
    neighbours[8*M+i][1] = 7*M+i;
    neighbours[8*M+i][2] = 9*M+i;
    
    neighbours[9*M+i][0] = 7*M+i;
    neighbours[9*M+i][1] = 8*M+i;
    neighbours[9*M+i][2] = 11*M+i;
    
    neighbours[10*M+i][0] = 6*M+i;
    neighbours[10*M+i][1] = 11*M+i;
    neighbours[10*M+i][2] = 11*M+(i-1+M)%M;
    
    neighbours[11*M+i][0] = 9*M+i;
    neighbours[11*M+i][1] = 10*M+i;
    neighbours[11*M+i][2] = 10*M+(i+1)%M;
  }

  return Graph(neighbours);
}


// smallest (?) polyhedron with only pentagons and heptagons
Graph example4(){
  const int M=7, N=28;
  neighbours_t neighbours(N,vector<node_t>(3));

  for(int i=0; i!=7; ++i){
    neighbours[i][0] = (i+1)%M;
    neighbours[i][1] = (i-1+M)%M;
    neighbours[i][2] = M+i;
    
    neighbours[1*M+i][0] =     i;
    neighbours[1*M+i][1] = 2*M+i;
    neighbours[1*M+i][2] = 2*M+(i-1+M)%M;
    
    neighbours[2*M+i][0] = 1*M+i%M;
    neighbours[2*M+i][1] = 1*M+(i+1)%M;
    neighbours[2*M+i][2] = 3*M+i;
    
    neighbours[3*M+i][0] = 2*M+i;
    neighbours[3*M+i][1] = 3*M+(i+1)%M;
    neighbours[3*M+i][2] = 3*M+(i-1+M)%M;
  }

  return Graph(neighbours);
}


Graph examples[4] = {example1(), example2(), example3(), example4()};

/* 
 Ex. 2 points from mathematica: {{1.58369, 4.85689, 1.94701}, {3.77283, 3.40739, 1.94783}, {1.40524, 
  2.27093, 1.94129}, {0.802452, 5.67691, 1.74995}, {4.88999, 3.67391, 
  1.74872}, {1.0775, 1.18288, 1.74496}, {1.53685, 5.81635, 
  1.58434}, {4.64367, 2.96501, 1.5879}, {0.593733, 1.74472, 
  1.5653}, {0.937614, 6.38161, 1.64838}, {5.44231, 3.1966, 
  1.65656}, {0.393613, 0.94514, 1.63707}, {0.10571, 5.92323, 
  1.09248}, {5.47139, 4.14713, 1.091}, {1.22062, 0.44689, 
  1.10692}, {0.628643, 6.85614, 1.06865}, {6.02303, 3.23402, 
  1.08254}, {0.145348, 0.418474, 1.07054}, {1.88329, 5.99409, 
  0.726333}, {4.64954, 2.57794, 0.730885}, {0.274608, 1.92851, 
  0.706214}, {0.0907154, 5.36617, 0.383971}, {5.00948, 4.44235, 
  0.377893}, {1.72661, 0.69798, 0.396817}, {0., 6.61621, 
  0.798106}, {6.1205, 3.90155, 0.807663}, {0.689221, 0., 
  0.824317}, {1.1756, 6.57309, 0.547627}, {5.50958, 2.89208, 
  0.572887}, {0.127775, 1.02774, 0.542233}, {2.5427, 5.26018, 
  0.194485}, {3.6849, 2.39241, 0.194494}, {0.606153, 2.85049, 
  0.183424}, {0.579674, 4.33679, 0.119217}, {3.86252, 4.54809, 
  0.119025}, {2.39719, 1.62552, 0.123872}, {0.0799958, 6.17694, 
  0.238099}, {5.71096, 4.04066, 0.240417}, {1.03746, 0.275945, 
  0.259574}, {0.844957, 6.07888, 0.0561757}, {5.25388, 3.42843, 
  0.065193}, {0.739878, 0.976791, 0.0532776}, {1.34965, 5.21765, 
  0.00712336}, {4.25376, 3.44237, 0.00467394}, {1.24434, 1.83824, 
  0.}, {2.25917, 3.51055, 2.03661}}
 */ 

int main(int ac, char **av)
{
  int Nex = ac>=2? strtol(av[1],0,0) : 1;

  vector<vector<int>> outer_faces{{1, 3, 4, 6, 17, 16, 15, 14},{1,4,13,22,34,33,21,9,3,46},{},{}};
  
  if(Nex < 1 || Nex > 4){cerr << "invalid graph chosen, aborting ..." << endl; abort();}

  string basename("gs-ex-"+to_string(Nex));

  ofstream output(("output/"+basename+".m").c_str());

  PlanarGraph g(examples[Nex-1]);
  output << "g = " << g << ";\n";
  output.flush();

  printf("Computing planar layout\n");
  g.outer_face = outer_faces[Nex-1]+(-1);
  g.layout2d = g.tutte_layout(g.outer_face);

  printf("Computing dual\n");
  PlanarGraph dg(g.dual_graph(8,true));

  printf("Computing planar layour of dual\n");
  dg.layout2d = dg.tutte_layout();

  output << "g2d  = " << g << ";\n";
  output << "dg2d = " << dg << ";\n";  
  output.flush();

  vector<face_t>  faces =  g.compute_faces(8,true);
  vector<face_t> dfaces = dg.compute_faces(3,true);

  output << "faces  = " << faces  << ";\n"
	 << "dfaces = " << dfaces << ";\n";
  output.flush();
  
  bool dual_is_triangulation = true, graph_is_cubic = true;
  for(auto &f: dfaces) if(f.size() != 3) dual_is_triangulation = false;
  for(auto &n: g.neighbours) if(n.size() != 3) graph_is_cubic = false;

  if(graph_is_cubic != dual_is_triangulation){
    fprintf(stderr,"Graph is %scubic, but dual is %striangulation.\n",
	    graph_is_cubic?"":"not ", dual_is_triangulation?"":"not ");
    abort();
  }
  
  vector<tri_t>  triangles(dfaces.begin(),dfaces.end());
  cout << "triangles = " << triangles << ";\n";

  Triangulation dt(dg,triangles);

  vector<int> spiral;
  jumplist_t jumplist;
  bool success = dt.get_spiral(spiral,jumplist);
  if(!success) cerr << "Canonical general spiral not found.\n";
//
  output << "spiral   = " << spiral << ";\n"
  	 << "jumplist = " << jumplist << ";\n";
  

  Polyhedron P0(g,g.zero_order_geometry(),10);
  Polyhedron dP0(dg,dg.zero_order_geometry(),3);

  output << "P0coordinates = " << P0.points << ";\n";

  {
    ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
    mol2 << P0.to_mol2();
    mol2.close();
  }

  {
    ofstream mol2(("output/"+basename+"-dP0.mol2").c_str());
    mol2 << dP0.to_mol2();
    mol2.close();
  }


  
  Polyhedron P(P0);

  printf("Optimizing polyhedron.\n");
  P.optimize();

  output << "Pcoordinates = "  << P.points << ";\n";
  output << "P = " << P << ";\n";
  
  Polyhedron D(P.dual(8,true));

  output << "PD = " << D << ";\n";
  
  output.close();


  {
    ofstream mol2(("output/"+basename+".mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();

    ofstream pov(("output/"+basename+".pov").c_str());
    pov << P.to_povray();
    pov.close();    
  }

  
  
  {
    ofstream xyz(("output/"+basename+".xyz").c_str());
    xyz << P.to_xyz();
    xyz.close();
  }

  {
    ofstream latex(("output/"+basename+".tex").c_str());
    latex << PlanarGraph(P).to_latex(10, 10, false, false, true);
    latex.close();
  }

  {
    ofstream mol2(("output/"+basename+"-dual.mol2").c_str());
    mol2 << D.to_mol2();
    mol2.close();
  }

  {
    ofstream latex(("output/"+basename+"-dual.tex").c_str());
    latex << PlanarGraph(D).to_latex(10, 10, false, false, true);
    latex.close();
  }

  return 0;
}
