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

int main(int ac, char **av)
{
  int Nex = ac>=2? strtol(av[1],0,0) : 1;

  if(Nex < 1 || Nex > 4){cerr << "invalid graph chosen, aborting ..." << endl; abort();}

  string basename("gs-ex-"+to_string(Nex));

  ofstream output(("output/"+basename+".m").c_str());

  PlanarGraph g(examples[Nex-1]);
  output << "g = " << g << ";\n";
  output.flush();

  printf("Computing planar layout\n");
  g.layout2d = g.tutte_layout();

  printf("Computing dual\n");
  PlanarGraph dg(g.dual_graph(8,true));

  printf("Computing planar layour of dual\n");
  dg.layout2d = dg.tutte_layout();

  output << "g = " << g << ";\n";
  output << "dg = " << dg << ";\n";  
  output.flush();

  vector<face_t> faces = dg.compute_faces_flat(3,true);
  vector<tri_t>  triangles(faces.begin(),faces.end());
  cout << "triangles = " << triangles << ";\n";

  Triangulation dt(dg,triangles);
  cout << "t" << endl;
  vector<int> spiral;
  Triangulation::jumplist_t jumplist;
//  bool success = dt.get_spiral(spiral,jumplist,true);
//  if(!success) cerr << "Canonical general spiral not found.\n";
//
//  output << "spiral   = " << spiral << ";\n"
//  	 << "jumplist = " << jumplist << ";\n";
  

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
