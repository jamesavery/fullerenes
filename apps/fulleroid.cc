#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

// int N_default = 60;
// int rspi_default[12] = {1,7,9,11,13,15,18,20,22,24,26,32};

// int N_default = 80;
// int rspi_default[12] = {1,8,10,12,14,16,28,30,32,34,36,42};

int N_default = 140;
int rspi_default[12] = {1,17,20,23,26,29,48,51,54,57,60,72};

vector<int> spiral = {5,5,5,5,5,5,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,7,5,5,5,7,7,7,5,5,5,7,7,7,5,5,5,7,7,7,5,5,5,7,7,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,5,5,7,7,5,7,7,5,5,5,5,5,5};


Triangulation from_rspi(int F, const vector<int>& rspi, const PlanarGraph::jumplist_t& jumps = PlanarGraph::jumplist_t())
{
  vector<int> spiral(F,6);
  for(int i=0;i<12;i++) spiral[rspi[i]] = 5;
  return Triangulation(spiral,jumps);
}

bool extrude(Triangulation &T, int p)
{
  fprintf(stderr,"Supposedly degree 5 node %d has %d neighbours\n",p,int(T.neighbours[p].size()));
  // Insert 5 new pentagons p0,...,p4 surrounding p

  int Hi[5], v2[5], v3[5], v4[5];
  for(int i=0;i<5;i++){
    Hi[i] = T.neighbours[p][i];
    if(T.neighbours[Hi[i]].size() != 6) return false; // p is not valid patch replacement site.

    const vector<node_t> &nHi(T.neighbours[Hi[i]]);

    int ip = 0; 
    for(;ip<nHi.size();ip++) if(nHi[ip] == p) break;

    v2[i] = nHi[(ip+2)%6];	// Assuming p's neighbours are all hexagons
    v3[i] = nHi[(ip+3)%6];
    v4[i] = nHi[(ip+4)%6];
  }

  for(int i=0;i<5;i++) {
    int pi = T.N + i, pim1 = T.N + ((i+4) % 5), pip1 = T.N + ((i+1)%5);

    T.neighbours.push_back(vector<int>({p,pim1,Hi[i],Hi[(i+1)%5],pip1})); // Neighbour list of pi
    T.neighbours[Hi[i]] = vector<int>({pi,pim1,Hi[(i+4)%5],v2[i],v3[i],v4[i],Hi[(i+1)%5]});  // Neighbour list of Hi
  }
  T.neighbours[p] = vector<node_t>({T.N,T.N+1,T.N+2,T.N+3,T.N+4}); // neighbours(p) = {p0,...,p4}
  T.N += 5;

  return true;
}


bool extrude(Polyhedron& P, int p) // Requires that P is a triangulation.
{
  assert(P.is_triangulation());

  int Hi[5], v2[5], v3[5], v4[5], pi[5];
  for(int i=0;i<5;i++){
    pi[i] = P.N+i;
    Hi[i] = P.neighbours[p][i];
    if(P.neighbours[Hi[i]].size() != 6) return false; // p is not valid patch replacement site.

    const vector<node_t> &nHi(P.neighbours[Hi[i]]);
    int ip = 0; 
    for(;ip<nHi.size();ip++) if(nHi[ip] == p) break;

    v2[i] = nHi[(ip+2)%6];	// Assuming p's neighbours are all hexagons
    v3[i] = nHi[(ip+3)%6];
    v4[i] = nHi[(ip+4)%6];
  }

  // Derive new planar layout 
  if(P.layout2d.size() == P.N){
    P.layout2d.resize(P.N+5);
    
    // coord2d cm;
    // for(int i=0;i<5;i++) cm += P.layout2d[Hi[i]]/5.;
    // P.layout2d[p] = cm;

    for(int i=0;i<5;i++) P.layout2d[pi[i]] = (P.layout2d[Hi[i]]+P.layout2d[Hi[(i+1)%5]] + P.layout2d[p])/3.0;
  }

  // Derive new approximate geometry
  if(P.points.size() == P.N){
    P.points.resize(P.N+5);
    
    // coord3d cm;
    // for(int i=0;i<5;i++) cm += P.points[Hi[i]]/5.;
    // P.points[p] = cm;

    for(int i=0;i<5;i++) P.points[pi[i]] = (P.points[Hi[i]]+P.points[Hi[(i+1)%5]] + P.points[p])/3.0;
  }

  // Update graph
  for(int i=0;i<5;i++) {
    P.neighbours.push_back(vector<int>({p,pi[(i+4)%5],Hi[i],Hi[(i+1)%5],pi[(i+1)%5]})); // Neighbour list of pi
    P.neighbours[Hi[i]] = vector<int>({pi[i],pi[(i+4)%5],Hi[(i+4)%5],v2[i],v3[i],v4[i],Hi[(i+1)%5]});  // Neighbour list of Hi
  }
  P.neighbours[p] = vector<node_t>({pi[0],pi[1],pi[2],pi[3],pi[4]}); 
  P.N += 5;

  P.update_from_neighbours();

  // TODO: Fix derived planar embedding instead
  //       My suspicion is that the problem arises from the outer face.
  P.layout2d = P.tutte_layout();

  // Update faces
  // TODO: Replace triangles directly instead of recomputing?
  P.faces    = P.compute_faces_flat(3,true);
  P.face_max = 3;


  return true;
}



int main(int ac, char **av)
{
  vector<int> rspi(12);
  int N = ac>1? strtol(av[1],0,0) : N_default;
  if(ac>12)
    for(int i=1;i<=12;i++) 
      rspi[i-1] = strtol(av[i+1],0,0)-1;
  else 
    for(int i=0;i<12;i++) 
      rspi[i] = rspi_default[i]-1;

  string basename("fulleroid-"+to_string(N));
  
  int F = N/2+2;
  Triangulation dual(from_rspi(F,rspi));
  dual.layout2d = dual.tutte_layout();

  Polyhedron PD(dual,dual.zero_order_geometry()*coord3d(1.4,1.4,1.4));

  {
    ofstream mol2(("output/"+basename+"-step0.mol2").c_str());
    mol2 << PD.to_mol2();
    mol2.close();
  }

  PD.optimize(4);

  {
    ofstream mol2(("output/"+basename+"-step0.5.mol2").c_str());
    mol2 << PD.to_mol2();
    mol2.close();
  }


  for(int p=0;p<12;p++)
    if(extrude(PD,rspi[p])){
      cerr << "Completed extrusion at pentagon " << p << endl;
      PD.update_from_neighbours();
      // {
      // 	ofstream mol2(("output/"+basename+"-step"+to_string(p+1)+".mol2").c_str());
      // 	mol2 << PD.to_mol2();
      // 	mol2.close();
      // }
      // PD.optimize(4);
      // PD.optimize(3);
      // {
      // 	ofstream mol2(("output/"+basename+"-step"+to_string(p+1)+".5.mol2").c_str());
      // 	mol2 << PD.to_mol2();
      // 	mol2.close();
      // }

    } else {
      cerr << "No extrusion patch replacement site at pentagon " << p << endl;
    }

  {
    ofstream mol2(("output/"+basename+"-step12.mol2").c_str());
    mol2 << PD.to_mol2();
    mol2.close();
  }
  PD.optimize(4);
  PD.optimize(3);
  {
    ofstream mol2(("output/"+basename+"-step12.5.mol2").c_str());
    mol2 << PD.to_mol2();
    mol2.close();
  }

  // PD.move_to_origin();
  // matrix3d If(PD.inertial_frame());
  // PD.points = If*PD.points;

  cout << "PD = " << PD << ";\n";
  cout << "dual = " << static_cast<PlanarGraph>(PD) << ";\n";

  assert(PD.is_triangulation());
  //  PD.orient_neighbours();
  //  PD.layout2d = PD.tutte_layout();

  fprintf(stderr,"Getting fulleroid as dual of dual, B.\n");

  Triangulation T(PD,true);
  vector<int> spiral;
  Triangulation::jumplist_t jumps;
  T.get_spiral(spiral,jumps,false,true,true);

  
  Polyhedron P  = PD.dual(3,true);

  cout << "G = " << static_cast<PlanarGraph>(PD) << ";\n";
  cout << "spiral = " << spiral << ";\n"
       << "jumps  = " << jumps << ";\n";
  {
    ofstream mol2(("output/"+basename+"-step-preopt.mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
  }

  P.optimize(4);
  P.optimize(3);

  
  
  {
    ofstream mol2(("output/"+basename+".mol2").c_str());
    mol2 << P.to_mol2();
    mol2.close();
     
    ofstream pov(("output/"+basename+".pov").c_str());
    pov << P.to_povray();
    pov.close();

    ofstream math(("output/"+basename+".m").c_str());
    math << "P = " << P;
    math.close();
  }   

  return 0;
}
