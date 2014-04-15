#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

// int N_default = 60;
// int rspi_default[12] = {1,7,9,11,13,15,18,20,22,24,26,32};

// int N_default = 80;
// int rspi_default[12] = {1,8,10,12,14,16,28,30,32,34,36,42};

int N_default = 140;
int rspi_default[12] = {1,17,20,23,26,29,48,51,54,57,60,72};



Triangulation from_rspi(int F, const vector<int>& rspi, const PlanarGraph::jumplist_t& jumps = PlanarGraph::jumplist_t())
{
  vector<int> spiral(F,6);
  for(int i=0;i<12;i++) spiral[rspi[i]] = 5;
  return Triangulation(spiral,jumps);
}

void extrude(Triangulation &T, int p)
{
  fprintf(stderr,"Supposedly degree 5 node %d has %d neighbours\n",p,int(T.neighbours[p].size()));
  // Insert 5 new pentagons p0,...,p4 surrounding p

  int Hi[5], v2[5], v3[5], v4[5];
  for(int i=0;i<5;i++){
    Hi[i] = T.neighbours[p][i];
    assert(T.neighbours[Hi[i]].size() == 6);

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
  
  int F = N/2+2;
  Triangulation dual(from_rspi(F,rspi));

  for(int i=0;i<12;i++) extrude(dual,rspi[i]);
  //  fprintf(stderr,"Edges before: %d\n",int(dual.edge_set.size()));
  //  extrude(dual,rspi[1]);
  dual.update_from_neighbours();
  //  fprintf(stderr,"Edges after: %d\n",int(dual.edge_set.size()));
  cout << "dual = " << dual << ";\n";

   fprintf(stderr,"-2\n");
   cout << dual.neighbours[0] << endl;
   //   dual.orient_neighbours();
   //   dual.update(true);
   fprintf(stderr,"-1\n");

   
   dual.layout2d = dual.tutte_layout();
   fprintf(stderr,"0\n");
   //   Polyhedron PD(dual,dual.zero_order_geometry());
   //   PD.optimize();

   string basename("fulleroid-"+to_string(N));
   // {
   //   ofstream mol2(("output/"+basename+"-dual.mol2").c_str());
   //   mol2 << PD.to_mol2();
   //   mol2.close();
   // }
   fprintf(stderr,"A\n");
   cout << "dual = " << dual << ";\n";
   PlanarGraph G  = PlanarGraph(dual,dual.layout2d).dual_graph(3,true);
   fprintf(stderr,"B\n");
   G.layout2d = G.tutte_layout();
   fprintf(stderr,"C\n");
   
   cout << "graph = " << G << ";\n";
   fprintf(stderr,"D\n");

   Polyhedron P(G,G.zero_order_geometry());
   fprintf(stderr,"E\n");
   {
     ofstream mol2(("output/"+basename+"-P0.mol2").c_str());
     mol2 << P.to_mol2();
     mol2.close();
   }
   fprintf(stderr,"F\n");
   P.optimize(4);

   P.move_to_origin();
   matrix3d If(P.inertial_frame());
   P.points = If*P.points;


   fprintf(stderr,"G\n");

   P.faces    = P.compute_faces_flat(7,true);
   P.face_max = 7;

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
