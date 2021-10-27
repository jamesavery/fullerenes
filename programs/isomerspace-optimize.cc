#include <limits.h>
#include <chrono>
#include <iostream>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

using namespace std;
using namespace std::chrono;

#include "fullerenes/gpu/isomerspace_forcefield.hh"

int face_size(const Graph &g, node_t u, node_t v)
{
  int d = 1;
  node_t u0 = u;
  while(v != u0){
    node_t w = v;
    v = g.next_on_face(u,v);
    u = w;
    d++;
  }
  return d;
}


int main(int ac, char **argv)
{
  
  if(ac<2){
    fprintf(stderr,"Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n",argv[0]);
    return -1;
  }
  int N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N

  string output_dir   = ac>=3? argv[2] : "output";    // Argument 2: directory to output files to
  int IPR             = ac>=4? strtol(argv[3],0,0):0; // Argument 3: Only generate IPR fullerenes?
  int only_nontrivial = ac>=5? strtol(argv[4],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?
  
  ofstream failures((output_dir+"/failures.txt").c_str()); // output/failures.txt contains list of any fullerenes that failed optimization

  size_t batch_size = IsomerspaceForcefield::get_batch_capacity(N);
  typedef IsomerspaceForcefield::device_real_t device_real_t;
  typedef IsomerspaceForcefield::device_node_t device_node_t;


  device_node_t   cubic_graph[batch_size*3*N], next_on_face[batch_size*3*N], prev_on_face[batch_size*3*N];
  uint8_t         face_right[batch_size*3*N]; // TODO: Reduce to 1 bit/arc 
  device_real_t            X[batch_size*3*N];
  device_real_t   bonds[batch_size*3*N], angles[batch_size*3*N], dihedrals[batch_size*3*N], bond_0[batch_size*3*N], angle_0[batch_size*3*N], dihedral_0[batch_size*3*N], gradients[batch_size*3*N];

  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial);  

  FullereneDual dualG;
  FullereneGraph G;
  G.N = N;
  G.neighbours = vector<vector<node_t>>(N,vector<node_t>(3));
  vector<coord3d> points(N);
  
  size_t I=0,			// Global isomer number at start of batch
         i=0;			// Isomer number within current batch
  bool more_to_do = true;
  auto T0 = system_clock::now();
  auto
    Tgen    = system_clock::now()-T0,
    Tupdate = system_clock::now()-T0,
    Tdual   = system_clock::now()-T0,    
    Ttutte  = system_clock::now()-T0,
    TX0     = system_clock::now()-T0,
    Tcopy   = system_clock::now()-T0,
    Topt    = system_clock::now()-T0,
    Tcheck  = system_clock::now()-T0;


  IsomerspaceForcefield isomer_batch = IsomerspaceForcefield(N);
  vector<Polyhedron> batch_polyhedra(batch_size);

  while(more_to_do){
    // Fill in a batch
    printf("Generating isomer bond graphs and initial geometries\n");
    for(i=0; (i<batch_size) && more_to_do; i++){
      //      printf("i=%ld, I=%ld, isomer_numer=%ld\n",i,I,I+i);
      auto t0 = system_clock::now();            
      more_to_do &= BuckyGen::next_fullerene(Q,dualG);
      if(!more_to_do) break;

      auto t1= system_clock::now(); Tgen += t1-t0;

      dualG.update();   		        // Update triangles
      auto t2= system_clock::now(); Tupdate += t2-t1;
      
      FullereneGraph   G = dualG.dual_graph();  // Construct fullerene graph
      auto t3= system_clock::now(); Tdual += t3-t2;
      G.layout2d         = G.tutte_layout();
      auto t4= system_clock::now(); Ttutte += t4-t3;
      vector<coord3d> X0 = G.zero_order_geometry(); // TODO: Faster, better X0
      auto t5= system_clock::now(); TX0    += t5-t4;
      isomer_batch.insert_isomer(G,X0);
      auto t6= system_clock::now(); Tcopy += t6-t5;
      Polyhedron P0(G,X0);
      string filename = output_dir+"/P0-C"+to_string(N)+"-"+to_string(I+i);
      Polyhedron::to_file(P0,filename+".mol2");   
      batch_polyhedra[i] = P0;   
    }

    size_t this_batch_size = i;
    printf("Optimizing %ld C%d fullerenes, isomer [%ld;%ld]\n",this_batch_size,N,I,I+this_batch_size-1);
    auto t0 = system_clock::now();
    isomer_batch.optimize_batch(N*5);
    Topt += system_clock::now()-t0;
    auto t6 = system_clock::now();
    isomer_batch.check_batch();
    Tcheck += system_clock::now()-t6;
    isomer_batch.get_cartesian_coordinates(X);
    isomer_batch.batch_statistics_to_file();
    isomer_batch.clear_batch();

    // Now do something with the optimized geometries

    for(size_t i=0;i<this_batch_size;i++){
      for(node_t u=0;u<N;u++){
      	size_t iu_idx = i*N*3 + u*3;
        batch_polyhedra[i].points[u] = {X[iu_idx], X[iu_idx+1], X[iu_idx+2]};
      }

      // Output molecular geometry files    
      string filename = output_dir+"/P-C"+to_string(N)+"-"+to_string(I+i);
      Polyhedron::to_file(batch_polyhedra[i],filename+".mol2");
    }

    // Output molecular geometry files
    
    I += this_batch_size;
    if (I > 10000){break;}
  }
  failures.close();
  cout << "Time spent on non:\n"
    "\tGenerating graphs = " << (Tgen/1ms)    << " ms\n"
    "\tUpdating metadata = " << (Tupdate/1ms) << " ms\n"
    "\tDualizing         = " << (Tdual/1ms)   << " ms\n"
    "\tTutte embedding   = " << (Ttutte/1ms)  << " ms\n"
    "\tInitial geometry  = " << (TX0/1ms)     << " ms\n"
    "\tCopying to buffer = " << (Tcopy/1ms)   << " ms\n"
    "\tFF Optimization   = " << (Topt/1ms)    << " ms\n"
    "\tFF Convergence Check   = " << (Tcheck/1ms)    << " ms\n";
  
  return 0;
}
