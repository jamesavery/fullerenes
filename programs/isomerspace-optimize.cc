#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>

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

  size_t batch_size = IsomerspaceForcefield::computeBatchSize(N);
  using IsomerspaceForcefield::device_real_t;
  using IsomerspaceForcefield::device_node_t;

  device_node_t   cubic_graph[batch_size*3*N], next_on_face[batch_size*3*N], prev_on_face[batch_size*3*N];
  uint8_t         face_right[batch_size*3*N]; // TODO: Reduce to 1 bit/arc 
  device_real_t            X[batch_size*3*N];
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial);  

  FullereneDual dualG;
  FullereneGraph G;
  G.N = N;
  G.neighbours = vector<vector<node_t>>(N,vector<node_t>(3));
  vector<coord3d> points(N);
  
  size_t I=0,			// Isomer number
         i=0;			// Isomer number within current batch
  bool more_to_do = true;
  while(more_to_do){
    // Fill in a batch
    for(i=0; (i<batch_size) && more_to_do; i++, I++){
      //      printf("i=%ld, I=%ld\n",i,I);
      more_to_do &= BuckyGen::next_fullerene(Q,dualG);
      if(!more_to_do) break;
      
      dualG.update();   		        // Update triangles
      FullereneGraph   G = dualG.dual_graph();  // Construct fullerene graph
      G.layout2d         = G.tutte_layout();      
      vector<coord3d> X0 = G.zero_order_geometry(); // TODO: Faster, better X0

      for(node_t u=0;u<N;u++)
      	for(int j=0;j<3;j++){
      	  node_t v  = G.neighbours[u][j];
      	  size_t arc_index = i*3*N + u*3 + j;
      	  cubic_graph [arc_index] = v;
      	  next_on_face[arc_index] = G.next_on_face(u,v);
      	  prev_on_face[arc_index] = G.prev_on_face(u,v);
      	  face_right  [arc_index] = face_size(G,u,v);
      	  X           [arc_index] = X0[u][j];	  
      	}
    }

    size_t this_batch_size = i;

    IsomerspaceForcefield::OptimizeBatch(X,cubic_graph, next_on_face, prev_on_face, face_right,
    					 N,this_batch_size);

    // Now do something with the optimized geometries
    for(size_t i=0;i<this_batch_size;i++){
      for(node_t u=0;u<N;u++){
	size_t node_ix = i*3*N + u*3;
	G.neighbours[u] = {cubic_graph[node_ix], cubic_graph[node_ix+1], cubic_graph[node_ix+2]};
	points[u]       = {X[node_ix], X[node_ix+1], X[node_ix+2]};
      }

      Polyhedron P(G,points);
      string filename = output_dir+"/C"+to_string(N)+"-"+to_string(I);        
      Polyhedron::to_file(P,filename+".mol2");
    }
    // Output molecular geometry files
    
  }
  failures.close();
  
  return 0;
}
