#include <limits.h>
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"

#include <iostream>


matrix<int> pentagon_distance(const Triangulation &G)
{
  vector<int> pentagon_indices(12);

  for(int u=0, i=0;u<G.N;u++) if(G.neighbours[u].size() == 5) pentagon_indices[i++] = u;

  return G.all_pairs_shortest_paths(pentagon_indices);
}

int max_length(const matrix<int>& P)
{
  int max_length = 0, A = -1, B = -1;
  for(int i=0;i<12;i++)
    for(int j=i+1;j<12;j++)
      if(P[i,j] > max_length)
	max_length = P[i,j];
  return max_length;
}

struct dendrogram_t {
  vector<int> class_tree, levels;
};

int bitcount16(uint16_t x)
{
  return 0;			// Tag fra bit-twiddling hacks 
}

dendrogram_t hierarchical_clustering(const matrix<int>& P, int K)
{
  matrix<int> dist = P;
  vector<int> class_tree({0,1,2,3,4,5,6,7,8,9,10,11}),
              levels(12);

  int k=0;
  while(k<K){
    int min_length = INT_MAX, A=-1,B=-1;

    // Find smallest distance between clusters
    for(int i=0;i<12;i++)
      for(int j=i;j<12;j++)
	if( dist[i,j] <  min_length)
	  min_length = dist[i,j], A = i, B = j;

    // Merge equivalence classes
    class_tree[max(A,B)] = min(A,B);
    levels[k] = dist[A,B];

    // Count equivalence classes
    uint16_t class_ids = 0;
    for(int i=0;i<12;i++) class_ids |= (1<<class_tree[i]);
    k = bitcount16(class_ids);
    
    // Update distance matrix
    for(int i=0;i<12;i++){
      dist[max(A,i),min(A,i)] = max(dist[A,i],dist[B,i]);
      dist[max(B,i),min(B,i)] = max(dist[A,i],dist[B,i]);
    }
  }
  return {class_tree,levels};
}

pair<bool,pair<int,int> > cluster_sizes(const matrix<int>& P, int d, int D)
{
  int mx_length = max_length(P);
  int Ad=0, AD=0, Bd=0, BD=0;

  for(int A=0;A<12;A++)
	for(int B=B+1;B<12;B++) if(P[A,B] == mx_length) {
		  for(int j=0;j<12;j++){
    			Ad += (P[A,j]<=d);		// Hvor mange femkanter er indenfor radius d fra A
 			AD += (P[A,j]<=D);		// Ditto D
    
    			Bd += (P[B,j]<=d);		// Osv.
    			BD += (P[B,j]<=D);        
  		}

  		if(Ad+BD == 12){ printf("A %d\n",mx_length); return {true,{Ad,BD}}; }
  		if(AD+Bd == 12){ printf("B %d\n",mx_length); return {true,{AD,Bd}}; }
 	}
 
  return {false,{0,0}};
}


int main(int ac, char **argv)
{
  int N                = strtol(argv[1],0,0);     // Argument 1: Number of vertices N
  
  if(ac<2 || N<20 || N==22 || N&1){
    fprintf(stderr,"Syntax: %s <N:int> [output_dir] [IPR:0|1] [only_nontrivial:0|1]\n",argv[0]);
    return -1;
  }

  string output_dir   = ac>=3? argv[2] : "output";    // Argument 2: directory to output files to
  int IPR             = ac>=4? strtol(argv[3],0,0):0; // Argument 3: Only generate IPR fullerenes?
  int only_nontrivial = ac>=5? strtol(argv[4],0,0):0; // Argument 4: Only generate fullerenes with nontrivial symmetry group?
  
  ofstream failures((output_dir+"/failures.txt").c_str()); // output/failures.txt contains list of any fullerenes that failed optimization

  int i=0;  
  Triangulation dualG;
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,IPR,only_nontrivial);  

  while(BuckyGen::next_fullerene(Q,dualG)){ // Generate all appropriate C_N isomer duals 
    i++;

    if(i%100000 == 0) cerr << "isomer "<< i << endl;

#if 0    
    auto C = cluster_sizes(pentagon_distance(dualG), 2, 4);
    if(C.first){
      cout << "C"<<i<<" = " << C.second<< ";\n";
      
    // if(0)
    //   if(max_length(pentagon_distance(dualG)) >= 10){
      //	cout << i << endl;
	
	dualG.update();		            // Update metadata
	FullereneGraph G = dualG.dual_graph();  // Construct fullerene graph
	vector<int> RSPI;
	jumplist_t jumps;
	G.get_rspi_from_fg(RSPI,jumps);
	cout << i << ", RSPI="<<(RSPI+1) << ", jumps="<<jumps<<endl;
      }
#endif
  }
  cout << i << " graphs.\n";
  failures.close();
  
  return 0;
}
