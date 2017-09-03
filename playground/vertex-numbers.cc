#include "libgraph/fullerenegraph.hh"
//#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"


// triangulations: we get the permutation of vertices. done.
// cubic graphs: get dual (how are the new vertex numbers derived from the old?), order faces, translate back (how are indices inherited?)
// leap-frog: (permutation of leap frog?), identify "


// elongated square bipyramid
PlanarGraph non_cubic_graph(){
  const int M=4, N=10;
  neighbours_t neighbours(N,vector<node_t>(4));

  for(int i=0; i!=M; ++i){
    neighbours[i][0] = 8;
    neighbours[i][1] = (i-1+M)%M;
    neighbours[i][2] = M+i;
    neighbours[i][3] = (i+1)%M;

    neighbours[M+i][0] = i;
    neighbours[M+i][1] = M+(i-1+M)%M;
    neighbours[M+i][2] = 9;
    neighbours[M+i][3] = M+(i+1)%M;
  }
  neighbours[8][0] = 0;
  neighbours[8][1] = 1;
  neighbours[8][2] = 2;
  neighbours[8][3] = 3;

  neighbours[9][0] = 7;
  neighbours[9][1] = 6;
  neighbours[9][2] = 5;
  neighbours[9][3] = 4;
  return PlanarGraph(neighbours,true);
}



int testN = 28;
int testRSPI[12] = {1,2,3,4,5,7,10,12,13,14,15,16};

int main(int ac, char **av)
{
  int N;
  Triangulation::jumplist_t jumps;
  vector<int> RSPI(12);
  if(ac==2){
    N = 0;
  } else if(ac<14){
    N = testN;
    for(int i=0;i<12;i++) RSPI[i] = testRSPI[i]-1;
  }
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  }
  if(ac>14){			// General RSPI: RSPI followed by jumps.
    for(int i=14;i<ac;i+=2)
      jumps.push_back(make_pair(strtol(av[i],0,0)-1,strtol(av[i+1],0,0)));
  }


  vector<int> spiral(N/2+2,6);
  for(int i=0;i<12;i++) spiral[RSPI[i]] = 5;

// vertex numbers in the triangulation
  cout << "*** trig" << endl;
  Triangulation T1(spiral,jumps);
  cout << T1 << endl;

  vector<int> S1;
  jumplist_t J1;
  vector<int> permutation1;
  T1.get_spiral(S1, J1, permutation1);
  cout << "vertex numbers of triangulation: " << permutation1 << endl;


  // vertex numbers in the cubic graph
  cout << "*** cub" << endl;
  FullereneGraph FG(N, RSPI,jumps);
  FG.layout2d = FG.tutte_layout();

  // get triangulation
  Triangulation T2(FG.dual_graph());
  vector<int> S2;
  jumplist_t J2;
  vector<int> perm2;
  // get spiral permutation (vertices in trig)
  T2.get_spiral(S2, J2, perm2);
  vector<int> vn2(FG.vertex_numbers(T2, perm2));
  cout << vn2 << endl;


// vertex numbers in a non cubic graph
  cout << "*** LF" << endl;
  PlanarGraph PG(non_cubic_graph());
  PG.is_oriented = true;

  Triangulation T3(PG.leapfrog_dual()); // and the first PG.N vertices are copied from PG with the order being preserved.  That's great!
  
  vector<int> S3;
  jumplist_t J3;
  vector<int> perm3;
  T3.get_spiral(S3, J3, perm3);

  vector<int> vn3(PG.vertex_numbers(T3, perm3));
  cout << vn3 << endl;

  return 0;
}

