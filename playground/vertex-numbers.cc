#include "libgraph/fullerenegraph.hh"
//#include "libgraph/polyhedron.hh"
#include "libgraph/triangulation.hh"


// triangulations: we get the permutation of vertices. done.
// cubic graphs: get dual (how are the new vertex numbers derived from the old?), order faces, translate back (how are indices inherited?)
// leap-frog: (permutation of leap frog?), identify "


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
  vector<int> permutation2;
  // get spiral permutation (vertices in trig)
  T2.get_spiral(S2, J2, permutation2);
  vector<int> vn(FG.vertex_numbers(T2, permutation2));
  cout << vn << endl;


// vertex numbers in a non cubic graph
  cout << "*** LF" << endl;
  vector<int> S3;
  jumplist_t J3;
  



  return 0;
}


