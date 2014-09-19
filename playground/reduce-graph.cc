#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"

int testN = 80;
int testRSPI[12] = {1, 2, 3, 4, 5, 6, 37, 38, 39, 40, 41, 42};

Triangulation reduce_triangulation(Triangulation T){
 
  //for(int k=T.N; k>12; --k){
  for( ; T.N>12; --T.N){
  
    // find a deg-6-node (if the vertices are sorted, then just take the last)
  
    // record its neighbours (vector hole = neighbours.last)
    vector<int> hole(T.neighbours.back());
  
    // remove node and neighbor entries
    for(vector<int>::iterator it=T.neighbours.back().begin(), to=T.neighbours.back().end(); it!=to; ++it){
      T.neighbours[*it].erase(std::remove(T.neighbours[*it].begin(), T.neighbours[*it].end(), *it), T.neighbours[*it].end());
      //alternatively sort & pop_back
    }
    T.neighbours.pop_back();
  
    // patch the hole:
    for(int i=0; i<hole.size(); ++i){
      if(T.neighbours[hole[i]].size() == 2){ // if hole[k] is deg-2,
        // then connect hole[k] with hole[k+2i%hole.size] and remove hole[k+1] from hole
        T.neighbours[hole[i]].push_back(hole[(i+2) % hole.size()]);
        T.neighbours[hole[(i+2) % hole.size()]].push_back(hole[i]);
        hole.erase(hole.begin()+i+1);
      }
    }
  
    // for the remaining patch:
    while(hole.size() > 3){
      // connect hole[k] with hole[k+2i%hole.size] and remove hole[k+1] from hole
      // (this is not a very good retriangulation but valid and delaunay on 12
      // points is fast anyway)
      T.neighbours[hole[0]].push_back(hole[2]);
      T.neighbours[hole[2]].push_back(hole[0]);
      hole.erase(hole.begin()+1);
    }
  }

  return T;
}


int main(int ac, char **av) {
  int N;
  Triangulation::jumplist_t jumps;
  vector<int> RSPI(12);
  bool from_file = false;
  if (ac == 2) {
    from_file = true;
    N = 0;
  } else if (ac < 14) {
    N = testN;
    for (int i = 0; i < 12; i++)
      RSPI[i] = testRSPI[i] - 1;
  }
  if (ac >= 14) {
    N = strtol(av[1], 0, 0);
    for (int i = 0; i < 12; i++)
      RSPI[i] = strtol(av[i + 2], 0, 0) - 1;
  }
  if (ac > 14) { // General RSPI: RSPI followed by jumps.
    for (int i = 14; i < ac; i += 2)
      jumps.push_back(
          make_pair(strtol(av[i], 0, 0) - 1, strtol(av[i + 1], 0, 0)));
  }

  // int isomer_number = ac>=3? strtol(av[2],0,0) : 1;
  // bool IPR = ac>=4? strtol(av[3],0,0) : false;

  // IsomerDB DB(N,IPR);

  // FullereneGraph
  // g(IsomerDB::makeIsomer(N,DB.getIsomer(N,isomer_number,IPR)));

  vector<int> spiral(N / 2 + 2, 6);
  for (int i = 0; i < 12; i++)
    spiral[RSPI[i]] = 5;

  Triangulation T1(spiral);

  PlanarGraph g = T1.dual_graph();
  g.layout2d = g.tutte_layout();

  Triangulation T(g.dual_graph());

  cout << "number vertices in T: " << T.N << endl;
  cout << "neighbours in T: " << T.neighbours << endl;

  Triangulation rT = reduce_triangulation(T);

  cout << "number vertices in rT: " << rT.N << endl;
  cout << "neighbours in rT: " << rT.neighbours << endl;

  return 0;
}

