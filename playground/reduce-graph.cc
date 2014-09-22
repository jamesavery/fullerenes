#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"


void insert_before(vector<int> &v, const int before_what, const int value){
  vector<int>::iterator pos = std::find(v.begin(), v.end(), before_what);
  v.insert(pos, value);
}


Triangulation reduce_triangulation(Triangulation T){

  for( int k=T.N; k>12; --k){
  //for( ; T.N>12; --T.N){

    cout << "n1" << T.neighbours << endl;
    cout << "n1" << T.neighbours << endl;
    cout << "N: " << T.N << endl;    
    cout << "T1=" << T << endl;
    cout << "deg N:" << T.neighbours[T.N-1].size() << endl;

    // find a deg-6-node (if the vertices are sorted, then just take the last)
  
    // save its neighbours
    vector<int> hole(T.neighbours.back());
 
    // remove node and neighbor entries
    for(vector<int>::iterator it=T.neighbours.back().begin(), to=T.neighbours.back().end(); it!=to; ++it){
//      cout << *it <<endl;
//      cout << T.neighbours[*it]<< endl;
      vector<int>::iterator to_erase = std::find(T.neighbours[*it].begin(), T.neighbours[*it].end(), T.N-1);
      //if(to_erase == T.neighbours[*it].end()) cout << "BAD" << endl;
      T.neighbours[*it].erase(to_erase);
    }
    T.neighbours.pop_back();
    --T.N;
  
    // patch the hole:
    // phase one: find nodes with degree-2 because they require a new connection
    for(int i=0; i<hole.size(); ++i){
      if(T.neighbours[hole[i]].size() == 2){ // if hole[k] is deg-2,
        cout << hole << "" << hole.size()<< endl;
        cout << "deg2 found, connecting " << hole[i] << " and " << hole[(i+2) % hole.size()] << endl;
        // then connect hole[k] with hole[k+2i%hole.size] and remove hole[k+1] from hole
        insert_before(T.neighbours[hole[i]],                 hole[(i-1+hole.size())%hole.size()], hole[(i+2)%hole.size()]);
        insert_before(T.neighbours[hole[(i+2)%hole.size()]], hole[(i+1)%hole.size()],             hole[i]);
        hole.erase(hole.begin() + ((i+1)%hole.size()));
      }
    }
  
    // phase two: triangulate the remaining hole in a fan-like manner:
    while(hole.size() > 3){

      int shift = 0;
      // check if 0 and 2 are connected already
      if(T.neighbours[hole[0]].end() != std::find(T.neighbours[hole[0]].begin(), T.neighbours[hole[0]].end(), hole[2])){
        shift = 1;
        // printf("%i and %i were connected already, connecting %i and %i instead.\n", hole[0], hole[2], hole[1], hole[3%hole.size()]);
      }

      // connect hole[k] with hole[k+2i%hole.size] and remove hole[k+1] from hole
      // (this is not a very good retriangulation but valid and delaunay on 12
      // points is fast anyway)
      insert_before(T.neighbours[hole[0+shift]],                 hole[(0+shift-1+hole.size())%hole.size()], hole[(0+shift+2)%hole.size()]);
      insert_before(T.neighbours[hole[(0+shift+2)%hole.size()]], hole[(0+shift+1)%hole.size()],             hole[0+shift]);
      hole.erase(hole.begin() + (0+shift+1)%hole.size());
    }
//    cout << "n3" << T.neighbours << endl;
//    cout << "T3=" << T << endl;
//  if(! T.is_consistently_oriented()) cout << "not consistently oriented" << endl;
//    cout << "------------------------" << endl;
  }

  return T;
}




int main(int ac, char **av) {
  int N;
  vector<int> RSPI(12);
  N = strtol(av[1], 0, 0);
  for (int i = 0; i < 12; i++)
    RSPI[i] = strtol(av[i + 2], 0, 0) - 1;

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

  Triangulation T(T1.sort_nodes());

  cout << "number vertices in T: " << T.N << endl;
  cout << "neighbours in T: " << T.neighbours << endl;

  Triangulation rT = reduce_triangulation(T);

  cout << "number vertices in rT: " << rT.N << endl;
  cout << "neighbours in rT: " << rT.neighbours << endl;
  cout << "rT=" << rT << endl;

  return 0;
}

