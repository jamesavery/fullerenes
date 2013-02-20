#include <stdlib.h>
#include <fstream>
#include <sstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/topological-indices.hh"

#include <vector>

using namespace std;

template <typename T> void printduplicate(string indextype, T value)
{
  cout << "Duplicate "<<indextype<<"-index: " << value << endl;
}

int main(int ac, char **av)
{
  FullereneGraph::jumplist_t jumps;
  char filename[0x100];

  if(ac<2) return -1;
  int N = strtol(av[1],0,0);

  sprintf(filename,"database/All/c%03d.rspi",N);

  //  cout << "Attempting to open " << filename << endl;
  ifstream rspifile(filename);
  string s;
  set<int>    wiener, hyperwiener, reversewiener, szeged;
  set<double> balaban;

  int linenumber =0;
  while(getline(rspifile,s)){
    vector<int> rspi(12);
    stringstream l(s);
    for(int i=0;i<12;i++){ l >> rspi[i]; rspi[i]--; if(l.fail()) abort(); }

    cout << (++linenumber) << ": " << rspi << endl;

    // Generate fullerene graph from spiral
    FullereneGraph g(N,rspi,jumps);
    // Are the topological indicators unique?
    TopologicalIndices T(g);
    
    int w = T.Wiener(), ww = T.hyperWiener(), rw = T.reverseWiener(), sz = T.Szeged();
    double b = T.Balaban();
    
    if(wiener.find(w) != wiener.end()) printduplicate("Wiener",w);
    else wiener.insert(w);

    if(hyperwiener.find(ww) != hyperwiener.end()) printduplicate("hyper-Wiener",ww);
    else hyperwiener.insert(ww);

    if(reversewiener.find(rw) != reversewiener.end()) printduplicate("reverse-Wiener",rw);
    else reversewiener.insert(rw);

    if(szeged.find(sz) != szeged.end()) printduplicate("Szeged",sz);
    else szeged.insert(sz);

    if(balaban.find(b) != balaban.end()) printduplicate("Balaban",b);
    else balaban.insert(b);
  }
  
  return 0;
}
