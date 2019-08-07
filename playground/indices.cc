#include <stdlib.h>
#include <fstream>
#include <sstream>
#include "libgraph/cubicgraph.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/topological-indices.hh"

#include <vector>

using namespace std;

template <typename T> void checkduplicates(const string& indextype, const map<T,int>& counts)
{
  cout << "Duplicate "<<indextype<<"-indices (value:count):\n";
  for(typename map<T,int>::const_iterator kv(counts.begin()); kv!=counts.end();kv++)
    if(kv->second>1){
      cout << "\t"<<kv->first <<": " << kv->second << endl;
    }
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
  map<int,int>    wiener, hyperwiener, reversewiener, szeged;
  map<double,int> balaban;

  int linenumber =0;
  while(getline(rspifile,s)){
    vector<int> rspi(12);
    stringstream l(s);
    for(int i=0;i<12;i++){ 
      l >> rspi[i]; rspi[i]--; 
      if(l.fail()){ cerr << "Malformed line: " << s << endl;  abort(); }
    }

    cerr << (++linenumber) << ": " << rspi << endl;

    // Generate fullerene graph from spiral
    FullereneGraph g(N,rspi,jumps);
    // Are the topological indicators unique?
    TopologicalIndices T(g);
    
    wiener[T.Wiener()]++;
    hyperwiener[T.hyperWiener()]++;
    reversewiener[T.reverseWiener()]++;
    szeged[T.Szeged()]++;
    balaban[T.Balaban()]++;
  }

  checkduplicates("Wiener",wiener);
  checkduplicates("hyperWiener",hyperwiener);
  checkduplicates("reverseWiener",reversewiener);
  checkduplicates("Szeged",szeged);
  checkduplicates("Balaban",balaban);
  
  return 0;
}
