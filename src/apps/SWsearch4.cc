#include <libgraph/auxiliary.hh>
#include <libgraph/triangulation.hh>
#include <libgraph/isomerdb.hh>
#include <libgraph/patches.hh>
#include <stack>

/*
 * Procedure:
 *  set<spiral>    already_seen;
 *  stack<spiral>  work_stack;
 *
 *  Begin by adding can_spiral of Isomer 1 to stack. Then, while the stack is nonempty:
 *  0. S = work_stack.pop();
 *  1. Windup S to triangulation G
 *  2. Identify SW sites in G
 *  3. For each SW site (p,i), let G' = SWtransform(p,i) and S' = canonical_spiral(G'). 
 *  3.1 If S' \notin already_seen, add S' to already_seen and to work_stack.
 *  3.2 add arc S -> S' annotated with whatever auxiliary information we want
 */

typedef vector<int> rspi_t;
typedef vector<int> spiral_t;

bool rspi_is_canonical(int N, const rspi_t& rspi)
{
  rspi_t rspic(12);
  FullereneDual g(N,rspi);

  bool success = g.get_rspi_standard(rspic,true,true);
  assert(success);
  if(rspi != rspic){
    cerr << "Non-canonical spiral: " << rspi << " != canonical " << rspic << "\n";
    FullereneDual gp(N,rspic);
    
    cerr << "g = " << g << ";\n"
	 << "gp = " << gp << ";\n";
  }
  return rspi == rspic;
}

template <typename S> class SparseMatrix: public Graph {
public:
  using Graph::Graph;
  
  map<dedge_t, S> values;
  
  const S& add_edge(const dedge_t& e, const S& v) {
    insert_edge(e,-1,-1);
    values[e] += v;
    return values[e];
  }

  friend ostream& operator<<(ostream& s, const SparseMatrix& A)
  {
    set<dedge_t> edges_set = A.directed_edges();
    vector<dedge_t> edges(edges_set.begin(),edges_set.end());
    vector<S> edge_values(edges.size());

    for(int i=0;i<edges.size();i++){
	auto x = A.values.find(edges[i]);
	edge_values[i] = x->second;
      }

    s << "{" << edges << ", " << edge_values  << "}";
    return s;
  }  
};


class SWsearch {
public:
  IDCounter<rspi_t>   isomer_list;

  int N, Nisomers; 
  vector<SparseMatrix<int>> As;

  SWsearch(int N) : N(N), Nisomers(IsomerDB::number_isomers(N)), As(4,SparseMatrix<int>(Nisomers)) {
    IsomerDB db = IsomerDB::readPDB(N);

    assert(db.entries.size() == Nisomers);
    
    for(int i=0;i<Nisomers;i++){
      rspi_t rspi = rspi_t(db.entries[i].RSPI, db.entries[i].RSPI+12);
      cerr << "Inserting isomer number " << i << ": " << rspi << endl;
      isomer_list.insert(rspi);
    }
    
    for(auto &iso: isomer_list){
      const rspi_t &rspi = iso.first;
      rspi_t rspip;
      int Sid = iso.second;
      FullereneDual g(N,rspi);
      Patches P(g);


      for(int d=0;d<=3;d++){
	vector<Patches::swsite_t> sites = P.swsites(d);
	cerr << "Stone-Wales level-" << d << " sites = " << sites << ";\n";
	
	As[d].add_edge({Sid,Sid},0); // Make sure that every node is included in graph.

	for(auto &site: sites){
	  FullereneDual gp = P.swtransform(site,d);
	  bool success = gp.get_rspi_standard(rspip,true,true);
	  assert(success);
	  
	  int Sidp = isomer_list(rspip);
	  if(!(Sidp >= 0)){
	    cerr << "RSPI = " << rspi << " -> RSPI' = " << rspip << endl;
	    abort();
	  }
	  
	  cerr << d<<"-edge " << vector<int>{Sid,Sidp} << endl;
	  As[d].add_edge({Sid,Sidp},1);
	}
      }
    }
  }
};

class SWsearch0 {
public:
  IDCounter<rspi_t>   isomer_list;

  int N, Nisomers; 
  SparseMatrix<int> A;

  SWsearch0(int N) : N(N), Nisomers(IsomerDB::number_isomers(N)), A(SparseMatrix<int>(Nisomers)) {
    IsomerDB db = IsomerDB::readPDB(N);

    assert(db.entries.size() == Nisomers);
    
    for(int i=0;i<Nisomers;i++){
      rspi_t rspi = rspi_t(db.entries[i].RSPI, db.entries[i].RSPI+12);
      cerr << "Inserting isomer number " << i << ": " << rspi << endl;
      isomer_list.insert(rspi);
    }
    
    for(auto &iso: isomer_list){
      cerr << "Isomer number " << iso.second << ": " << iso.first << endl;
      const rspi_t &rspi = iso.first;
      rspi_t rspip;
      int Sid = iso.second;
      FullereneDual g(N,rspi);
      Patches P(g);


      vector<Patches::swsite_t> sites = P.swsites();
      cerr << "Stone-Wales level-" << 0 << " sites = " << sites << ";\n";
	
      A.add_edge({Sid,Sid},0); // Make sure that every node is included in graph.

      for(auto &site: sites){
	FullereneDual gp = P.swtransform(site);
	bool success = gp.get_rspi_standard(rspip,true,true);
	assert(success);
	  
	int Sidp = isomer_list(rspip);
	if(!(Sidp >= 0)){
	  cerr << "RSPI = " << rspi << " -> RSPI' = " << rspip << endl;
	  abort();
	}
	  
	cerr << 0<<"-edge " << vector<int>{Sid,Sidp} << endl;
	A.add_edge({Sid,Sidp},1);
      }
    }
  }
};



int main(int ac, char **av)
{
  if(ac<=1) return -1;

  int N = strtol(av[1],0,0);

  SWsearch0 SWConnect0(N);
  cout << "A0p0 = " << SWConnect0.A << ";\n";

  SWsearch SWConnect(N);
  cout << "A0p = " << SWConnect.As[0] << ";\n"
       << "A1p = " << SWConnect.As[1] << ";\n"
       << "A2p = " << SWConnect.As[2] << ";\n"
       << "A3p = " << SWConnect.As[3] << ";\n"
       << "Aps = {A0p,A1p,A2p,A3p};\n"
       << "As = Table[SparseArray[Table[(Aps[[d,1,i]]+1)->Aps[[d,2,i]],{i,Length[Aps[[d,1]]]}]],{d,Length[Aps]}];\n";


  return 0;
}
