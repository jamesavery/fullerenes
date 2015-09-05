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

class SWsearch {
public:
  IDCounter<rspi_t> already_seen;
  stack<int>          work_stack;

  int N; 
  Graph               G;

  SWsearch(int N, const rspi_t& S0) : N(N), G(IsomerDB::number_isomers(N)) {
    int Sid = already_seen.insert(S0);
    work_stack.push(Sid);

    while(!work_stack.empty()){
      int           Sid = work_stack.top(); work_stack.pop();
      rspi_t        rspi = already_seen.invert(Sid);
      cerr << "Popping RSPI " << Sid << ": " << rspi << endl;
      FullereneDual g(N,rspi);
      Patches       P(g);
      vector<Patches::swsite_t> sites = P.swsites();
      cerr << "Stone-Wales sites: " << sites << endl;

      for(auto &site: sites){
	FullereneDual gp = P.swtransform(site);
	rspi_t        rspi;
	Triangulation::jumplist_t    j;
	bool success = gp.get_rspi(rspi,j,true,false);
	assert(success);

	int Spid = already_seen(rspi);
	if(Spid < 0){ 
	  Spid = already_seen.insert(rspi);
	  work_stack.push(Spid);
	  cerr << "Inserting new RSPI number " << Spid << ": " << rspi << endl;
	}
	G.insert_edge(Sid,Spid);
      }
    }
  }
};


int main(int ac, char **av)
{
  if(ac<=13) return -1;

  int N = strtol(av[1],0,0);
  rspi_t rspi;
  for(int i=0;i<12;i++) rspi.push_back(strtol(av[i+2],0,0));

  {
    rspi_t rspip(12);
    jumplist_t j;
    FullereneDual g(N,rspi);
    bool success = g.get_rspi(rspip,j,true,false);
    assert(success);
    cout << "RSPI = " << rspi << "; RSPI' = " << rspip << endl;
  }

  SWsearch IsomerSearch(N,rspi);
  cout << "G = " << IsomerSearch.G << endl;
  return 0;
}
