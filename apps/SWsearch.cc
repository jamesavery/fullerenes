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
	rspi_t        rspip;
	bool success = gp.get_rspi_standard(rspip,true,true);
	assert(success);

	int Sidp = already_seen(rspip);
	if(Sidp < 0){ 
	  Sidp = already_seen.insert(rspip);
	  work_stack.push(Sidp);
	  cerr << "Inserting new RSPI number " << Sidp << ": " << rspip << endl;
	  assert(rspi_is_canonical(N,rspip));
	}
	G.insert_edge({Sid,Sidp});
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

  assert(rspi_is_canonical(N,rspi));

  SWsearch IsomerSearch(N,rspi);
  cout << "G = " << IsomerSearch.G << endl;
  return 0;
}
