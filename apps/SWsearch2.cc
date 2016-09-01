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
  IDCounter<rspi_t>   isomer_list;

  int N, Nisomers; 
  Graph               G0, G1, G2, G3;

  SWsearch(int N) : N(N), Nisomers(IsomerDB::number_isomers(N)), G0(Nisomers), G1(Nisomers), G2(Nisomers), G3(Nisomers) {
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

      vector<Patches::swsite_t> sites0 = P.swsites(), sites1 = P.swsites1(), sites2 = P.swsites2(), sites3 = P.swsites3();
      
      for(auto &site: sites0){
	FullereneDual gp = P.swtransform(site);
	bool success = gp.get_rspi_standard(rspip,true,true);
	assert(success);

	int Sidp = isomer_list(rspip);
	assert(Sidp >= 0);

	cerr << "0-edge " << vector<int>{Sid,Sidp} << endl;
	G0.insert_edge({Sid,Sidp},-1,-1);
      }

      for(auto &site: sites1){
      	FullereneDual gp = P.swtransform1(site);
      	bool success = gp.get_rspi_standard(rspip,true,true);
      	assert(success);

      	int Sidp = isomer_list(rspip);
      	assert(Sidp >= 0);

      	cerr << "1-edge " << vector<int>{Sid,Sidp} << endl;
      	G1.insert_edge({Sid,Sidp},-1,-1);
      }

      for(auto &site: sites2){
      	FullereneDual gp = P.swtransform2(site);
      	bool success = gp.get_rspi_standard(rspip,true,true);
      	assert(success);

      	int Sidp = isomer_list(rspip);
      	assert(Sidp >= 0);

      	cerr << "2-edge " << vector<int>{Sid,Sidp} << endl;
      	G2.insert_edge({Sid,Sidp},-1,-1);
      }

      for(auto &site: sites3){
      	FullereneDual gp = P.swtransform3(site);
      	bool success = gp.get_rspi_standard(rspip,true,true);
      	assert(success);

      	int Sidp = isomer_list(rspip);
      	assert(Sidp >= 0);

      	cerr << "3-edge " << vector<int>{Sid,Sidp} << endl;
      	G3.insert_edge({Sid,Sidp},-1,-1);
      }
    }
  }
};



int main(int ac, char **av)
{
  if(ac<=1) return -1;

  int N = strtol(av[1],0,0);

  SWsearch SWConnect(N);
  cout << "G0 = " << SWConnect.G0 << ";\n"
       << "G1 = " << SWConnect.G1 << ";\n"
       << "G2 = " << SWConnect.G2 << ";\n"
       << "G3 = " << SWConnect.G3 << ";\n";

  return 0;
}
