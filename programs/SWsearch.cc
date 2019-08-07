#include <libgraph/auxiliary.hh>
#include <libgraph/triangulation.hh>
#include <libgraph/isomerdb.hh>
#include <libgraph/patches.hh>
#include <contrib/buckygen-wrapper.hh>

// TODO:
//  1. Check SW-site detection. Why is it not symmetric?
//  2. Generate Schlegel-diagram flips

Triangulation reverse(const Triangulation& G)
{
  Triangulation Greversed = G;
  for(node_t u=0;u<G.N;u++) reverse(Greversed.neighbours[u].begin(),Greversed.neighbours[u].end());

  return Greversed;
}

template <typename S> class SparseMatrix: public Graph {
public:
  using Graph::Graph;
  
  map<dedge_t, S> values;
  
  const S& add_edge(const dedge_t& e, const S& v) {
    insert_edge(e);
    values[e] += v;
    return values[e];
  }

  friend ostream& operator<<(ostream& s, const SparseMatrix& A)
  {
    vector<dedge_t> edges_set = A.directed_edges();
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
  int N;
  IDCounter<general_spiral>  isomer_list;
  vector<SparseMatrix<int>> As;

  static vector<general_spiral> generate_all_isomers(int N) {
    vector<general_spiral> isomers;

    size_t i = 0;
    Triangulation G;    
    BuckyGen::buckygen_queue Q = BuckyGen::start(N,false,false);

    while(BuckyGen::next_fullerene(Q,G)){
      if(i%10000 == 0) cerr << "Isomer " << i<< endl;
      isomers.push_back( FullereneDual(G).get_rspi() );
    }

    BuckyGen::stop(Q);
    cerr << " done.\n";
    
    return isomers;
  }
  
  SWsearch(int N, int SWmaxlevel) : N(N), isomer_list(generate_all_isomers(N)),
				    As(SWmaxlevel+1,SparseMatrix<int>(isomer_list.size()))
  {
    cerr << "Generated list of "<<isomer_list.size()<<" isomers...\n";
    
    size_t i=0;
    for(auto &iso: isomer_list){
      const general_spiral &rspi = iso.first;
      const size_t           Sid = iso.second;
      FullereneDual G(N,rspi);
      Patches P(G);
      
      cerr << "Processing isomer " << (++i) << ", rspi="<< rspi << "... ";
      size_t Nsites = 0;
      for(int d=0;d<=SWmaxlevel;d++){
	vector<Patches::swsite_t> sites = P.swsites(d);
	Nsites += sites.size();
	//	cerr << "Stone-Wales level-" << d << " sites = " << sites << ";\n";
	
	//	As[d].add_edge({Sid,Sid},0); // Make sure that every node is included in graph.

	for(auto &site: sites){
	  FullereneDual Gp = P.swtransform(site,d);
	  general_spiral rspip = Gp.get_rspi();
	  
	  size_t Sidp = isomer_list(rspip);
	  if(!(Sidp >= 0)){
	    cerr << "RSPI = " << rspi << " -> RSPI' = " << rspip << endl;
	    abort();
	  }
	  
	  //	  cerr << d<<"-edge " << vector<size_t>{Sid,Sidp} << endl;
	  As[d].add_edge({Sid,Sidp},1);
	}
      }
      cerr << "Total number of sites: " << Nsites << endl;
    }
  }
};


int main(int ac, char **av)
{
  if(ac<=1) return -1;

  int N          = strtol(av[1],0,0);
  int SWmaxlevel = ac>2? strtol(av[2],0,0) : 0;

  SWsearch SWConnect(N,SWmaxlevel);
  for(int i=0;i<=SWmaxlevel;i++)
    cout << "A"<<i<<" = " << SWConnect.As[i] << ";\n";

  return 0;
}
