#include <libgraph/triangulation.hh>
#include <libgraph/spiral.hh>
#include <contrib/buckygen-wrapper.hh>

struct topindex {
  size_t Wiener, Wi_min, Wi_max;
  double erho, rho;

  topindex() : Wiener(0), Wi_min(INT_MAX), Wi_max(0), erho(0), rho(0) {}
};

pair<topindex,topindex> calculate_topindex(const FullereneDual& Gd)
{
  topindex T, Td;

  // Topological indices for dual
  matrix<int> Dd = Gd.all_pairs_shortest_paths();
  //  cout << "Dd = " << Dd << endl;
  
  for(node_t u=0;u<Gd.N;u++){
    size_t Wi = 0;
    for(node_t v=0;v<Gd.N;v++) Wi += Dd(u,v);

    Td.Wi_max = max(Td.Wi_max,Wi);
    Td.Wi_min = min(Td.Wi_min,Wi);
    Td.Wiener += Wi;
  }

  Td.Wiener /= 2;
  Td.erho = Td.Wi_max / (double) Td.Wi_min;
  Td.rho =  2*Td.Wiener / (double) (Gd.N*Td.Wi_min);

  // Topological indices for fullerene graph
  PlanarGraph G = Gd.dual_graph();
  matrix<int> D = G.all_pairs_shortest_paths();
  //cout << "G = " << G << "\n\n";
  //  cout << "D = " << D << "\n\n";
  
  for(node_t u=0;u<G.N;u++){
    size_t Wi = 0;
    for(node_t v=0;v<G.N;v++) Wi += D(u,v);

    T.Wi_max = max(T.Wi_max,Wi);
    T.Wi_min = min(T.Wi_min,Wi);
    T.Wiener += Wi;
  }
  T.Wiener /= 2;
  T.erho = T.Wi_max / (double) T.Wi_min;
  T.rho =  2*T.Wiener / (double) (G.N*T.Wi_min);

  return {T,Td};
}

int main(int ac, char **av)
{
  int N          = ac>1? strtol(av[1],0,0) : 60;

  Graph g;
  Triangulation G;    
  BuckyGen::buckygen_queue Q = BuckyGen::start(N,false,false);
  
  while(BuckyGen::next_fullerene(Q,g)){
    FullereneDual G(g);
    pair<topindex,topindex> Ts = calculate_topindex(G);
    spiral_nomenclature name(G, spiral_nomenclature::FULLERENE,spiral_nomenclature::TRIANGULATION);

    printf("%s %ld %ld %ld %g %g %ld %ld %ld %g %g\n",
	   name.to_string().c_str(),
	   Ts.first.Wiener, Ts.first.Wi_min, Ts.first.Wi_max, Ts.first.erho, Ts.first.rho,
	   Ts.second.Wiener, Ts.second.Wi_min, Ts.second.Wi_max, Ts.second.rho, Ts.second.erho);
  }
  
  BuckyGen::stop(Q);
  
}
