
#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

#include <iostream>
#include <ctime>

using namespace std;
typedef Triangulation::jumplist_t jumplist_t;

int main(int ac, char **av) {
  int N = strtol(av[1], 0, 0);
  int maxcount = strtol(av[2], 0, 0);
  bool pentagon_start = strtol(av[3], 0, 0);

  FullereneDual g;
  switch (N) {
  case 380: {
    g = FullereneDual(N, {44, 69, 70, 81, 82, 109, 118, 119, 143, 183, 184, 191},
		         {{109, 2}});
    break;
  }

  case 384: {
    g = FullereneDual(N, {28, 29, 30, 48, 144, 145, 169, 170, 189, 190, 191, 193},
		          {{48, 1}});
    break;
  }

  case 440: {
    vector<int> input_rspi();
    g = FullereneDual(N, {38, 39, 40, 61, 169, 170, 196, 197, 217, 218, 219, 221},
		          {{61, 1}});
    break;
  }

  case 672: {
    g = FullereneDual(N, {50, 52, 108, 110, 219, 251, 287, 301, 303, 305, 319, 337},
		          {{141, 1}});
    break;
  }

  default:
    return 1;
  }
  g.update(false);

  cout << "(* Order: {{Norig,Ngc,k,l},{time_RSPI,rate_RSPI},{time_GC,rate_GC}, jumps,rspi} *)\n";
  cout << "result = {\n";
  for (int k = 1; k < 10; k++) {
    for (int l = 1; l < 10; l++) {
      int Ngc = N*(k*k+k*l+l*l);
      if (Ngc > maxcount) continue;

      //      fprintf(stderr,"(N,k,l,Ngc) = (%d,%d,%d,%d)\n", N,k,l,Ngc);

      auto gc_start = clock();
      FullereneDual gkl = g.GCtransform(k,l);
      auto gc_end   = clock();
      
      //      cerr << "GC done " << endl;


      vector<int> rspi(12, 0);
      Triangulation::jumplist_t jumps;

      auto rspi_start = clock();
      //      for(int i=0;i<10;i++)
	gkl.get_rspi(rspi, jumps, true, true, pentagon_start);
      auto rspi_end   = clock();

      for(auto &j: jumps) j.first++;
      
      double time_rspi = (rspi_end-rspi_start) / (1.0 * CLOCKS_PER_SEC);
      double time_gc   = (gc_end-gc_start) / (1.0  * CLOCKS_PER_SEC);
      double rate_rspi = time_rspi/Ngc;
      double rate_gc   = time_gc/Ngc;

      std::cout.precision(10);

      cout << std::fixed
	   << "{" << (vector<int>{{N,Ngc,k,l}}) << ","
	   << "{" << time_rspi << "," << rate_rspi << "},"
	   << "{" << time_gc   << "," << rate_gc << "},"
	   << jumps << ","
	   << (rspi+1) << "},\n";
    }
  }
  cout << "{}};\n";


  return 0;
}
