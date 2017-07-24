
#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
typedef Triangulation::jumplist_t jumplist_t;

FullereneDual FullereneDual_from_rspi(int N, const vector<int> &rspi,
				      const jumplist_t jumps = jumplist_t())
{
  vector<int> spiral_string(N/2+2,6);
  for(auto r: rspi) spiral_string[r] = 5;

  return Triangulation(spiral_string,jumps);
}


int main(int ac, char **av) {
  int N = strtol(av[1], 0, 0);
  int k = strtol(av[2], 0, 0);
  int l = strtol(av[3], 0, 0);
  bool pentagon_start = strtol(av[4], 0, 0);
  //  cerr << N << ", " << k << ", " << l << endl;

  FullereneDual g;
  switch (N) {
  case 140: {
    g = FullereneDual_from_rspi(N, {0, 16, 19, 22, 25, 28, 47, 50, 53, 56, 59, 71});
    break;
  }

  case 380: {    
    g = FullereneDual_from_rspi(N, {44, 69, 70, 81, 82, 109, 118, 119, 143, 183, 184, 191},jumplist_t{{109, 2}});
    break;
  }

  case 384: {
    g = FullereneDual_from_rspi(N, {28, 29, 30, 48, 144, 145, 169, 170, 189, 190, 191, 193},jumplist_t{{48, 1}});
    break;
  }

  case 440: {
    g = FullereneDual_from_rspi(N,{38, 39, 40, 61, 169, 170, 196, 197, 217, 218, 219, 221},jumplist_t{{61, 1}});
    break;
  }

  case 672: {
    g = FullereneDual_from_rspi(N,{50, 52, 108, 110, 219, 251, 287, 301, 303, 305, 319, 337},
				jumplist_t{{141, 1}});
    break;
  }

  default:
    return 1;
  }

  
  FullereneDual gkl;
  auto gc_start = std::clock();
  gkl = FullereneDual(g.GCtransform(k, l));
  auto gc_end   = std::clock();
  cerr << "GC done " << endl;

  vector<int> rspi(12, 0);
  jumplist_t jumps;
  
  auto grspi_start = std::clock();
  gkl.get_rspi(rspi, jumps, true, true, pentagon_start);
  auto grspi_end   = std::clock();
  cerr << "Got RSPI: " << rspi << endl;
  
  // ofstream output(("spiral-" + to_string(N) + "-" + to_string(k) + "-" +
  //                  to_string(l) + "-" + to_string(pentagon_start)).c_str());

  for (int i = 0; i < 12; i++) rspi[i]++;
  for (auto it = jumps.begin(); it != jumps.end(); it++) it->first++;
  auto t4 = Clock::now();
  double grspi_time = (grspi_end-grspi_start) * 1.0/CLOCKS_PER_SEC;
  double gc_time = (gc_end-gc_start) * 1.0/CLOCKS_PER_SEC;
  
  cout << "Ngc       = " << N *(k * k + k * l + l * l) << ";\n"
       << "{n,k,l}   = " << vector<int>{N,k,l} << ";\n"
       << "jumplist  = " << jumps << ";\n"
       << "spiral    = " << rspi << ";\n"
       << "gctime    = " << gc_time << ";\n"
       << "grspitime = " << grspi_time << ";\n";

  
  
  //  output.close();

  return 0;
}
