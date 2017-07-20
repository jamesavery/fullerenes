
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
  cout << N << ", " << k << ", " << l << endl;

  FullereneDual g;
  switch (N) {
  case 140: {

    break;
  }

  case 380: {
    g = FullereneDual_from_rspi(N, {44, 69, 70, 81, 82, 109, 118, 119, 143, 183, 184, 191});
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

  cout << "g = " << g << ";\n";

  auto t0 = Clock::now();
  cout << "Don't calculate layout, stay with triangulation" << endl;

  auto t1 = Clock::now();
  FullereneDual gkl;
  gkl = FullereneDual(g.GCtransform(k, l));
  cout << "gc done " << endl;

  auto t2 = Clock::now();
  cout << "Still don't calculate layout" << endl;

  auto t3 = Clock::now();
  vector<int> rspi(12, 0);
  jumplist_t jumps;
  gkl.get_rspi(rspi, jumps, true, true, pentagon_start);
  cout << "Got RSPI: " << rspi << endl;
  
  ofstream output(("spiral-" + to_string(N) + "-" + to_string(k) + "-" +
                   to_string(l) + "-" + to_string(pentagon_start)).c_str());
  for (int i = 0; i < 12; i++) rspi[i]++;
  for (auto it = jumps.begin(); it != jumps.end(); it++) it->first++;
  auto t4 = Clock::now();
  output << "N = " << N *(k * k + k * l + l * l) << ";\n"
         << "jumplist = " << jumps << ";\n"
         << "spiral   = " << rspi << ";\n";

  output << "Delta t1-t0 (layout): " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0) .count() << " ms" << std::endl;
  output << "Delta t2-t1 (gc): " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1) .count() << " ms" << std::endl;
  output << "Delta t3-t2 (layout): " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2) .count() << " ms" << std::endl;
  output << "Delta t4-t3 (spiral): " << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3) .count() << " ms" << std::endl;
  output.close();

  return 0;
}
