
#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

int main(int ac, char **av) {
  int N = strtol(av[1], 0, 0);
  int k = strtol(av[2], 0, 0);
  int l = strtol(av[3], 0, 0);
  bool pentagon_start = strtol(av[4], 0, 0);
  cout << N << ", " << k << ", " << l << endl;

  FullereneGraph g;
  switch (N) {
  case 140: {
    vector<int> input_rspi({0, 16, 19, 22, 25, 28, 47, 50, 53, 56, 59, 71});
    Triangulation::jumplist_t input_jumps = Triangulation::jumplist_t();
    g = FullereneGraph(N, input_rspi, input_jumps);
    break;
  }

  case 380: {
    vector<int> input_rspi({44, 69, 70, 81, 82, 109, 118, 119, 143, 183, 184, 191});
    Triangulation::jumplist_t input_jumps({{109, 2}});
    g = FullereneGraph(N, input_rspi, input_jumps);
    break;
  }

  case 384: {
    vector<int> input_rspi({28, 29, 30, 48, 144, 145, 169, 170, 189, 190, 191, 193});
    Triangulation::jumplist_t input_jumps({{48, 1}});
    g = FullereneGraph(N, input_rspi, input_jumps);
    break;
  }

  case 440: {
    vector<int> input_rspi({38, 39, 40, 61, 169, 170, 196, 197, 217, 218, 219, 221});
    Triangulation::jumplist_t input_jumps({{61, 1}});
    g = FullereneGraph(N, input_rspi, input_jumps);
    break;
  }

  case 672: {
    vector<int> input_rspi({50, 52, 108, 110, 219, 251, 287, 301, 303, 305, 319, 337});
    Triangulation::jumplist_t input_jumps({{141, 1}});
    g = FullereneGraph(N, input_rspi, input_jumps);
    break;
  }

  default:
    return 1;
  }

  cout << "g = " << g << ";\n";

  auto t0 = Clock::now();
  g.layout2d = g.tutte_layout();
  g.layout_is_spherical = false;
  g.orient_neighbours(); 
  cout << "Tutte done" << endl;

  auto t1 = Clock::now();
  FullereneGraph gkl;
  // if (l == 0) {
  //   gkl = g.halma_fullerene(k - 1, false);
  // } else {
  gkl = g.GCtransform(k, l, false);
  //  }
  cout << "gc done " << endl;

  auto t2 = Clock::now();
  gkl.layout2d = gkl.tutte_layout();
  gkl.layout_is_spherical = false;
  cout << "Tutte done" << endl;

  auto t3 = Clock::now();
  vector<int> rspi(12, 0);
  Triangulation::jumplist_t jumps;
  gkl.get_rspi_from_fg(rspi, jumps, true, true, pentagon_start);

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
