
#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

int main(int ac, char **av) {
  int N;
  Triangulation::jumplist_t input_jumps;
  vector<int> RSPI(12);
  if(ac>=14) {
    N = strtol(av[1],0,0);
    for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;
  } else{assert(false);} 
  
  if(ac>14){			// General RSPI: RSPI followed by jumps.
    for(int i=14;i<ac;i+=2)
      input_jumps.push_back(make_pair(strtol(av[i],0,0)-1,strtol(av[i+1],0,0)));
  }

  cout << "N = " << N  << ";\n"
       << "jumplist = " << input_jumps << ";\n"
       << "spiral   = " << RSPI << ";\n";

  auto t0 = Clock::now();
  FullereneGraph g(N, RSPI, input_jumps);
  auto t1 = Clock::now();
  cout << "Delta t1-t0 (gen graph): " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0) .count() << " ms" << std::endl;

  cout << "g = " << g << ";\n";

  auto t2 = Clock::now();
  g.layout2d = g.tutte_layout();
  g.layout_is_spherical = false;
  cout << "Tutte done" << endl;

  auto t3 = Clock::now();
  cout << "Delta t3-t2 (layout): " << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2) .count() << " ms" << std::endl;
  auto t4 = Clock::now();

  vector<int> rspi(12, 0);
  Triangulation::jumplist_t jumps;
  bool pentagon_start = true;

  for (int foo=0; foo<10 ; foo++){
    t3 = Clock::now();
    g.get_rspi_from_fg(rspi, jumps, true, pentagon_start);
    t4 = Clock::now();
    cout << "Delta t4-t3 (spiral): " << std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3) .count() << " ms" << std::endl;
  }
    for (int i = 0; i < 12; i++) rspi[i]++;
    for (auto it = jumps.begin(); it != jumps.end(); it++) it->first++;
  cout << "N = " << N  << ";\n"
       << "jumplist = " << jumps << ";\n"
       << "spiral   = " << rspi << ";\n";

  system("cat /proc/cpuinfo | grep 'MHz'");

  return 0;
}
