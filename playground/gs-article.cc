
#include "libgraph/fullerenegraph.hh"
#include "libgraph/triangulation.hh"
#include "libgraph/polyhedron.hh"

int main(int ac, char **av) {
  int N = strtol(av[1], 0, 0);
  int k = strtol(av[2], 0, 0);
  int l = strtol(av[3], 0, 0);
  cout << N << ", " << k << ", " << l << endl;

  FullereneGraph g;
  switch (N) {
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

  g.layout2d = g.tutte_layout();
  g.layout_is_spherical = false;
  cout << "Tutte done" << endl;


      FullereneGraph gkl(g.GCtransform(k, l, true));
      cout << "gc done " << endl;

      gkl.layout2d = gkl.tutte_layout();
      gkl.layout_is_spherical = false;
      cout << "Tutte done" << endl;

      vector<int> rspi(12, 0);
      Triangulation::jumplist_t jumps;
      bool pentagon_start = true;
      gkl.get_rspi_from_fg(rspi, jumps, true, true, pentagon_start);

      ofstream output(("spiral-" + to_string(N) + "-" + to_string(k) + "-" + to_string(l)).c_str());
      for (int i = 0; i < 12; i++) rspi[i]++;
      for (auto it = jumps.begin(); it != jumps.end(); it++) it->first++;
      output << "N = " << N*(k * k + k * l + l * l) << ";\n"
             << "jumplist = " << jumps << ";\n"
             << "spiral   = " << rspi << ";\n";
      output.close();

  return 0;
}
