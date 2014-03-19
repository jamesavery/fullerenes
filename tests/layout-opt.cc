#include "libgraph/cubicgraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"
#include "libgraph/triangulation.hh"



int main(int ac, char **av)
{
  vector<int> rspi(12);
  FullereneGraph::jumplist_t jumps;

  if(ac != 14) return -1;

  int N = atol(av[1]);
  cout << N << endl;
  for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;

  cout << "Attempting to create graph from spiral indices " << rspi << endl;
 
  FullereneGraph g(N, rspi, jumps);
  cout << "done" << endl;

  cout << g << endl;

  g.layout2d = g.tutte_layout();
  cout << g.layout2d << endl;
  {
    ofstream latex(("output/test-layout-1.tex"));
    latex << g.to_latex(10, 10, false, true, true);
    latex.close();
  }

  g.optimize_layout(10, 10, 10);
  cout << g.layout2d << endl;
  {
    ofstream latex(("output/test-layout-10-10-10.tex"));
    latex << g.to_latex(10, 10, false, false, true);
    latex.close();
  }

  g.optimize_layout(10, 10, 20);
  cout << g.layout2d << endl;
  {
    ofstream latex(("output/test-layout-10-10-20.tex"));
    latex << g.to_latex(10, 10, false, false, true);
    latex.close();
  }

  g.optimize_layout(10, 10, 30);
  cout << g.layout2d << endl;
  {
    ofstream latex(("output/test-layout-10-10-30.tex"));
    latex << g.to_latex(10, 10, false, false, true);
    latex.close();
  }

  g.optimize_layout(10, 10, 50);
  cout << g.layout2d << endl;
  {
    ofstream latex(("output/test-layout-10-10-50.tex"));
    latex << g.to_latex(10, 10, false, false, true);
    latex.close();
  }

  g.optimize_layout(10, 10, 100);
  cout << g.layout2d << endl;
  {
    ofstream latex(("output/test-layout-10-10-100.tex"));
    latex << g.to_latex(10, 10, false, false, true);
    latex.close();
  }

  return 0;
}
