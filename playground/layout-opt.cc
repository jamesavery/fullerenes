#include "libgraph/cubicgraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"
#include "libgraph/triangulation.hh"



int main(int ac, char **av)
{
  Polyhedron P;
  PlanarGraph g;
  bool from_file = ac<14;
  int N = 0;
  if(from_file){
    P = Polyhedron(av[1]);
    N = P.N;
    g = P;
    face_t outer_face;
    for(int i=2;i<ac;i++) outer_face.push_back(strtol(av[i],0,0));
    cout << "Reading in graph from file '"<<av[1]<<"'; outer face set to " << outer_face << endl;
    g.layout2d = g.tutte_layout(outer_face);//face_t({0,15,30,8,10,35,2})
  } else {
    vector<int> rspi(12);
    FullereneGraph::jumplist_t jumps;
    
    N = atol(av[1]);
    for(int i=0;i<12;i++) rspi[i] = strtol(av[2+i],0,0)-1;

    cout << "Attempting to create graph from spiral indices " << rspi << endl;
    
    g = FullereneGraph(N, rspi, jumps);
    cout << "done" << endl;
    g.layout2d = g.tutte_layout();
  }

  cout << g << endl;


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
