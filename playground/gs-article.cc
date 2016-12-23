
// fullerene_graph_ptr goldberg_coxeter_(const fullerene_graph_ptr *g, const int *k, const int *l)
// {
//   FullereneGraph fg(**g);
//   fg.layout2d = fg.tutte_layout(); // FIXME remove, and pass argument to GCtransform?
//   fg.layout_is_spherical = false;
//   return new FullereneGraph(fg.GCtransform(*k,*l));
// }


#include "libgraph/fullerenegraph.hh"
#include "libgraph/polyhedron.hh"
#include "libgraph/isomerdb.hh"
#include "libgraph/triangulation.hh"


int main(int ac, char **av)
{
  int k = strtol(av[1],0,0);
  int l = strtol(av[2],0,0);
  cout << k << ", " << l << endl;

//  const int N=380;
//  vector<int> input_rspi({44,69,70,81,82,109,118,119,143,183,184,191});
//  Triangulation::jumplist_t input_jumps({{109,2}});

//  const int N=384;
//  vector<int> input_rspi({28,29,30,48,144,145,169,170,189,190,191,193});
//  Triangulation::jumplist_t input_jumps({{48,1}});

//  const int N=440;
//  vector<int> input_rspi({38,39,40,61,169,170,196,197,217,218,219,221});
//  Triangulation::jumplist_t input_jumps({{61,1}});

  const int N=672;
  vector<int> input_rspi({50,52,108,110,219,251,287,301,303,305,319,337});
  Triangulation::jumplist_t input_jumps({{141,1}});


  FullereneGraph g(N,input_rspi,input_jumps);
  cout << "g = " << g << ";\n";

  g.layout2d = g.tutte_layout();
  g.layout_is_spherical = false;
cout << "Tutte done" << endl;

  //const int k=3,l=5;
  FullereneGraph gkl(g.GCtransform(k,l));
cout << "gc done " << endl;

  gkl.layout2d = gkl.tutte_layout();
  gkl.layout_is_spherical = false;
cout << "Tutte done" << endl;

  vector<int> rspi(12,0);
  Triangulation::jumplist_t jumps;
  gkl.get_rspi_from_fg(rspi, jumps, true, true);

  cout << "spiral   = " << rspi << ";\n"
  	 << "jumplist = " << jumps << ";\n";
  
  return 0;
}
