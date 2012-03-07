#include "libgraph/fullerenegraph.hh"

coord3d surface(const coord2d& polar_angle)
{
  double theta = polar_angle.first, phi = polar_angle.second;
  return coord3d(cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi));
}

vector<coord3d> map_surface(const vector<coord2d>& angles)
{
  vector<coord3d> S(angles.size());
  for(unsigned int i=0;i<angles.size();i++) S[i] = surface(angles[i]);
  return S;
}


int main()
{
  FullereneGraph c20(FullereneGraph::C20());
  c20.layout2d = c20.tutte_layout();
  
#if 0
  cout << "c20 = " << c20 << endl;

  FullereneGraph 
    halma1(c20.halma_fullerene(1)),     // C80
    halma11(halma1.halma_fullerene(1)), // C320
    halma2(c20.halma_fullerene(2)),	// C180
    halma3(c20.halma_fullerene(3)),	// C320
    halma4(c20.halma_fullerene(4)),	// C500
    halma5(c20.halma_fullerene(5)),	// C720
    halma6(c20.halma_fullerene(6)),	// C980
    halma7(c20.halma_fullerene(7)),	// C1280
    halma111(halma11.halma_fullerene(1)); // C1280


  cout << "halma1 = "   << halma1 << endl; 
  cout << "halma11 = "  << halma11 << endl; 
  cout << "halma111 = " << halma111 << endl;

  cout << "halma2 = " << halma2 << endl; 
  cout << "halma3 = " << halma3 << endl; 
  cout << "halma4 = " << halma4 << endl; 
  cout << "halma5 = " << halma5 << endl; 
  cout << "halma6 = " << halma6 << endl; 
  cout << "halma7 = " << halma7 << endl; 
#endif

  FullereneGraph frog1(c20.leapfrog_fullerene());
  FullereneGraph frog2(frog1.leapfrog_fullerene());
  vector<coord3d> surface(map_surface(frog1.spherical_projection(frog1.layout2d)));

  
  cout << frog1.to_latex(frog1.layout2d,true) << endl;
  //  cout << frog1.to_latex(surface,false,false,true) << endl;
  // cout << "frog1 = " << frog1 << endl;
  return 0;
}
