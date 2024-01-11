#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/spiral.hh"
#include <fstream>

vector<int> rspi_c240_ih = {0,32,36,40,44,48,74,78,82,86,90,121};
string spiralcode_c240_ih = "[GS: 1,33,37,41,45,49,75,79,83,87,91,122]-fullerene";
string spiralcode_c60_ih  = "[GS: 1,7,9,11,13,15,18,20,22,24,26,32]-fullerene";

polygon C240_unfolding_outline{vector<Eisenstein>{{
  {0,14},   {0,16},   {2,16},   {4,16},  
  {6,14},   {8,12},   {10,12},  {8,14},  {8,16},
  {10,16},  {12,14},  {12,16},  {14,16}, 
  {16,14},  {16,12},  {16,10}, 		   
  {14,10},  {12,10},  {12, 8},  {14, 8}, {16,6},
  {16, 4},  {14, 4},  {16, 2},  {16, 0},
  {14, 0},  {12, 2},  {10, 4}, 		   
  {10, 6},  {10, 8},  {8, 10},  {8, 8},{6, 8},
  {4, 10},  {4, 12},  {2, 12}
}}};

polygon C60_unfolding_outline{vector<Eisenstein>{{
  {0,7},   {0,8},   {1,8},   {2,8},  
  {3,7},   {4,6},   {5,6},  {4,7},  {4,8},
  {5,8},  {6,7},  {6,8},  {7,8}, 
  {8,7},  {8,6},  {8,5}, 		   
  {7,5},  {6,5},  {6, 4},  {7, 4}, {8,3},
  {8, 2},  {7, 2},  {8, 1},  {8, 0},
  {7, 0},  {6, 1},  {5, 2}, 		   
  {5, 3},  {5, 4},  {4, 5},  {4, 4},{3, 4},
  {2, 5},  {2, 6},  {1, 6}
}}};


int main(int ac, char **av)
{
  //  spiral_nomenclature spiral(spiralcode_c240_ih);
  spiral_nomenclature spiral(spiralcode_c60_ih);
  Triangulation G = Triangulation(spiral).sort_nodes();
  //  G.triangles = G.compute_faces();
  
  FullereneGraph cubic = G.dual_graph();

  cubic.layout2d = cubic.tutte_layout();
  Polyhedron P(cubic,cubic.zero_order_geometry(),6);
  P.points = cubic.optimized_geometry(P.points);

  tri_t T0 = {0,12,G.next(0,12)};
  
  //  Unfolding U(G,C240_unfolding_outline,T0);
  Unfolding U(G,C60_unfolding_outline,T0);


  
  auto keys = get_keys(U.arc_coords);
  auto values = get_values(U.arc_coords);
  for(int i=0;i<keys.size();i++){
    cout << "arc " << keys[i] << " at " << values[i] << endl;
  }
  
  cout << "dual_neighours = " << G.neighbours << ";\n\n";
  cout << "cubic_neighbours = " << G.dual_graph().neighbours << ";\n\n";
  cout << "triangles = " << G.compute_faces() << ";\n\n";  
  cout << "dual_arc_list = " << get_keys(U.arc_coords) << ";\n\n";
  cout << "dual_arc_coordinates = " << get_values(U.arc_coords) << ";\n\n";  
  cout << "points_opt = " << P.points << ";\n\n";
  return 0;
}
