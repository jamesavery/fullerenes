#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/auxiliary.hh"
#include <fstream>


Polyhedron fullerene_dual_polyhedron(const Triangulation& dg)
{
  PlanarGraph pg(dg.dual_graph());
  cout << "pg = " << pg << endl;

  FullereneGraph g(pg);
  g.layout2d = g.tutte_layout();

  vector<coord3d> points = g.zero_order_geometry();
  points = g.optimized_geometry(points);

  vector<coord3d> dual_points(dg.N);

  vector<face_t> faces(dg.N);
  for(int i=0;i<dg.triangles.size();i++)
    for(int j=0;j<3;j++)
      faces[dg.triangles[i][j]].push_back(i);

  for(int i=0;i<faces.size();i++)
    dual_points[i] = faces[i].centroid(points);

  return Polyhedron(dg, dual_points);
}

vector<pair<tri_t,vector<Eisenstein>>> EisenTris(const Triangulation& dual)
{
  vector<pair<tri_t,vector<Eisenstein>>> result;
  switch(2*dual.N-4){
  case 120: {
    Eisenstein i{1,0}, j{0,1};
    vector<Eisenstein> 
      A{{1,-4},{2,-1},{-3,5}},
      B{{1,3}, {-2,1},{1,-4}},
      C{{-1,-1},{3,0},{-2,1}},
      D{{3,-3},{0,3},{-3,0}};
    
      result = {{{21,0,1},  A},
		{{1,20,21}, A*-i},
		{{20,1,2},  A},
		{{2,19,20}, A*-i},
		{{19,2,3},  A},
		{{3,18,19}, A*-i},
		{{3,10,18}, B},
		{{10,3,8},  B*-i},
		{{10,8,9},  A},
		{{3,4,8},   C},
		{{5,6,4},   C*(j-i)},
		{{7,8,6},   C*-j},
		{{4,6,8},   D},
		{{18,11,16},D*j},
		{{17,18,16},C*j},
		{{10,11,18},C*-i},
		{{15,16,11},C*(i-j)},
		{{15,13,14},B*(i-j)},
		{{13,15,11},B*(j-i)},
		{{13,11,12},A*(i-j)}};
  }
    break;
  default:
    cerr << "Eisentris not implemented yet.\n";
    break;
  }
  return result;
}

vector<tri_t> Tris(const Triangulation& dual)
{
  vector<tri_t> result;		// TODO: Vi kan bruge de rigtige knudenavne (u,v,w) -- den orienterede kant fastlægger nummeret i outline. 
  switch(2*dual.N-4){
  case 120: 
    result = {{0,1,21},{1,2,20},{2,3,19},{3,4,5},{3,5,10},{5,6,10},{10,11,3},
	      {20,21,1},{19,20,2},{18,19,3},{3,11,18},{11,13,18},{11,12,13},
	      {6,7,8},{6,10,8},{8,9,10},{17,18,13},{13,14,15},{15,16,17},{13,15,17}};
    break;
  case 102:
    result = {{0,1,2},{2,3,4},{4,5,6},{6,7,8},{4,6,8},{2,4,8},{2,8,9},{0,2,9},{0,9,10},
	      {0,10,20},{0,20,21},{10,11,13},{11,12,13},{13,14,15},{10,13,19},{10,19,20},
	      {15,16,19},{16,17,18}};
    break;
  default:
    cerr << "Eisentris not implemented yet.\n";
    break;
  }
  return result;
}

string latextris(const Unfolding& UF, const vector<tri_t>& tris)
{
  string s;
  ostringstream latex(s);
  for(auto t: tris)
    latex << "\\draw[overlay] ("<<t[0]<<") -- ("<<t[1]<<") -- ("<<t[2]<<") -- ("<<t[0]<<") -- cycle;\n";

  return latex.str();
}

coord3d barycentric(const coord2d& x1,const coord2d& x2,const coord2d& x3, const coord2d& x){
  double det     = ((x2.second-x3.second)*(x1.first-x3.first) +
		    (x3.first-x2.first)*(x1.second-x3.second));
  double s = ((x2.second-x3.second)*(x.first-x3.first) +
		    (x3.first-x2.first)*(x.second-x3.second))/det;
  double t = ((x3.second-x1.second)*(x.first-x3.first) + 
		    (x1.first-x3.first)*(x.second-x3.second))/det;
  double r = 1-t-s;

  return coord3d{s,t,r};
}

vector<pair<node_t,coord3d>> barycentric_coordinates(const tri_t& t, const Folding& F)

{  
  vector<Eisenstein> t_IJ = {F.outline[t[0]].first,F.outline[t[1]].first,F.outline[t[2]].first};
  vector<coord2d>    t_x  = {t_IJ[0].coord(), t_IJ[1].coord(), t_IJ[2].coord()};

  vector<Eisenstein> allpoints   = polygon(t_IJ).allpoints();

  vector<pair<node_t,coord3d>> bs(allpoints.size());
  for(int i=0;i<allpoints.size();i++){
    bs[i].first  = F.grid(allpoints[i]);
    bs[i].second = barycentric(t_x[0],t_x[1],t_x[2], allpoints[i].coord());
  }
  //TODO: Polygon::allpoints() returns overapproximation to polygon. Filter out points outside triangle.

  return bs;
}


vector< pair<node_t,coord3d> > interpolate_triangle(const tri_t& t, const Folding& F, const vector<coord3d>& points)
{
  vector<node_t>  t_nodes = {F.outline[t[0]].second,F.outline[t[1]].second,F.outline[t[2]].second};
  vector<coord3d> t_xs    = {points[t_nodes[0]],points[t_nodes[1]],points[t_nodes[2]]};
  vector<pair<node_t,coord3d>> lambdas = barycentric_coordinates(t,F);
  vector<pair<node_t,coord3d>> xs(lambdas.size());

  cerr << "t_nodes = " << t_nodes << ";\n"
       << "t_xs    = " << t_xs    << ";\n";

  for(int i=0;i<lambdas.size();i++){
    const coord3d &lambda(lambdas[i].second);
    xs[i].first  = lambdas[i].first;
    xs[i].second = t_xs[0]*lambda[0] + t_xs[1]*lambda[1] + t_xs[2]*lambda[2];
  }

  return xs;
}


int main(int ac, char **av)
{
  assert(ac >= 13);

  int N = strtol(av[1],0,0);
  vector<int> RSPI(12);
  for(int i=0;i<12;i++) RSPI[i] = strtol(av[i+2],0,0)-1;

  vector<int> spiral(N/2+2, 6);
  for (int i = 0; i < 12; i++){
    spiral[RSPI[i]] = 5;
  }

  Triangulation dual(spiral);
  dual = dual.sort_nodes();
  
  dual.layout2d = dual.tutte_layout();
  Unfolding uf(dual);

  ofstream output("output/C"+to_string(N)+"-unfold.m");
  output << "dual    = " << dual << ";\n"
	 << "outline = " << uf.outline << ";\n"
	 << "arcs   = " << get_keys(uf.arc_coords) << ";\n"
	 << "arcpos = " << get_values(uf.arc_coords) << ";\n";
  output.close();

  vector<tri_t> triangles = Tris(dual);

  Unfolding UF = uf.straighten_lines();

  Folding F(UF*Eisenstein{1,0});
  Triangulation f(F.fold());

  output.open("output/C"+to_string(N)+"-star.m");

  output << "dual    = " << dual << ";\n"
	 << "outline = " << UF.outline << ";\n"
	 << "arcs   = " << get_keys(UF.arc_coords) << ";\n"
	 << "arcpos = " << get_values(UF.arc_coords) << ";\n"
	 << "folded   = " << f << ";\n";

  output.close();

  ofstream latex_output("output/C"+to_string(N)+"-unfold.tex");
  latex_output << "\\newcommand{\\figunfolding}{"<<uf.to_latex(1,0,2) << "}\n\n";
  latex_output << "\\newcommand{\\figstar}{"<<UF.to_latex(1,0,2) << "\n" << "\\begin{pgfonlayer}{bg}\n"<<latextris(UF,triangles) << "\\end{pgfonlayer}}\n\n";
  latex_output.close();

  Polyhedron P(fullerene_dual_polyhedron(dual));

  vector<coord3d> pentagon_points(&P.points[0],&P.points[12]);
  
  // Triangulation C12d(FullereneGraph::C20().dual_graph());
  // Polyhedron CH(C12d,pentagon_points);
  // CH = CH.convex_hull();

  // output.open("output/C"+to_string(N)+"-polyhedron.m");
  // output << "P  = " << P << ";\n"
  //   	 << "CH = " << CH << ";\n";
  // output.close();


  vector<vector<pair<node_t,coord3d>>> barycentrics(triangles.size());
  vector<vector<pair<node_t,coord3d>>> interpolated_points(triangles.size());
  vector<coord3d> new_points(f.N);

  for(int i=0;i<triangles.size();i++) {
    barycentrics[i] = barycentric_coordinates(triangles[i],F);
    interpolated_points[i] = interpolate_triangle(triangles[i],F,pentagon_points);

    cout << "barycentrics = " << barycentrics[i] << endl;
    cout << "interpolated_points = " << interpolated_points[i] << endl;


    for(auto &p: interpolated_points[i]){
      if((new_points[p.first]-coord3d{0,0,0}).norm() < 1e-8){
	cerr << "new_points["<<p.first<<"] = " << new_points[p.first] << " or " << p.second << "\n";
      }
      new_points[p.first] = p.second;
    }
  }


  
  Polyhedron Pfolded(f,new_points);

  string basename("unfold+3D-"+to_string(N));  
  Polyhedron::to_file(P,"output/"+basename+".mol2");
  Polyhedron::to_file(P,"output/"+basename+"-dual.mol2");

  return 0;
}
