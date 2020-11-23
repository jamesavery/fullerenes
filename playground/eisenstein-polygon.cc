#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/geometry.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/triangulation.hh"
#include <fstream>

vector< pair<Eisenstein,
node_t> > outline{
  {{0,14},  1 }, {{0,16},  13}, {{2,16},  5 }, {{4,16},   12},
  {{6,14},  3 }, {{8,12},  14}, {{10,12}, 2 }, {{8,14},   14}, {{8,16}, 3},
  {{10,16}, 12}, {{12,14}, 4 }, {{12,16}, 12}, {{14,16},  5 },
  {{16,14}, 13}, {{16,12}, 9 }, {{16,10}, 15},		   
  {{14,10}, 7 }, {{12,10}, 16}, {{12, 8}, 6 }, {{14, 8}, 16},{{16,6},7},
  {{16, 4}, 15}, {{14, 4}, 8 }, {{16, 2}, 15}, {{16, 0}, 9},
  {{14, 0}, 13}, {{12, 2}, 1 }, {{10, 4}, 17},		   
  {{10, 6}, 11}, {{10, 8}, 18}, {{8, 10}, 10}, {{8, 8}, 18},{{6, 8},11},
  {{4, 10}, 17}, {{4, 12}, 0 }, {{2, 12}, 17}
};

void draw_nodes_outline(const vector<pair<Eisenstein,node_t>> &outline, matrix<int> &grid, node_t& last_node)
{
  int x,y;
  
  for(auto o: outline){
    tie(x,y)  = o.first;
    node_t v  = o.second;
    last_node = max(last_node,v);

    grid(x,y) = v;
  }
}

void draw_nodes_polygon(const polygon &p, const Eisenstein w,
			  matrix<int> &grid, node_t& last_node)
{
  auto scan = (p*w).scanConvert(); // Coordinate-transformed polygon  
  Eisenstein iw = w.invertn();	   // Inverse coordinate transformation

  int X,Y;			   // Un-transformed coordinates
  
  for(int i=0;i<scan.xs.size();i++){
    auto &xs = scan.xs[i];    
    int y = scan.minY + i;

    for(int j=0;j<xs.size()/2;j++){    
      int x0 = xs[2*j];      
      int x1 = xs[2*j+1];

      for(int x=x0;x<=x1;x++){
	tie(X,Y) = Eisenstein{x,y}*iw;
	assert(X < grid.m && Y < grid.n);
	
	if(grid(X,Y) < 0) 	// Grid-point not yet initialized, new node
	  grid(X,Y) = ++last_node;
      }
    }
  }
}

matrix<int> node_grid_from_outline(const vector<pair<Eisenstein,node_t>> &outline,
				   node_t& last_node /* return value */)
{
  polygon P = get_keys(outline);

  int minx=INT_MAX, maxx=INT_MIN, miny=INT_MAX, maxy=INT_MIN, x,y;
  for(auto xy: P.reduced_outline){
    tie(x,y) = xy;

    minx = min(minx,x);
    maxx = max(maxx,x);
    miny = min(miny,y);
    maxy = max(maxy,y);
  }

  assert(minx == 0 && miny == 0); // TODO: Make it so
  
  matrix<int> node_grid(maxx+1,maxy+1,-1);
  last_node = -1;

  draw_nodes_outline(outline,   node_grid, last_node);
  draw_nodes_polygon(P, {1,0},  node_grid, last_node);
  draw_nodes_polygon(P, {0,1},  node_grid, last_node);
  draw_nodes_polygon(P, {1,-1}, node_grid, last_node);
  
  return node_grid;
}

vector<Eisenstein> cubic_inner_faces(polygon P, matrix<int> node_grid)
{
  vector<Eisenstein> inner_faces;
  
  auto scan = P.scanConvert();
  int X, Y;
	
  for(int i=0;i+1<scan.xs.size();i++){
    auto &xs = scan.xs[i];
    int y = scan.minY + i;

    for(int j=0;j<xs.size()/2;j++){    
      int x0 = xs[2*j];      
      int x1 = xs[2*j+1];

      for(int x=x0+1;x<x1;x++){
	Eisenstein xy = {x,y}, omega_n = {1,0}, omega = {0,1};
	bool inner_face = (node_grid(x,y)>=0);
	
	for(int i=0;i<6;i++, omega_n *= omega){
	  tie(X,Y) = xy + omega_n; // ith neighbour to (x,y)
	  inner_face &= (node_grid(X,Y) >= 0);
	}

	if(inner_face) inner_faces.push_back(xy);
      }
    }
  }
  return inner_faces;
}

// matrix<int> cubic_edges_from_outline(const vector<pair<Eisenstein,node_t>> &outline)
// {
//   int N;
//   matrix<int> node_grid = node_grid_from_outline(outline,N);
//   matrix<node_t> neighbours(N,3);

//   polygon P = get_keys(outline);
//   connect_nodes_polygon(P,0,{1,0}, node_grid, neighbours);
//   connect_nodes_polygon(P,1,{0,1}, node_grid, neighbours);
//   connect_nodes_polygon(P,2,{1,-1},node_grid, neighbours);
// }

void latex_scanconversion(const char *name, Eisenstein w, polygon p)
{
  Eisenstein iw = w.invertn();
  auto scan = (p*w).scanConvert();

  printf("\\newcommand{\\controlpoints%s}{",name);
  for(int i=0;i<scan.xs.size();i++){
    int y = scan.minY + i;
    for(int j=0;j<scan.xs[i].size();j++){
      int x = scan.xs[i][j];
      auto xy = Eisenstein{x,y}*iw;
      printf("%d/%d%s",xy.first,xy.second,(i+1==scan.xs.size() && j+1==scan.xs[i].size())?"":",");
    }
  }
  printf("}\n");

  printf("\\newcommand{\\dualedges%s}{",name);
  for(int i=0;i<scan.xs.size();i++){
    int y = scan.minY + i;
    for(int j=0;j<scan.xs[i].size()/2;j++){
      int x0 = scan.xs[i][2*j];      
      int x1 = scan.xs[i][2*j+1];
      for(int x=x0;x<x1;x++){
	auto u = Eisenstein{x,y}*iw;
	auto v = Eisenstein{x+1,y}*iw;
	printf("%d/%d/%d/%d%s",u.first,u.second,v.first,v.second,((i+1==scan.xs.size()) && (j+1==scan.xs[i].size()/2) && (x+1==x1))?"":",");
      }
    }
  }
  printf("}\n");
}

void outline_to_latex(vector< pair<Eisenstein,node_t> > outline)
{
  polygon p(get_keys(outline));
  
  printf("\\newcommand{\\outline}{");
  for(int i=0;i<outline.size();i++){
    auto o = outline[i];
    printf("%d/%d/%d%s",o.second,o.first.first,o.first.second,i+1==outline.size()?"":",");
  }
  printf("}\n\n");
  printf("\\newcommand{\\xzero}{0}\n"
	 "\\newcommand{\\yzero}{0}\n"
	 "\\newcommand{\\nx}{16}\n"
	 "\\newcommand{\\ny}{16}\n\n");
  
  latex_scanconversion("A",{1, 0},p);
  latex_scanconversion("B",{0, 1},p);
  latex_scanconversion("C",{1,-1},p);    

  int N = -1;
  matrix<int> node_grid          = node_grid_from_outline(outline,N);
  vector<Eisenstein> inner_faces = cubic_inner_faces(p, node_grid);

  printf("\\newcommand{\\innernodes}{");
  for(int i=0;i<inner_faces.size();i++)
    printf("%d/%d%s",inner_faces[i].first,inner_faces[i].second,i+1<inner_faces.size()?",":"");
  printf("}\n");

  printf("\\newcommand{\\innerhexagons}{");
  for(int i=0;i<inner_faces.size();i++){
    Eisenstein xy = inner_faces[i], omega = {0,1}, omega_n = {1,0};
    double cx, cy;

    for(int j=0;j<6;j++, omega_n *= omega){
      Eisenstein C = (xy) + (xy+omega_n) + (xy+omega_n*omega);

      printf("%.2f/%.2f%s",C.first/3.0, C.second/3.0, j+1<6?"/":"");
    }
    if(i+1<inner_faces.size()) printf(",\n");
  }
  printf("}\n");  

  
}


int main(int ac, char **av)
{
  Unfolding U(outline);
  Folding   F(U);

  return 0;
  ofstream latex_file("/tmp/test.tex");
  latex_file << U.to_latex(1,0,3,1,1) << endl;
  latex_file.close();
  
  Triangulation T = F.fold();

  cout << "T.spiral = " << T.get_general_spiral() << endl;
  
  return 0;
}
