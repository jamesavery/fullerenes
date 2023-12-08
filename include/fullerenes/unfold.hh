#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/eisenstein.hh"
#include <stack>


class Unfolding {
public:
  typedef pair<Eisenstein,Eisenstein> arccoord_t;
  
  typedef arccoord_t arc_coord_t;  

  Triangulation graph;

  map<arc_t,arc_coord_t> arc_coords; 
  map<arc_t,size_t>      arc_to_tri_id;	     // Unique triangle containing arc
  
  vector< pair<Eisenstein,node_t> > outline; // Polygon outline in the Eisenstein plane. This is always initialized.
  vector<int> degrees;

  vector<vector<Eisenstein>> tri_coords()
  {
    vector<tri_t> triangles = graph.compute_faces_oriented();
    vector<vector<Eisenstein>> coords(triangles.size());

    
    for(int i=0;i<triangles.size();i++){
      vector<Eisenstein> x(6);
      tri_t t = triangles[i];
      for(int j=0;j<3;j++) tie(x[2*j],x[2*j+1]) = arc_coords.at({t[j],t[(j+1)%3]});
      for(int j=0;j<2;j++)
	if(x[2*j+1] != x[2*(j+1)]){
	  cerr << "Broken triangle: " << t << " has coordinates " << x << ". :'-(\n";
	  abort();
	}

      coords[i] = {{x[0],x[2],x[4]}};
    }
    return coords;
  }
  
  // There are 3 ways to create an "unfolding":
  // 
  // 1. Provide a triangulation of the sphere, i.e., the dual of a planar cubic graph.
  Unfolding(const Triangulation& graph, arc_t first_arc = {0,0}) : graph(graph), degrees(graph.N) 
  {
    // Store degrees of each node   
    for(int u=0;u<graph.N;u++) degrees[u] = graph.degree(u);
    
    unfold(graph, first_arc);
    outline = get_outline(arc_coords); 
  }

  // 2. Provide a triangulation of the sphere 'G', the final Eisenstein-coordinate cutout 'outline' (in CW order),
  //    and a starting triangle 'T0' (in CW order)
  //    Assumption: Only full triangles in outline, and T0=[u,v,w] is such that the coordinates of:
  //    u is outline[0], v is outline[0]+(1,0), and w is outline[0]+(0,1).
  Unfolding(const Triangulation& G, const polygon& outline_polygon, const tri_t T0) : graph(G), degrees(G.N) 
  {   
    for(int u=0;u<G.N;u++) degrees[u] = G.neighbours[u].size();     // Store degrees of each node
    arc_coords = unfold(G,outline_polygon,T0);
    //    outline    = get_outline(arc_coords);
  }
  
  // 2. Provide a simple polygon in the Eisenstein plane, with polygon vertices annotated with
  //    the corresponding node numbers.
  Unfolding(const vector< pair<Eisenstein,node_t> > &outline) : outline(outline) {
    // Calculate degrees of each node and store directed edge coordinates.
    Eisenstein x, y;
    node_t u,v;

    node_t N_outline = 0;
    for(auto o: outline) N_outline = max(N_outline,o.second);
    degrees = vector<int>(N_outline,0);
    
    polygon P(get_keys(outline));
    for(int i=0;i<outline.size();i++){
      tie(x,u) = outline[i];
      tie(y,v) = outline[(i+1)%outline.size()];


      Eisenstein unit(1,0);
      for(int j=0;j<6;j++,unit = unit.nextCW())
	if(P.point_inside(x+unit)) degrees[u]++;
	
      arc_coords[{u,v}] = {x,y};
    }
  }

  // Transform unfolding into a polygon with only straight lines between nodes of
  // degree != 6; i.e., degree 6 nodes can only be interior or lie on a straight line.
  Unfolding straighten_lines() const;

  // This function unfolds a triangulation and lays it out on an equilateral
  // triangular grid, such that if one were to cut along the outline
  // and glue together the nodes with the same labels, one would obtain
  // again the original fullerene dual.
  //
  // Preconditions: Triangles are oriented consistently, i.e. CW or CCW.
  void unfold(const Triangulation& G, const arc_t first_arc={0,0});
  map<arc_t,arccoord_t> unfold(const Triangulation& G, const polygon& outline, const tri_t T0);
  
  // Compute outline in CW order of the map returned from unfold().
  static vector< pair<Eisenstein,node_t> > get_outline(const map<arc_t,arccoord_t>& arc_coords);

  // Simple transformations in the Eisenstein plane.
  Unfolding& operator *= (const Eisenstein& y){ for(int i=0;i<outline.size();i++) outline[i].first *= y; return *this;  }
  Unfolding operator*(const Eisenstein& y) const { Unfolding U(outline); return (U *= y); }

  Unfolding& operator += (const Eisenstein& y){ for(int i=0;i<outline.size();i++) outline[i].first += y; return *this;  }
  Unfolding operator+(const Eisenstein& y) const { Unfolding U(outline); return (U += y); }

  // Unfolding& operator /= (const Eisenstein& y){ for(int i=0;i<outline.size();i++) outline[i].first /= y; return *this;  }
  // Unfolding operator/(const Eisenstein& y) const { Unfolding U(outline); return (U /= y); }

  static void transform_line(const arccoord_t& l1, const arccoord_t& l2, Eisenstein& x0, Eisenstein& x0p, Eisenstein& w);



  // Output
  string to_latex(int K=1, int L=0,int label_vertices=1, bool draw_equilaterally=true, bool include_headers=false) const;
  string to_mathematica() const;
};




class Folding {
public:
  typedef Unfolding::arccoord_t arccoord_t;

  const polygon P;

  IDCounter<Eisenstein>             grid;
  map<Eisenstein, node_t>           final_grid; // 
  vector<vector<Eisenstein>>        node_pos;	 // Move grid, node_pos, and final_grid to Unfolding class?
  vector<node_t>                    same_as;
  vector< pair<Eisenstein,node_t> > outline;
  
  // Debug flags:
  enum {WRITE_FILE=1, DONT_ROTATE=2, DONT_CONNECT_POLYGON = 4, DONT_CONNECT_ACROSS = 8, DONT_IDENTIFY_NODES = 16, DO_NONPLANAR_LAYOUT = 32};
  int debug_flags;
  ostream &debug_file;

  Folding(const Unfolding& U, int debug_flags=0, ostream& debug_file = std::cerr) : 
    outline(U.outline), P(get_keys(U.outline)), debug_flags(debug_flags), debug_file(debug_file)
  {
    node_t u,v;
    Eisenstein xu, xv;

    // First transfer grid from unfolding's arc-coordinate grid
    for(auto kv: U.arc_coords){
      tie(u,v)   = kv.first;
      tie(xu,xv) = kv.second;

      grid[xu]         = u;
      grid.nextid = grid.nextid>=u? grid.nextid : (u+1);
    }

    // Build node id -> Eisenstein coordinate lookup table
    grid.reverse = vector<Eisenstein>(grid.nextid,0);
    for(auto kv: grid){
      tie(xu,u) = kv;
      grid.reverse[u] = xu;
    }
  
    
    
    //    cout << "grid.reverse = " << grid.reverse << endl;
    
    // Next, fill in uninitialized grid points
    // (in case only outline arcs are filled in)
    const vector<Eisenstein> allpoints(P.allpoints());
    for(auto &x: allpoints) grid.insert(x);

    cout << "grid_keys = "   << get_keys(grid) << endl;
    cout << "grid_values = " << get_values(grid) << endl;    

    // Now reduce redundant grid to the real node names
    // if(!(debug_flags & DONT_IDENTIFY_NODES)){
    node_t N = 0;
    same_as = identify_nodes(grid, outline);
    for(auto &kv: grid){
      tie(xu,u) = kv;
      final_grid[xu] = same_as[u];
      N = max(N,same_as[u]+1);
    }
    // cout << "grid:\n";
    // for(auto k: get_keys(grid))
    //   cout << "\t" << k << ": " << grid[k] << endl;

    // cout << "\nfinal_grid:\n";
    // for(auto k: get_keys(grid))
    //   cout << "\t" << k << ": " << final_grid[k] << endl;      

    node_pos = vector<vector<Eisenstein>>(N); // All coordinates of each node    
    node_pos[same_as[u]].push_back(xu);    

    cout << "final_grid_values = " << get_values(final_grid) << endl;    
  }

  void connect(neighbours_t &neighbours);    

  void connect_polygon(int i_omega, neighbours_t& neighbours);  // Connect scanlines of polygon rotated by w.
  void connect_cross(int i_omega, neighbours_t& neighbours);    // Connect horizontal lines across edges or polygon rotated by w.
  void connect(int i_omega, neighbours_t& neighbours);	        // Both of the above operations.

  // Collect nodes in unfolding that will correspond to the same nodes in the folded graph.
  vector<node_t> identify_nodes(const IDCounter<Eisenstein>& grid, const vector< pair<Eisenstein,node_t>>& outline) const;
  vector<node_t> outline_nodes() const;

  Triangulation fold();

  // Output
  string to_latex(int K=1, int L=0,int label_vertices=1, bool draw_equilaterally=true, bool include_headers=false) const;
  string to_mathematica() const;  
};
