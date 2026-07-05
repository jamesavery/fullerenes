#pragma once
#include "fullerenes/triangulation.hh"
#include "fullerenes/eisenstein.hh"
#include "fullerenes/permutation.hh"
#include <array>
#include <sstream>
#include <stack>
#include <stdexcept>


class Unfolding {
public:
  typedef pair<Eisenstein,Eisenstein> arccoord_t;

  typedef arccoord_t arc_coord_t;

  enum { VERBOSE = 1 };
  int debug_flags = 0;

  Triangulation graph;

  map<arc_t,arc_coord_t> arc_coords; 
  map<arc_t,size_t>      arc_to_tri_id;	     // Unique triangle containing arc
  
  vector< pair<Eisenstein,node_t> > outline; // Polygon outline in the Eisenstein plane. This is always initialized.
  vector<int> degrees;

  // Cones-first relabelling applied by the Triangulation constructor:
  // cone_perm[orig]=new (non-hexagon cones get ids 0..11, hexagons after),
  // cone_perm.inverse()[new]=orig recovers the original labels (Folding::fold()
  // does this restore). Empty for the outline-based constructors (no relabel).
  Permutation cone_perm;

  vector<vector<Eisenstein>> tri_coords() const
  {
    vector<tri_t> triangles = graph.compute_faces_oriented();
    vector<vector<Eisenstein>> coords(triangles.size());

    
    for(int i=0;i<triangles.size();i++){
      vector<Eisenstein> x(6);
      tri_t t = triangles[i];
      for(int j=0;j<3;j++) tie(x[2*j],x[2*j+1]) = arc_coords.at({t[j],t[(j+1)%3]});
      for(int j=0;j<2;j++)
	if(x[2*j+1] != x[2*(j+1)]){
	  std::ostringstream msg;
	  msg << "Unfolding::tri_coords: broken triangle " << t << " coords " << x;
	  throw std::logic_error(msg.str());
	}

      coords[i] = {{x[0],x[2],x[4]}};
    }
    return coords;
  }
  
  // There are 3 ways to create an "unfolding":
  // 
  // 1. Provide a triangulation of the sphere, i.e., the dual of a planar cubic graph.
  Unfolding(const Triangulation& G, arc_t first_arc = {0,0}) : degrees(G.N)
  {
    // Relabel cones-first (sort_flat_last: non-hexagon cones get ids 0..11, hexagons
    // after) so straighten_lines' cone-indexed tables are valid and cone iteration is
    // 0..11. cone_perm[orig]=new; cone_perm.inverse()[new]=orig recovers the original
    // labels, which Folding::fold() restores on the folded graph.
    cone_perm = G.sort_flat_last();
    graph = G;
    graph.apply_permutation(cone_perm);

    // Store degrees of each (relabelled) node
    for(int u=0;u<graph.N;u++) degrees[u] = graph.degree(u);

    // A caller-supplied first_arc is given in original labels; map it cones-first.
    // The {0,0} default is a sentinel (handled inside unfold), not a real arc.
    if(!(first_arc==arc_t{0,0}))
      first_arc = {cone_perm[first_arc.first], cone_perm[first_arc.second]};

    unfold(graph, first_arc);
    outline = get_outline(arc_coords);
  }

  // 2. Provide a triangulation of the sphere 'G', the final Eisenstein-coordinate cutout 'outline' (in CW order),
  //    and a starting triangle 'T0' (in CW order)
  //    Assumption: Only full triangles in outline, and T0=[u,v,w] is such that the coordinates of:
  //    u is outline[0], v is outline[0]+(1,0), and w is outline[0]+(0,1).
  Unfolding(const Triangulation& G, const polygon& outline_polygon, const tri_t T0) : graph(G), degrees(G.N) 
  {   
    for(int u=0;u<G.N;u++) degrees[u] = G.degree(u);     // Store degrees of each node
    arc_coords = unfold(G,outline_polygon,T0);
    //    outline    = get_outline(arc_coords);
  }
  
  // 2. Provide a simple polygon in the Eisenstein plane, with polygon vertices annotated with
  //    the corresponding node numbers (a cone-star: the outline corners are the cones).
  //
  // @pre cone-star: every labelled outline vertex is a cone of degree <= 5. A cone is
  //      never internal, so the outline carries all of them; degree 6 lies on an edge,
  //      and degree >= 7 cannot develop into a simple polygon -- hence Ncones <= 12.
  Unfolding(const vector< pair<Eisenstein,node_t> > &outline_in) {
    // Relabel cones-first so straighten_lines' precondition holds. Build a dense
    // permutation of [0,Vmax]: the Ncones distinct outline labels -> 0..Ncones-1 (first-
    // appearance order), the rest -> Ncones..Vmax in any order. cone_perm[orig]=new;
    // cone_perm.inverse() recovers the original labels (Folding::fold() does the restore).
    node_t Vmax = 0;
    for(const auto& o: outline_in) Vmax = max(Vmax, o.second);
    const int n = Vmax + 1;
    cone_perm = Permutation(n);
    vector<bool> is_cone(n, false);
    int next_id = 0;
    for(const auto& o: outline_in)
      if(!is_cone[o.second]){ is_cone[o.second]=true; cone_perm[o.second]=next_id++; }  // cones -> 0..Ncones-1
    for(int id=0; id<n; id++)
      if(!is_cone[id]) cone_perm[id]=next_id++;                                          // rest -> Ncones..Vmax

    outline.reserve(outline_in.size());
    for(const auto& o: outline_in) outline.push_back({o.first, cone_perm[o.second]});    // relabelled outline

    // Calculate degrees of each node and store directed edge coordinates.
    Eisenstein x, y;
    node_t u,v;
    degrees = vector<int>(n, 0);

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
  //
  // @pre cones-first: the (12, for a fullerene) non-hexagon cones are labelled 0..11,
  //      so they index the Ncones-wide cactus tables directly. The Triangulation
  //      constructor ESTABLISHES this (via sort_flat_last) and stores cone_perm to
  //      recover the original labels; the outline-based constructors do NOT, so a
  //      caller using those must ensure cones-first itself. Violation throws
  //      std::logic_error (see src/c++/unfold.cc). The result carries cone_perm.
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
  // Note: both outline AND arc_coords must be transformed consistently.
  Unfolding& operator *= (const Eisenstein& y){
    for(int i=0;i<outline.size();i++) outline[i].first *= y;
    for(auto& kv : arc_coords) { kv.second.first *= y; kv.second.second *= y; }
    return *this;
  }
  Unfolding operator*(const Eisenstein& y) const { Unfolding U(*this); return (U *= y); }

  Unfolding& operator += (const Eisenstein& y){
    for(int i=0;i<outline.size();i++) outline[i].first += y;
    for(auto& kv : arc_coords) { kv.second.first += y; kv.second.second += y; }
    return *this;
  }
  Unfolding operator+(const Eisenstein& y) const { Unfolding U(*this); return (U += y); }

  // Unfolding& operator /= (const Eisenstein& y){ for(int i=0;i<outline.size();i++) outline[i].first /= y; return *this;  }
  // Unfolding operator/(const Eisenstein& y) const { Unfolding U(outline); return (U /= y); }

  static void transform_line(const arccoord_t& l1, const arccoord_t& l2, Eisenstein& x0, Eisenstein& x0p, Eisenstein& w);

  // Divide each outline segment by the Eisenstein GCD of all segments,
  // producing a reduced outline suitable for folding back to the original graph.
  static vector<pair<Eisenstein,node_t>> GCDreduce(const vector<pair<Eisenstein,node_t>>& outline);



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
  map<node_t, node_t>               compact;     // canonical grid id -> contiguous compacted id (0..M-1)
  vector< pair<Eisenstein,node_t> > outline;
  map<arc_t,arccoord_t>             arc_coords;  // copied from U; for the cone wedge-assembly
  vector<int>                       degrees;     // per compacted node id (cones: degree != 6)

  // Label-restore state (from the Unfolding): cone_perm[orig]=new (empty => no relabel
  // known; may be larger than N_orig, e.g. a cone-star's dense [0,Vmax] permutation),
  // N_orig = the tracked-vertex boundary (dev labels < N_orig are tracked originals).
  // fold() maps the folded graph back to original labels (tracked originals via
  // cone_perm.inverse(), every other vertex to a remaining id of [0,N)).
  Permutation cone_perm;
  node_t      N_orig = 0;
  
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
      if((size_t)(u+1) > grid.nextid) grid.nextid = u+1;
    }

    // Original-vertex count, captured BEFORE the GC-subdivision points below bump
    // nextid. For a Triangulation-derived Unfolding this equals cone_perm.size().
    // N_orig = the tracked-vertex boundary: development labels < N_orig are tracked
    // originals. grid.nextid here (after the arc-loop, before GC points bump it) is
    // exactly that count -- G.N for a Triangulation-derived unfolding (arc_coords
    // covers every vertex), Ncones for a cone-star outline (only the cones).
    cone_perm = U.cone_perm;
    N_orig    = (node_t)grid.nextid;

    // Build node id -> Eisenstein coordinate lookup table
    grid.reverse = vector<Eisenstein>(grid.nextid,0);
    for(auto kv: grid){
      tie(xu,u) = kv;
      grid.reverse[u] = xu;
    }

    // Next, fill in uninitialized grid points
    // (in case only outline arcs are filled in)
    const vector<Eisenstein> allpoints(P.allpoints());
    for(auto &x: allpoints) grid.insert(x);

    if(debug_flags & WRITE_FILE){
      cout << "grid_keys = "   << get_keys(grid) << endl;
      cout << "grid_values = " << get_values(grid) << endl;
    }

    // Now reduce redundant grid to the real node names
    same_as = identify_nodes(grid, outline);

    // Compact canonical IDs to be contiguous (0..N-1).
    // Merges and unused IDs can leave gaps in the canonical numbering.
    node_t N = 0;
    for(auto &kv: grid){
      node_t canonical = same_as[kv.second];
      if(compact.find(canonical) == compact.end())
        compact[canonical] = N++;
    }

    for(auto &kv: grid){
      tie(xu,u) = kv;
      final_grid[xu] = compact[same_as[u]];
    }

    node_pos = vector<vector<Eisenstein>>(N); // All coordinates of each node
    for(const auto &kv: grid){
      node_t cid = compact[same_as[kv.second]];
      node_pos[cid].push_back(kv.first);
    }

    // Keep arc_coords and per-(compacted)-node degree for the cone assembly.
    arc_coords = U.arc_coords;
    degrees = vector<int>(N, 6);
    for(const auto &kv: grid){
      node_t orig = kv.second;
      if(orig >= 0 && orig < (node_t)U.degrees.size())
        degrees[compact[same_as[orig]]] = U.degrees[orig];
    }

    if(debug_flags & WRITE_FILE) cout << "final_grid_values = " << get_values(final_grid) << endl;
  }

  // Build the folded graph's oriented neighbour slots by scan-converting the
  // outline polygon in each of the 3 lattice edge-directions and drawing the
  // horizontal edges (interior + across the seam). neighbours[u][d] is u's
  // neighbour in Eisenstein unit-direction d; slots i_omega and i_omega+3 are
  // the +/- of rotation i_omega.
  void connect(vector<array<node_t,6>>& neighbours);
  void connect_polygon(int i_omega, vector<array<node_t,6>>& neighbours);  // Scan-convert interior edges (rotated by omega).
  void connect_cross(int i_omega, vector<array<node_t,6>>& neighbours);    // Connect edges across the seam.
  void connect(int i_omega, vector<array<node_t,6>>& neighbours);          // Both of the above, for one rotation.

  // Assemble a cone's CCW neighbour cycle: chain its wedges via the outline
  // seam-segment matches, glue each wedge into one frame with transform_line,
  // read the unit-distance neighbours from final_grid in each wedge's sector,
  // and round about U. Only for the cones (degree != 6). seg_index maps an
  // outline segment's (srcLabel,tgtLabel) to its outline index (for finding the
  // reverse-label mate that names the adjacent wedge).
  vector<node_t> assemble_cone(node_t cid, const map<pair<node_t,node_t>,size_t>& seg_index) const;

  // Collect nodes in unfolding that will correspond to the same nodes in the folded graph.
  vector<node_t> identify_nodes(const IDCounter<Eisenstein>& grid, const vector< pair<Eisenstein,node_t>>& outline) const;
  vector<node_t> outline_nodes() const;

  // Relabel the folded graph T (and the id-keyed Folding state: final_grid, node_pos,
  // degrees) from compacted ids to ORIGINAL labels. Each PRESENT tracked original (its
  // class's canonical grid id < N_orig) is restored to its original id
  // cone_perm.inverse()[canonical]; every other vertex -- untracked originals OR
  // transform-introduced ones (e.g. Goldberg-Coxeter subdivision) -- takes one of the
  // remaining ids of [0,N). Correct for any PARTIAL subset of originals present, which
  // is the usual case after a transform. Single-shot: finalizes the id space (fold()
  // calls it once).
  // @pre  cone_perm.size() >= N_orig (covers the tracked range) or empty (identity restore)
  // @post the relabel is a permutation of [0,N) -- T stays a valid oriented triangulation
  //       -- and every present tracked original keeps its original id.
  void restore_original_labels(Triangulation& T);

  Triangulation fold();

  // Output
  string to_latex(int K=1, int L=0,int label_vertices=1, bool draw_equilaterally=true, bool include_headers=false) const;
  string to_mathematica() const;  
};
