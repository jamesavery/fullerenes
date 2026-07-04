#pragma once
// unfold.hh -- BYTE-FAITHFUL namespaced copy of the production library's
// Unfolding + Folding classes, from include/fullerenes/unfold.hh.
//
// Baseline copy for the fold-unfold sub-project: NO algorithm changes. Its only
// purpose is to reproduce production behaviour exactly so the bug can be
// instrumented here without touching the production library.
//
// Wiring:
//   - Everything is wrapped in `namespace foldunfold`. This is REQUIRED: the
//     production libfullerenes.so (which we keep linking) already exports global
//     Unfolding/Folding/polygon symbols; without a namespace we'd get ODR /
//     duplicate-symbol clashes.
//   - Unqualified `polygon` resolves to foldunfold::polygon (the local copy,
//     declared in fu_geometry.hh included below). The local polygon also clashes
//     by-name with the lib's global ::polygon that arrives transitively via
//     triangulation.hh, so the namespace is what disambiguates.
//   - Triangulation, Graph, Eisenstein, spiral_nomenclature, node_t, arc_t,
//     tri_t, edge_t, IDCounter, Deque, get_keys, get_values, sgn, wedge,
//     GRAPH_DMAX and the vector-scalar operators resolve to the GLOBAL library
//     types/templates (header-only / linked from libfullerenes) -- NOT copied.
//   - Eisenstein::unit[7] is the one Eisenstein symbol that is NOT header-only
//     (defined in the lib's eisenstein.cc); it stays linked from libfullerenes.
//   - polygon::draw_line is a STATIC member, called as polygon::draw_line(...)
//     in fold.cc; it stays a static member of foldunfold::polygon.

#include "fullerenes/triangulation.hh"
#include "fullerenes/eisenstein.hh"
#include "fullerenes/permutation.hh"   // Permutation (cones-first relabelling)
#include "fullerenes/geometry.hh"   // local foldunfold::polygon
#include "fullerenes/seam_step.hh"     // SeamAtlas, ConeRing, seam_rings (fold Steps 2+3)
#include <array>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>


// One recorded step of straighten_lines(), for step-by-step illustration. A
// non-null trace passed to straighten_lines is filled with one of these after
// the init and after each Step-1 cancellation / Step-2 bend-straighten.
struct StraightenStep {
  std::string kind;                              // "init" | "step1" | "step2"
  std::vector<arc_t>                  cancelled; // step1: the U<->V pairs removed this round
  int U=-1, V=-1, W=-1;                          // step2: the bend U->V->W (W = cycle successor)
  Eisenstein Wx{0,0};                            // step2: straightened position inserted for W
  int ins=-1;                                    // step2: O-index where W was inserted
  std::vector<std::pair<Eisenstein,node_t>> O;   // straightened outline AFTER this step
  std::vector<arc_t>                  active;    // workset (remaining arcs) AFTER this step
};

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

  // The whole fold/unfold runs in CONES-FIRST labelling: the (<=12, for a
  // fullerene) non-hexagon cones are relabelled to ids 0..11 (so straighten_lines'
  // 12x12 adjacency / 0..11 loops are valid, cone iteration is 0..11, and the
  // cone-only subgraph is the first 12 vertices). The relabelling is the canonical
  // TriangulationView::sort_flat_last. cone_perm[orig] = sort_flat_last id;
  // cone_perm.inverse()[id] recovers the original label.
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
  Unfolding(const Triangulation& G0, arc_t first_arc = {0,0}) : degrees(G0.N)
  {
    // Relabel cones-first via the canonical sort_flat_last (non-hexagon cones
    // 0..11, hexagons after) and work in that labelling throughout.
    cone_perm = G0.sort_flat_last();
    graph = G0;
    graph.apply_permutation(cone_perm);
    for(int u=0;u<graph.N;u++) degrees[u] = graph.degree(u);
    if(!(first_arc==arc_t{0,0})) first_arc = {cone_perm[first_arc.first], cone_perm[first_arc.second]};

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
  //    the corresponding node numbers.
  //
  // Degrees computed here are PROVISIONAL (whole incident unit triangles per
  // corner copy, exact centroid-in-3x-polygon test). A corner wedge can be
  // cut FRACTIONALLY by the seams (a sector triangle sliced at arbitrary
  // angle contributes its parts at two mate copies), so no local count is
  // exact -- the Folding constructor recomputes every corner label's degree
  // exactly via the glue-chained cone holonomy walk (seam_step.hh).
  // (The pre-LLM unit-NEIGHBOUR count was further off: the star unfolding's
  // source cone came out 16 and a leaf cone 3.)
  Unfolding(const vector< pair<Eisenstein,node_t> > &outline) : outline(outline) {
    Eisenstein x, y;
    node_t u,v;

    node_t N_outline = 0;
    for(auto o: outline) N_outline = max(N_outline,o.second);
    degrees = vector<int>(N_outline + 1, 0);  // +1: indices 0..N_outline inclusive

    vector<Eisenstein> corners3;
    for(auto o: outline) corners3.push_back(o.first * 3);
    polygon P3(corners3);
    for(size_t i=0;i<outline.size();i++){
      tie(x,u) = outline[i];
      tie(y,v) = outline[(i+1)%outline.size()];

      Eisenstein ua(1,0);
      for(int j=0;j<6;j++,ua = ua.nextCW()){
        const Eisenstein ub = ua.nextCW();
        if(P3.point_inside(x*3 + ua + ub)) degrees[u]++;   // sector triangle {x, x+ua, x+ub}
      }

      arc_coords[{u,v}] = {x,y};
    }
  }

  // Transform unfolding into a polygon with only straight lines between nodes of
  // degree != 6; i.e., degree 6 nodes can only be interior or lie on a straight line.
  // If trace is non-null, it is filled with the step-by-step record (no effect
  // on the returned Unfolding).
  Unfolding straighten_lines(std::vector<StraightenStep>* trace = nullptr) const;

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


// Complete a partial id assignment to a bijection of [0,n): pin[i] is item i's desired
// id, kept iff it is a distinct, in-range id; everything else takes the remaining ids in
// increasing order. pin[i] = -1 means "no preference". Shared by Folding's label-
// preserving compaction and by restore_original_labels (same retain-then-fill pattern).
inline vector<node_t> fill_free_ids(const vector<node_t>& pin){
  const node_t n = (node_t)pin.size();
  vector<node_t> out(n, -1);
  vector<bool>   used(n, false);
  vector<node_t> rest;
  for(node_t i=0;i<n;i++){
    const node_t t = pin[i];
    if(t >= 0 && t < n && !used[t]){ out[i]=t; used[t]=true; }
    else rest.push_back(i);
  }
  node_t next = 0;
  for(node_t i : rest){ while(used[next]) ++next; out[i]=next; used[next]=true; }
  return out;
}


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
    // nextid. N_orig = the tracked-vertex boundary: development labels < N_orig are
    // tracked originals (G.N for a Triangulation-derived unfolding; Ncones for a
    // cone-star). cone_perm carries the cones-first relabel so fold() can restore.
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

    // Compact canonical IDs to contiguous [0,N), but PRESERVE the labels that arrived
    // pre-specified (canonical < N_orig, i.e. set from arc_coords); mint fresh ids only
    // for the vertices folding created (canonical >= N_orig, the subdivision points).
    // Same retain-then-fill as restore_original_labels. Sorted iteration is deterministic.
    std::set<node_t> canon;
    for(const auto &kv: grid) canon.insert(same_as[kv.second]);
    const vector<node_t> canon_v(canon.begin(), canon.end());        // sorted ascending
    vector<node_t> pin(canon_v.size());
    for(size_t i=0;i<canon_v.size();i++)
      pin[i] = (canon_v[i] < N_orig) ? canon_v[i] : node_t(-1);       // keep pre-specified labels
    const vector<node_t> assign = fill_free_ids(pin);
    for(size_t i=0;i<canon_v.size();i++) compact[canon_v[i]] = assign[i];
    const node_t N = (node_t)canon_v.size();

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

    // Exact corner degrees + cycles (the ctor-3 degrees above are provisional).
    assemble_corners();
  }

  // Fold Step 2, exact: per outline-corner label, the plane-CCW neighbour
  // cycle (corner_cycles) and degree via the glue-chained cone holonomy walk
  // (seam_step primitives; defined in fold.cc). Called by the constructor.
  // corner_rings keeps each entry's development data (copy, direction, chart,
  // cell) -- fold() writes the RECIPROCAL slots at the cones' flat neighbours
  // from it.
  void assemble_corners();
  map<node_t, vector<node_t>> corner_cycles;
  map<node_t, ConeRing>       corner_rings;

  // Build the folded graph's oriented neighbour slots by scan-converting the
  // outline polygon in each of the 3 lattice edge-directions and drawing the
  // horizontal interior edges. neighbours[u][d] is u's neighbour in Eisenstein
  // unit-direction d; slots i_omega and i_omega+3 are the +/- of rotation i_omega.
  void connect(vector<array<node_t,6>>& neighbours);
  void connect_polygon(int i_omega, vector<array<node_t,6>>& neighbours);  // Scan-convert interior edges (rotated by omega).
  void connect(int i_omega, vector<array<node_t,6>>& neighbours);          // connect_polygon for one rotation.

  // Write arc p0--p1 (rotated back by omega_inv) into both endpoints' oriented
  // slots i_omega / i_omega+3. Idempotent -- a shared edge written twice is harmless.
  void connect_arc(vector<array<node_t,6>>& nb, Eisenstein p0, Eisenstein p1,
                   Eisenstein omega_inv, int i_omega) const;

  // Collect nodes in unfolding that will correspond to the same nodes in the folded graph.
  vector<node_t> identify_nodes(const IDCounter<Eisenstein>& grid, const vector< pair<Eisenstein,node_t>>& outline) const;
  // Exact gcd-geodesic alternative, for A/B validation only (see fold.cc).
  vector<node_t> identify_nodes_exact(const IDCounter<Eisenstein>& grid, const vector< pair<Eisenstein,node_t>>& outline) const;
  vector<node_t> outline_nodes() const;

  Triangulation fold();

  // Map the folded graph's compacted ids back to original labels: each tracked original
  // (canonical < N_orig) to cone_perm.inverse()[canonical], every other vertex to a
  // remaining id of [0,N). Applied to T and to the id-keyed state (final_grid, node_pos,
  // degrees) so the whole Folding shares one id space.
  void restore_original_labels(Triangulation& T);

  // Output
  string to_latex(int K=1, int L=0,int label_vertices=1, bool draw_equilaterally=true, bool include_headers=false) const;
  string to_mathematica() const;
};

