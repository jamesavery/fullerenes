#include "fullerenes/unfold.hh"
#include <list>
#include <iterator>
/************************************************************/ 
/************************************************************/ 
/*                UNFOLDING IMPLEMENTATION                  */
/************************************************************/ 
/************************************************************/ 


// This function unfolds a triangulation and lays it out on an equilateral
// triangular grid, such that if one were to cut along the outline
// and glue together the nodes with the same labels, one would obtain
// again the original fullerene dual.
//
// Preconditions: Triangles are oriented consistently, i.e. CW or CCW.
// TODO: A more intricate directed edge selection scheme could lead to
// more compact unfoldings (i.e. with more interior points).
void Unfolding::unfold(const Triangulation& T, arc_t first_arc)
{
  size_t Nf = T.N, N = 2*(Nf-2);

  vector<arc_t> arc_queue_buf(6*N);
  Deque<arc_t> arc_queue(arc_queue_buf);
  // auto stay_close_to_center =
  //   [&](arc_t uv, arc_t st) {
  //     Eisenstein ux, vx, sx, tx;
  //     tie(ux,vx) = arc_coords.at({uv.second,uv.first});
  //     tie(sx,tx) = arc_coords.at({st.second,st.first});

  //     if(ux.norm2() + vx.norm2() < sx.norm2() + tx.norm2()) return true;
  //     if(ux.norm2() + vx.norm2() > sx.norm2() + tx.norm2()) return false;
  //     return uv < st;
  //   };
  // std::set<arc_t, decltype(stay_close_to_center)> arc_queue(stay_close_to_center);

  vector<bool> tri_done(N);
  vector<bool> face_done(Nf);

  Eisenstein zero(0,0), veci(1,0), vecj(0,1);

  // Initialize helper structures
  auto tris = T.triangles();
  for(int U=0;U<N;U++){
    const tri_t &t = tris[U];

    for(int i=0;i<3;i++){
      node_t u = t[i], v = t[(i+1)%3];
      arc_to_tri_id[{u,v}] = U;
    }
  }  
  
  auto place_triangle = [&](arc_t uv, Eisenstein ux, Eisenstein vx) {
    node_t u=uv.first, v = uv.second,  w = T.next_on_face(u,v);
    if(w<0){
      cerr << tri_t{u,v} << " is not part of a triangle in " << T << endl;
      abort();
    }
    
    Eisenstein wx = ux + (vx-ux).nextCCW();

    //    cout << "# Attempting to place " << vector<node_t>{{u,v,w}} << " at " << vector<Eisenstein>{{ux,vx,wx}} << endl;
    size_t uv_tri_id = arc_to_tri_id.at(arc_t{u,v});
    
    if(!tri_done[uv_tri_id]){
	tri_done[uv_tri_id] = true;
	
	arc_coords[{u,v}] = {ux,vx};
	arc_coords[{v,w}] = {vx,wx};
	arc_coords[{w,u}] = {wx,ux};
	//	cout << "   Success.\n";//, arc_coords = " << get_keys(arc_coords) << " |-> " << get_values(arc_coords) <<".\n";

	if(!tri_done[arc_to_tri_id.at({v,u})]) arc_queue.push_back({v,u});//arc_queue.insert({v,u});//
	if(!tri_done[arc_to_tri_id.at({u,w})]) arc_queue.push_back({u,w});//arc_queue.insert({u,w});//
	if(!tri_done[arc_to_tri_id.at({w,v})]) arc_queue.push_back({w,v});//arc_queue.insert({w,v});//
	//	cout << "   Queueing " << tri_t{w,v,u} << " -> " << arc_queue << endl;
    } else {
      //      cout << "   ALREADY PLACED, MAN!\n";
    }
  };

  
  // 1. Place first triangle.
  // If no first arc is given, take first arc emanating from node 0
  if(first_arc==arc_t{0,0}) first_arc = {0,T.nbrs(0)[0]};

  place_triangle(first_arc,zero,veci);

  // 2. Keep placing every unplaced triangle that is connected to the boundary
  //    of the already placed triangles until we are done.
  int i=0;
  while(!arc_queue.empty()){
    //    cout << "Next element from queue " << arc_queue << " is ";
    arc_t uv = arc_queue.pop_front(); // Next placeable directed edge
    //arc_t uv = *arc_queue.begin(); arc_queue.erase(uv); // Next placeable directed edge
    //    cout << uv << endl;
    
    node_t u(uv.first), v(uv.second), w(T.next_on_face(u,v));

    Eisenstein ux, vx;
    tie(vx,ux) = arc_coords.at({v,u});

    place_triangle(uv,ux,vx);
    // cout << "arcs"<<i << " = array(" << get_keys(arc_coords) << ");\n";
    // cout << "arcpos"<<(i++) << " = array(" << get_values(arc_coords) << ");\n";    
  }
  if(arc_coords.size() != N*3){
    cerr << "Number of arcs placed is " << arc_coords.size() << " != " << (3*N) << ": Incomplete unfolding.\n";
    abort();
  }
}


map<arc_t,Unfolding::arccoord_t> Unfolding::unfold(const Triangulation& G, const polygon& outline, const tri_t T0)
{
  map<arc_t, arccoord_t > arc_coords;
  std::stack<arc_t> workset;
  map<arc_t,bool> seen;	// TODO: flat array, neighbours-shape

  node_t u,v;
  Eisenstein xu, xv, xw, direction;

  // TODO: Use scanconverted polygon to test in constant time for whether a point is inside the polygon
  //    polygon::scanline S[3] = {outline.scanConvert(), (outline*CCW).scanConvert(), (outline*CW).scanConvert()};
   

  // The idea is simply to fill out the polygon from known positions. I.e., start with one triangle placed
  // explicitly, then grow the unfolding from the perimeter, but only placing triangles that are on the interior
  // of the polygon.
  // This is done in the following way:
  //  A workset 'W' of keeps track of arcs on the perimeter where the corresponding triangle can be placed.
  //  A rasterization 'S' of the polygon lets us look up in constant time whether an arc is on the inside of the polygon
  //  A boolean map 'seen' keeps track of which arcs we have already processed.
  //
  // An arc u->v is ready to be processed and placed into the work-set iff
  //  1. We have placed the reverse arc v->u (so its position is well-defined)
  //  2. !seen[{u->v}], i.e., we have not already placed v->u
  //  3. v->u is an internal edge to the polygon

  auto triangle_is_internal = [&](node_t u, node_t v) -> bool {
    //    printf("triangle_is_internal(%d,%d)",u,v);
    Eisenstein xu, xv, xw, direction;      
    tie(xv,xu) = arc_coords[{v,u}];      // Arc v->u has already been placed if u->v is in workset, giving us two coordinates
    direction = xv - xu;                 // Direction of arc u->v
    xw        = xu + direction.nextCW(); // The final coordinate of the triangle is found by turning one step CW

    //    cout << " - " << vector<Eisenstein>{{xu,xv,xw}} << " - ";
    bool internal= outline.point_included(xu) & outline.point_included(xv) & outline.point_included(xw);
    //    printf(" = %d\n",internal);
    return internal;
  };
    
  auto arc_can_be_processed = [&](node_t u, node_t v) -> bool {
    return !seen[{u,v}] && triangle_is_internal(u,v);
  };
    
  auto place_triangle = [&](const tri_t &T, const Eisenstein position[3]) {

    if(debug_flags & VERBOSE) cout << "Placing triangle " << T << " at " << vector<Eisenstein>{{position[0],position[1],position[2]}} << endl;
    for(int i=0;i<3;i++) seen[{T[i],T[(i+1)%3]}] = true;
    
    for(int i=0;i<3;i++){
      node_t      u = T[i],         v = T[(i+1)%3];
      Eisenstein xu = position[i], xv = position[(i+1)%3]; // Fix coordinates of each arc 

      arc_coords[{u,v}] = {xu,xv}; // Fix coordinates of each arc 

      if(arc_can_be_processed(v,u)) workset.push({v,u});             // ...and stack reverse arc on workset if it is to be processed
    }
  };
    

  // Now the algorithm becomes simply:
    // 1. Place the first triangle T0
    Eisenstein pos[3] = {outline.reduced_outline[0],outline.reduced_outline[0]+Eisenstein{0,1}, outline.reduced_outline[0]+Eisenstein{1,0}};
    place_triangle(T0,pos);

    // 2. Place the rest of the triangles
    while(!workset.empty()){
      tie(u,v) = workset.top(); workset.pop();
      if(arc_can_be_processed(u,v)){
	node_t w = G.next(u,v);
	
	tie(xv,xu) = arc_coords[{v,u}]; // Arc v->u has already been placed if u->v is in workset
	
	Eisenstein direction = xv - xu; // Direction of arc u->v
	assert(direction.norm2() == 1);		      
	
	Eisenstein xuvw[3] = {xu, xv, xu + direction.nextCW()};
	place_triangle({u,v,w}, xuvw);
      }
    }
    //    cout << "Done!\n";
    return arc_coords;    
}

// Given the output of unfold(), this function efficiently computes the polygon outlining
// the unfolded triangulation and returns it in clockwise order.
vector< pair<Eisenstein,node_t> > Unfolding::get_outline(const map<arc_t,Unfolding::arccoord_t>& arc_coords)
{
  map<Eisenstein,node_t>    label;
  map<Eisenstein,Eisenstein> next;

  // Collect the directed edges u->v whose positions do not coincide with the reverse edge v->u. 
  // These form the outline of the polygon. 
  for(map<arc_t,arccoord_t>::const_iterator i(arc_coords.begin()); i!= arc_coords.end(); i++){
    const arc_t &uv(i->first), vu(uv.second,uv.first);
    const arccoord_t &uvpos(i->second), vupos(arc_coords.find(vu)->second);

    if(uvpos != make_pair(vupos.second,vupos.first)){
      next[uvpos.first]   = uvpos.second; 
      label[uvpos.first]  = uv.first;
    }
  }
  
  // Now we're ready to find a CW walk through the polygon outline coordinates 
  // and to assign the corresponding nodes in the triangulation graph as labels.
  vector< pair<Eisenstein,node_t> > outline(next.size());
  Eisenstein nextpos = next.begin()->first;

  for(int i=0;i<outline.size();i++, nextpos = next[nextpos]){
      outline[i] = make_pair(nextpos,label[nextpos]);
  }

  // If the outline doesn't form a closed loop, something is wrong.
  assert(nextpos == next.begin()->first);

  // CCW to CW order
  reverse(outline.begin(),outline.end());
  return outline;
}

void Unfolding::transform_line(const Unfolding::arccoord_t& l1, const Unfolding::arccoord_t& l2,
			       Eisenstein& x0, Eisenstein& x0p, Eisenstein& w)
{
  Eisenstein Duv(l1.second-l1.first), Dvu(l2.first-l2.second), Tuvvu((Duv.complex_conj()*Dvu)/Dvu.norm2());

  x0  = l1.first;
  x0p = l2.second;
  w   = Tuvvu;
}

vector<pair<Eisenstein,node_t>> Unfolding::GCDreduce(const vector<pair<Eisenstein,node_t>>& outline)
{
  vector<Eisenstein> segments(outline.size());
  for(size_t i=0;i<outline.size();i++)
    segments[i] = outline[(i+1)%outline.size()].first - outline[i].first;

  Eisenstein d(Eisenstein::gcd(segments));
  for(size_t i=0;i<segments.size();i++) segments[i] = segments[i].div(d);

  vector<pair<Eisenstein,node_t>> new_outline(outline.size());
  new_outline[0] = outline[0];
  for(size_t i=0;i+1<outline.size();i++)
    new_outline[i+1] = make_pair(new_outline[i].first+segments[i], outline[i+1].second);

  return new_outline;
}

namespace {
  // @post  for every seam segment U->V on the boundary, succ[U*Ncones+V] is the label
  //        W of the segment V->W that follows it WITHIN their common cactus cycle
  //        (pred is the inverse); entries with no segment stay -1.
  // @inv   stk holds the cones on the current open path; pos[c] is c's stack index.
  // Decompose the closed cone-sequence into its simple cycles (unique, since the cone
  // graph is a cactus) by a stack walk: push cones; when the current cone is already on
  // the stack, the run back to it is one cycle. This is the starshape proof's "successor
  // in C_1", which differs from the outline order exactly at cut-vertices.
  void cone_cycle_links(const vector<int>& seq, int Ncones,
                        vector<int>& succ, vector<int>& pred){
    const int m = (int)seq.size();
    vector<int> stk, pos(Ncones,-1);                 // pos[label] = stack index, or -1
    auto set_cycle = [&](const vector<int>& cyc){
      const int k = (int)cyc.size();
      for(int t=0;t<k;t++){
        int U=cyc[t], V=cyc[(t+1)%k], W=cyc[(t+2)%k];
        succ[U*Ncones+V]=W; pred[V*Ncones+W]=U;
      }
    };
    for(int step=0; step<=m; step++){
      const int c = seq[step%m];
      if(pos[c]>=0){                                 // closed a cycle back to c
        const int p = pos[c];                        // capture: the pop loop clears pos[c]
        set_cycle(vector<int>(stk.begin()+p, stk.end()));
        while((int)stk.size()>p){ pos[stk.back()]=-1; stk.pop_back(); }
        if(step<m){ pos[c]=(int)stk.size(); stk.push_back(c); }  // c stays as the cut-vertex
      } else { pos[c]=(int)stk.size(); stk.push_back(c); }
    }
  }

  typedef list<pair<Eisenstein,node_t>>::iterator Onode_t;

  // The developed cone-segment graph awaiting straightening: the outline O (cone copies,
  // the eventual output), the cactus cycle structure (succ/pred), per-segment geometry
  // (XUV = endpoints, XUv/XVu = the two boundary arms), a stable handle UVindex into O
  // per segment, and the set of segments still to reduce.
  struct ConeSegments {
    list<pair<Eisenstein,node_t>>     O;
    vector<int>                       succ, pred;
    map<arc_t,Unfolding::arccoord_t>  XUV, XUv, XVu;
    map<arc_t,Onode_t>                UVindex;
    set<arc_t>                        workset;
  };

  // @pre    cones_first: every non-hexagon outline label is < Ncones
  // @post   result.O are the cone copies in outline order; result.{succ,pred,XUV,XUv,
  //         XVu,UVindex,workset} annotate each cone-to-cone seam segment.
  // @throws std::logic_error when the cones-first pre-condition is violated.
  // Develop the cone-only outline and its per-segment annotations (the starshape
  // proof's graph G, carrying geometry); the reduction below consumes this.
  ConeSegments develop_cone_segments(const vector<pair<Eisenstein,node_t>>& outline,
                                     const vector<int>& degrees, int Ncones){
    ConeSegments g;
    g.succ.assign(Ncones*Ncones,-1);
    g.pred.assign(Ncones*Ncones,-1);

    // Cone copies in outline order, with a stable handle (Onode) and outline index
    // (Oindex) for each -- the index locates the incident boundary arms below.
    vector<int>     Oindex;
    vector<Onode_t> Onode;
    for(int i=0;i<(int)outline.size();i++)
      if(degrees[outline[i].second] != 6){
        if(outline[i].second >= Ncones)
          throw std::logic_error("straighten_lines: cone label >= "+std::to_string(Ncones)+
                                 " -- the Unfolding must be cones-first (sort_flat_last-relabelled)");
        Oindex.push_back(i);
        g.O.push_back(outline[i]);
        Onode.push_back(std::prev(g.O.end()));
      }

    for(int i=0;i<(int)Onode.size();i++){
      int j = (i+1)%(int)Onode.size();
      int i1 = (Oindex[i]+1)%outline.size(), j1 = (Oindex[j]-1+outline.size())%outline.size();
      arc_t UV(Onode[i]->second,Onode[j]->second);
      g.workset.insert(UV);
      Eisenstein Ux(Onode[i]->first), vx(outline[i1].first), Vx(Onode[j]->first), ux(outline[j1].first);
      g.UVindex[UV] = Onode[i];
      g.XUV[UV] = Unfolding::arccoord_t(Ux,Vx);
      g.XUv[UV] = Unfolding::arccoord_t(Ux,vx);
      g.XVu[UV] = Unfolding::arccoord_t(Vx,ux);
    }

    vector<int> seq; for(auto& n: Onode) seq.push_back(n->second);
    cone_cycle_links(seq, Ncones, g.succ, g.pred);
    return g;
  }

  // @inv    g.workset holds exactly the seam segments not yet straightened.
  // @post   g.O is a simple cone-star polygon: straight cone-to-cone segments, the
  //         hexagon detours removed.
  // @throws std::logic_error when a reduction invariant is violated.
  // Straighten every bend by the starshape reduction: retire finished leaf 2-cycles,
  // then straighten one bend U->V->W onto a single segment U->W, where W is the CYCLE
  // successor of U->V. The repositioning is transform_line; V survives in g.O.
  void reduce_cone_star(ConeSegments& g, int Ncones){
    while(!g.workset.empty()){
      // Retire each 2-cycle {U->V, V->U} -- a finished leaf edge (cones stay in O).
      for(node_t U=0;U<Ncones;U++)
        for(node_t V=0;V<Ncones;V++)
          if(g.succ[U*Ncones+V]==U){
            g.succ[U*Ncones+V]=-1; g.pred[U*Ncones+V]=-1;
            g.succ[V*Ncones+U]=-1; g.pred[V*Ncones+U]=-1;
            g.workset.erase(arc_t(U,V)); g.workset.erase(arc_t(V,U));
          }
      if(g.workset.empty()) break;

      arc_t UV(*g.workset.begin());
      node_t U(UV.first), V(UV.second);
      int Wi = g.succ[U*Ncones+V];
      if(Wi<0) throw std::logic_error("straighten_lines: active segment has no cycle successor");
      node_t W = (node_t)Wi;
      if(U == W) throw std::logic_error("straighten_lines: 2-cycle survived leaf removal");

      arc_t VW(V,W), UW(U,W);
      Eisenstein x0, x0p, omega;
      Unfolding::transform_line(g.XUv[VW], reverse(g.XVu[UV]), x0, x0p, omega);
      Eisenstein Wxp = g.XUV[VW].second, Wx((Wxp-x0)*omega+x0p);

      // Geometry of the new straight segment U->W and its arms.
      Unfolding::arccoord_t Wuxp(g.XVu[VW]);
      g.XUV[UW] = Unfolding::arccoord_t(g.XUV[UV].first, Wx);
      g.XUv[UW] = g.XUv[UV];
      g.XVu[UW] = Unfolding::arccoord_t((Wuxp.first-x0)*omega+x0p, (Wuxp.second-x0)*omega+x0p);

      // Insert W right after U's copy (O(1); list handles stay valid).
      Onode_t it_U = g.UVindex[UV];
      g.O.insert(std::next(it_U), make_pair(Wx,W));
      g.UVindex[UW] = it_U;

      // ...T->U->V->W->Y... becomes ...T->U->W->Y... (V dropped from the cycle).
      int Y = g.succ[V*Ncones+W];
      int T = g.pred[U*Ncones+V];
      if(Y<0 || T<0) throw std::logic_error("straighten_lines: cycle neighbour missing");
      g.succ[U*Ncones+W]=Y;  g.pred[U*Ncones+W]=T;
      g.pred[W*Ncones+Y]=U;  g.succ[T*Ncones+U]=W;
      g.succ[U*Ncones+V]=-1; g.pred[U*Ncones+V]=-1;
      g.succ[V*Ncones+W]=-1; g.pred[V*Ncones+W]=-1;
      g.workset.erase(UV); g.workset.erase(VW);
      g.workset.insert(UW);
    }
  }
}

// Reduce the cone outline to a simple "cone-star" polygon -- straight cone-to-cone seam
// segments -- by the starshape proof's Phase-1 reduction (starshape.tex). The geometry
// is in transform_line; the combinatorial core is taking, for each bend U->V, the CYCLE
// successor V->W (the continuation of U->V around its cactus cycle).
//
// @pre    cones_first: every non-hexagon (cone) outline label is < 12
// @post   result.outline is a simple polygon whose vertices are the cones
// @throws std::logic_error when the cones-first pre-condition or a reduction invariant
//         is violated.
// TODO:   Preserve the graph; handle negative-curvature graphs.
Unfolding Unfolding::straighten_lines() const
{
  const int Ncones = 12;   // a fullerene has exactly 12 pentagon (degree-5) cones
  ConeSegments g = develop_cone_segments(outline, degrees, Ncones);
  reduce_cone_star(g, Ncones);
  return Unfolding(vector<pair<Eisenstein,node_t>>(g.O.begin(), g.O.end()));
}

string Unfolding::to_latex(int K, int L, int label_vertices,  bool draw_equilaterally, bool include_headers) const 
{
  string result;
  ostringstream latexfile(result);

  if(include_headers)
    latexfile << 
"\\documentclass{standalone}\n\
\\usepackage{tikz}\n\
\\begin{document}\n\
\\definecolor{darkgreen}{rgb}{.4, .7, .2}\n\
\\tikzstyle{outline}=[draw=black, ultra thick, fill opacity=.2, fill=darkgreen]\n\
\\tikzstyle{vertex}=[circle, draw, inner sep=0, fill=white, minimum width=4.00000mm]\n\
\n\
\\pgfmathsetmacro{\\xcoord}{cos(60)}\n\
\\pgfmathsetmacro{\\ycoord}{sin(60)}\n\
\n\
\\newcommand{\\drawEGrid}{\n\
    \\path[clip,preaction = {draw=black}] (0,0) -- (\\last,0) -- (\\cols,\\rows) --(\\first,\\rows) -- cycle;\n\
    \\draw (\\first,0) grid (\\last,\\rows);\n\
    \\foreach \\x in {1,2,...,\\total}\n\
        \\draw (-\\x,\\x*2) -- (\\x,0);\n\
}\n\
\\newcommand{\\drawRGrid}{\n\
    \\path[clip,preaction = {draw=black}] (0,0) -- (\\last,0) -- (\\cols,\\rows) --(0,\\rows) -- cycle;\n\
    \\draw (0,0) grid (\\cols,\\rows);\n\
    \\foreach \\x in {1,2,...,\\total}\n\
        \\draw (-\\x,\\x*2) -- (\\x,0);\n\
}\n\
\\begin{tikzpicture}\n\
";

  vector<Eisenstein> outline_gc(outline.size());
  for(int i=0;i<outline.size();i++) outline_gc[i] = outline[i].first * Eisenstein{K,L};

  // Extract (I,J)-bounds
  int imin = INT_MAX, imax = INT_MIN, jmin = INT_MAX, jmax = INT_MIN;
  for(int i=0;i<outline_gc.size();i++){
    Eisenstein x(outline_gc[i]);

    if(draw_equilaterally){
      coord2d X(x.coord());
      x.first  = floor(X.first);
      x.second = floor(X.second);
    }

    if(x.first < imin) imin = x.first;
    if(x.first > imax) imax = x.first;
    if(x.second < jmin) jmin = x.second;
    if(x.second > jmax) jmax = x.second;

  }
  Eisenstein gcmin(imin-1,jmin-1), gcmax(imax+1,jmax+1);
  if(draw_equilaterally){
    gcmin = Eisenstein(coord2d(gcmin.first-1,gcmin.second));
    gcmax = Eisenstein(coord2d(gcmax.first,gcmax.second+1));
  }

  
  // Define bounds
  Eisenstein D(gcmax-gcmin);
  latexfile << "\\newcommand*{\\cols}{"<<D.first<<"}\n"
	    << "\\newcommand*{\\rows}{"<<D.second<<"}\n"
	    << "\\newcommand*{\\total}{"<<(D.first+D.second)<<"}\n"
	    << "\\newcommand*{\\first}{-\\rows/2}\n"
	    << "\\newcommand*{\\last}{\\cols+\\rows/2}\n";

  // Draw E-grid
  latexfile << "\\bgroup\n";
  if(draw_equilaterally){
    // -- as equilateral triangles
    latexfile << "\\pgftransformcm{1}{0}{\\xcoord}{\\ycoord}{\\pgfpointorigin}\n"
	      << "\\drawEGrid{}\n";
  } else // -- or as a regular grid
    latexfile << "\\drawRGrid{}\n";

  // Place vertex labels according to scheme chosen in parameter 'label_vertices':
  latexfile << "\\foreach \\place/\\name/\\lbl in {";
  switch(label_vertices){
  case 0: break; // Don't label vertices at all.
  case 1:        // Only label non-hexagon vertices on polygon outline.
    if(label_vertices == 1)	
      for(int i=0;i<outline.size();i++) 
	if(degrees[outline[i].second] != 6) {
	const Eisenstein &IJ(outline_gc[i]-gcmin);

	const node_t &u(outline[i].second);
	latexfile << "{(" << IJ.first << "," << IJ.second << ")/"<<i<<"/"<<u<<(i+1<outline.size()?"},":"}");
      }
    break;
  case 2:        // Only label vertices on polygon outline.
    if(label_vertices == 2)	
      for(int i=0;i<outline.size();i++){
	const Eisenstein &IJ(outline_gc[i]-gcmin);

	const node_t &u(outline[i].second);
	latexfile << "{(" << IJ.first << "," << IJ.second << ")/"<<i<<"/"<<u<<(i+1<outline.size()?"},":"}");
      }
    break;
  case 3: // Label all original vertices, including internal ones
    {
      int i=0;
      for(const auto& uv_ij: arc_coords){
	const arc_t& uv(uv_ij.first);

	const arccoord_t& ij(uv_ij.second);

	node_t u = uv.first;
	Eisenstein I = ij.first;
	Eisenstein IJ(I*Eisenstein{K,L}-gcmin);

	//	latexfile << "{(" << IJ.first << "," << IJ.second << ")/"<<i<<"/"<<u<<(++it != arc_coords.end()? "},":"}");
      }
    }
    break;
  case 4: // Label ALL vertices, both old and new.
    break;
  }
  latexfile << "}\n"
	    << "\t \\node[vertex] (\\name) at \\place {\\lbl};\n\n";

  
  // Draw outline polygon
  latexfile << "\\begin{pgfonlayer}{bg}\n"
	    << "\\draw[outline] (0.center) \\foreach \\i in {1,...,"<<(outline.size()-1)<<"}{ -- (\\i.center) } -- cycle;\n\n";
  
  latexfile << "\\end{pgfonlayer}\n"
	    << "\\egroup\n\n";

  if(include_headers) 
    latexfile << "\\end{tikzpicture}\n"
	      << "\\end{document}\n";
  
  return latexfile.str();
}



// static vector<Unfolding> generate_all_unfoldings(const Triangulation& graph)
// {
//   vector<arc_t> workset;

//   return 
// }


struct unfolding_parent_state {
 map<arc_t,bool> arc_seen;
 vector<arc_t>   arc_boundary;
 int level;
};

unfolding_parent_state place_triangle(const Triangulation& G, const arc_t arc, const unfolding_parent_state &S)
{
  unfolding_parent_state Snext;
  node_t u = arc.first, v = arc.second, w = G.next_on_face(u,v);

  Snext.level    = S.level + 1;
  Snext.arc_seen = S.arc_seen;
  Snext.arc_seen[{u,v}] = true;
  Snext.arc_seen[{v,w}] = true;
  Snext.arc_seen[{w,u}] = true;

  Snext.arc_boundary.reserve(S.arc_boundary.size()); // Memory need for new boundary will be abound the same
  for(const arc_t& arc: S.arc_boundary){	     // Update boundary
    // 3 cases:
    // 1: Neither u->w or v->w are on the boundary => u->v is replaced by u->w->v
    // 2: u->v->w is a boundary segment            => u->v->w is replaced by u->w
    // 3: u->v and v->w exist separately on the B  => u->v is replaced by u->w
    
  }
  return Snext;
}

void generate_unfolding_subtree(const Triangulation& G, const unfolding_parent_state &S, vector<vector<arc_t>> arc_boundaries)
{
  if(S.level == G.N){ // If we have placed all N triangles, the unfolding is complete, and we can add it to our results
    arc_boundaries.push_back(S.arc_boundary);
    return;
  } else	      // Otherwise, recursively proceed through each unproccessed arc on the outline:
    for(arc_t arc: S.arc_boundary) if(!S.arc_seen.at(arc)) {
      unfolding_parent_state Snext = place_triangle(G,arc,S); 
      generate_unfolding_subtree(G,Snext,arc_boundaries);
    }
}

