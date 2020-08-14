#include "fullerenes/unfold.hh"
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
void Unfolding::unfold(const Triangulation& T, dedge_t first_arc)
{
  size_t Nf = T.N, N = 2*(Nf-2);

  Deque<arc_t> arc_queue(6*N);  
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
  for(int U=0;U<N;U++){
    const tri_t &t = T.triangles[U];

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
  if(first_arc==arc_t{0,0}) first_arc = {0,T.neighbours[0][0]};

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


map<dedge_t,Unfolding::dedgecoord_t> Unfolding::unfold(const Triangulation& G, const polygon& outline, const tri_t T0)
{
  map<arc_t, dedgecoord_t > arc_coords;
  std::stack<arc_t> workset;
  map<arc_t,bool> seen;	// TODO: flat array, neighbours-shape

  node_t u,v;
  Eisenstein xu, xv, xw, direction;

  // TODO: Use scanconverted polygon to test in constant time for whether a point is inside the polygon
  //    polygon::scanline S[3] = {outline.scanConvert(), (outline*CCW).scanConvert(), (outline*CW).scanConvert()};
   

  // The idea is simply to fill out the polygon from known positions. I.e., start with one triangle placed
  // explicitly, then grow the nfolding from the perimeter, but only placing triangles that are on the interior
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

    cout << "Placing triangle " << T << " at " << vector<Eisenstein>{{position[0],position[1],position[2]}} << endl;
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
vector< pair<Eisenstein,node_t> > Unfolding::get_outline(const map<arc_t,Unfolding::dedgecoord_t>& arc_coords)
{
  map<Eisenstein,node_t>    label;
  map<Eisenstein,Eisenstein> next;

  // Collect the directed edges u->v whose positions do not coincide with the reverse edge v->u. 
  // These form the outline of the polygon. 
  for(map<arc_t,dedgecoord_t>::const_iterator i(arc_coords.begin()); i!= arc_coords.end(); i++){
    const arc_t &uv(i->first), vu(uv.second,uv.first);
    const dedgecoord_t &uvpos(i->second), vupos(arc_coords.find(vu)->second);

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

void Unfolding::transform_line(const Unfolding::dedgecoord_t& l1, const Unfolding::dedgecoord_t& l2,
			       Eisenstein& x0, Eisenstein& x0p, Eisenstein& w)
{
  Eisenstein Duv(l1.second-l1.first), Dvu(l2.first-l2.second), Tuvvu((Duv.invertn()*Dvu)/Dvu.norm2());

  x0  = l1.first;
  x0p = l2.second;
  w   = Tuvvu;
}

// TODO: Preserve the graph!
#include <unistd.h>
Unfolding Unfolding::straighten_lines() const 
// ASSUMES: that nodes 1-12 are pentagons. 
// ASSUMES: that this is a triangulation of a fullerene
// TODO:    Work out how to do this for negative-curvature graphs.
{
  vector< int > Oindex;
  vector< pair<Eisenstein,node_t> > O;  // Straightened outline

  for(int i=0;i<outline.size();i++)	// Find non-hexagon node outline
    if(degrees[outline[i].second] != 6){
      Oindex.push_back(i);
      O.push_back(outline[i]);
    }
  
  // Directed graph defined by non-hexagon node outline
  vector<bool> A(12*12,false);	// Always at most 12 nodes of degree 5 or less
  set<arc_t> workset;

  // Arc annotations
  map<arc_t,pair<int,int> > UVindex;
  map<arc_t,dedgecoord_t> XUV;
  map<arc_t,dedgecoord_t> XUv;
  map<arc_t,dedgecoord_t> XVu;

  for(int i=0;i<O.size();i++){
    int j = (i+1)%O.size();
    int i1 = (Oindex[i]+1)%outline.size(), j1 = (Oindex[j]-1+outline.size())%outline.size();

    arc_t UV(O[i].second,O[j].second);

    workset.insert(UV);
    A[UV.first*12+UV.second] = true;

    Eisenstein Ux(O[i].first), vx(outline[i1].first), Vx(O[j].first), ux(outline[j1].first);

    UVindex[UV]  = make_pair(i,j);
    XUV[UV] = dedgecoord_t(Ux,Vx);
    XUv[UV] = dedgecoord_t(Ux,vx);
    XVu[UV] = dedgecoord_t(Vx,ux);
  }

  // Now repeatedly eliminate arcs by the following rules:
  while(!workset.empty()){
    
    fprintf(stderr,"Step 1\n");
    cerr << "workset = " << workset << ";\n";
    //  1. If u->v and v->u are both in the digraph, u->v matches up with v->u as
    //     desired, and we can remove the cycle u<->v from the digraph.
    for(node_t U=0;U<12;U++)
      for(node_t V=U+1;V<12;V++)
	if(A[U*12+V] && A[V*12+U]){
	  fprintf(stderr,"Found %d->%d and %d->%d, removing both\n",U,V,V,U);
	  A[U*12+V] = false;
	  A[V*12+U] = false;

	  workset.erase(arc_t(U,V));
	  workset.erase(arc_t(V,U));
	}
    
    fprintf(stderr,"\nStep 2\n");
    cerr << "workset = " << workset << ";\n";

    // 2. When this step is reached, edges in workset are part of cycles of length >=3 and must be reduced. 
    // 2.1 Find first length-3 segment U->V->W
    arc_t UV(*workset.begin());
    node_t U(UV.first), V(UV.second), W;
    for(W=0;W<12;W++) if(A[V*12+W]) break; 
    if(W==12){
      if(!workset.empty())
	fprintf(stderr,"straighten_lines: workset not empty, but no arcs to process.\n");
	
      break;
    }

    // 2.2 Transform W
    if(U != W){    
      arc_t VW(V,W);
      fprintf(stderr,"%d->%d->%d at ",U,V,W); cerr << UVindex[UV] << " and " << UVindex[VW] << endl;
      Eisenstein x0, x0p, omega;
      
      transform_line(XUv[VW], reverse(XVu[UV]), x0, x0p, omega);
      cerr << "XUV = " << XUV[UV] << "; XWv = " << XUv[VW] << "; XUv = " << XVu[UV] << ";\n";
      Eisenstein Wxp = XUV[VW].second,  Wx((Wxp-x0)*omega+x0p);
      cerr << "Wxp = " << Wxp << "; Wx = " << Wx << endl;


      // 2.3 Create annotation for new U->W arc
      arc_t UW(U,W);
      dedgecoord_t Wuxp(XVu[VW]), Wux(dedgecoord_t((Wuxp.first-x0)*omega+x0p, (Wuxp.second-x0)*omega+x0p));
      XUV[UW] = dedgecoord_t(XUV[UV].first,Wx);
      XUv[UW] = XUv[UV];
      XVu[UW] = Wux;

      // 2.4 Replace U->V by U->W->V in new outline O
      int Uindex(UVindex[UV].first);
      O.insert(O.begin()+Uindex+1, make_pair(Wx,W));

      // 2.5 Remove U->V and V->W from workset
      fprintf(stderr,"Removing %d->%d and %d->%d\n",U,V,V,W);
      A[U*12+V] = false;
      A[V*12+W] = false;
      workset.erase(UV);
      workset.erase(VW);

      // 2.6 Add U->W to workset

      fprintf(stderr,"Adding %d->%d\n",U,W);
      A[U*12+W] = true;
      workset.insert(UW);
    } else {
      fprintf(stderr,"Not implemented yet: %d->%d->%d\n",U,V,W);
      A[U*12+V] = false;
      A[V*12+W] = false;
      workset.erase({U,V});
      workset.erase({V,W});      
    }
  }
  
  return Unfolding(O);
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
	const dedgecoord_t& ij(uv_ij.second);

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


