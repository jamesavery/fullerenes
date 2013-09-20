#include "unfold.hh"
/************************************************************/ 
/************************************************************/ 
/*                UNFOLDING IMPLEMENTATION                  */
/************************************************************/ 
/************************************************************/ 


struct dedge_sort : public std::binary_function<dedge_t, dedge_t, bool>
{
    bool operator()(const dedge_t &x, const dedge_t &y) const
    {   
      int maxx = max(x.first,x.second), maxy = max(y.first,y.second);
      return maxx < maxy || (maxx == maxy && min(x.first,x.second) < min(y.first,y.second));
    }
};


// This function unfolds a triangulation and lays it out on an equilateral
// triangular grid, such that if one were to cut along the outline
// and glue together the nodes with the same labels, one would obtain
// again the original fullerene dual.
//
// Preconditions: Triangles are oriented consistently, i.e. CW or CCW.
// TODO: A more intricate directed edge selection scheme could lead to
// more compact unfoldings (i.e. with more interior points).
map<dedge_t,Unfolding::dedgecoord_t> Unfolding::unfold(const vector<tri_t> &triangulation)
{
#define set_dedge(u,v,ux,vx) {	          \
  dedge_t uv(u,v), vu(v,u);               \
  dedge_done[uv] = true;                  \
  workset.erase(uv);                      \
  dedge_position[vu] = make_pair(vx,ux);  \
  if(!dedge_done[vu])                     \
    workset.insert(vu);			  \
}

  // A single directed edge uniquely defines the third node in the oriented triangle
  map<dedge_t,node_t> nextNode;
  for(int i=0;i<triangulation.size();i++){
    const tri_t &t(triangulation[i]);
    for(int j=0;j<3;j++)
      nextNode[dedge_t(t[j],t[(j+1)%3])] = t[(j+2)%3];
  }

  map<dedge_t,bool> dedge_done;
  set<dedge_t, dedge_sort>   workset;
  map<dedge_t, dedgecoord_t > dedge_position
;
  map<Eisenstein,node_t> grid;
  Eisenstein zero(0,0), veci(1,0), vecj(0,1);

  // 1. Place first triangle. 
  tri_t t(triangulation[0]);

  set_dedge(t[0],t[1],zero,veci);
  set_dedge(t[1],t[2],veci,veci-vecj);
  set_dedge(t[2],t[0],veci-vecj,zero);

  // 2. Keep placing every unplaced triangle that is connected to the boundary
  //    of the already placed triangles until we are done.
  while(!workset.empty()){
    dedge_t uv(*workset.rbegin()); // Next placeable directed edge 
    // set_triangle(uv)
    node_t u(uv.first), v(uv.second), w(nextNode[uv]);

    dedgecoord_t uvpos(dedge_position[uv]);
    Eisenstein ux(uvpos.first), vx(uvpos.second), wx(ux+(vx-ux).nextCW());

    set_dedge(u,v,ux,vx);
    set_dedge(v,w,vx,wx);
    set_dedge(w,u,wx,ux);
  }
  return dedge_position;
}


// Given the output of unfold(), this function efficiently computes the polygon outlining
// the unfolded triangulation and returns it in clockwise order.
vector< pair<Eisenstein,node_t> > Unfolding::get_outline(const map<dedge_t,Unfolding::dedgecoord_t>& edgecoords)
{
  map<Eisenstein,node_t>    label;
  map<Eisenstein,Eisenstein> next;

  // Collect the directed edges u->v whose positions do not coincide with the reverse edge v->u. 
  // These form the outline of the polygon. 
  for(map<dedge_t,dedgecoord_t>::const_iterator i(edgecoords.begin()); i!= edgecoords.end(); i++){
    const dedge_t &uv(i->first), vu(uv.second,uv.first);
    const dedgecoord_t &uvpos(i->second), vupos(edgecoords.find(vu)->second);

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

void Unfolding::transform_line(const Unfolding::dedgecoord_t& l1, const Unfolding::dedgecoord_t& l2, Eisenstein& x0, Eisenstein& x0p, Eisenstein& w)
{
  Eisenstein Duv(l1.second-l1.first), Dvu(l2.first-l2.second), Tuvvu((Duv.invertn()*Dvu)/Dvu.norm2());

  x0  = l1.first;
  x0p = l2.second;
  w   = Tuvvu;
}

#include <unistd.h>
Unfolding Unfolding::straighten_lines() const 
{
  vector< int > Oindex;
  vector< pair<Eisenstein,node_t> > O;

  for(int i=0;i<outline.size();i++)	// Find non-hexagon node outline
    if((degrees.find(outline[i].second))->second != 6){
      Oindex.push_back(i);
      O.push_back(outline[i]);
    }
  
  // Directed graph defined by non-hexagon node outline
  vector<bool> A(12*12);	// Always at most 12 nodes of degree 5 or less
  set<dedge_t> workset;

  // Arc annotations
  map<dedge_t,pair<int,int> > UVindex;
  map<dedge_t,dedgecoord_t> XUV;
  map<dedge_t,dedgecoord_t> XUv;
  map<dedge_t,dedgecoord_t> XVu;

  for(int i=0;i<O.size();i++){
    int j = (i+1)%O.size();
    int i1 = (Oindex[i]+1)%outline.size(), j1 = (Oindex[j]-1+outline.size())%outline.size();

    dedge_t UV(O[i].second,O[j].second);

    workset.insert(UV);
    A[UV.first*12+UV.second] = true;

    Eisenstein Ux(O[i].first), vx(outline[i1].first), Vx(O[j].first), ux(outline[j1].first);

    UVindex[UV]  = make_pair(i,j);
    XUV[UV] = dedgecoord_t(Ux,Vx);
    XUv[UV] = dedgecoord_t(Ux,vx);
    XVu[UV] = dedgecoord_t(Vx,ux);
  }

  // Now repeatedly eliminate dedges by the following rules:
  cerr << "workset = " << workset << ";\n";
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

	  workset.erase(dedge_t(U,V));
	  workset.erase(dedge_t(V,U));
	}
    
    fprintf(stderr,"\nStep 2\n");
    cerr << "workset = " << workset << ";\n";

    // 2. When this step is reached, edges in workset are part of cycles of length >=3 and must be reduced. 
    // 2.1 Find first length-3 segment U->V->W
    dedge_t UV(*workset.begin());
    node_t U(UV.first), V(UV.second), W;
    for(W=0;W<12;W++) if(A[V*12+W]) break; 
    if(W==12){
      assert(workset.empty());
      break;
    }

    // 2.2 Transform W
    dedge_t VW(V,W);
    fprintf(stderr,"%d->%d->%d at ",U,V,W); cerr << UVindex[UV] << " and " << UVindex[VW] << endl;
    Eisenstein x0, x0p, omega;

    transform_line(XUv[VW], reverse(XVu[UV]), x0, x0p, omega);
    cout << "XUV = " << XUV[UV] << "; XWv = " << XUv[VW] << "; XUv = " << XVu[UV] << ";\n";
    Eisenstein Wxp = XUV[VW].second,  Wx((Wxp-x0)*omega+x0p);
    cout << "Wxp = " << Wxp << "; Wx = " << Wx << endl;


    // 2.3 Create annotation for new U->W arc
    dedge_t UW(U,W);
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
  }
  
  return Unfolding(O);
}
