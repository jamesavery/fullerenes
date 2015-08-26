#include "unfold.hh"

/************************************************************/ 
/************************************************************/ 
/*                 FOLDING IMPLEMENTATION                   */
/************************************************************/ 
/************************************************************/ 


vector<edge_t> Folding::connect_cross(const Eisenstein& w)
{
  vector<edge_t> edges;

  // Register reverse arcs
  map<dedge_t,dedgecoord_t> reverse_arc;
  vector< pair<Eisenstein,Eisenstein> > EC, ECp;

  for(int i=0;i<outline.size();i++){
    const Eisenstein xu(outline[i].first), xv(outline[(i+1)%outline.size()].first);
    const node_t u(outline[i].second), v(outline[(i+1)%outline.size()].second);
    reverse_arc[dedge_t(v,u)] = dedgecoord_t(xu,xv);
  }

  Eisenstein iw(w.invertn());
  for(int i=0;i<outline.size();i++){
    Eisenstein X0(outline[i].first*w), X1(outline[(i+1)%outline.size()].first*w);
    node_t     U(outline[i].second), V(outline[(i+1)%outline.size()].second);
    dedgecoord_t Xuv(X0,X1), Xvu(reverse_arc[dedge_t(U,V)]);
    Xvu.first  = Xvu.first *w;
    Xvu.second = Xvu.second*w;

    Eisenstein x0,x0p,T;
    Unfolding::transform_line(Xuv,Xvu, x0,x0p, T);

    vector<Eisenstein> segment(polygon::draw_line(X0,X1)), 
      revsegment(polygon::draw_line(X1,X0));
    reverse(revsegment.begin(),revsegment.end());
    assert(segment.size() == revsegment.size());

    for(int j=0;j<segment.size();j++){
      const Eisenstein& x(segment[j]), y(revsegment[j]);
      if(x != y){
        Eisenstein xp((x-x0)*T+x0p);
        Eisenstein yp((y-x0)*T+x0p);
        node_t u = grid(x*iw), v = grid(yp*iw);
        assert(u>=0 && v>=0);
        edges.push_back(edge_t(u,v));
        
        EC.push_back(make_pair(x,y));
        ECp.push_back(make_pair(xp,yp));
      }
      
    }
  }

  polygon P(w*get_keys(outline));
  if(debug_flags & WRITE_FILE)
    debug_file
      << "points = " << P.allpoints() << ";\n"
      << "outlinecoords = " << P.outline << ";\n" 
      << "EC = " << EC << ";\n"
      << "ECp = " << ECp << ";\n"
      ;

  return edges;
}

vector<edge_t> Folding::connect_polygon(const Eisenstein& w)
{
  Eisenstein iw(w.invertn());
  vector<Eisenstein> outline_coords(w*get_keys(outline));
  polygon P(outline_coords);
  polygon::scanline S(P.scanConvert());
  
  vector<edge_t> edges;
  for(int i=0;i<S.xs.size();i++){
    for(int j=0;j<S.xs[i].size()/2;j++){
      int start = S.xs[i][2*j], end = S.xs[i][2*j+1];
      for(int x=start;x<end;x++){
        node_t u = grid(Eisenstein(x,i+S.minY)*iw), v = grid(Eisenstein(x+1,i+S.minY)*iw);
        if(u>=0 && v>=0) edges.push_back(edge_t(u,v));
      }
    }
  }
  return edges;
}

vector<edge_t> Folding::connect(const Eisenstein& w)
{
  vector<edge_t> 
    polygon_edges(connect_polygon(w)), 
    cross_edges(connect_cross(w)),
    edges;
  
  if(!(debug_flags & DONT_CONNECT_POLYGON)) copy(polygon_edges.begin(),polygon_edges.end(),back_inserter(edges));
  if(!(debug_flags & DONT_CONNECT_ACROSS))  copy(cross_edges.begin(),cross_edges.end(),back_inserter(edges));

  return edges;
}


vector<edge_t> Folding::connect()
{
  vector<edge_t> 
    e(connect(Eisenstein(1,0))),
    eCW(connect(Eisenstein(1,-1))),
    eCCW(connect(Eisenstein(0,1))), edges;
  
  copy(e.begin(),   e.end(),back_inserter(edges));
  if(!(debug_flags & DONT_ROTATE)){
    copy(eCW.begin(), eCW.end(),back_inserter(edges));
    copy(eCCW.begin(),eCCW.end(),back_inserter(edges));
  }

  return edges;
}


vector<int> Folding::identify_nodes()
{
  set<edge_t> same_as;

  // Register reverse arcs
  map<dedge_t,dedgecoord_t> reverse_arc;

  for(int i=0;i<outline.size();i++){
    const Eisenstein xu(outline[i].first), xv(outline[(i+1)%outline.size()].first);
    const node_t u(outline[i].second), v(outline[(i+1)%outline.size()].second);
    reverse_arc[dedge_t(v,u)] = dedgecoord_t(xu,xv);
  }

  for(int i=0;i<outline.size();i++){
    Eisenstein X0(outline[i].first), X1(outline[(i+1)%outline.size()].first);
    node_t     U(outline[i].second), V(outline[(i+1)%outline.size()].second);
    dedgecoord_t Xuv(X0,X1), Xvu(reverse_arc[dedge_t(U,V)]);
    Xvu.first  = Xvu.first ;
    Xvu.second = Xvu.second;

    Eisenstein x0,x0p,T;
    Unfolding::transform_line(Xuv,Xvu, x0,x0p, T);

    //TODO: Handle horizontal lines.
    vector<Eisenstein> segment(polygon::draw_line(X0,X1)), 
      revsegment(polygon::draw_line(X1,X0));
    reverse(revsegment.begin(),revsegment.end());
    assert(segment.size() == revsegment.size());

    for(int j=0;j<segment.size();j++){
      const Eisenstein& x(segment[j]), y(revsegment[j]);
      if(x == y){
        Eisenstein xp((x-x0)*T+x0p);
        node_t u = grid(x), v = grid(xp);
        assert(u>=0 && v>=0);
        same_as.insert(edge_t(u,v));
      }
    }
  }


  // Find connected components
  vector<int> same(grid.size());
  for(int i=0;i<grid.size();i++) same[i] = i;

  Graph S(same_as);
  list< list<node_t> > components(S.connected_components());

  for(list< list<node_t> >::const_iterator s(components.begin()); s!=components.end(); s++){
    node_t canonical = *s->begin();

    for(list<node_t>::const_iterator t(s->begin()); t!=s->end(); t++){
      same[*t] = canonical;
    }
  }

  if(debug_flags&WRITE_FILE) 
    debug_file << "samecomponents = " << components << ";\n";

  return same;
}

PlanarGraph Folding::fold()
{
  vector<int> same_as;
  vector<edge_t> edge_list(connect());

  if(!(debug_flags & DONT_IDENTIFY_NODES)){
    same_as = identify_nodes();
    for(int i=0;i<edge_list.size();i++){
      node_t u = edge_list[i].first, v = edge_list[i].second;
      edge_list[i].first = same_as[u];
      edge_list[i].second = same_as[v];
    }
  }
  
  // Compactify node names
  set<edge_t> edges;
  IDCounter<int> new_names;

  if(!(debug_flags & DONT_IDENTIFY_NODES)) {
    for(int i=0;i<edge_list.size();i++){
      edge_t e = edge_list[i];
      node_t u = e.first, v = e.second;

      int newu = new_names.insert(u);
      int newv = new_names.insert(v);

      edges.insert(edge_t(newu,newv));
    }
  } else {
    edges = set<edge_t>(edge_list.begin(),edge_list.end());
    for(int i=0;i<grid.size();i++) new_names.insert(i);
  }

  int N = new_names.size();
  vector<coord2d> layout2d;

  if(debug_flags & DO_NONPLANAR_LAYOUT){
    layout2d.resize(N);
    for(IDCounter<Eisenstein>::const_iterator xi(grid.begin()); xi!=grid.end(); xi++){
      Eisenstein xy(xi->first);
      node_t      u(new_names(xi->second));

      if(u>=0) layout2d[u] = coord2d(xy.first,xy.second);
    }
  }

  PlanarGraph G(edges,layout2d);

  if(debug_flags & WRITE_FILE) debug_file << "G = " << G << ";\n";

  return G;
}
