class Patches {
public:
  Triangulation g; 

  typedef pair<node_t,int> swsite_t;

  Patches(const Triangulation& g) : g(g) {}
 
  //  Before:       After: 	  
  //  v(6)---q(5)   v(5)---q(6) 
  //  |  \    |	    |    /  |	  
  //  |    \  |	    |  /    |	  
  //  p(5)---w(6)   p(6)---w(5) 
  Triangulation swtransform(const swsite_t& s) const {
    Triangulation gp(g);

    int    i = s.second;
    node_t p = s.first;
    node_t v = g.neighbours[p][i];
    node_t w = g.next(p,v);
    node_t q = g.prev(v,w);

    gp.remove_edge({w,v});
    gp.insert_edge({p,q},w,v);

    return gp;
  }


  //  Before:              After: 	  
  //  v(6)---y(?)---q(5)   v(5)---y(?)---q(6)  
  //  |  \    |  \    |    |     / |    / |   	  
  //  |    \  |    \  |    |  /    |  /   |   	  
  //  p(5)---x(?)---w(6)   p(6)---x(?)---w(5)  
  Triangulation swtransform1(const swsite_t& s) const {
    Triangulation gp(g);

    int    i = s.second;
    node_t p = s.first;
    node_t v = g.neighbours[p][i];
    node_t x = g.next(p,v);
    node_t y = g.prev(v,x);
    node_t w = g.next(x,y);
    node_t q = g.prev(y,w);

    gp.remove_edge({v,x});
    gp.insert_edge({p,y},x,v);

    gp.remove_edge({y,w});
    gp.insert_edge({x,q},w,y);

    assert(gp.degree(p) == 6 && gp.degree(v) == 5 && 
	   gp.degree(w) == 5 && gp.degree(q) == 6);    
    return gp;
  }

  //  Before:              After: 	  
  //  v(6)---y1(?)---y2(?)---q(5)  v(6)---y1(?)---y2(?)---q(5)   
  //  |  \    |   \    |  \    |   |     / |     /  |     / |      
  //  |    \  |     \  |    \  |   |  /    |  /     |  /    |      
  //  p(5)---x1(?)---x2(?)---w(6)  p(5)---x1(?)---x2(?)---w(6)   
  Triangulation swtransform2(const swsite_t& s) const {
    Triangulation gp(g);

    int    i = s.second;
    node_t p = s.first;
    node_t v = g.neighbours[p][i];
    node_t x1 = g.next(p,v);
    node_t y1 = g.prev(v,x1);
    node_t x2 = g.next(x1,y1);
    node_t y2 = g.prev(y1,x2);
    node_t w  = g.next(x2,y2);
    node_t q  = g.prev(y2,w);

    gp.remove_edge({v,x1});
    gp.insert_edge({p,y1},x1,v);

    gp.remove_edge({y1,x2});
    gp.insert_edge({x1,y2},x2,y1);

    gp.remove_edge({y2,w});
    gp.insert_edge({x2,q},w,y2);

    assert(gp.degree(p) == 6 && gp.degree(v) == 5 && 
	   gp.degree(w) == 5 && gp.degree(q) == 6);    
    return gp;
  }

  //  Before:               After:
  //  v---y1---y2---y3---q  v---y1---y2--y3---q 
  //  | \ |  \ |  \  | \ |  | /  | /  | / | / | 
  //  p---x1---x2---x3---w  p---x1---x2--x3---w 
  Triangulation swtransform3(const swsite_t& s) const {
    Triangulation gp(g);

    int    i = s.second;
    node_t p = s.first;
    node_t v = g.neighbours[p][i];
    node_t x1 = g.next(p,v); 	    // deg ?
    node_t y1 = g.prev(v,x1); 	    // deg ?
    node_t x2 = g.next(x1,y1); 	    // deg ?
    node_t y2 = g.prev(y1,x2); 	    // deg ?
    node_t x3 = g.next(x2,y2); 	    // deg ?
    node_t y3 = g.prev(y2,x3); 	    // deg ?
    node_t w  = g.next(x3,y3);	    // deg 6
    node_t q  = g.prev(y3,w);	    // deg 5

    gp.remove_edge({v,x1});
    gp.insert_edge({p,y1},x1,v);

    gp.remove_edge({y1,x2});
    gp.insert_edge({x1,y2},x2,y1);

    gp.remove_edge({y2,x3});
    gp.insert_edge({x2,y3},x3,y2);

    gp.remove_edge({y3,w});
    gp.insert_edge({x3,q},w,y3);

    assert(gp.degree(p) == 6 && gp.degree(v) == 5 && 
	   gp.degree(w) == 5 && gp.degree(q) == 6);    
    return gp;
  }

  //  Before:               After:
  //  v---y1-- ... --yd---q  v---y1-- ... --yd---q 
  //  | \ |  \    \  | \ |  | /  | /         | / | 
  //  p---x1-- ... --xd---w  p---x1-- ... --xd---w 
  Triangulation swtransform(const swsite_t& s, int d) const {
    Triangulation gp(g);

    node_t p = s.first;
    node_t v = s.second;
    node_t x=p, y = v, xp, yp;

    for(int i=0;i<=d;i++, x=xp, y=yp){
      xp = g.next(x,y);
      yp = g.next(xp,y);

      gp.remove_edge({y,xp});
      gp.insert_edge({x,yp},xp,y);
    }

    node_t w  = x;	    // deg 6
    node_t q  = y;	    // deg 5

    assert(gp.degree(p) == 6 && gp.degree(v) == 5 && 
	   gp.degree(w) == 5 && gp.degree(q) == 6);    
    return gp;
  }

  
  
  bool swmatch(node_t p, node_t v) const {
    if(g.degree(p) != 5 || g.degree(v) != 6) return false;

    node_t w = g.next(p,v); 
    node_t q = g.prev(v,w);

    return g.degree(w) == 6 && g.degree(q) == 5;
  }

  //  Match:
  //  v---y---q
  //  | \ | \ |
  //  p---x---w 
  bool swmatch1(node_t p, node_t v) const {
    if(g.degree(p) != 5 || g.degree(v) != 6) return false;

    node_t x = g.next(p,v); 	   // deg ?
    node_t y = g.prev(v,x); 	   // deg ?
    node_t w = g.next(x,y);	   // deg 6
    node_t q = g.prev(y,w);	   // deg 5

    return g.degree(w) == 6 && g.degree(q) == 5;
  }

  //  Match:
  //  v---y1---y2---q
  //  | \ |  \ |  \ |
  //  p---x1---x2---w 
  bool swmatch2(node_t p, node_t v) const {
    if(g.degree(p) != 5 || g.degree(v) != 6) return false;

    node_t x1 = g.next(p,v); 	    // deg ?
    node_t y1 = g.prev(v,x1); 	    // deg ?
    node_t x2 = g.next(x1,y1); 	    // deg ?
    node_t y2 = g.prev(y1,x2); 	    // deg ?
    node_t w  = g.next(x2,y2);	    // deg 6
    node_t q  = g.prev(y2,w);	    // deg 5
    
    return g.degree(w) == 6 && g.degree(q) == 5;
  }


  //  Match:
  //  v---y1---y2---y3---q
  //  | \ |  \ |  \  | \ |
  //  p---x1---x2---x3---w 
  bool swmatch3(node_t p, node_t v) const {
    if(g.degree(p) != 5 || g.degree(v) != 6) return false;

    node_t x1 = g.next(p,v); 	    // deg ?
    node_t y1 = g.prev(v,x1); 	    // deg ?
    node_t x2 = g.next(x1,y1); 	    // deg ?
    node_t y2 = g.prev(y1,x2); 	    // deg ?
    node_t x3 = g.next(x2,y2); 	    // deg ?
    node_t y3 = g.prev(y2,x3); 	    // deg ?
    node_t w  = g.next(x3,y3);	    // deg 6
    node_t q  = g.prev(y3,w);	    // deg 5
    
    return g.degree(w) == 6 && g.degree(q) == 5;
  }

  //  Match:
  //  v---y1-- ... --yd---q
  //  | \ |          |  \ |
  //  p---x1-- ... --xd---w 
  bool swmatch(node_t p, node_t v, int d) const { // TODO: Test for inverse Stone-Wales as well!
    if(g.degree(p) != 5 || g.degree(v) != 6) return false;

    // TODO: For d>=5 we risk going around in circles! Check that we never visit a vertex twice.
    vector<bool> seen(g.N,false);
    seen[p] = true;
    seen[v] = true;

    node_t x=p, y = v, xp, yp;
    for(int i=0;i<d;i++, x=xp, y=yp){
      xp = g.next(x,y);  if(seen[xp]) return false; else seen[xp] = true;
      yp = g.next(xp,y); if(seen[yp]) return false; else seen[yp] = true;
    }
    node_t w  = g.next(x,y);	   // deg 6
    node_t q  = g.prev(y,w);	   // deg 5
    
    return g.degree(w) == 6 && g.degree(q) == 5 && !seen[w] && !seen[q];
  }

  vector< swsite_t > swsites() const {
    vector<swsite_t> sites;
    for(node_t p=0;p<g.N;p++){
      if(g.degree(p) == 5)
	for(int i=0;i<5;i++){
	  node_t v = g.neighbours[p][i];
	  if(swmatch(p,v)) sites.push_back({p,v});
	}
    }
    return sites;
  }

  vector< swsite_t > swsites1() const {
    vector<swsite_t> sites;
    for(node_t p=0;p<g.N;p++){
      if(g.degree(p) == 5)
	for(int i=0;i<5;i++){
	  node_t v = g.neighbours[p][i];
	  if(swmatch1(p,v)) sites.push_back({p,v});
	}
    }
    return sites;
  }

  vector< swsite_t > swsites2() const {
    vector<swsite_t> sites;
    for(node_t p=0;p<g.N;p++){
      if(g.degree(p) == 5)
	for(int i=0;i<5;i++){
	  node_t v = g.neighbours[p][i];
	  if(swmatch2(p,v)) sites.push_back({p,v});
	}
    }
    return sites;
  }

  vector< swsite_t > swsites3() const {
    vector<swsite_t> sites;
    for(node_t p=0;p<g.N;p++){
      if(g.degree(p) == 5)
	for(int i=0;i<5;i++) {
	  node_t v = g.neighbours[p][i];
	  if(swmatch3(p,v)) sites.push_back({p,v});
	}
    }
    return sites;
  }

  vector< swsite_t > swsites(int d) const {
    vector<swsite_t> sites;
    for(node_t p=0;p<g.N;p++){
      if(g.degree(p) == 5)
	for(int i=0;i<5;i++){
	  node_t v = g.neighbours[p][i];
	  if(swmatch(p,v,d)) sites.push_back({p,v});
	}
    }
    return sites;
  }

  PlanarGraph swlayout(node_t p, node_t v, int d) const {
    assert(g.degree(p) == 5 && g.degree(v) == 6 && g.next(p,v) != -1);
    matrix2d hex_ccw = matrix2d::rotation(2*M_PI/6), penta_ccw = matrix2d::rotation(2*M_PI/5),
             hex_cw  = matrix2d::rotation(-2*M_PI/6), penta_cw = matrix2d::rotation(-2*M_PI/5);
    matrix2d rotate_ccw[2] = {hex_ccw,penta_ccw};
    matrix2d rotate_cw[2]  = {hex_cw,penta_cw};

    vector<coord2d> coords(g.N);
    vector<bool>    placed(g.N);
    coords[p] = {0,0}; placed[p] = true;
    coords[v] = {1,0}; placed[v] = true;
    coord2d direction;

    vector<node_t> workstack = {p,v};
    node_t x=p, y = v, xp, yp;
    for(int i=0;i<d;i++, x=xp, y=yp){
      xp = g.next(x,y);  direction = rotate_ccw[g.degree(x)==5]*(coords[y]-coords[x]); coords[xp] = coords[x]+direction; placed[xp] = true;
      yp = g.prev(y,xp); direction = rotate_cw[g.degree(x)==5]*(coords[xp]-coords[y]); coords[yp] = coords[y]+direction; placed[yp] = true;

      workstack.push_back(xp);
      workstack.push_back(yp);
    }

    while(!workstack.empty()){ // Place everything around known nodes
      node_t  u  = workstack.back(); workstack.pop_back();
      coord2d ux = coords[u];
      
      const vector<node_t> &nu(g.neighbours[u]);
      matrix2d rot = rotate_ccw[nu.size()==5];

      // Place unplaced vertices with placed predecessors
      // TODO: Make sure to place all vertices?
      for(int i=0;i<nu.size();i++){
	node_t v = nu[i], w = g.prev(u,v);
	if(!placed[v] && placed[w]){
	  direction = rot*(coords[w]-coords[u]); coords[v] = coords[u]+direction; placed[v] = true;
	}
      }
    }

    set<edge_t> placed_edges;
    for(node_t u=0;u<g.N;u++)
      for(node_t v: g.neighbours[u])
	if(placed[u] && placed[v]) placed_edges.insert(edge_t{u,v});

    return PlanarGraph(Graph(placed_edges), coords);
  }

};


