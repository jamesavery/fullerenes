#include "fullerenes/triangulation.hh"
#include "fullerenes/fullerenegraph.hh"
#include <stack>
#include <queue>

struct sort_by_radius : public std::binary_function<arc_t, arc_t, bool>
{
  const vector<coord2d> &node_layout;
  
  sort_by_radius(const vector<coord2d> &node_layout) : node_layout(node_layout) {  }
  
  bool operator()(const arc_t &uv, const arc_t &st) const
  {
    node_t u,v, s,t;
    tie(u,v) = uv;
    tie(s,t) = st;
    
    coord2d xu = node_layout[u], xv = node_layout[v], xs = node_layout[s], xt = node_layout[t];
    double  ru = xu.norm2(),     rv = xv.norm2(),     rs = xs.norm2(),     rt = xt.norm2();
    
    return min(ru,rv) < min(rs,rt);
  }
};


double circumscribed_circle_radius(const vector<coord2d> &xs, const coord2d &xc)
{
  double radius_sqr = 0;
  for(auto x: xs) radius_sqr = max(radius_sqr, (x-xc).norm2());

  return radius_sqr;
}


bool all_points_outside(const vector<coord2d> &new_face,
		   const vector<coord2d> &old_face_centers,
		   const vector<double>  &old_face_radii_sqr)
{
  for(int i=0;i<old_face_centers.size();i++){
    coord2d xc    = old_face_centers[i];
    double  r_sqr = old_face_radii_sqr[i];
    
    for(auto x: new_face)
      if((x-xc).norm2() <= r_sqr) return false;
  }
  return true;
}

bool circles_disjoint(coord2d x1, double r1, coord2d x2, double r2)
{
  return (x2-x1).norm2() > (r1+r2)*(r1+r2);
}

class CubicUnfolding {
  typedef arc_t arc_t;
  typedef pair<coord2d,coord2d> arccoord_t;

  FullereneGraph graph;
  map<arc_t,arccoord_t> edgecoords;

  map<arc_t, int> arc_to_face_id;	       // Each arc belongs to a unique face. This maps each arc to the face_id of this face.

  vector<int>     placed_faces;                // Face id's in order of placement
  vector<coord2d> placed_face_centers;         // Face centers in order of placement
  vector<vector<coord2d>> face_coords;         // For convenience. In order of placement sequence, not in order of id
  vector< pair<coord2d,node_t> > outline;    

  
  CubicUnfolding(const FullereneGraph& graph, int root_face_id=0,int Fmax=6) : graph(graph)
  {
    size_t Nf = graph.N/2+2;	        // Number of faces in a cubic polyhedral graph is always N/2+2
    vector<bool>   face_placed(Nf,false);	// Is the face already placed? Indexed by face_id
    queue<arc_t>   work_queue;		// Next arc to process
    vector<face_t> faces = graph.compute_faces(Fmax);

    double   cosine[Fmax-2], sine[Fmax-2];
    double   circumscribed_radius[Fmax-2];
    double   inscribed_radius[Fmax-2];
    
    // Step 0: Initialize helper data structures
    //  r: distance from center to vertices
    //  l: distance from center to edges
    for(int d=3;d<Fmax;d++){
      double s/*in*/, c/*os*/, t/*an*/, r, l;
      sincos(M_PI/d,&s,&c);
      t = s/c;
      r = 1/(2*s);  /* circumscribed radius */
      l = c*r;      /* inscribed_radius */

      cosine[d-3] = c;
      sine[d-3]   = s;
      circumscribed_radius[d-3] = r;
      inscribed_radius[d-3]     = l;
    }
    
    for(int face_id=0;face_id<Nf;face_id++){
      const face_t &f = faces[face_id];
      size_t        d = f.size();
      
      for(int i=0;i<d;i++){
	node_t u = f[i], v = f[(i+1)%d];
	arc_to_face_id[{u,v}] = face_id;
      }
    }

    
    auto place_face = [&](int face_id, const coord2d& xc, const coord2d& xu, const coord2d& xv) {
      const face_t& face = faces[face_id];      
      int d              = face.size();                	// face is regular d-gon
      double radius      = circumscribed_radius[d-3];	// radius of d-gon with unit edge lengths
      double length      = inscribed_radius[d-3];	// distance from center to edges
      double c = cosine[d-3], s = sine[d-3], t = s/c;

      face_placed[face_id]  = true;

      // Compute transformation matrix A that takes {c,-s} to xu-xc and {c,s} to xv-xc
      matrix2d Ci = matrix2d(t,-1,t,1); /* A^T = CTinverse * ({{xu,yu},{xv,yv}}-xc) => xc+ A*{{c,c},{-s,s}}*r = {{xu,xv},{yu,yv}} */
      matrix2d A  = matrix2d(xu.first, xu.second, xv.first, xv.second)*Ci;
      
      vector<coord2d> face_coords(d);
      complex<double> x{c,-s}, edge_angle = complex<double>{c,s}; // Start in (x,y)=(c,-s). In each step, add 2pi/d to the angle
	
      for(int i=0;i<d;i++, x*= edge_angle){
	coord2d x_transformed = A*coord2d(radius*x.real(),radius*x.imag());
	face_coords[i] = x_transformed + xc;
      }
	
      for(int i=0;i<d;i++){
	node_t u = face[i], v = face[(i+1)%d];
	edgecoords[{u,v}] = {face_coords[i],face_coords[(i+1)%d]};
	work_queue.push({v,u});
      }

      placed_faces.push_back(face_id);
      placed_face_centers.push_back(xc);
    };

    // /*placeable: face not seen, no corners inside circumcircle of previously placed face*/    
    // auto placeable = [&](const vector<coord2d> &face_coords) -> bool {
    // 		       return !face_
    // };
		     
    

    // Step 1: Place root face
    double c, s;
    face_t root_face = faces[root_face_id];
    int d         = root_face.size();                   	// root_face is regular d-gon
    double radius = circumscribed_radius[d-3];			// radius of d-gon with unit edge lengths

    node_t u,v;
    coord2d Xu = {radius,0}, Xv = {radius*c,radius*s};
    place_face(root_face_id, {0,0}, Xu,Xv);   
    
    // Step 2: Place the rest of the fucking faces
    while(!work_queue.empty()){
      tie(u,v)   = work_queue.front(); work_queue.pop();

      tie(Xv,Xu) = edgecoords[{v,u}];

      int face_id = arc_to_face_id.at({u,v});
      face_t f    = faces[face_id];

      if(face_placed[face_id]) continue; // Face is already placed, don't place it again.
      
      int d = f.size();
      double s, c;
      sincos(M_PI/d,&s,&c); 	// Cos and sin of half-angle
      // Radius:                             1/2 = r_d * sin(pi/d) <=> r_d = 1/(2 sin(pi/d))
      // Distance from face center to edge:  l_d = r_d * cos(pi/d)
      double rd = 1/(2*s);
      double ld = c*rd; 
      vector<coord2d> corners(d);
      
      coord2d Xuv = Xv-Xu, X_mid = (Xu+Xv)/2;    // u->v direction vector and midpoint between u and v
      coord2d normal = {Xuv.second, -Xuv.first}; // Clockwise normal
      normal *= 1/normal.norm();
      coord2d Xc = X_mid + normal*ld;          // Center of new d-gon is l_d along the normal from the midpoint
      
      if(1
	 /*placeable: face not seen, no corners inside circumcircle of previously placed face*/
      ){	  

	
      } else work_queue.push({u,v});
    }

  }
};


int main()
{
  return 0;
}
