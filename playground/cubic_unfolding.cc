#include "fullerenes/triangulation.hh"
#include <stack>
#include <queue>

struct sort_by_radius : public std::binary_function<dedge_t, dedge_t, bool>
{
  const vector<coord2d> &node_layout;
  
  sort_by_radius(const vector<coord2d> &node_layout) : node_layout(node_layout) {  }
  
  bool operator()(const dedge_t &uv, const dedge_t &st) const
  {
    node_t u,v, s,t;
    tie(u,v) = uv;
    tie(s,t) = st;
    
    coord2d xu = node_layout[u], xv = node_layout[v], xr = node_layout[r], xt = node_layout[t];
    double  ru = xu.norm2(),     rv = xv.norm2(),     rr = xr.norm2(),     rt = xt.norm2();
    
    return min(ru,rv) < min(rr,rt);
  }
};


double circumscribed_circle_radius(const vector<coord2d> &xs, const coord2d &xc)
{
  double radius_sqr = 0;
  for(auto x: xs) radius = max(radius, (x-xc)).norm2();

  return make_pair(xc,radius_sqr);
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

class CubicUnfolding {
  typedef dedge_t arc_t;
  typedef pair<coord2d,coord2d> dedgecoord_t;

  CubicGraph graph;
  map<arc_t,dedgecoord_t> edgecoords;

  vector<arc_t>   face_id_to_rep;	       // faces[i] has canonical arc-representation face_id_to_rep[i]
  map<arc_t, int> face_rep_to_id;	       // face_id_to_rep[face_rep_to_id[faces[i]]] = i

  vector<int>     placed_faces;                // Face id's in order of placement
  vector<vector<coord2d>> face_coords;         // For convenience. In order of placement sequence, not in order of id
  vector< pair<coord2d,node_t> > outline;    

  CubicUnfolding(const Triangulation& graph, int root_face_id=0,int Fmax=6) : graph(graph)
  {
    map<arc_t,bool> face_seen;
    queue<arc_t> work_queue;
    vector<face_t>  faces = graph.compute_faces_oriented(Fmax);

    vector<coord2d> placed_face_centers;
    vector<double>  placed_face_radii_sqr;

    // Step 1: Place root face
    face_t root_face = faces[root_face_id];
    int d         = root_face.size();                   	// root_face is regular d-gon
    double radius = 1/(2*sin(M_PI/d));				// radius of d-gon with unit edge lengths

    node_t u,v;
    coord2d Xu, Xv;

    placed_faces.push_back(root_face_id);
    placed_face_centers.push_back({0,0});
    placed_face_radii_sqr.push_back(radius);

    auto place_face = [&](int face_id, const coord2d& face_center, const coord2d& xu, const coord2d& xv) {
      const face_t& face = faces[face_id];      
      int d              = face.size();                   	// face is regular d-gon
      double radius      = 1/(2*sin(M_PI/d));		// radius of d-gon with unit edge lengths

      
      vector<coord2d> face_coords(d);
      for(int i=0;i<d;i++)
	face_coords[i] = coord2d{radius*cos(i*2*M_PI/d),radius*sin(i*2*M_PI/d)} - face_center;

      for(int i=0;i<d;i++){
	u = face[i], v = face[(i+1)%d];
	edgecoords[{u,v}] = {face_coords[i],face_coords[(i+1)%d]};
	face_seen[{u,v}]  = true;
	work_queue.push({v,u});
      }
    };
    
    // Step 2: Place the rest of the fucking faces
    while(!work_queue.empty()){
      tie(u,v)   = work_queue.front(); work_queue.pop();
      
      tie(Xv,Xu) = edgecoords[{v,u}];

      int face_id = face_rep_to_id.at({u,v});
      face_t f    = faces[face_id];
      
      int d = f.size();
      double rd = 1/(2*sin(2*M_PI/d));
      double ld = cot(2*M_PI/d)/2;
      vector<coord2d> corners(d);
      
      coord2d Xuv = Xv-Xu, X_mid = (Xu+Xv)/2;    // u->v direction vector and midpoint between u and v
      coord2d normal = {Xuv.second, -Xuv.first}; // Clockwise normal
      normal /= normal.norm();
      coord2d Xc = Xuv_mid + normal*ld;          // Center of new d-gon is l_d along the normal from the midpoint
      
      if(/*placeable: face not seen, no corners inside circumcircle of previously placed face*/
      ){	  

	// place face
	
	for(int i=0;i<d;i++){
	  u = root_face[i], v = root_face[(i+1)%d];
	  edgecoords[{u,v}] = {radius*cos(i*2*M_PI/d),radius*sin(i*2*M_PI/d)};
	  face_seen[{u,v}]  = true;
	  work_queue.push({v,u});
	}

	
	
      } else work_queue.push({u,v});
    }

  }
};
