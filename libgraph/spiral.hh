#pragma once
//TODO: The plan is to gather the spiral stuff here, so it's not scattered around the place.

typedef pair<int,int>  jump_t;
typedef vector<jump_t> jumplist_t;

// This general spiral type can represent any polyhedral graph.
struct general_spiral {
  typedef enum { CUBIC, TRIANGULATION, GENERAL } spiral_type;

  vector<int> spiral;  
  jumplist_t  jumps;
  spiral_type type;

  general_spiral(const vector<int>& spiral,
		 const jumplist_t& jumps = jumplist_t(),
		 const polyhedron_type& type) : spiral(spiral),
						jumps(jumps),
						type(type) {}

  general_spiral(const PlanarGraph &g, bool compatibility=false);
  general_spiral(const string& name);
  
  bool operator<(const general_spiral &s) const;
  friend ostream &operator<<(ostream &s, const general_spiral &GS);
};


struct spiral_nomenclature {
  bool is_a_fullerene;
  string atom;
  string point_group;
  
  PlanarGraph   graph;
  Triangulation triangulation;
  bool compatibility;

  general_spiral GS;
  Permutation permutation;

  spiral_nomenclature(const PlanarGraph &g, const string& atom="", bool compatibility=false);
  spiral_nomenclature(const string& name);
  
  friend ostream& operator<<(ostream &s, const spiral_nomenclature &n);


  
  
};


#endif
