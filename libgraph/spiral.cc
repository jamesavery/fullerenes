#include "spiral.cc"

bool operator<(const general_spiral &s) const
{
  return jumps.size() < s.jumps.size() ||
	(jumps.size() == s.jumps.size() && jumps < s.jumps) ||
        (jumps == s.jumps && spiral < s.spiral);
}


ostream &operator<<(ostream &s, const general_spiral &GS)
{
  return s << make_pair(GS.jumps,GS.spiral); 
}


spiral_nomenclature::spiral_nomenclature(const PlanarGraph &g, const string& atom="",
					 bool compatibility=false) : atom(atom), graph(g),cs(compatibility) {
  if(g.is_triangulation()){
    graph_type     = TRIANGULATION;
    is_a_fullerene = g.dual_graph().is_a_fullerene();
    triangulation  = g;
  } else if(g.is_cubic()){
    graph_type     = CUBIC;
    triangulation  = g.dual_graph();
    is_a_fullerene = g.is_a_fullerene();
  } else {
    graph_type     = GENERAL;
    triangulation  = g.leapfrog_dual();
    is_a_fullerene = false;
  }

  permutation.resize(triangulation.N);
  bool spiral_success = triangulation.get_spiral(GS.spiral,GS.jumps,!cs);
  assert(spiral_success);

  // Find spiral order with respect to triangulation vertices
  permutation.resize(triangulation.N);
  spiral_success = triangulation.get_spiral(GS.spiral,GS.jumps,permutation,!cs);
}

template <typename T> vector<T> read_list
spiral_nomenclature::spiral_nomenclature(const string& name)
{
  // Grammar for nomenclature as described in spiral paper:
  // 
  // MOLECULE_NAME   ::= (POINTGROUP '-')? POLYHEDRON_NAME '-' FORMULA '-' SUFFIX
  // POLYHEDRON_NAME ::= '[' METHOD_SPEC? (JUMP+';')? SPIRAL ']'
  //                 |   '[' METHOD_SPEC? (JUMP+';')? PI     ']'
  //                 |   '[' METHOD_SPEC? (JUMP+';')? FI     "]-(" FACE_SIZES ")" BASE_FACE
  //
  // METHOD_SPEC     ::= ("D"|"LF"|"CS"|"GS")+ ':' 
  // SUFFIX          ::= 'cage' | 'fullerene' | 'fulleroid'
  // SPIRAL          ::= NUMBER_LIST
  // PI              ::= NUMBER_LIST
  // FI              ::= NUMBER_LISTS
  // NUMBER_LIST     ::= number (',' number)*
  // NUMBER_LISTS    ::= NUMBER_LIST (';' NUMBER_LIST)*
  // FORMULA         ::= string
  // POINTGROUP      ::= string

  // Simplified grammar which accepts the above + some more:
  //
  // SEPARATOR  ::= " []-"
  // POLYHEDRON_NAME ::= '[' METHOD_SPEC? NUMBER_LISTS ']'
  //                 |   '[' METHOD_SPEC? NUMBER_LISTS "]-(" NUMBER_LIST ")" number
  //
  // FI              ::= NUMBER_LIST (';' NUMBER_LIST)*
  // FORMULA         ::= string
  // POINTGROUP      ::= string
  
  // The grammar is so simple that we'll just brutally parse by hand.

  for(int i=0;i<name.size();i++){
  }

}

friend ostream& operator<<(ostream &s, const spiral_nomenclature &n){
  string graph_type_string[3] = {"","D,","LF,"};
  s << "["<<graph_type_string[n.graph_type] <<(n.compatibility?"CS":"GS") << ": "
    << (n.GS.jumps.empty()? "": (jumps_to_string(n.GS.jumps)+"; "))
    << (n.is_a_fullerene? spiral_to_rspi_string(n.GS.spiral) : spiral_to_string(n.GS.spiral))
    << "]-" <<n.atom<< n.graph.N <<"-" << (n.is_a_fullerene? "fullerene" : "cage");
  return s;
}
