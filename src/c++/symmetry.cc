#include "fullerenes/symmetry.hh"
#include "fullerenes/auxiliary.hh"
using namespace std;

//////////////////////////////////////////////////////////////////////
//	      POINT GROUPS FOR FULLEROIDS OF DEGREE <= 6
//////////////////////////////////////////////////////////////////////
// Reference: Deza-2009, Theorem 2.1 (iii)
// For future reference: Theorem 2.2 lists symmetry groups for all triangulations
// of degree <= 6, arranged by signature (p3,p4,p5). 
PointGroup PointGroup::FullereneSymmetries[28] = {
  {C,1},      {C,2},       {C,REF_I},   {C,REF_S},
  {C,3},      {D,2},       {S,4},       {C,2,REF_V},
  {C,2,REF_H},{D,3},       {S,6},       {C,3,REF_V},
  {C,3,REF_H},{D,2,REF_H}, {D,2,REF_D}, {D,5},  
  {D,6},      {D,3,REF_H}, {D,3,REF_D}, {T},  
  {D,5,REF_H},{D,5,REF_D}, {D,6,REF_H}, {D,6,REF_D}, 
  {T,REF_D},  {T,REF_H},   {I},         {I,REF_H}
};

PointGroup::PointGroup(const string& name_) : sym_type(UNKNOWN), n(0), sym_reflection(NONE)
{
  string name(name_);
  const char ts[7] = {'?','C','D','T','S','O','I'};
  const char rs[6] = {' ','v','h','d','i','s'};

  // trim name - TODO: fix
  while(name[0] == ' ') name = name.substr(1);

  for(int st=6;st>0;st--) 
    if(name[0] == ts[st]) sym_type = symmetry_type(st);
  
  if(name.size()  == 1) return;

  size_t idx = 0;
  if(name[1] >= '0' && name[1] <='9')
    n = stol(name.substr(1),&idx,0);
  
  if(name.size() <= idx+1) return;
  
  for(int sr=5;sr>0;sr--)
    if(name[idx+1] == rs[sr]) sym_reflection = symmetry_reflection(sr);
}


string PointGroup::to_string() const {
  const char ts[7] = {'?','C','D','T','S','O','I'};
  const char rs[6] = {' ','v','h','d','i','s'};

  string result;
  result += ts[sym_type];
  if(n>0) result += std::to_string(n);
  if(sym_reflection != NONE) result += rs[sym_reflection];
  return result;
}



//////////////////////////////////////////////////////////////////////
//			PERMUTATION DATA TYPE
//////////////////////////////////////////////////////////////////////
Permutation Permutation::inverse() const {
  const vector<int> &p(*this);
  vector<int> ip(size());

  for(int i=0;i<size();i++) ip[p[i]] = i;
  return Permutation(ip);
}


Permutation Permutation::identity(int N)
{
  Permutation Id(N);
  for(int i=0;i<N;i++) Id[i] = i;
  return Id;
}

int Permutation::order() const {
  Permutation Id(identity(size())), power(*this), p(*this);

  int order=1;
  while(power != Id){ power = power*p; order++; }
  return order;
}

Permutation Permutation::operator*(const Permutation& q) const {
  assert(size() == q.size());

  const vector<int> &p(*this);
  vector<int> r(size());
  for(int i=0;i<size();i++) r[i] = q[p[i]];

  return Permutation(r);
}

// NB: Removed. This seems identical to vector '=='-operator?
// bool Permutation::operator==(const Permutation& q) const {
//   const vector<int> &p(*this);

//   if(size() != q.size()) return false;
//   for(int i=0;i<size();i++) if(p[i] != q[i]) return false;
//   return true;
// }


//////////////////////////////////////////////////////////////////////
//		  SYMMETRY-DETECTION IMPLEMENTATION
//////////////////////////////////////////////////////////////////////
vector<Permutation> Symmetry::tri_permutation(const vector<Permutation>& Gf) const {
  auto tris = triangles();
  assert(tris.size() == (N-2)*2); // Triangulation is cubic dual
  vector<Permutation> Gtri(Gf.size(),Permutation(tris.size()));
  IDCounter<tri_t> tri_id;

  for(int i=0;i<tris.size();i++) tri_id.insert(tris[i].sorted());

  for(int j=0;j<Gf.size();j++){
    const Permutation& pi = Gf[j];
    for(int i=0;i<tris.size();i++){
      const tri_t &t  = tris[i];
      const tri_t &tp = {pi[t[0]], pi[t[1]], pi[t[2]]}; 
      
      int tp_id = tri_id(tp.sorted());

      Gtri[j][i] = tp_id;

      if(tp_id < 0){
	cout << "SYMMETRY OPERATION DOES NOT MAP TRIANGLE TO EXISTING TRIANGLE.\n";
	cout << "tp_id = " << tp_id << endl;
	cout << "pi["<<j<<"] = " << pi << endl;
	cout << "t = " << t <<endl;
	cout << "tp = " << tp << endl;
	cout << "tp.sorted() = " << tp.sorted() << endl;

	cout << "triangles = " << tris << endl;
	
	assert(tp_id >= 0);
      }
    }
  }
  return Gtri;
}

vector<Permutation> Symmetry::edge_permutation(const vector<Permutation>& Gf) const
{
  vector<Permutation> Gedge(Gf.size(),Permutation(edge_id.size()));
    
  for(int j=0;j<Gf.size();j++){
    for(const auto &ei: edge_id){
      edge_t e = ei.first;
      int    i = ei.second;
      
      Gedge[j][i] = edge_id({Gf[j][e.first],Gf[j][e.second]});
    }
  }
  return Gedge;
}

vector<Permutation> Symmetry::arc_permutation(const vector<Permutation>& Gf) const {

  vector<Permutation> Gedge(Gf.size(),Permutation(arc_id.size()));
    
  for(int j=0;j<Gf.size();j++){
    for(const auto &ei: arc_id){
      arc_t e = ei.first;
      int     i = ei.second;
      
      Gedge[j][i] = arc_id({Gf[j][e.first],Gf[j][e.second]});
    }
  }
  return Gedge;
}

vector<Permutation> Symmetry::permutation_representation() const
{
  // get_spiral_implementation returns position→node maps: perm[i] = node at spiral position i.
  // When the graph was constructed from a spiral, node indices match spiral positions (P0 = identity),
  // so these maps are directly valid node→node automorphisms. But when the graph comes from an
  // external source (e.g. buckygen), node indices differ from spiral positions, and we must convert:
  //   automorphism[P0[i]] = P[i]  ⟺  automorphism = P0⁻¹ ∘ P
  // where P0 is the position→node map for the canonical spiral and P is for each symmetric start.
  //
  // Collect all position→node maps that produce spiral S0 with jumps J0.
  // get_spiral_implementation checks S0 degree-by-degree but does NOT check J0,
  // so we must verify jumps == J0 ourselves. A map with different jumps has
  // different position-space connectivity (boundary rotations differ), making
  // P0⁻¹∘P invalid.
  vector<Permutation> numberings;

  for(node_t u=0;u<N;u++){
    if(degree(u) == S0[0])
      for(const node_t &v: nbrs(u)){
	if(degree(v) == S0[1]){
	  vector<int> spiral,permutation;
	  jumplist_t  jumps;

	  node_t wCCW = next(u,v), wCW = prev(u,v);

	  if(degree(wCCW) == S0[2] && get_spiral_implementation(u,v,wCCW,spiral,jumps,permutation,true,S0,J0) && jumps == J0)
	    numberings.push_back(permutation);
	  if(degree(wCW)  == S0[2] && get_spiral_implementation(u,v,wCW,spiral,jumps,permutation,true,S0,J0) && jumps == J0)
	    numberings.push_back(permutation);
	}
      }
  }

  if(numberings.empty()) return {};

  // Use the first found map as P0 (the canonical reference). P0⁻¹∘P0 = identity.
  Permutation P0_inv = numberings[0].inverse();

  vector<Permutation> pi;
  for(auto& P : numberings)
    pi.push_back(P0_inv * P);

  return pi;
}

vector<int> Symmetry::site_symmetry_counts(const vector<Permutation>& pi) const
{
  vector<int> m(12); 		// Probably needs to be changed for more general point groups
  int order = pi.size(), M = pi[0].size();

  vector<bool> seen(M,false);
  for(int i=0;i<M;i++){	// Calculate length and site-symmetry group order of every orbit. Iterate through all sites, but skip seen one to only do each orbit once.
    if(seen[i]) continue;
    seen[i] = true;
      
    int orbit_length = 1;
    for(int j=0;j<order;j++){
      int I = pi[j][i];
      assert(I<M); // TODO: No asserts/aborts allowed. Eliminate failure modes instead, throw exception on should-never-happen-failures, incorporate remaining-failures-that-can-happen into return type or status param.
      assert(I>=0);
      if(seen[I]) continue;
      seen[I] = true;
      orbit_length++;
    }	
    int site_order = order/orbit_length; 
    assert(site_order <= 12); // TODO: Only holds for fullerenes. Remove or make general.
    m[site_order-1]++;
  }
  return m;
}

vector<int> Symmetry::involutions() const // Returns the involutions *except* from the identity
{
  Permutation Id(Permutation::identity(N));
  vector<int> result;
  for(int i=0;i<G.size();i++) 
    if(G[i]*G[i] == Id && G[i] != Id) result.push_back(i);
  return result;
}

vector< vector<int> > Symmetry::multiplication_table() const 
{
  IDCounter<Permutation> pid;

  for(int i=0;i<G.size();i++) pid.insert(G[i]);

  vector< vector<int> > table(G.size(), vector<int>(G.size()));

  for(int i=0;i<G.size();i++)
    for(int j=0;j<G.size();j++)
      table[i][j] = pid(G[i]*G[j]);
    
  return table;
}

vector<int> Symmetry::fixpoints(const Permutation& pi) const {
  vector<int> fp;
  for(int i=0;i<N;i++) if(pi[i] == i) fp.push_back(i);
  return fp;
}


vector<int> Symmetry::group_fixpoints(const vector<Permutation>& G) const { 
  vector<int> fp;
  for(int i=0;i<N;i++){
    bool fixed = true;
    for(int j=0;j<G.size();j++) if(G[j][i] != i) fixed = false;
    if(fixed) fp.push_back(i);
  }
  return fp;
}


bool Symmetry::reverses_orientation(const Permutation& pi) const
{
  // Compare how pi maps the CCW neighbor ordering of vertex 0
  // against the CCW neighbor ordering of vertex pi[0].
  //
  // If pi preserves orientation, the images pi[nb0[0]], pi[nb0[1]], ...
  // appear in the same cyclic order around pi[0] in the original graph.
  // If pi reverses orientation, they appear in the opposite cyclic order.

  auto nb0 = (*this)[0];         // neighbors of vertex 0 in CCW order
  auto nb_pi0 = (*this)[pi[0]];  // neighbors of pi[0] in CCW order
  int deg = nb_pi0.size();

  // Find where pi[nb0[0]] sits in pi[0]'s neighbor list.
  int pos0 = -1;
  for (int i = 0; i < deg; i++)
    if (nb_pi0[i] == pi[nb0[0]]) { pos0 = i; break; }
  assert(pos0 >= 0);

  // Check if pi[nb0[1]] follows pi[nb0[0]] in CCW or CW direction.
  node_t pi_nb1 = pi[nb0[1]];
  if (nb_pi0[(pos0 + 1) % deg] == pi_nb1) return false;  // same cyclic order
  if (nb_pi0[(pos0 + deg - 1) % deg] == pi_nb1) return true;  // reversed

  fprintf(stderr, "reverses_orientation: pi is not a valid automorphism\n");
  abort();
  return false;
}


PointGroup Symmetry::point_group() const
{
  vector<int> 
    mF = site_symmetry_counts(G),
    mV = site_symmetry_counts(Gtri),
    mE = site_symmetry_counts(Gedge);

  vector<int> mS(13,0);
  for(int i=0;i<12;i++) mS[i+1] = mF[i] + mV[i] + mE[i];
  
  switch(G.size()){
  case 1:
    return PointGroup("C1");

  case 2:
    switch(mS[2]){
    case 0: 
      return PointGroup("Ci");
    case 2:
      return PointGroup("C2");
    default:
      if(mS[2]>2) 
	return PointGroup("Cs");
    }
    break;
  case 3:
    return PointGroup("C3");

  case 4: 
    switch(mS[4]){
    case 0:
      switch(mS[2]){
      case 1:
	return PointGroup("S4");
      case 3: 
	return PointGroup("D2");
      default:
	if(mS[2]>3) return PointGroup("C2h");
      }
      break;
    case 2:
      return PointGroup("C2v");
    default:
      break;
    }
    break;
  case 5:  // No fullerene groups of order 5 -- fill out for fulleroids
    break;

  case 6:
    switch(mS[6]){
    case 0:
      switch(mS[2]){
      case 0:
	return PointGroup("S6");
      case 2:
	return PointGroup("D3");
      default:
	if(mS[2]>2) return PointGroup("C3h");	
      }
      break;
    case 2:
      return PointGroup("C3v");
    default:
      break;
    }
    break;

  case 7:  // No fullerene groups of order 7 -- fill out for fulleroids
    break;

  case 8:
    switch(mS[4]){
    case 1: 
      return PointGroup("D2d");
    case 3:
      return PointGroup("D2h");
    default:
      break;
    }
    break;

  case 10:
    return PointGroup("D5");

  case 12:
    switch(mS[6]){
    case 0:
      return PointGroup("T");
    case 1:
      switch(mS[4]){
      case 0:
	switch(mS[2]){
	case 2:
	  return PointGroup("D6");
	default:
	  if(mS[2]>2) return PointGroup("D3d");
	}
	break;
      case 2:
	return PointGroup("D3h");
      default:
	break;
      }     
    default:
      break;
    }
    break;

  case 20:
    switch(mS[4]){
    case 0:
      return PointGroup("D5d");
    case 2:
      return PointGroup("D5h");
    default:
      break;
    }
    break;

  case 24:
    switch(mS[12]){
    case 0:
      switch(mS[6]){
      case 0:
	return PointGroup("Th");
      case 2:
	return PointGroup("Td");
      default:
	break;
      }
      break;
    case 1:
      switch(mS[4]){
      case 0:
	return PointGroup("D6d");
      case 2:
	return PointGroup("D6h");
      default:
	break;
      }
    default:
      break;
    }
    break;
    
  case 60:
    return PointGroup("I");

  case 120:
    return PointGroup("Ih");

  default:
    break;
  }

  return PointGroup();
}


vector< pair<int,int> > Symmetry::NMR_pattern() const
{
  vector<int>  mV = site_symmetry_counts(Gtri);
  vector< pair<int,int> > NMR;
  
  // F&M
  int order = G.size();
  for(int K=6;K>=1;K--)
    if(mV[K-1] != 0) NMR.push_back(make_pair(mV[K-1],order/K));
  
  return NMR;
}


vector<vector<node_t>> Symmetry::equivalence_classes(const vector<Permutation>& G) const {
  size_t N = G[0].size();
  Graph E(N);

  for(auto &pi: G)
    for(node_t u=0;u<N;u++) E.insert_edge({u,pi[u]});

  return E.connected_components();
}

// ============================================================================
//    3D representation of the point group (rotation/reflection matrices)
// ============================================================================

namespace {

// --- Elementary 3D operations ---

// Rotation about z by angle theta.
matrix3d Rz(double theta) {
  double c = cos(theta), s = sin(theta);
  return matrix3d(c, -s, 0,
                  s,  c, 0,
                  0,  0, 1);
}

// Rotation about x by angle theta.
matrix3d Rx(double theta) {
  double c = cos(theta), s = sin(theta);
  return matrix3d(1, 0,  0,
                  0, c, -s,
                  0, s,  c);
}

// Rotation about an arbitrary unit axis by angle theta (Rodrigues).
matrix3d Raxis(const coord3d& axis, double theta) {
  double c = cos(theta), s = sin(theta), t = 1 - c;
  double x = axis[0], y = axis[1], z = axis[2];
  return matrix3d(t*x*x + c,   t*x*y - s*z, t*x*z + s*y,
                  t*x*y + s*z, t*y*y + c,   t*y*z - s*x,
                  t*x*z - s*y, t*y*z + s*x, t*z*z + c);
}

const matrix3d I3 = matrix3d::unit_matrix();
const matrix3d neg_I3 = matrix3d(-1,0,0, 0,-1,0, 0,0,-1);  // inversion
const matrix3d sigma_h = matrix3d(1,0,0, 0,1,0, 0,0,-1);   // xy-plane reflection
const matrix3d sigma_v = matrix3d(1,0,0, 0,-1,0, 0,0,1);   // xz-plane reflection
const matrix3d C2x = matrix3d(1,0,0, 0,-1,0, 0,0,-1);      // 180 deg about x
const matrix3d C2z = Rz(M_PI);

// Reflection in the plane containing the z-axis and the vector (cos(phi), sin(phi), 0).
matrix3d sigma_v_rotated(double phi) {
  // Rotate sigma_v by phi about z: Rz(phi) * sigma_v * Rz(-phi)
  double c2 = cos(2*phi), s2 = sin(2*phi);
  return matrix3d(c2,  s2, 0,
                  s2, -c2, 0,
                   0,   0, 1);
}

// Generate group closure from a set of generators.
// Returns all unique matrices (using Frobenius norm of difference for equality).
vector<matrix3d> generate_group(const vector<matrix3d>& gens, int max_order = 240) {
  vector<matrix3d> group = {I3};

  auto contains = [](const vector<matrix3d>& v, const matrix3d& m) {
    for (auto& g : v)
      if ((g - m).norm() < 1e-10) return true;
    return false;
  };

  bool changed = true;
  while (changed && (int)group.size() < max_order) {
    changed = false;
    size_t n = group.size();
    for (size_t i = 0; i < n; i++) {
      for (auto& gen : gens) {
        matrix3d product = group[i] * gen;
        if (!contains(group, product)) {
          group.push_back(product);
          changed = true;
          if ((int)group.size() >= max_order) break;
        }
        // Also try gen * group[i] (non-abelian groups)
        product = gen * group[i];
        if (!contains(group, product)) {
          group.push_back(product);
          changed = true;
          if ((int)group.size() >= max_order) break;
        }
      }
      if ((int)group.size() >= max_order) break;
    }
  }
  return group;
}

// Matrix order: smallest k >= 1 such that M^k ≈ I.
int matrix_order(const matrix3d& M) {
  matrix3d P = M;
  for (int k = 1; k <= 120; k++) {
    if ((P - I3).norm() < 1e-8) return k;
    P = P * M;
  }
  return -1;
}

// Determinant of 3x3 matrix.
double det3(const matrix3d& M) {
  return M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
       - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
       + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
}

// Build the multiplication table for a set of matrices.
// table[i][j] = index k such that matrices[i]*matrices[j] ≈ matrices[k].
// Returns empty if any product is missing (shouldn't happen for a closed group).
vector<vector<int>> matrix_mult_table(const vector<matrix3d>& matrices) {
  int n = matrices.size();
  vector<vector<int>> table(n, vector<int>(n, -1));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      matrix3d P = matrices[i] * matrices[j];
      for (int k = 0; k < n; k++)
        if ((matrices[k] - P).norm() < 1e-8) {
          table[i][j] = k;
          break;
        }
    }
  return table;
}

// Find an isomorphism between two groups given by their multiplication tables.
// Returns perm[i] = index in target group corresponding to source group element i.
// Empty return means no isomorphism found.
// Uses backtracking with constraint propagation.
vector<int> find_isomorphism(const vector<vector<int>>& src_table,
                             const vector<vector<int>>& tgt_table,
                             const vector<int>& src_orders,
                             const vector<int>& tgt_orders,
                             const vector<int>& src_chars,  // +1 proper, -1 improper
                             const vector<int>& tgt_chars)
{
  int n = src_table.size();
  if ((int)tgt_table.size() != n) return {};

  // The identity (index 0 in both) must map to identity.
  vector<int> mapping(n, -1);
  vector<bool> used(n, false);
  mapping[0] = 0;
  used[0] = true;

  // Recursive backtracking.
  function<bool(int)> solve = [&](int pos) -> bool {
    if (pos == n) return true;
    if (mapping[pos] >= 0) return solve(pos + 1);

    for (int t = 0; t < n; t++) {
      if (used[t]) continue;
      if (src_orders[pos] != tgt_orders[t]) continue;
      if (src_chars[pos] != tgt_chars[t]) continue;

      // Check consistency with all already-assigned elements.
      bool consistent = true;
      for (int i = 0; i < n && consistent; i++) {
        if (mapping[i] < 0) continue;
        // src[pos]*src[i] should map to tgt[t]*tgt[mapping[i]]
        int src_prod = src_table[pos][i];
        int tgt_prod = tgt_table[t][mapping[i]];
        if (mapping[src_prod] >= 0 && mapping[src_prod] != tgt_prod) consistent = false;

        // src[i]*src[pos] should map to tgt[mapping[i]]*tgt[t]
        src_prod = src_table[i][pos];
        tgt_prod = tgt_table[mapping[i]][t];
        if (mapping[src_prod] >= 0 && mapping[src_prod] != tgt_prod) consistent = false;
      }
      if (!consistent) continue;

      mapping[pos] = t;
      used[t] = true;

      // Propagate: fill in any forced assignments from multiplication table.
      vector<pair<int,int>> forced;
      bool ok = true;
      for (int i = 0; i < n && ok; i++) {
        if (mapping[i] < 0) continue;
        // pos*i -> product must be assigned
        int sp = src_table[pos][i], tp = tgt_table[t][mapping[i]];
        if (mapping[sp] < 0 && !used[tp]) {
          if (src_orders[sp] != tgt_orders[tp] || src_chars[sp] != tgt_chars[tp]) { ok = false; break; }
          mapping[sp] = tp; used[tp] = true; forced.push_back({sp, tp});
        } else if (mapping[sp] >= 0 && mapping[sp] != tp) { ok = false; }

        sp = src_table[i][pos]; tp = tgt_table[mapping[i]][t];
        if (mapping[sp] < 0 && !used[tp]) {
          if (src_orders[sp] != tgt_orders[tp] || src_chars[sp] != tgt_chars[tp]) { ok = false; break; }
          mapping[sp] = tp; used[tp] = true; forced.push_back({sp, tp});
        } else if (mapping[sp] >= 0 && mapping[sp] != tp) { ok = false; }
      }

      if (ok && solve(pos + 1)) return true;

      // Undo
      mapping[pos] = -1;
      used[t] = false;
      for (auto& [s, tt] : forced) { mapping[s] = -1; used[tt] = false; }
    }
    return false;
  };

  if (solve(1)) return mapping;
  return {};
}

// Golden ratio for icosahedral geometry.
static const double PHI = (1.0 + sqrt(5.0)) / 2.0;

// Standard icosahedral C5 rotation axis (along a vertex of the icosahedron).
// Standard orientation: C5 along z (0,0,1).
// C3 along (1,1,1)/sqrt(3) -- but for icosahedron we need a specific convention.
//
// Icosahedron vertices (unit sphere): take the axis through (0, 0, 1) as C5.
// The 5 neighbors of the north pole are at latitude atan(1/2) ≈ 26.57 deg.
// C3 axis through center of a face adjacent to north pole.

// Generate icosahedral rotation group I (60 elements).
vector<matrix3d> generate_icosahedral_group() {
  // Use two generators:
  // C5 about z-axis: rotation by 2*pi/5
  // C3 about an axis through a face center adjacent to the "north pole"
  //
  // The standard icosahedron has vertices at:
  //   (0, 0, +-1) and (cos(2*pi*k/5)*sin(a), sin(2*pi*k/5)*sin(a), +-cos(a))
  //   where a = atan(2)
  // A face connects the north pole (0,0,1) to two adjacent equatorial vertices.
  // The C3 axis goes through the centroid of this face.

  double a = atan(2.0);  // colatitude of upper ring
  coord3d north(0, 0, 1);
  coord3d v1(sin(a), 0, cos(a));
  coord3d v2(sin(a)*cos(2*M_PI/5), sin(a)*sin(2*M_PI/5), cos(a));

  coord3d face_center = north + v1 + v2;
  double fn = face_center.norm();
  coord3d c3_axis = face_center / fn;

  matrix3d gen_c5 = Rz(2*M_PI/5);
  matrix3d gen_c3 = Raxis(c3_axis, 2*M_PI/3);

  return generate_group({gen_c5, gen_c3}, 60);
}

// Generate tetrahedral rotation group T (12 elements).
vector<matrix3d> generate_tetrahedral_group() {
  // Generators: C3 about (1,1,1)/sqrt(3), C2 about z
  coord3d c3_axis(1, 1, 1);
  c3_axis /= c3_axis.norm();
  return generate_group({Raxis(c3_axis, 2*M_PI/3), C2z}, 12);
}

// Generate the standard 3D matrices for a given point group.
vector<matrix3d> standard_matrices(const PointGroup& pg) {
  string name = pg.to_string();

  // --- Order 1 ---
  if (name == "C1") return {I3};

  // --- Order 2 ---
  if (name == "C2") return generate_group({C2z});
  if (name == "Ci") return generate_group({neg_I3});
  if (name == "Cs") return generate_group({sigma_h});

  // --- Order 3 ---
  if (name == "C3") return generate_group({Rz(2*M_PI/3)});

  // --- Order 4 ---
  if (name == "S4") return generate_group({Rz(M_PI/2) * sigma_h});  // S4 = C4 * sigma_h
  if (name == "C2v") return generate_group({C2z, sigma_v});
  if (name == "C2h") return generate_group({C2z, sigma_h});
  if (name == "D2") return generate_group({C2z, C2x});

  // --- Order 6 ---
  if (name == "S6") return generate_group({Rz(M_PI/3) * sigma_h}); // S6 = C6*sigma_h = C3*i
  if (name == "C3v") return generate_group({Rz(2*M_PI/3), sigma_v});
  if (name == "C3h") return generate_group({Rz(2*M_PI/3), sigma_h});
  if (name == "D3") return generate_group({Rz(2*M_PI/3), C2x});

  // --- Order 8 ---
  if (name == "D2h") return generate_group({C2z, C2x, sigma_h});
  if (name == "D2d") return generate_group({Rz(M_PI/2) * sigma_h, C2x});  // S4z, C2x

  // --- Order 10 ---
  if (name == "D5") return generate_group({Rz(2*M_PI/5), C2x});

  // --- Order 12 ---
  if (name == "D6") return generate_group({Rz(M_PI/3), C2x});
  if (name == "D3h") return generate_group({Rz(2*M_PI/3), C2x, sigma_h});
  if (name == "D3d") return generate_group({Rz(M_PI/3) * sigma_h, C2x}); // S6z, C2x
  if (name == "T") return generate_tetrahedral_group();

  // --- Order 20 ---
  if (name == "D5h") return generate_group({Rz(2*M_PI/5), C2x, sigma_h});
  if (name == "D5d") return generate_group({Rz(M_PI/5) * sigma_h, C2x}); // S10z, C2x

  // --- Order 24 ---
  if (name == "D6h") return generate_group({Rz(M_PI/3), C2x, sigma_h});
  if (name == "D6d") return generate_group({Rz(M_PI/6) * sigma_h, C2x}); // S12z, C2x
  if (name == "Td") {
    coord3d c3(1,1,1); c3 /= c3.norm();
    return generate_group({Raxis(c3, 2*M_PI/3), Rz(M_PI/2) * sigma_h}); // C3(111), S4z
  }
  if (name == "Th") {
    coord3d c3(1,1,1); c3 /= c3.norm();
    return generate_group({Raxis(c3, 2*M_PI/3), C2z, neg_I3});
  }

  // --- Order 60 ---
  if (name == "I") return generate_icosahedral_group();

  // --- Order 120 ---
  if (name == "Ih") {
    auto I_group = generate_icosahedral_group();
    // Add inversion to generators
    vector<matrix3d> gens;
    for (auto& m : I_group) gens.push_back(m);
    gens.push_back(neg_I3);
    return generate_group(gens, 120);
  }

  // Fallback: identity only
  return {I3};
}

}  // anonymous namespace


Representation3D Symmetry::representation_3d() const
{
  Representation3D rep;
  int n = G.size();
  if (n == 0) return rep;

  // Special case: C1 (trivial group)
  if (n == 1) {
    rep.R = {I3};
    return rep;
  }

  // Step 1: Generate standard 3D matrices for this point group.
  PointGroup pg = point_group();
  vector<matrix3d> std_mats = standard_matrices(pg);

  // Verify we generated the right number of elements.
  if ((int)std_mats.size() != n) {
    fprintf(stderr, "representation_3d: group %s expected %d elements, got %zu\n",
            pg.to_string().c_str(), n, std_mats.size());
    // Fallback: identity for everything
    rep.R.assign(n, I3);
    return rep;
  }

  // Step 2: Build multiplication tables for both representations.

  // Permutation multiplication table.
  // perm_table[i][j] = k such that G[i]*G[j] == G[k].
  IDCounter<Permutation> perm_id;
  for (int i = 0; i < n; i++) perm_id.insert(G[i]);
  vector<vector<int>> perm_table(n, vector<int>(n));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      perm_table[i][j] = perm_id(G[i] * G[j]);

  // Matrix multiplication table.
  vector<vector<int>> mat_table = matrix_mult_table(std_mats);

  // Step 3: Classify elements by (order, orientation character).
  vector<int> perm_orders(n), perm_chars(n);
  for (int i = 0; i < n; i++) {
    perm_orders[i] = G[i].order();
    perm_chars[i] = reverses_orientation(G[i]) ? -1 : +1;
  }

  vector<int> mat_orders(n), mat_chars(n);
  for (int i = 0; i < n; i++) {
    mat_orders[i] = matrix_order(std_mats[i]);
    mat_chars[i] = det3(std_mats[i]) > 0 ? +1 : -1;
  }

  // Step 4: Find isomorphism between the two multiplication tables.
  // mapping[i] = j means permutation G[i] corresponds to matrix std_mats[j].
  vector<int> mapping = find_isomorphism(perm_table, mat_table,
                                         perm_orders, mat_orders,
                                         perm_chars, mat_chars);

  if (mapping.empty()) {
    fprintf(stderr, "representation_3d: failed to find isomorphism for group %s\n",
            pg.to_string().c_str());
    rep.R.assign(n, I3);
    return rep;
  }

  // Step 5: Build the result: R[i] = std_mats[mapping[i]].
  rep.R.resize(n);
  for (int i = 0; i < n; i++)
    rep.R[i] = std_mats[mapping[i]];

  return rep;
}
