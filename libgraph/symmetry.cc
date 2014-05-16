#include "symmetry.hh"
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


PointGroup::PointGroup(const Triangulation& T)
{
  // TODO: 
  // 1) Port the fortran symmetry detection routine, or
  // 2) Make the fortran routine callable from C++, or
  // 3) Implement the more efficient and general symmetry detection based on Brinkman's method.

  // Orders of fullerene point groups
  // 1: C1
  // 2: C2, Ci, Cs
  // 3: C3
  // 4: S4, D2, C2h, C2v
  // 6: S6, D3, C3h, C3v
  // 8: D2d, D2h
  // 10: D5
  // 12: T, D6, D3d, D3h
  // 20: D5d, D5h
  // 24: Th, Td, D6d, D6h
  // 60: I
  //120: Ih

}

string PointGroup::to_string() const {
  const char ts[7] = {'?','C','D','T','S','O','I'};
  const char rs[6] = {' ','v','h','d','i','s'};
  char result[4]   = {0,0,0,0};
  result[0] = ts[sym_type];
  result[1] = n>0? '0'+n : ' ';
  result[2] = rs[sym_reflection];
  return string(result);
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

bool Permutation::operator==(const Permutation& q) const {
  const vector<int> &p(*this);

  if(size() != q.size()) return false;
  for(int i=0;i<size();i++) if(p[i] != q[i]) return false;
  return true;
}


//////////////////////////////////////////////////////////////////////
//		  SYMMETRY-DETECTION IMPLEMENTATION
//////////////////////////////////////////////////////////////////////
vector<Permutation> Symmetry::tri_permutation(const vector<Permutation>& Gf) const {
  vector<Permutation> Gtri(Gf.size(),Permutation(triangles.size()));
  IDCounter<tri_t> tri_id;
    
  for(int i=0;i<triangles.size();i++) tri_id.insert(triangles[i].sorted());
    
  for(int j=0;j<Gf.size();j++)
    for(int i=0;i<triangles.size();i++){
      const tri_t &t(triangles[i]);
      tri_t tp(Gf[j][t[0]],Gf[j][t[1]],Gf[j][t[2]]);
      Gtri[j][i] = tri_id(tp.sorted());
    }
  return Gtri;
}

vector<Permutation> Symmetry::edge_permutation(const vector<Permutation>& Gf) const {
  vector<Permutation> Gedge(Gf.size(),Permutation(edge_set.size()));
  IDCounter<edge_t> edge_id;
    
  for(auto e(edge_set.begin());e!=edge_set.end();e++) edge_id.insert(*e);
    
  for(int j=0;j<Gf.size();j++){
    int i=0;
    for(auto e=edge_set.begin();e!=edge_set.end();e++,i++)
      Gedge[j][i] = edge_id(edge_t(Gf[j][e->first],Gf[j][e->second]));
  }
  return Gedge;
}

vector<Permutation> Symmetry::permutation_representation() const
{
  vector<Permutation> pi;

  for(node_t u=0;u<N;u++){
    const vector<node_t>& nu(neighbours[u]);
    if(nu.size() == S0[0]) // u has same degree as vertex 1: possible spiral start
      for(int i=0;i<nu.size();i++){
	node_t v = nu[i];
	const vector<node_t>& nv(neighbours[v]);
	if(nv.size() == S0[1]){ // v has same degree as vertex 2: still possible spiral start
	  vector<int> spiral,permutation;
	  jumplist_t  jumps;

	  node_t wCCW = nextCCW(dedge_t(u,v)), wCW = nextCW(dedge_t(u,v));

	  const vector<node_t> &nwCCW(neighbours[wCCW]), &nwCW(neighbours[wCW]);

	  if(nwCCW.size() == S0[2] && get_spiral_implementation(u,v,wCCW,spiral,jumps,permutation,false,S0))
	    pi.push_back(permutation);
	  if(nwCW.size() == S0[2] && get_spiral_implementation(u,v,wCW,spiral,jumps,permutation,false,S0))
	    pi.push_back(Permutation(permutation));
	}
      }
  }	    
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
    for(int j=1;j<order;j++){ 
      int I = pi[j][i];
      assert(I<M);
      if(seen[I]) continue;
      seen[I] = true;
      orbit_length++;
    }	
    int site_order = order/orbit_length; 
    assert(site_order <= 12);
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
  Triangulation piG(neighbours,true);

  for(node_t u=0;u<N;u++){
    const vector<node_t>& nu(neighbours[u]);
    for(int i=0;i<nu.size();i++) piG.neighbours[pi[u]][i] = pi[nu[i]];
  }
  if(piG.nextCCW(dedge_t(0,1)) == 2) return false;
  if(piG.nextCW(dedge_t(0,1))!=2){
    fprintf(stderr,"G.next(0,1) == {%d,%d} (CW,CCW)\n",
	    nextCW(dedge_t(0,1)),
	    nextCCW(dedge_t(0,1)));
    fprintf(stderr,"pi(G).next(0,1) == {%d,%d} (CW,CCW)\n",
	    piG.nextCW(dedge_t(0,1)),
	    piG.nextCCW(dedge_t(0,1)));
    cout << "pi(G).neighbours[0] = " << piG.neighbours[0] << ";\n"
	 << "pi(G).neighbours[1] = " << piG.neighbours[1] << ";\n";
    abort();
  }
  return true;
}

