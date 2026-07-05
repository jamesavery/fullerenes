// Unfolding::straighten_lines must reduce the cone outline to a SIMPLE polygon (the
// "cone-star") for every fullerene. This checks it over the full buckygen enumeration
// of a few sizes; the exhaustive C20-C110 sweep (4.03M isomers) lives in the
// fold-unfold sub-project's validate_conestar.

#include "fullerenes/triangulation.hh"
#include "fullerenes/unfold.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/eisenstein.hh"

#include <gtest/gtest.h>

using namespace std;

static int sgnw(Eisenstein u, Eisenstein v){ long w=wedge(u,v); return (w>0)-(w<0); }
// Do the open segments a-b and c-d properly cross?
static bool seg_cross(Eisenstein a, Eisenstein b, Eisenstein c, Eisenstein d){
  int o1=sgnw(b-a,c-a), o2=sgnw(b-a,d-a), o3=sgnw(d-c,a-c), o4=sgnw(d-c,b-c);
  return o1!=o2 && o3!=o4 && o1 && o2 && o3 && o4;
}
// Is the cyclic outline a simple polygon (no coincident vertices, no crossing edges)?
static bool is_simple(const vector<pair<Eisenstein,node_t>>& O){
  int n=O.size();
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){
    if(O[i].first==O[j].first) return false;                       // coincident vertices
    if(j==i+1 || (i==0 && j==n-1)) continue;                       // adjacent edges share a vertex
    if(seg_cross(O[i].first,O[(i+1)%n].first,O[j].first,O[(j+1)%n].first)) return false;
  }
  return true;
}

// Straighten every C_N isomer and require a simple cone-star each time. The
// Unfolding(Triangulation) ctor now establishes the cones-first precondition itself
// (via sort_flat_last), so no manual relabelling is needed at the call site.
static void check_all(int N){
  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Graph G; long count=0;
  while(BuckyGen::next_fullerene(Q,G)){
    Triangulation T(G);
    Unfolding S = Unfolding(T).straighten_lines();
    ASSERT_TRUE(is_simple(S.outline)) << "C" << N << " isomer " << count << ": cone-star not simple";
    count++;
  }
  BuckyGen::stop(Q);
  ASSERT_GT(count, 0) << "no C" << N << " isomers enumerated";
}

TEST(StraightenLines, C20){ check_all(20); }
TEST(StraightenLines, C28){ check_all(28); }
TEST(StraightenLines, C60){ check_all(60); }
