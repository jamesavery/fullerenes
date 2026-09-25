#include <gtest/gtest.h>
#include <string>
#include <algorithm>
#include "fullerenes/triangulation.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/spiral.hh"

using namespace std;

// C28 dual spiral
static const vector<int> C28_spiral = {5,5,5,6,5,6,5,6,5,5,5,5,5,5,5,6};

// Representative GC(k,l) pairs covering small to medium sizes.
// N_carbon = 2*(V0-2)*T = 28*T where T = k^2 + kl + l^2, V0 = 16.
static const pair<int,int> test_cases[] = {
  {1, 0},   // T=1,    N=28
  {2, 0},   // T=4,    N=112
  {2, 1},   // T=7,    N=196
  {3, 0},   // T=9,    N=252
  {3, 1},   // T=13,   N=364
  {3, 2},   // T=19,   N=532
  {4, 0},   // T=16,   N=448
  {4, 1},   // T=21,   N=588
  {4, 3},   // T=37,   N=1036
  {5, 0},   // T=25,   N=700
  {5, 2},   // T=39,   N=1092
  {6, 0},   // T=36,   N=1008
  {6, 1},   // T=43,   N=1204
  {7, 0},   // T=49,   N=1372
  {7, 3},   // T=67,   N=1876
  {8, 0},   // T=64,   N=1792
  {10, 0},  // T=100,  N=2800
  {10, 3},  // T=139,  N=3892
  {12, 5},  // T=199,  N=5572
  {15, 0},  // T=225,  N=6300
  {15, 7},  // T=274,  N=7672
  {20, 0},  // T=400,  N=11200
  {25, 0},  // T=625,  N=17500
  {26, 0},  // T=676,  N=18928
};
static const int n_cases = sizeof(test_cases) / sizeof(test_cases[0]);

// Validate that a Triangulation is a valid fullerene dual
static void check_fullerene_dual(const Triangulation& T) {
  EXPECT_TRUE(T.is_consistently_oriented()) << "Orientation is inconsistent";

  int Nf = T.N;
  int expected_triangles = 2 * (Nf - 2);
  EXPECT_EQ((int)T.triangles().size(), expected_triangles)
    << "Triangle count violates Euler formula";

  int deg5 = 0, deg6 = 0, other = 0;
  for(int u = 0; u < Nf; u++) {
    int d = T.degree(u);
    if(d == 5) deg5++;
    else if(d == 6) deg6++;
    else other++;
  }
  EXPECT_EQ(deg5, 12) << "Must have exactly 12 degree-5 nodes";
  EXPECT_EQ(other, 0) << "All nodes should be degree 5 or 6";
}

// Check that nb_list is a cyclic rotation of expected_list
static bool is_cyclic_rotation(const vector<node_t>& a, const vector<node_t>& b) {
  if(a.size() != b.size()) return false;
  int n = a.size();
  for(int offset = 0; offset < n; offset++){
    bool match = true;
    for(int j = 0; j < n; j++){
      if(a[j] != b[(j + offset) % n]){ match = false; break; }
    }
    if(match) return true;
  }
  return false;
}

// Test: graph -> spiral -> graph round-trip using the permutation.
// get_spiral returns perm where perm[i] = G-node at spiral position i.
// G' (wound-up) has node i at spiral position i.
// So G'.neighbours[i] should be the relabeled version of G.neighbours[perm[i]],
// i.e. G' = inv_perm(G).
TEST(SpiralRoundtrip, PermutationIsomorphism) {
  Triangulation C28dual(C28_spiral);

  for(int i = 0; i < n_cases; i++) {
    auto [k, l] = test_cases[i];
    int T = k*k + k*l + l*l;
    int N_carbon = 28 * T;
    SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ") N=" + to_string(N_carbon));

    // Generate graph via GC transform
    Triangulation G = C28dual.GCtransform(k, l);
    int N = G.N;


    // Compute canonical spiral with permutation
    vector<int> spiral_code;
    jumplist_t jumps;
    vector<vector<node_t>> permutations;
    G.get_spiral(spiral_code, jumps, permutations, true, true);

    ASSERT_FALSE(permutations.empty()) << "No spiral permutation found";
    const vector<node_t>& perm = permutations[0];

    ASSERT_EQ((int)perm.size(), N) << "Permutation size mismatch";

    // Build inverse permutation: inv_perm[perm[i]] = i
    vector<node_t> inv_perm(N);
    for(int u = 0; u < N; u++)
      inv_perm[perm[u]] = u;

    // Wind up the spiral to get G'
    Triangulation Gprime(spiral_code, jumps, true);

    ASSERT_EQ(G.N, Gprime.N)
      << "Node count mismatch: G has " << G.N << " vs G' has " << Gprime.N;

    // Verify G' = inv_perm(G):
    // For each node i in G', its neighbour list should be
    // {inv_perm[v] : v in G.neighbours[perm[i]]} in cyclic order.
    for(int u = 0; u < N; u++) {
      auto gprime_nb = Gprime.nbrs(u);
      auto g_nb = G.nbrs(perm[u]);

      ASSERT_EQ(gprime_nb.size(), g_nb.size())
        << "Degree mismatch at G'-node " << u << " (G-node " << perm[u] << ")";

      // Relabel G's neighbours through inv_perm
      vector<node_t> relabeled(g_nb.size());
      for(size_t j = 0; j < g_nb.size(); j++)
        relabeled[j] = inv_perm[g_nb[j]];

      // The canonical spiral may use a CW or CCW starting triple;
      // the windup always builds CW. So the round-trip may invert orientation.
      vector<node_t> relabeled_rev(relabeled.rbegin(), relabeled.rend());

      EXPECT_TRUE(is_cyclic_rotation(vector<node_t>(gprime_nb.begin(), gprime_nb.end()), relabeled) ||
                  is_cyclic_rotation(vector<node_t>(gprime_nb.begin(), gprime_nb.end()), relabeled_rev))
        << "Neighbour list mismatch at G'-node " << u << " (G-node " << perm[u] << ")";
    }
  }
}

// Test: graph -> spiral -> graph round-trip with CW_only=true.
// When restricted to CW starting triples, orientation is preserved through
// the round-trip, so we can check exact cyclic match (no reversal needed).
TEST(SpiralRoundtrip, PermutationIsomorphismCW) {
  Triangulation C28dual(C28_spiral);

  for(int i = 0; i < n_cases; i++) {
    auto [k, l] = test_cases[i];
    int T = k*k + k*l + l*l;
    int N_carbon = 28 * T;
    SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ") N=" + to_string(N_carbon));

    // Generate graph via GC transform
    Triangulation G = C28dual.GCtransform(k, l);
    int N = G.N;

    // Compute canonical CW-only spiral with permutation
    vector<int> spiral_code;
    jumplist_t jumps;
    vector<vector<node_t>> permutations;
    G.get_spiral(spiral_code, jumps, permutations, true, true, true);  // CW_only=true

    ASSERT_FALSE(permutations.empty()) << "No spiral permutation found";
    const vector<node_t>& perm = permutations[0];

    ASSERT_EQ((int)perm.size(), N) << "Permutation size mismatch";

    // Build inverse permutation: inv_perm[perm[i]] = i
    vector<node_t> inv_perm(N);
    for(int u = 0; u < N; u++)
      inv_perm[perm[u]] = u;

    // Wind up the spiral to get G'
    Triangulation Gprime(spiral_code, jumps, true);

    ASSERT_EQ(G.N, Gprime.N)
      << "Node count mismatch: G has " << G.N << " vs G' has " << Gprime.N;

    // Verify G' = inv_perm(G) with EXACT cyclic match (orientation preserved):
    for(int u = 0; u < N; u++) {
      auto gprime_nb = Gprime.nbrs(u);
      auto g_nb = G.nbrs(perm[u]);

      ASSERT_EQ(gprime_nb.size(), g_nb.size())
        << "Degree mismatch at G'-node " << u << " (G-node " << perm[u] << ")";

      // Relabel G's neighbours through inv_perm
      vector<node_t> relabeled(g_nb.size());
      for(size_t j = 0; j < g_nb.size(); j++)
        relabeled[j] = inv_perm[g_nb[j]];

      // CW_only ensures orientation is preserved — exact cyclic match required
      EXPECT_TRUE(is_cyclic_rotation(vector<node_t>(gprime_nb.begin(), gprime_nb.end()), relabeled))
        << "Neighbour list mismatch at G'-node " << u << " (G-node " << perm[u] << ")";
    }
  }
}

// Test: spiral -> graph -> spiral round-trip.
// Wind up a spiral, compute canonical spiral of the result, verify it matches.
TEST(SpiralRoundtrip, SpiralToGraphToSpiral) {
  Triangulation C28dual(C28_spiral);

  for(int i = 0; i < n_cases; i++) {
    auto [k, l] = test_cases[i];
    int T = k*k + k*l + l*l;
    int N_carbon = 28 * T;
    SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ") N=" + to_string(N_carbon));

    // Step 1: Generate graph and compute its canonical spiral
    Triangulation graph = C28dual.GCtransform(k, l);
    spiral_nomenclature sn1(graph, spiral_nomenclature::FULLERENE,
                            spiral_nomenclature::TRIANGULATION, true);

    // Step 2: Wind up from the spiral
    Triangulation from_spiral(sn1);

    // Step 3: Compute canonical spiral again from wound-up graph
    spiral_nomenclature sn2(from_spiral, spiral_nomenclature::FULLERENE,
                            spiral_nomenclature::TRIANGULATION, true);

    // Verify: canonical spiral data (jumps + spiral_code) and the name match.
    EXPECT_EQ(sn1.spiral, sn2.spiral)
      << "Spiral data differs after spiral->graph->spiral round-trip";
    EXPECT_EQ(sn1.to_string(), sn2.to_string());
  }
}

// The seven-vertex triangulation of the torus (K7): vertex i has the
// neighbours i+1, i+3, i+2, i+6, i+4, i+5 (mod 7) in rotation order, which is
// a consistent orientation since i -> i+1 is an orientation-preserving
// automorphism.  Every vertex has degree 6, and no spiral of a sphere
// triangulation describes a torus, so the general-spiral search must report
// failure by throwing SpiralSearchFailed -- the only way it can fail.
static Triangulation torus_K7() {
  return Triangulation(Graph{{1,3,2,6,4,5}, {2,4,3,0,5,6}, {3,5,4,1,6,0},
                             {4,6,5,2,0,1}, {5,0,6,3,1,2}, {6,1,0,4,2,3},
                             {0,2,1,5,3,4}});
}

// Rarest-special starts on an all-hexagonal triangulation: there is no
// vertex of degree other than 6, hence no starting triple.
TEST(SpiralSearchFailed, NoStartingTriple) {
  Triangulation T = torus_K7();
  try {
    T.get_general_spiral(/*rarest_start=*/true);
    FAIL() << "general spiral search returned on the torus";
  } catch(const SpiralSearchFailed& e) {
    EXPECT_EQ(e.reason, SpiralSearchFailed::NO_SPIRAL);
    EXPECT_EQ(e.N, 7);
    EXPECT_EQ(e.start[0], -1);
    EXPECT_TRUE(e.only_rarest_special);
    EXPECT_FALSE(e.CW_only);
  }
}

// All starts: every starting triple exists, and the general spiral from the
// first one tried does not close.
TEST(SpiralSearchFailed, NoClosingSpiral) {
  Triangulation T = torus_K7();
  try {
    T.get_general_spiral(/*rarest_start=*/false);
    FAIL() << "general spiral search returned on the torus";
  } catch(const SpiralSearchFailed& e) {
    EXPECT_EQ(e.reason, SpiralSearchFailed::NO_SPIRAL);
    EXPECT_EQ(e.N, 7);
    EXPECT_GE(e.start[0], 0);
    EXPECT_FALSE(e.only_rarest_special);
  }
}

// The nomenclature propagates the exception unchanged.
TEST(SpiralSearchFailed, NomenclaturePropagates) {
  EXPECT_THROW(spiral_nomenclature(torus_K7(), spiral_nomenclature::FULLERENE,
                                   spiral_nomenclature::TRIANGULATION, true),
               SpiralSearchFailed);
}

// ---------------------------------------------------------------------------
// Names: which graph a name describes, and when the search scheme is written.
// ---------------------------------------------------------------------------

// The triangulation whose vertex spiral has degree-5 vertices at the 1-based
// positions `pentagons` (and jumps as given, 0-based positions).
static Triangulation dual_from_rspi(int N_carbon, const vector<int>& pentagons,
                                    const jumplist_t& jumps = jumplist_t()) {
  vector<int> spiral(N_carbon/2 + 2, 6);
  for(int p: pentagons) spiral[p-1] = 5;
  return Triangulation(spiral, jumps);
}

static const vector<int> C60_Ih_rspi = {1,7,9,11,13,15,18,20,22,24,26,32};
static const string C60_Ih_indices = "1,7,9,11,13,15,18,20,22,24,26,32";

// Same vertex count and the same rotation at every vertex (a rotation is a
// cyclic order: where a row starts is not part of the graph).
static void expect_same_graph(const PlanarGraph& A, const PlanarGraph& B) {
  ASSERT_EQ(A.N, B.N);
  for(node_t u = 0; u < A.N; u++){
    auto a = A.nbrs(u), b = B.nbrs(u);
    EXPECT_TRUE(is_cyclic_rotation(vector<node_t>(a.begin(), a.end()),
                                   vector<node_t>(b.begin(), b.end())))
      << "rotation differs at vertex " << u;
  }
}

// The construction-scheme member records the graph that was named: a
// triangulation named as such carries "T", the cubic graph does not, and the
// fullerene named from its dual (the canonical fullerene name) does not.
// The test also checks that the constructors initialise the member: each
// round builds the "T" name in fresh heap storage, where an uninitialised
// member would hold whatever the allocator hands back instead of the scheme
// requested.
TEST(SpiralNomenclature, NameRecordsTheNamedGraphInAnInitialisedMember) {
  const Triangulation T = dual_from_rspi(60, C60_Ih_rspi);
  const PlanarGraph   C = T.dual_graph();

  for(int rep = 0; rep < 3; rep++){
    auto *sT = new spiral_nomenclature(T, spiral_nomenclature::FULLERENE,
                                       spiral_nomenclature::TRIANGULATION);
    EXPECT_EQ(sT->construction_scheme, spiral_nomenclature::TRIANGULATION);
    EXPECT_EQ(sT->to_string(), "[T:" + C60_Ih_indices + "]-fullerene");
    delete sT;

    spiral_nomenclature sC(C, spiral_nomenclature::FULLERENE, spiral_nomenclature::CUBIC);
    EXPECT_EQ(sC.construction_scheme, spiral_nomenclature::CUBIC);
    EXPECT_EQ(sC.to_string(), "[" + C60_Ih_indices + "]-fullerene");

    spiral_nomenclature sF = spiral_nomenclature::fullerene_from_dual(T);
    EXPECT_EQ(sF.construction_scheme, spiral_nomenclature::CUBIC);
    EXPECT_EQ(sF.naming_scheme, spiral_nomenclature::FULLERENE);
    EXPECT_EQ(sF.to_string(), "[" + C60_Ih_indices + "]-fullerene");
    EXPECT_EQ(FullereneDual(T).name().to_string(), sF.to_string());
  }
}

// CS_NONE: the graph's shape decides, and the member holds the decision.
TEST(SpiralNomenclature, SchemeNoneRecordsTheDecision) {
  const Triangulation T = dual_from_rspi(60, C60_Ih_rspi);
  spiral_nomenclature sT(T, spiral_nomenclature::FULLERENE);
  EXPECT_EQ(sT.construction_scheme, spiral_nomenclature::TRIANGULATION);
  EXPECT_EQ(sT.to_string(), "[T:" + C60_Ih_indices + "]-fullerene");

  spiral_nomenclature sC(PlanarGraph(T.dual_graph()), spiral_nomenclature::FULLERENE);
  EXPECT_EQ(sC.construction_scheme, spiral_nomenclature::CUBIC);
  EXPECT_EQ(sC.to_string(), "[" + C60_Ih_indices + "]-fullerene");

  // The square pyramid: neither cubic nor a triangulation -> its leapfrog dual
  // (whose apex has degree 8: leapfrog_dual sizes its rows for that).
  const PlanarGraph pyramid(Graph{{1,4,3}, {2,4,0}, {3,4,1}, {0,4,2}, {0,1,2,3}});
  spiral_nomenclature sL(pyramid);
  EXPECT_EQ(sL.construction_scheme, spiral_nomenclature::LEAPFROG);
  EXPECT_EQ(sL.to_string().substr(0, 3), "[LF");
  EXPECT_EQ(spiral_nomenclature(sL.to_string()).construction_scheme,
            spiral_nomenclature::LEAPFROG);
}

// Every form ever written parses to the same spiral, and each builds the right
// graph: Triangulation(sn) is the dual for a cubic name and the named graph for
// a "T" name; PlanarGraph(sn) is the named graph.
TEST(SpiralNomenclature, OldAndNewFormsParseToTheSameGraph) {
  const Triangulation T = dual_from_rspi(60, C60_Ih_rspi);
  const FullereneGraph F = T.dual_graph();
  for(const string& name : {"[" + C60_Ih_indices + "]-fullerene",
                            "[GS:" + C60_Ih_indices + "]-fullerene",
                            "[T:" + C60_Ih_indices + "]-fullerene",
                            "[T,GS:" + C60_Ih_indices + "]-fullerene"}){
    SCOPED_TRACE(name);
    const spiral_nomenclature sn(name);
    const bool triangulation_form = name.find("T") != string::npos;
    EXPECT_EQ(sn.construction_scheme, triangulation_form ? spiral_nomenclature::TRIANGULATION
                                                         : spiral_nomenclature::CUBIC);
    expect_same_graph(Triangulation(sn), T);
    expect_same_graph(FullereneGraph(sn), F);
    EXPECT_EQ(PlanarGraph(sn).N, triangulation_form ? T.N : F.N);
    // With or without "GS", the search scheme is the default one.
    EXPECT_EQ(sn.search_scheme, spiral_nomenclature::CANONICAL_GENERALIZED_SPIRAL);
    // Re-written in the current form: the default search scheme is not written.
    EXPECT_EQ(sn.to_string(), triangulation_form ? "[T:" + C60_Ih_indices + "]-fullerene"
                                                 : "[" + C60_Ih_indices + "]-fullerene");
  }
}

// A C380 dual wound up from a general spiral (the one in
// programs/c380-title.cc); its canonical general spiral needs jumps.  The
// default search scheme is not written even with jumps, the jumps being
// recognised by the "<jumps>; <numbers>" form, and a name with jumps
// round-trips through the graph in both the cubic and the "T" form.
TEST(SpiralNomenclature, JumpSpiralOmitsTheDefaultSearchScheme) {
  const vector<int> rspi = {45,70,71,82,83,110,119,120,144,184,185,192};
  const Triangulation T = dual_from_rspi(380, rspi, jumplist_t{{109,2}});
  check_fullerene_dual(T);
  ASSERT_EQ(T.N, 380/2 + 2);

  const spiral_nomenclature sG = spiral_nomenclature::fullerene_from_dual(T);
  ASSERT_FALSE(sG.spiral.jumps.empty());
  const string gs = sG.to_string();
  string jumps;
  for(const auto& [position, rotations] : sG.spiral.jumps)
    jumps += (jumps.empty() ? "" : ",") + to_string(position + 1) + "," + to_string(rotations);
  EXPECT_EQ(gs.substr(0, jumps.size() + 3), "[" + jumps + "; ");
  EXPECT_EQ(gs.find(':'), string::npos);

  // name -> triangulation -> name.
  const spiral_nomenclature parsed(gs);
  EXPECT_EQ(parsed.search_scheme, spiral_nomenclature::CANONICAL_GENERALIZED_SPIRAL);
  EXPECT_EQ(parsed.construction_scheme, spiral_nomenclature::CUBIC);
  EXPECT_EQ(parsed.spiral.jumps, sG.spiral.jumps);
  const Triangulation Tp(parsed);
  check_fullerene_dual(Tp);
  EXPECT_EQ(spiral_nomenclature::fullerene_from_dual(Tp).to_string(), gs);
  EXPECT_EQ(parsed.to_string(), gs);

  // The legacy "GS" tag names the same spiral and is dropped on re-writing.
  const spiral_nomenclature legacy("[GS:" + gs.substr(1));
  EXPECT_EQ(legacy.search_scheme, spiral_nomenclature::CANONICAL_GENERALIZED_SPIRAL);
  EXPECT_EQ(legacy.spiral.jumps, sG.spiral.jumps);
  EXPECT_EQ(legacy.to_string(), gs);

  // The "T" form of the same spiral names the triangulation Tp itself.
  const string tgs = "[T:" + gs.substr(1);
  const spiral_nomenclature parsedT(tgs);
  EXPECT_EQ(parsedT.construction_scheme, spiral_nomenclature::TRIANGULATION);
  EXPECT_EQ(parsedT.search_scheme, spiral_nomenclature::CANONICAL_GENERALIZED_SPIRAL);
  expect_same_graph(Triangulation(parsedT), Tp);
  EXPECT_EQ(PlanarGraph(parsedT).N, Tp.N);
  EXPECT_EQ(PlanarGraph(parsed).N, 380);
  EXPECT_EQ(parsedT.to_string(), tgs);
  EXPECT_EQ(spiral_nomenclature(Tp, spiral_nomenclature::FULLERENE,
                                spiral_nomenclature::TRIANGULATION).to_string(), tgs);

  // The compatibility search is always written, and its name round-trips to
  // the compatibility scheme and the same graph.
  const spiral_nomenclature sC = spiral_nomenclature::fullerene_from_dual(T, /*rarest=*/false);
  const string cs = sC.to_string();
  EXPECT_EQ(cs.substr(0, 4), "[CS:");
  const spiral_nomenclature parsedC(cs);
  EXPECT_EQ(parsedC.search_scheme, spiral_nomenclature::COMPATIBILITY_CANONICAL_SPIRAL);
  EXPECT_EQ(parsedC.spiral.jumps, sC.spiral.jumps);
  EXPECT_EQ(parsedC.to_string(), cs);
  // The graph it builds is the same fullerene: both canonical names agree.
  const Triangulation Tc(parsedC);
  check_fullerene_dual(Tc);
  EXPECT_EQ(spiral_nomenclature::fullerene_from_dual(Tc, /*rarest=*/false).to_string(), cs);
  EXPECT_EQ(spiral_nomenclature::fullerene_from_dual(Tc).to_string(), gs);
}

// A jump-free compatibility spiral also carries its "CS" tag and round-trips
// to the compatibility scheme.
TEST(SpiralNomenclature, JumpFreeCompatibilitySpiralIsTagged) {
  const Triangulation T = dual_from_rspi(60, C60_Ih_rspi);
  const spiral_nomenclature sC = spiral_nomenclature::fullerene_from_dual(T, /*rarest=*/false);
  ASSERT_TRUE(sC.spiral.jumps.empty());
  const string cs = sC.to_string();
  EXPECT_EQ(cs, "[CS:" + C60_Ih_indices + "]-fullerene");
  const spiral_nomenclature parsed(cs);
  EXPECT_EQ(parsed.search_scheme, spiral_nomenclature::COMPATIBILITY_CANONICAL_SPIRAL);
  EXPECT_EQ(parsed.construction_scheme, spiral_nomenclature::CUBIC);
  EXPECT_EQ(parsed.to_string(), cs);
  expect_same_graph(Triangulation(parsed), T);

  const spiral_nomenclature parsedT("[T,CS:" + C60_Ih_indices + "]-fullerene");
  EXPECT_EQ(parsedT.search_scheme, spiral_nomenclature::COMPATIBILITY_CANONICAL_SPIRAL);
  EXPECT_EQ(parsedT.construction_scheme, spiral_nomenclature::TRIANGULATION);
  EXPECT_EQ(parsedT.to_string(), "[T,CS:" + C60_Ih_indices + "]-fullerene");
}

// A name without a search tag parses to the default search scheme, for every
// graph type.
TEST(SpiralNomenclature, UntaggedNameParsesToTheDefaultSearch) {
  for(const string& name : {"[" + C60_Ih_indices + "]-fullerene",
                            "[T:" + C60_Ih_indices + "]-fullerene",
                            string("[4,4,4,4,4,4]-cage")}){
    SCOPED_TRACE(name);
    EXPECT_EQ(spiral_nomenclature(name).search_scheme, spiral_nomenclature::CANONICAL_GENERALIZED_SPIRAL);
    EXPECT_EQ(spiral_nomenclature(name).to_string(), name);
  }
}

// For a fullerene the jump-free canonical spirals of the two searches coincide
// (spiral.hh, 2.), which is why the untagged name is unambiguous.
TEST(SpiralNomenclature, JumpFreeSearchesCoincideOnFullerenes) {
  Triangulation C28dual(C28_spiral);
  for(int i = 0; i < n_cases; i++) {
    auto [k, l] = test_cases[i];
    if(28*(k*k + k*l + l*l) > 1100) continue;   // the all-starts search is quadratic
    SCOPED_TRACE("GC(" + to_string(k) + "," + to_string(l) + ")");
    Triangulation G = C28dual.GCtransform(k, l);
    const spiral_nomenclature gs = spiral_nomenclature::fullerene_from_dual(G, true);
    const spiral_nomenclature cs = spiral_nomenclature::fullerene_from_dual(G, false);
    if(gs.spiral.jumps.empty()){
      EXPECT_EQ(gs.spiral, cs.spiral);
      EXPECT_EQ("[CS:" + gs.to_string().substr(1), cs.to_string());
    }
  }
}

// ---------------------------------------------------------------------------
// The name parser refuses every malformed name with the defect named
// (spiral.hh, the name grammar, 4.).
// ---------------------------------------------------------------------------

static void expect_refused(const string& name, SpiralNameMalformed::defect_t defect) {
  SCOPED_TRACE(name);
  try {
    const spiral_nomenclature sn(name);
    ADD_FAILURE() << "accepted as " << sn.to_string(true);
  } catch(const SpiralNameMalformed& e) {
    EXPECT_EQ(e.defect, defect) << e.what();
    EXPECT_EQ(e.name, name);
    EXPECT_NE(string(e.what()).find(name), string::npos) << e.what();
  }
}

// C60-Ih's pentagon list with its first entry replaced.
static string c60_with_first(const string& first) {
  return "[" + first + C60_Ih_indices.substr(1) + "]-fullerene";
}

TEST(SpiralNameMalformed, Brackets) {
  using M = SpiralNameMalformed;
  expect_refused(C60_Ih_indices + "]-fullerene", M::BRACKETS);
  expect_refused("[" + C60_Ih_indices + "-fullerene", M::BRACKETS);
  expect_refused("]" + C60_Ih_indices + "[-fullerene", M::BRACKETS);
  expect_refused("[[" + C60_Ih_indices + "]]-fullerene", M::BRACKETS);
  expect_refused("[" + C60_Ih_indices + "]]-fullerene", M::BRACKETS);
  expect_refused("", M::BRACKETS);
}

TEST(SpiralNameMalformed, SchemeTag) {
  using M = SpiralNameMalformed;
  expect_refused("[X:" + C60_Ih_indices + "]-fullerene", M::SCHEME_TAG);
  expect_refused("[T,LF:" + C60_Ih_indices + "]-fullerene", M::SCHEME_TAG);
  expect_refused("[GS,CS:" + C60_Ih_indices + "]-fullerene", M::SCHEME_TAG);
  expect_refused("[T,T:" + C60_Ih_indices + "]-fullerene", M::SCHEME_TAG);
  expect_refused("[:" + C60_Ih_indices + "]-fullerene", M::SCHEME_TAG);
  expect_refused("[T,:" + C60_Ih_indices + "]-fullerene", M::SCHEME_TAG);
  expect_refused("[GS:T:" + C60_Ih_indices + "]-fullerene", M::SCHEME_TAG);
}

TEST(SpiralNameMalformed, Number) {
  using M = SpiralNameMalformed;
  expect_refused(c60_with_first("a"), M::NUMBER);
  expect_refused(c60_with_first("-1"), M::NUMBER);
  expect_refused(c60_with_first("+1"), M::NUMBER);
  expect_refused(c60_with_first("1.5"), M::NUMBER);
  expect_refused(c60_with_first("99999999999"), M::NUMBER);
  expect_refused("[1,," + C60_Ih_indices.substr(2) + "]-fullerene", M::NUMBER);
  expect_refused("[" + C60_Ih_indices + ",]-fullerene", M::NUMBER);
}

TEST(SpiralNameMalformed, SegmentCount) {
  using M = SpiralNameMalformed;
  expect_refused("[]-fullerene", M::SEGMENT_COUNT);
  expect_refused("[GS:]-fullerene", M::SEGMENT_COUNT);
  expect_refused("[5,1;7,1;" + C60_Ih_indices + "]-fullerene", M::SEGMENT_COUNT);
  expect_refused("[4,4,4,4,4,4;4,4,4,4,4,4;4,4,4,4,4,4]-cage", M::SEGMENT_COUNT);
  expect_refused("[1,2,3,4,5,6]-(4,5)-fulleroid", M::SEGMENT_COUNT);
}

TEST(SpiralNameMalformed, Jumps) {
  using M = SpiralNameMalformed;
  expect_refused("[GS:110;" + C60_Ih_indices + "]-fullerene", M::JUMPS);       // odd
  expect_refused("[GS:2,1;" + C60_Ih_indices + "]-fullerene", M::JUMPS);       // position < 3
  expect_refused("[GS:10,0;" + C60_Ih_indices + "]-fullerene", M::JUMPS);      // no rotation
  expect_refused("[GS:10,1,10,1;" + C60_Ih_indices + "]-fullerene", M::JUMPS); // not increasing
  expect_refused("[GS:10,1,9,1;" + C60_Ih_indices + "]-fullerene", M::JUMPS);
  expect_refused("[GS:;" + C60_Ih_indices + "]-fullerene", M::JUMPS);          // empty
  expect_refused("[GS:7,1;4,4,4,4,4,4]-cage", M::JUMPS);                        // past the spiral
}

TEST(SpiralNameMalformed, IndexRange) {
  using M = SpiralNameMalformed;
  expect_refused(c60_with_first("0"), M::INDEX_RANGE);
  expect_refused("[" + C60_Ih_indices.substr(0, C60_Ih_indices.size()-2) + "1500000000]-fullerene",
                 M::INDEX_RANGE);
  expect_refused("[0,1,2,3,4,5]-(4)-fulleroid", M::INDEX_RANGE);
}

TEST(SpiralNameMalformed, IndexOrder) {
  using M = SpiralNameMalformed;
  expect_refused("[7,1,9,11,13,15,18,20,22,24,26,32]-fullerene", M::INDEX_ORDER);
  expect_refused("[1,7,7,11,13,15,18,20,22,24,26,32]-fullerene", M::INDEX_ORDER);
  // One face listed as a square and as a pentagon.
  expect_refused("[1,2,3;3,5,6,7,8,9,10,11,12]-(4,5)-fulleroid", M::INDEX_ORDER);
}

TEST(SpiralNameMalformed, IndexCount) {
  using M = SpiralNameMalformed;
  expect_refused("[1,2,3]-fullerene", M::INDEX_COUNT);                       // 3 pentagons
  expect_refused("[" + C60_Ih_indices + ",40]-fullerene", M::INDEX_COUNT);   // 13
  expect_refused("[1,2,3,4,5]-(4)-fulleroid", M::INDEX_COUNT);               // Euler: 10 != 12
  expect_refused("[1,2;3,4,5,6,7,8,9,10,11,12]-(4,5)-fulleroid", M::INDEX_COUNT);  // 4+10
}

TEST(SpiralNameMalformed, Degrees) {
  using M = SpiralNameMalformed;
  expect_refused("[2,4,4,4,4,4,4,4]-cage", M::DEGREES);                     // degree below 3
  expect_refused("[3,3,3]-cage", M::DEGREES);                               // Euler: 9 != 12
  expect_refused("[4,4,4,4,4,4,4]-cage", M::DEGREES);                       // Euler: 14 != 12
  expect_refused("[1,2,3,4,5,6]-(2)-fulleroid", M::DEGREES);
  expect_refused("[1,2,3,4,5,6]-(6)-fulleroid", M::DEGREES);
  expect_refused("[1,2,3;4,5,6]-(4,4)-fulleroid", M::DEGREES);
  expect_refused("[1,2,3,4,5,6]-()-fulleroid", M::DEGREES);
}

TEST(SpiralNameMalformed, Suffix) {
  using M = SpiralNameMalformed;
  expect_refused("[" + C60_Ih_indices + "]-fullerine", M::SUFFIX);
  expect_refused("[" + C60_Ih_indices + "]-", M::SUFFIX);
  expect_refused("[" + C60_Ih_indices + "]--fullerene", M::SUFFIX);
  expect_refused("[" + C60_Ih_indices + "]-C60-C61-fullerene", M::SUFFIX);
  expect_refused("[" + C60_Ih_indices + "]-(5)-fullerene", M::SUFFIX);
  expect_refused("[1,2,3,4,5,6]-fulleroid", M::SUFFIX);
  expect_refused("[1,2,3,4,5,6]-(4)-(5)-fulleroid", M::SUFFIX);
  expect_refused("[1,2,3,4,5,6]-(4-fulleroid", M::SUFFIX);
  expect_refused("[1,2,3,4,5,6]-(4)_7-fulleroid", M::SUFFIX);
}

// The refusal is an std::invalid_argument, as the parser's refusals always were.
TEST(SpiralNameMalformed, IsAnInvalidArgument) {
  EXPECT_THROW(spiral_nomenclature("[1,2,3]-fullerene"), std::invalid_argument);
}

// Well-formed names in every grammatical shape still parse: the fulleroid
// face-degree group before or after the formula, a base-degree tail, the
// point-group prefix, and a cage.
TEST(SpiralNameMalformed, WellFormedNamesParse) {
  // The octahedron: six degree-4 vertices, the dual of the cube.
  const Triangulation O(vector<int>{4,4,4,4,4,4});
  const spiral_nomenclature sO(PlanarGraph(O.dual_graph()), spiral_nomenclature::FULLEROID,
                               spiral_nomenclature::CUBIC);
  const string written = sO.to_string();
  EXPECT_EQ(written, "[1,2,3,4,5,6]-(4)-fulleroid");

  for(const string& name : {written, string("[1,2,3,4,5,6]-(4)-C8-fulleroid"),
                            string("[1,2,3,4,5,6]-C8-(4)-fulleroid"),
                            string("[1,2,3,4,5,6]-(4)_6-fulleroid"),
                            string("Oh-[1,2,3,4,5,6]-(4)6-C8-fulleroid")}){
    SCOPED_TRACE(name);
    const spiral_nomenclature sn(name);
    EXPECT_EQ(sn.naming_scheme, spiral_nomenclature::FULLEROID);
    EXPECT_EQ(sn.face_degrees, vector<int>{4});
    EXPECT_EQ(sn.chemical_formula, name.find("C8") != string::npos ? "C8" : "");
    expect_same_graph(Triangulation(sn), O);
    EXPECT_EQ(PlanarGraph(sn).N, 8);                  // the cube
  }
  EXPECT_EQ(spiral_nomenclature("Oh-[1,2,3,4,5,6]-(4)6-C8-fulleroid").point_group, "Oh");

  const spiral_nomenclature cage("[T:4,4,4,4,4,4]-cage");
  EXPECT_EQ(cage.naming_scheme, spiral_nomenclature::CAGE);
  expect_same_graph(Triangulation(cage), O);
  expect_same_graph(Triangulation(spiral_nomenclature("[T:4,4,4,4,4,4]")), O);
}

// ---------------------------------------------------------------------------
// leapfrog_dual: rows sized for every vertex and face, and the round trip
// G -> leapfrog_dual(G) -> inverse_leapfrog_dual -> G.
// ---------------------------------------------------------------------------

// The wheel W_n: an n-cycle 0..n-1 and a hub n adjacent to all of it,
// oriented.  The hub has degree n and the rim face is an n-gon, so the
// leapfrog dual has a vertex of degree 2n and a face-centre of degree n.
static PlanarGraph wheel(int n) {
  vector<vector<node_t>> nbrs(n+1);
  for(int i = 0; i < n; i++) nbrs[i] = {(i+1)%n, n, (i+n-1)%n};
  for(int i = 0; i < n; i++) nbrs[n].push_back(i);
  return PlanarGraph(Graph(nbrs));
}

TEST(LeapfrogDual, RoundTrip) {
  for(int n : {4, 5, 7, 12}){   // W_3 is the tetrahedron, a triangulation
    SCOPED_TRACE("W_" + to_string(n));
    const PlanarGraph G = wheel(n);
    ASSERT_TRUE(G.is_consistently_oriented());
    const size_t F = G.compute_faces_oriented().size();

    const PlanarGraph L = G.leapfrog_dual();
    ASSERT_EQ(size_t(L.N), G.N + F);
    EXPECT_TRUE(L.is_consistently_oriented());
    EXPECT_TRUE(L.is_triangulation());
    for(node_t u = 0; u < G.N; u++) EXPECT_EQ(L.degree(u), 2*G.degree(u));
    EXPECT_EQ(L.degree(n), 2*n);                       // the hub

    // The original vertices keep their labels 0..N-1.
    const Triangulation T(L);
    expect_same_graph(T.inverse_leapfrog_dual(), G);

    // Through the name: the "LF" name builds a graph isomorphic to G (the
    // spiral relabels it; the degree sequence and size are compared).  The
    // spiral search handles vertex degrees up to SpiralSearchFailed::degree_limit,
    // so the hub of W_n qualifies for 2n <= 16 and is refused by name above.
    if(2*n > SpiralSearchFailed::degree_limit){
      try {
        const spiral_nomenclature sn(G);
        ADD_FAILURE() << "named as " << sn.to_string();
      } catch(const SpiralSearchFailed& e) {
        EXPECT_EQ(e.reason, SpiralSearchFailed::DEGREE_LIMIT) << e.what();
        EXPECT_EQ(e.max_degree, 2*n);
      }
      continue;
    }
    const spiral_nomenclature sn(G);
    ASSERT_EQ(sn.construction_scheme, spiral_nomenclature::LEAPFROG);
    const PlanarGraph H(spiral_nomenclature(sn.to_string()));
    ASSERT_EQ(H.N, G.N);
    EXPECT_TRUE(H.is_consistently_oriented());
    vector<int> dG, dH;
    for(node_t u = 0; u < G.N; u++){ dG.push_back(G.degree(u)); dH.push_back(H.degree(u)); }
    sort(dG.begin(), dG.end());
    sort(dH.begin(), dH.end());
    EXPECT_EQ(dH, dG);
  }
}

// The n-gonal bipyramid: an n-cycle 0..n-1 and two apexes n, n+1 each
// adjacent to all of it, oriented.  It is a triangulation of the sphere whose
// apexes have degree n.
static Triangulation bipyramid(int n) {
  vector<vector<node_t>> nbrs(n+2);
  for(int i = 0; i < n; i++) nbrs[i] = {(i+1)%n, n, (i+n-1)%n, n+1};
  for(int i = 0; i < n; i++) nbrs[n].push_back(i);
  for(int i = n-1; i >= 0; i--) nbrs[n+1].push_back(i);
  return Triangulation(Graph(nbrs));
}

// A vertex degree above the search's limit is refused by name, stating the
// degree found, from every entry point of the search.
TEST(SpiralSearchFailed, DegreeLimit) {
  const int n = SpiralSearchFailed::degree_limit + 1;
  const Triangulation T = bipyramid(n);
  ASSERT_TRUE(T.is_consistently_oriented());
  ASSERT_TRUE(T.is_triangulation());
  for(bool rarest : {true, false}){
    SCOPED_TRACE(rarest ? "rarest-special starts" : "all starts");
    try {
      T.get_general_spiral(rarest);
      ADD_FAILURE() << "general spiral search returned on a vertex of degree " << n;
    } catch(const SpiralSearchFailed& e) {
      EXPECT_EQ(e.reason, SpiralSearchFailed::DEGREE_LIMIT) << e.what();
      EXPECT_EQ(e.N, n + 2);
      EXPECT_EQ(e.max_degree, n);
      EXPECT_NE(string(e.what()).find("degree " + to_string(n)), string::npos) << e.what();
    }
  }
  EXPECT_THROW(spiral_nomenclature(T, spiral_nomenclature::CAGE,
                                   spiral_nomenclature::TRIANGULATION), SpiralSearchFailed);

  // At the limit the search runs and the name round-trips.
  const Triangulation B = bipyramid(SpiralSearchFailed::degree_limit);
  const spiral_nomenclature sn(B, spiral_nomenclature::CAGE, spiral_nomenclature::TRIANGULATION);
  EXPECT_EQ(Triangulation(spiral_nomenclature(sn.to_string())).N, B.N);
}
