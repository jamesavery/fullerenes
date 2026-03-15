#include <gtest/gtest.h>
#include <cstdlib>
#include <fstream>
#include <map>
#include <array>
#include <string>
#include "fullerenes/triangulation.hh"
#include "fullerenes/symmetry.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/planargraph.hh"
#include "fullerenes/buckygen-wrapper.hh"
#include "fullerenes/auxiliary.hh"  // pad_string

using namespace std;

static bool tally_only = false;

// Read Nmax from environment variable, clamped to [20, hard_max].
static int get_nmax(const char* env_var, int default_val, int hard_max) {
  if(const char* s = getenv(env_var)){
    int v = atoi(s);
    if(v >= 20 && v % 2 == 0) return min(v, hard_max);
  }
  return default_val;
}

// Check if the database file exists for a given N.
static bool database_exists(int N) {
  string filename = IsomerDB::database_path + "/All/c"
    + pad_string(to_string(N), 3, '0') + "all.database";
  return ifstream(filename).good();
}

// Build RSPI->group map from the All database for a given N.
// Returns only nontrivial entries (group != "C1").
static map<array<int,12>, string> load_nontrivial_rspi(int N) {
  IsomerDB db = IsomerDB::readPDB(N, false);
  map<array<int,12>, string> result;
  for(auto& e : db.entries){
    string group = PointGroup(string(e.group, 3)).to_string();
    if(group == "C1") continue;
    array<int,12> rspi;
    for(int i = 0; i < 12; i++) rspi[i] = e.RSPI[i] - 1;  // 1-indexed -> 0-indexed
    result[rspi] = group;
  }
  return result;
}

// Build reference group->count from the embedded symmetry data in IsomerDB.
// Returns empty map if N is out of range.
static map<string, int> load_embedded_counts(int N) {
  map<string, int> counts;
  vector<string> groups = IsomerDB::symmetries(N, false);
  for(auto& g : groups){
    string name = PointGroup(g).to_string();
    if(name == "C1") continue;
    int64_t c = IsomerDB::number_isomers(N, g, false);
    if(c > 0) counts[name] = c;
  }
  return counts;
}

// Per-isomer validation for one N: check each isomer's point group against DB.
static void check_per_isomer(int N, int& total) {
  auto ref = load_nontrivial_rspi(N);

  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Graph g;
  int count = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);
    Symmetry S(T);
    string computed = S.point_group().to_string();

    vector<int> rspi_vec;
    jumplist_t jumps;
    FullereneDual(T).get_rspi(rspi_vec, jumps);
    array<int,12> rspi;
    copy(rspi_vec.begin(), rspi_vec.end(), rspi.begin());

    auto it = ref.find(rspi);
    string expected = (it != ref.end()) ? it->second : "C1";

    EXPECT_EQ(computed, expected)
      << "N=" << N << " isomer #" << count;
  }
  total += count;
  fprintf(stderr, "N=%3d: %6d isomers, per-isomer check (%zu nontrivial in ref)\n",
          N, count, ref.size());
}

// Tally validation for one N: check group counts against embedded data.
static void check_tally(int N, int& total) {
  auto ref_counts = load_embedded_counts(N);

  BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
  Graph g;
  map<string, int> computed_counts;
  int count = 0;
  while(BuckyGen::next_fullerene(Q, g)){
    count++;
    Triangulation T(g);
    Symmetry S(T);
    string group = S.point_group().to_string();
    if(group != "C1")
      computed_counts[group]++;
  }

  for(auto& [group, expected] : ref_counts){
    EXPECT_EQ(computed_counts[group], expected)
      << "N=" << N << " group=" << group;
  }
  for(auto& [group, cnt] : computed_counts){
    EXPECT_TRUE(ref_counts.count(group) > 0)
      << "N=" << N << " unexpected group " << group << " (count=" << cnt << ")";
  }

  total += count;
  fprintf(stderr, "N=%3d: %6d isomers, tally check (%zu nontrivial groups)\n",
          N, count, computed_counts.size());
}

// Unified test: for each N, prefer per-isomer check (if DB file exists),
// otherwise fall back to tally check (using embedded count data).
// With --test_type=tally, always use tally check.
TEST(SymmetryValidation, PointGroups) {
  if(const char* s = getenv("FULLERENE_DB"))
    IsomerDB::database_path = s;

  int Nmax = get_nmax("NMAX", 60, 200);
  int total = 0, n_per_isomer = 0, n_tally = 0, n_skipped = 0;

  for(int N = 20; N <= Nmax; N += 2){
    if(N == 22) continue;

    if(!tally_only && database_exists(N)){
      check_per_isomer(N, total);
      n_per_isomer++;
    } else if(!load_embedded_counts(N).empty()){
      check_tally(N, total);
      n_tally++;
    } else {
      fprintf(stderr, "N=%3d: skipped (no reference data)\n", N);
      n_skipped++;
    }
  }

  fprintf(stderr, "\nTotal: %d isomers, %d per-isomer, %d tally, %d skipped\n",
          total, n_per_isomer, n_tally, n_skipped);
}

// ============================================================================
// Representation3D validation
// ============================================================================

// Determinant of 3x3 matrix.
static double det3(const matrix3d& M) {
  return M(0,0)*(M(1,1)*M(2,2) - M(1,2)*M(2,1))
       - M(0,1)*(M(1,0)*M(2,2) - M(1,2)*M(2,0))
       + M(0,2)*(M(1,0)*M(2,1) - M(1,1)*M(2,0));
}

// Verify Representation3D for a given symmetry object.
static void verify_rep3d(const Symmetry& S, const string& label) {
  Representation3D rep = S.representation_3d();
  int n = S.G.size();

  ASSERT_EQ((int)rep.R.size(), n) << label << ": wrong number of matrices";

  matrix3d I3 = matrix3d::unit_matrix();

  // 1. All matrices are orthogonal: R^T * R == I
  for (int i = 0; i < n; i++) {
    matrix3d RtR = rep.R[i].transpose() * rep.R[i];
    EXPECT_LT((RtR - I3).norm(), 1e-10)
      << label << ": R[" << i << "] not orthogonal";
  }

  // 2. det matches orientation character
  for (int i = 0; i < n; i++) {
    double d = det3(rep.R[i]);
    bool reverses = S.reverses_orientation(S.G[i]);
    if (reverses)
      EXPECT_NEAR(d, -1.0, 1e-10) << label << ": R[" << i << "] should be improper";
    else
      EXPECT_NEAR(d, +1.0, 1e-10) << label << ": R[" << i << "] should be proper";
  }

  // 3. Multiplication table consistency: R[i]*R[j] == R[k] when G[i]*G[j] == G[k]
  IDCounter<Permutation> pid;
  for (int i = 0; i < n; i++) pid.insert(S.G[i]);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      int k = pid(S.G[i] * S.G[j]);
      ASSERT_GE(k, 0) << label << ": perm table broken at " << i << "," << j;
      matrix3d prod = rep.R[i] * rep.R[j];
      EXPECT_LT((prod - rep.R[k]).norm(), 1e-8)
        << label << ": R[" << i << "]*R[" << j << "] != R[" << k << "]";
    }

  // 4. R[0] == identity (G[0] is always the identity permutation)
  EXPECT_LT((rep.R[0] - I3).norm(), 1e-10) << label << ": R[0] != identity";
}

// Test representation_3d on all C20..C60 fullerenes with nontrivial symmetry.
TEST(Representation3D, AllNontrivialUpToC60) {
  int Nmax = 60;
  int tested = 0, skipped_C1 = 0;
  map<string, int> group_counts;

  for (int N = 20; N <= Nmax; N += 2) {
    if (N == 22) continue;

    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Graph g;
    int idx = 0;
    while (BuckyGen::next_fullerene(Q, g)) {
      idx++;
      Triangulation T(g);
      Symmetry S(T);
      string pg = S.point_group().to_string();

      if (pg == "C1") { skipped_C1++; continue; }

      string label = "C" + to_string(N) + " #" + to_string(idx) + " " + pg;
      verify_rep3d(S, label);
      group_counts[pg]++;
      tested++;
    }
  }

  fprintf(stderr, "\nRepresentation3D: tested %d nontrivial isomers (skipped %d C1)\n",
          tested, skipped_C1);
  fprintf(stderr, "Groups found: ");
  for (auto& [pg, cnt] : group_counts)
    fprintf(stderr, "%s(%d) ", pg.c_str(), cnt);
  fprintf(stderr, "\n");
}

// Build a fullerene dual triangulation from N and isomer index (via IsomerDB).
static Triangulation make_sym_dual(int N, int idx, bool IPR = false) {
  IsomerDB db = IsomerDB::readPDB(N, IPR);
  FullereneGraph G = IsomerDB::makeIsomer(N, db.entries[idx]);
  PlanarGraph PG(G);
  return Triangulation(PG.dual_graph());
}

// Focused test on specific high-symmetry fullerenes.
TEST(Representation3D, HighSymmetry) {
  // C20 (Ih)
  {
    Triangulation T = make_sym_dual(20, 0);
    Symmetry S(T);
    EXPECT_EQ(S.point_group().to_string(), "Ih");
    verify_rep3d(S, "C20 Ih");
    fprintf(stderr, "C20 Ih: |G|=%zu, representation_3d OK\n", S.G.size());
  }

  // C60 Ih (IPR #1)
  {
    Triangulation T = make_sym_dual(60, 0, true);
    Symmetry S(T);
    EXPECT_EQ(S.point_group().to_string(), "Ih");
    verify_rep3d(S, "C60 Ih");
    fprintf(stderr, "C60 Ih: |G|=%zu, representation_3d OK\n", S.G.size());
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  // Check for --test_type=tally among remaining args.
  for(int i = 1; i < argc; i++){
    string arg = argv[i];
    if(arg == "--test_type=tally" || arg == "tally")
      tally_only = true;
  }

  return RUN_ALL_TESTS();
}
