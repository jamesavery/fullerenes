#include <gtest/gtest.h>
#include <cstdlib>
#include <map>
#include <array>
#include <string>
#include "fullerenes/triangulation.hh"
#include "fullerenes/symmetry.hh"
#include "fullerenes/isomerdb.hh"
#include "fullerenes/buckygen-wrapper.hh"

using namespace std;

// Override database path from environment if set.
static void configure_database_path() {
  if(const char* s = getenv("FULLERENE_DB"))
    IsomerDB::database_path = s;
}

// Read Nmax from environment variable, clamped to [20, hard_max].
static int get_nmax(const char* env_var, int default_val, int hard_max) {
  if(const char* s = getenv(env_var)){
    int v = atoi(s);
    if(v >= 20 && v % 2 == 0) return min(v, hard_max);
  }
  return default_val;
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

// Build group->count map from the All database for a given N.
// Returns only nontrivial entries.
static map<string, int> load_nontrivial_counts(int N) {
  IsomerDB db = IsomerDB::readPDB(N, false);
  map<string, int> counts;
  for(auto& e : db.entries){
    string group = PointGroup(string(e.group, 3)).to_string();
    if(group == "C1") continue;
    counts[group]++;
  }
  return counts;
}

// Test 1: Per-isomer point group validation.
// For each buckygen isomer, compute symmetry and compare against database RSPI->group.
TEST(SymmetryValidation, PointGroupPerIsomer) {
  configure_database_path();
  int Nmax = get_nmax("NMAX_FULL", 70, 70);
  int total = 0, checked = 0;

  for(int N = 20; N <= Nmax; N += 2){
    if(N == 22) continue;

    auto ref = load_nontrivial_rspi(N);
    if(ref.empty() && N > 20) continue;  // Skip if no database

    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Graph g;
    int count = 0;
    while(BuckyGen::next_fullerene(Q, g)){
      count++;
      Triangulation T(g);
      Symmetry S(T);
      string computed = S.point_group().to_string();

      // Get canonical RSPI to look up in reference
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
    checked++;
    fprintf(stderr, "N=%3d: %6d isomers checked (%zu nontrivial in ref)\n",
            N, count, ref.size());
  }
  fprintf(stderr, "\nTotal: %d isomers across %d sizes\n", total, checked);
}

// Test 2: Group count validation.
// For each N, tally symmetry groups from buckygen and compare against database counts.
TEST(SymmetryValidation, GroupCounts) {
  configure_database_path();
  int Nmax = get_nmax("NMAX_TALLY", 80, 130);
  int checked = 0;

  for(int N = 20; N <= Nmax; N += 2){
    if(N == 22) continue;

    auto ref_counts = load_nontrivial_counts(N);

    BuckyGen::buckygen_queue Q = BuckyGen::start(N, false, false);
    Graph g;
    map<string, int> computed_counts;
    int total_isomers = 0;
    while(BuckyGen::next_fullerene(Q, g)){
      total_isomers++;
      Triangulation T(g);
      Symmetry S(T);
      string group = S.point_group().to_string();
      if(group != "C1")
        computed_counts[group]++;
    }

    // Check each reference group appears with the right count
    for(auto& [group, expected] : ref_counts){
      EXPECT_EQ(computed_counts[group], expected)
        << "N=" << N << " group=" << group;
    }
    // Check no extra groups computed that aren't in reference
    for(auto& [group, count] : computed_counts){
      EXPECT_TRUE(ref_counts.count(group) > 0)
        << "N=" << N << " unexpected group " << group << " (count=" << count << ")";
    }

    checked++;
    fprintf(stderr, "N=%3d: %6d isomers, %zu nontrivial groups\n",
            N, total_isomers, computed_counts.size());
  }
  fprintf(stderr, "\n%d sizes checked\n", checked);
}
