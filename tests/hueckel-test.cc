// Correctness tests for the Hueckel module (include/fullerenes/hueckel.hh).
//
// The isomer database is the oracle: its NeHOMO / NedgeHOMO / HLgap fields
// are the legacy Fortran program's own Hueckel output for every isomer, so
// full-set agreement validates the port end to end. Spectral-moment
// identities pin the eigensolver wiring independently of the database.

#include "fullerene-test-main.hh"

#include "fullerenes/hueckel.hh"
#include "fullerenes/isomerdb.hh"

#include "fullerenes/graph.hh"

#include <array>
#include <cmath>
#include <cstdio>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

TEST(Hueckel, CubicSpectralInvariants) {
  if (!IsomerDB::is_installed(20)) GTEST_SKIP() << "C20 database not installed";
  IsomerDB db = IsomerDB::readPDB(20, false);
  FullereneGraph g = IsomerDB::makeIsomer(20, db.entries[0]);
  hueckel::Analysis H = hueckel::analyze(g);
  const double N = g.N;

  const std::array<double,16> M = spectral_moments(H.x);
  ASSERT_EQ(H.x.size(), size_t(g.N));
  EXPECT_NEAR(H.x[0], 3.0, 1e-10);  // Perron root of a connected cubic graph
  EXPECT_NEAR(M[0], N, 1e-9);
  EXPECT_NEAR(M[1], 0.0, 1e-9);      // trace A = 0
  EXPECT_NEAR(M[2], 3 * N, 1e-8);    // 2|E|, cubic
  EXPECT_NEAR(M[3], 0.0, 1e-8);      // triangle-free
}

TEST(Hueckel, MatchesIsomerDatabase) {
  ASSERT_FALSE(fullerene_test::sizes().empty()) << "--sizes= would make this sweep assert nothing";
  // A partial install narrows the sweep rather than failing it; a corpus that
  // holds none of the requested sizes skips.  The n_compared guard below is
  // what refuses a sweep that compared nothing.
  std::vector<int> installed;
  for (int N : fullerene_test::sizes()) if (IsomerDB::is_installed(N)) installed.push_back(N);
  if (installed.empty()) GTEST_SKIP() << "none of the requested sizes is installed";
  size_t n_compared = 0;
  for (int N : installed) {
    IsomerDB db = IsomerDB::readPDB(N, false);
    ASSERT_GT(db.entries.size(), 0u) << "empty database for N=" << N;

    for (size_t i = 0; i < db.entries.size(); i++) {
      const IsomerDB::Entry& e = db.entries[i];
      FullereneGraph g = IsomerDB::makeIsomer(N, e);
      hueckel::Analysis H = hueckel::analyze(g);

      // The spectral identities hold for every cubic triangle-free graph,
      // so the sweep is where they belong, not one C20 isomer.
      // Universal for every fullerene: cubic (M2), triangle- and square-free
      // (M3, M4), and exactly 12 pentagons (M5 = 10*12, M7 = 140*12 -- both
      // independent of N, so they pin the topology and not just the size).
      // Tolerance scaled by N and 3^k, ~100x the measured deviation.
      const auto mtol = [&](int k){ return 1e-12 * N * std::pow(3.0, k); };
      const std::array<double,16> M = spectral_moments(H.x);
      ASSERT_EQ(H.x.size(), size_t(N)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(M[1], 0.0, mtol(1)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(M[2], 3.0 * N, mtol(2)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(M[3], 0.0, mtol(3)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(M[4], 15.0 * N, mtol(4)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(M[5], 120.0, mtol(5)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(M[6], 93.0 * N - 120.0, mtol(6)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(M[7], 1680.0, mtol(7)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(H.x[0], 3.0, 1e-9) << "N=" << N << " isomer " << i + 1;
      // The class invariants (hueckel.hh @inv), on every isomer.
      int64_t deg = 0, ne = 0;
      for (const auto& L : H.levels) { deg += L.degeneracy; ne += L.electrons; }
      EXPECT_EQ(deg, int64_t(H.x.size())) << "N=" << N << " isomer " << i + 1;
      EXPECT_EQ(ne, int64_t(H.n_electrons)) << "N=" << N << " isomer " << i + 1;
      if (H.shell != hueckel::Shell::open)
        EXPECT_GE(H.gap, hueckel::degeneracy_tol) << "N=" << N << " isomer " << i + 1;
      EXPECT_EQ(H.homo_electrons(), int(e.NeHOMO)) << "N=" << N << " isomer " << i + 1;
      EXPECT_EQ(H.homo_degeneracy(), int(e.NedgeHOMO)) << "N=" << N << " isomer " << i + 1;
      EXPECT_NEAR(H.gap, e.HLgap, IsomerDB::Entry::HLgap_column_tol) << "N=" << N << " isomer " << i + 1;
      if (H.shell == hueckel::Shell::open)
        EXPECT_EQ(H.gap, 0.0) << "N=" << N << " isomer " << i + 1;
      n_compared++;
    }
  }
  printf("[          ] %zu isomers compared with the database\n", n_compared);
  ASSERT_GT(n_compared, 0u);
}

// The spectrum constructor is public and reusable, so its edges are part of
// the contract: the electron count must fit, and a completely filled set of
// levels has no LUMO to classify or measure a gap against.
// A simple undirected graph from an edge list; these cases need no database.
static Graph graph_of(int N, const std::vector<std::pair<int,int>>& edges) {
  Graph g(size_t(N), GRAPH_DMAX);
  for (auto [u, v] : edges) g.insert_edge({u, v});
  return g;
}

// The vocabulary's own guards, unmasked: through the constructor, the
// filled-shell check fires first and hides aufbau_fill's capacity check.
TEST(Hueckel, VocabularyGuards) {
  const std::vector<double> x{2.0, 1.0, 1.0, -1.0};
  const std::vector<hueckel::Level> L = hueckel::degenerate_levels(x);
  ASSERT_EQ(L.size(), 3u);                                  // 2 | 1,1 | -1
  EXPECT_EQ(L[1].degeneracy, 2);
  for (const auto& l : L) EXPECT_EQ(l.electrons, 0);        // @post unoccupied

  EXPECT_THROW(hueckel::aufbau_fill(L, 9), std::invalid_argument);   // capacity: 9 > 2*4
  EXPECT_THROW(hueckel::aufbau_fill(L, 0), std::invalid_argument);   // positivity
  EXPECT_THROW(hueckel::aufbau_fill(L, -1), std::invalid_argument);

  // Filling every orbital is fine for the aufbau itself -- and leaves a HOMO
  // with no LUMO, which BOTH readers of levels[homo+1] must say rather than
  // read past the end.
  const std::vector<hueckel::Level> full = hueckel::aufbau_fill(L, 8);
  EXPECT_EQ(hueckel::homo_index(full), 2);
  EXPECT_EQ(full[2].electrons, 2);
  int64_t ne = 0; for (const auto& l : full) ne += l.electrons;
  EXPECT_EQ(ne, 8);                                          // @post conserved
  EXPECT_THROW(hueckel::classify_shell(full, 2), std::logic_error);
  EXPECT_THROW(hueckel::homo_lumo_gap(full, 2), std::logic_error);

  // The fill is a function of its argument, not of the argument's history:
  // filling the already-filled levels again gives the smaller occupation.
  const std::vector<hueckel::Level> refilled = hueckel::aufbau_fill(full, 2);
  ne = 0; for (const auto& l : refilled) ne += l.electrons;
  EXPECT_EQ(ne, 2);
  EXPECT_EQ(hueckel::homo_index(refilled), 0);
  EXPECT_THROW(hueckel::homo_index(L), std::invalid_argument);   // nothing occupied
}

TEST(Hueckel, SpectrumConstructorBoundaries) {
  const std::vector<double> x{2.0, 1.0, 1.0, -1.0};   // 4 orbitals, a degenerate pair
  EXPECT_NO_THROW(hueckel::analyze(x, 1));
  const hueckel::Analysis last_level_open = hueckel::analyze(x, 7);   // open shell AT the last level: legal
  EXPECT_EQ(last_level_open.shell, hueckel::Shell::open);
  EXPECT_EQ(last_level_open.gap, 0.0);
  EXPECT_THROW(hueckel::analyze(x, 0), std::invalid_argument);
  EXPECT_THROW(hueckel::analyze(x, -2), std::invalid_argument);
  EXPECT_THROW(hueckel::analyze(x, 9), std::invalid_argument);            // 9 > 2*4
  EXPECT_THROW(hueckel::analyze(x, 8), std::invalid_argument);            // filled: no LUMO
  EXPECT_THROW(hueckel::analyze({1.0, 2.0}, 2), std::invalid_argument);   // not descending
  EXPECT_THROW(hueckel::analyze({}, 1), std::invalid_argument);

  // Aufbau and level structure on a hand-checked case: 1 electron in the
  // top orbital is an open shell with zero gap; 2 fills it and the gap runs
  // to the degenerate pair below.
  const hueckel::Analysis open = hueckel::analyze(x, 1);
  EXPECT_EQ(open.shell, hueckel::Shell::open);
  EXPECT_EQ(open.gap, 0.0);
  EXPECT_EQ(open.homo_electrons(), 1);
  EXPECT_EQ(open.homo_degeneracy(), 1);
  const hueckel::Analysis closed = hueckel::analyze(x, 2);
  EXPECT_NE(closed.shell, hueckel::Shell::open);
  EXPECT_NEAR(closed.gap, 1.0, 1e-12);
  EXPECT_EQ(closed.levels.size(), 3u);              // 2 | 1,1 | -1
  EXPECT_EQ(closed.levels[1].degeneracy, 2);
  EXPECT_NEAR(closed.E_pi, 4.0, 1e-12);             // 2 electrons at x = 2
}

// bipartivity == 1 exactly for a bipartite graph (its spectrum is symmetric,
// so every odd moment vanishes), and below 1 otherwise.
// Exact spectra, no database: the cube and K3,3 are bipartite (bipartivity
// exactly 1, every odd moment 0), the Petersen graph is not, and all three
// pin quantities the fullerene sweep has no oracle for.
TEST(Hueckel, KnownGraphs) {
  const hueckel::Analysis Q3 = hueckel::analyze(graph_of(8, {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},
                                        {0,4},{1,5},{2,6},{3,7}}));   // spectrum {3,1,1,1,-1,-1,-1,-3}
  EXPECT_EQ(bipartivity(Q3.x), 1.0);
  EXPECT_NEAR(estrada_index(Q3.x), 29.393807800447, 1e-9);          // 2cosh3 + 6cosh1
  EXPECT_EQ(Q3.levels.size(), 4u);
  EXPECT_EQ(Q3.homo_level, 1);
  EXPECT_EQ(Q3.homo_degeneracy(), 3);
  EXPECT_EQ(Q3.homo_electrons(), 6);
  EXPECT_EQ(Q3.shell, hueckel::Shell::properly_closed);
  EXPECT_NEAR(Q3.gap, 2.0, 1e-12);
  EXPECT_NEAR(Q3.E_pi, 12.0, 1e-12);
  EXPECT_NEAR(Q3.E_resonance, 4.0, 1e-12);
  EXPECT_NEAR(Q3.TRE, hueckel::babic_TRE(12.0, 8), 1e-12);
  EXPECT_NEAR(spectral_moments(Q3.x)[4], 168.0, 1e-9);         // 15N + 8*(6 four-cycles)
  for (size_t k = 1; k < 16; k += 2) EXPECT_NEAR(spectral_moments(Q3.x)[k], 0.0, 1e-6) << "odd moment " << k;

  const hueckel::Analysis K33 = hueckel::analyze(graph_of(6, {{0,3},{0,4},{0,5},{1,3},{1,4},{1,5},{2,3},{2,4},{2,5}}));
  EXPECT_EQ(bipartivity(K33.x), 1.0);                          // spectrum {3,0,0,0,0,-3}
  EXPECT_NEAR(estrada_index(K33.x), 24.135323991556, 1e-9);          // 2cosh3 + 4
  EXPECT_EQ(K33.shell, hueckel::Shell::open);       // a 4-fold level, half filled
  EXPECT_EQ(K33.gap, 0.0);
  EXPECT_EQ(K33.homo_degeneracy(), 4);
  EXPECT_EQ(K33.homo_electrons(), 4);
  EXPECT_NEAR(K33.E_pi, 6.0, 1e-12);

  // Petersen: cubic, girth 5, spectrum {3, 1x5, -2x4} -- not bipartite, and
  // it carries a fullerene's universal moments without being one.
  const hueckel::Analysis P = hueckel::analyze(graph_of(10, {{0,1},{1,2},{2,3},{3,4},{4,0},{0,5},{1,6},{2,7},{3,8},{4,9},
                                        {5,7},{7,9},{9,6},{6,8},{8,5}}));
  EXPECT_NEAR(bipartivity(P.x), 0.959482505474, 1e-12);
  EXPECT_NEAR(estrada_index(P.x), 34.218287198429, 1e-9);
  EXPECT_NEAR(spectral_moments(P.x)[4], 150.0, 1e-9);           // 15N, girth 5
  EXPECT_NEAR(spectral_moments(P.x)[5], 120.0, 1e-9);           // 10 * 12 pentagons
  EXPECT_NEAR(spectral_moments(P.x)[7], 1680.0, 1e-8);          // 140 * 12 pentagons
}

TEST(Hueckel, BipartivityDetectsBipartite) {
  const std::vector<double> symmetric{2.0, 0.5, -0.5, -2.0};
  EXPECT_NEAR(bipartivity(symmetric), 1.0, 1e-12);
  const std::array<double,16> M = spectral_moments(symmetric);
  for (size_t k = 1; k < M.size(); k += 2) EXPECT_NEAR(M[k], 0.0, 1e-9) << "odd moment " << k;
  // A fullerene has pentagons, hence odd cycles, hence bipartivity < 1.
  if (!IsomerDB::is_installed(60, true)) GTEST_SKIP() << "IPR C60 database not installed";
  IsomerDB db = IsomerDB::readPDB(60, true);
  const hueckel::Analysis H = hueckel::analyze(IsomerDB::makeIsomer(60, db.entries[0]));
  EXPECT_GT(bipartivity(H.x), 0.5);
  EXPECT_LT(bipartivity(H.x), 1.0);
}

TEST(Hueckel, C60IhReference) {
  if (!IsomerDB::is_installed(60, true)) GTEST_SKIP() << "IPR C60 database not installed";
  IsomerDB db = IsomerDB::readPDB(60, true);  // the single IPR C60 isomer = Ih
  ASSERT_EQ(db.entries.size(), 1u);
  FullereneGraph g = IsomerDB::makeIsomer(60, db.entries[0]);
  hueckel::Analysis H = hueckel::analyze(g);

  // Babic's TRE C60 reference constant is defined by this very isomer.
  // TRE_C60 is hueckel.f:133's own constant, so this pins E_pi through it.
  EXPECT_NEAR(H.dTRE_C60, 0.0, 1e-11);
  // Golden values from the legacy program's own C60-Ih transcript.
  EXPECT_NEAR(H.E_pi, 93.161603794370, 1e-9);
  EXPECT_NEAR(H.E_resonance, 33.161603794370, 1e-9);
  EXPECT_NEAR(H.excitation_eV, 2.2853040, 1e-6);
  EXPECT_NEAR(estrada_index(H.x), 197.450227717, 1e-6);
  EXPECT_NEAR(bipartivity(H.x), 0.992966035976, 1e-9);
  const double M_C60[8] = {60, 0, 180, 0, 900, 120, 5460, 1680};
  for (int k = 0; k < 8; k++) EXPECT_NEAR(spectral_moments(H.x)[k], M_C60[k], 1e-8) << "moment " << k;

  // Buckminsterfullerene: properly closed shell, 5-fold degenerate h_u HOMO.
  EXPECT_EQ(H.shell, hueckel::Shell::properly_closed);
  EXPECT_EQ(H.homo_degeneracy(), 5);
  EXPECT_EQ(H.homo_electrons(), 10);
  EXPECT_NEAR(H.gap, db.entries[0].HLgap, IsomerDB::Entry::HLgap_column_tol);

  // Pentagons make the graph non-bipartite.
  EXPECT_GT(estrada_index(H.x), 0.0);
  EXPECT_GT(bipartivity(H.x), 0.0);
  EXPECT_LT(bipartivity(H.x), 1.0);
}

// The default sweep is quick and complete; --sizes N,N,... widens it.
int main(int argc, char** argv) {
  return fullerene_test::run(argc, argv, {20, 24, 26, 28, 30, 32, 34, 36, 38, 40}, "hueckel");
}
