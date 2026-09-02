#pragma once

#include "fullerenes/graphview.hh"

#include <span>
#include <vector>

// =====================================================================
// hueckel -- simple Hückel (tight-binding) analysis of a carbon cage's
// pi system.
//
//   E. Hückel, "Quantentheoretische Beiträge zum Benzolproblem",
//   Z. Physik 70, 204-286 (1931).
//
// Every atom of a carbon cage contributes one pi electron, and the Hückel
// Hamiltonian is H = alpha*I + beta*A with A the adjacency matrix of the
// cage.  Its eigenvectors are the molecular orbitals and its eigenvalues
// are E_i = alpha + x_i*beta, so the whole analysis is a function of the
// adjacency spectrum x_1 >= ... >= x_N alone (GraphView::adjacency_spectrum;
// x_i in [-3,3] for a cubic cage).  Since beta < 0, the occupied orbitals
// are the ones with the LARGEST x, and a quantity "in units of beta" grows
// more binding as it grows more positive.
//
// The analysis is a pipeline of named steps, each of which is a function
// here and each of which a proof can cite:
//
//   x  --degenerate_levels-->  levels  --aufbau_fill-->  occupied levels
//      --{homo_index, classify_shell, homo_lumo_gap, pi_energy, babic_TRE}-->
//      the quantities, gathered by hueckel::analyze into an Analysis.
//
// The spectral graph invariants the legacy program reported alongside these
// (Estrada index, bipartivity, spectral moments) have no electron in them
// and are functions of the graph: see graphview.hh.
//
// This is a port of two legacy Fortran routines (Fullerene 4.5,
// Schwerdtfeger/Wirz/Avery), whose conventions it reproduces so its output
// stays comparable with the isomer database they produced; each function
// below cites the lines it ports.  Three DELIBERATE DEVIATIONS, all in
// quantities the database does not store (of the three it does store,
// NeHOMO and NedgeHOMO are integers and reproduced exactly, and HLgap to
// within the F7.5 column it is stored in -- IsomerDB::Entry::HLgap_column_tol):
//
//   * The eigensolver is the in-house cyclic Jacobi rather than the
//     legacy tred2l/tqlil.  The spectra agree to ~1e-12, but degeneracy
//     grouping at a fixed tolerance is a DISCONTINUOUS function of the
//     spectrum, so the level structure is validated against the database
//     (every isomer C20-C70), not proved.
//   * The Estrada index and the bipartivity are functions of the raw
//     spectrum here; hueckel.f:124-125 sums degeneracy-weighted level
//     representatives, so the legacy values carry the grouping tolerance
//     where these exact spectral definitions must not.
//   * resonance_energy measures against the ELECTRON count (the localized
//     reference is one unit of beta per electron) where hueckel.f:141
//     subtracts the atom count.  The two coincide for the neutral cage
//     this was ported for, and the electron form stays right when the
//     analysis is handed another filling.  babic_TRE keeps the atom count,
//     because that quantity is defined per atom.
//
// INHERITED IMPRECISION, not a deviation: pi_energy sums level
// representatives exactly as hueckel.f:122 does, so it carries the grouping
// tolerance (a level spans up to (degeneracy-1)*tol under chaining).  The
// exact per-orbital sum is the honest reading; changing it would change a
// number the legacy program published, so it is a decision, not a cleanup.
// =====================================================================

namespace hueckel {

// --- Parameters.  alpha and beta are the legacy program's DFT-calibrated
// values (hueckel.f:69), NOT textbook constants, and are kept verbatim so
// reported energies stay comparable with the database.
inline constexpr double alpha = -0.21;       // au
inline constexpr double beta  = -0.111;      // au; beta < 0, hence "largest x is lowest E"
inline constexpr double au_to_eV   = 27.2117;      // hueckel.f:157 (sic: the program's constant)
inline constexpr double au_to_kcal = 627.509541;   // hueckel.f:134
// Two tolerances that share a value in the Fortran (hueckel.f:69 vs the
// literal at :109) while answering different questions.
inline constexpr double degeneracy_tol = 1e-5;  // when two eigenvalues are one level
inline constexpr double lumo_zero_tol  = 1e-5;  // when the LUMO still counts as bonding
// Babic's topological-resonance-energy regression and its C60 reference.
//   D. Babic, "Topological resonance energy of fullerenes", ... hueckel.f:131,:133.
inline constexpr double TRE_a = 1.024296, TRE_b = 1.562211;
inline constexpr double TRE_C60 = 2.82066353359331501e-2;

// How the highest occupied level sits relative to the lowest empty one.
// (hueckel.f:109, dispatched at :152-155 to formats 1009/1010.)
enum class Shell {
  open,             // the HOMO is fractionally occupied
  properly_closed,  // filled HOMO, antibonding LUMO (x_LUMO <= -lumo_zero_tol)
  pseudo_closed     // filled HOMO, bonding LUMO    (x_LUMO >  -lumo_zero_tol)
};
const char* name(Shell s);

// One degenerate level of the spectrum, with its aufbau occupation.
struct Level {
  double x_rep;    // the level's representative eigenvalue: the first and, the
                   // spectrum being descending, the largest of its group --
                   // the convention hueckel.f:73-88 keeps
  int degeneracy;  // eigenvalues coincident within degeneracy_tol
  int electrons;   // 0 .. 2*degeneracy; 0 until aufbau_fill
};

// Every orbital of the level doubly occupied.
inline bool filled(const Level& L)  { return L.electrons == 2 * L.degeneracy; }
// The level is bonding: its orbitals lower the energy of an electron in them.
inline bool bonding(const Level& L) { return L.x_rep > -lumo_zero_tol; }

// --- The pipeline ----------------------------------------------------

// The spectrum grouped into degenerate levels: an eigenvalue within tol of
// its predecessor joins that level (chained, so a level may span more than
// tol), and the level keeps the first of its group as representative.
// @anchor hueckel-degenerate-levels
// @ref    hueckel.f:73-88 = spiral.f DualAnalyze:1911-1926
// @pre    descending: is_sorted_descending(x)
// @post   partitions: sum_of(result, &Level::degeneracy) == int64_t(x.size())
// @post   unoccupied: all_of(result, [](const Level& L){ return L.electrons == 0; })
std::vector<Level> degenerate_levels(std::span<const double> x, double tol = degeneracy_tol);

// The same levels with n_electrons filled in from the top, two per orbital.
// @anchor hueckel-aufbau-fill
// @ref    hueckel.f:101-119
// @pre    fits: 0 < n_electrons && n_electrons <= 2*sum_of(levels, &Level::degeneracy)
// @post   conserved: sum_of(result, &Level::electrons) == int64_t(n_electrons)
std::vector<Level> aufbau_fill(std::vector<Level> levels, int n_electrons);

// The highest occupied level.
// @anchor hueckel-homo-index
// @pre    occupied: any_of(levels, [](const Level& L){ return L.electrons > 0; })
// @post   in_range: 0 <= result && result < int64_t(levels.size())
int homo_index(std::span<const Level> levels);

// The shell the occupation makes.
// @anchor hueckel-classify-shell
// @ref    hueckel.f:109
// @pre    homo_in_range: 0 <= homo && homo < int64_t(levels.size())
// @pre    has_lumo: !filled(levels[homo]) || homo + 1 < int64_t(levels.size())
Shell classify_shell(std::span<const Level> levels, int homo);

// x_HOMO - x_LUMO, the spectral gap, in units of beta.  It is a gap between
// LEVELS, so it is at least the grouping tolerance; an open shell has no gap
// in this sense and the caller decides what to report for one (the database
// stores 0 -- spiral.f:292-296, its writer, not DualAnalyze which computes
// the gap for open shells too).
// @anchor hueckel-homo-lumo-gap
// @pre    has_lumo: homo + 1 < int64_t(levels.size())
// @post   at_least_tol: result >= degeneracy_tol
double homo_lumo_gap(std::span<const Level> levels, int homo);

// E_pi = sum_j n_j x_j over the levels, in units of beta.
// @anchor hueckel-pi-energy
// @ref    hueckel.f:122
double pi_energy(std::span<const Level> levels);

// The Hückel resonance energy: the pi energy less the localized reference of
// one unit of beta per electron (see the deviation note above).
// @anchor hueckel-resonance-energy
// @ref    hueckel.f:141
inline double resonance_energy(double E_pi, int n_electrons) { return E_pi - n_electrons; }

// Babic's topological resonance energy, per atom, in units of beta.
// @anchor hueckel-babic-tre
// @ref    hueckel.f:131
// @pre    positive: n_atoms > 0
inline double babic_TRE(double E_pi, int n_atoms) { return TRE_a * E_pi / n_atoms - TRE_b; }

// A quantity in units of beta, in physical units.  beta < 0, and an
// excitation is reported positive (hueckel.f:157).
inline double excitation_in_eV(double gap_in_beta) { return -gap_in_beta * beta * au_to_eV; }
inline double in_kcal(double x_in_beta)            { return  x_in_beta * beta * au_to_kcal; }

// --- The report ------------------------------------------------------

// What the pipeline computes about one pi system, gathered.  Every field is
// one of the functions above applied to `levels`; nothing here is computed
// that is not named there.
// @anchor hueckel-analysis
// @inv levels_partition_spectrum: sum_of(levels, &Level::degeneracy) == int64_t(x.size())
// @inv occupation_conserved:      sum_of(levels, &Level::electrons)  == int64_t(n_electrons)
// @inv gap_iff_open: iff(shell == Shell::open, gap == 0) &&
//          implies(shell != Shell::open, gap >= degeneracy_tol)
struct Analysis {
  std::vector<double> x;       // the adjacency spectrum, descending
  std::vector<Level> levels;   // grouped and filled
  int n_electrons;

  int    homo_level;    // index into levels of the highest occupied level
  Shell  shell;
  double gap;           // the spectral gap, or 0 for an open shell -- the
                        // database's convention (spiral.f:292-296)
  double excitation_eV;
  double E_pi;
  double E_resonance;
  double TRE;           // Babic, per atom
  double dTRE_C60;      // TRE - TRE(C60)
  double dTRE_C60_kcal;

  int homo_degeneracy() const { return levels[homo_level].degeneracy; }
  int homo_electrons()  const { return levels[homo_level].electrons; }
};

// The analysis of a spectrum with n_electrons electrons in it.  The electron
// count must leave the top level unfilled: a completely filled spectrum has
// no LUMO to classify the shell against or to measure a gap to.
// @anchor hueckel-analyze-spectrum
// @pre  descending: is_sorted_descending(x)
// @pre  electrons_fit: 0 < n_electrons && n_electrons < 2*int64_t(x.size())
Analysis analyze(std::vector<double> x, int n_electrons);

// The analysis of a neutral pi system: one pi electron per atom, which is the
// carbon cage's convention and nothing else's.
// @anchor hueckel-analyze-graph
// @pre  nonempty: g.N > 0
inline Analysis analyze(const GraphView& g) { return analyze(g.adjacency_spectrum(), g.N); }

}  // namespace hueckel
