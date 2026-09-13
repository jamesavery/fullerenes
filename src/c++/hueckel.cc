#include "fullerenes/hueckel.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// The Hückel pipeline (contracts on the declarations in hueckel.hh).

namespace hueckel {

const char* name(Shell s)
{
  switch(s){
    case Shell::open:            return "open";
    case Shell::properly_closed: return "properly closed";
    case Shell::pseudo_closed:   return "pseudo closed";
  }
  return "?";
}

std::vector<Level> degenerate_levels(std::span<const double> x, double tol)
{
  std::vector<Level> levels;
  levels.reserve(x.size());
  for(size_t i = 0; i < x.size(); i++){
    const bool joins = i > 0 && std::abs(x[i-1] - x[i]) < tol;
    if(joins) levels.back().degeneracy++;
    else      levels.push_back({x[i], 1, 0});
  }
  return levels;
}

std::vector<Level> aufbau_fill(std::vector<Level> levels, int n_electrons)
{
  int capacity = 0;
  for(const Level& L: levels) capacity += 2 * L.degeneracy;
  if(n_electrons <= 0 || n_electrons > capacity)
    throw std::invalid_argument("hueckel::aufbau_fill: " + std::to_string(n_electrons) +
                                " electrons must be positive and fit " + std::to_string(capacity/2) + " orbitals");

  int remaining = n_electrons;
  for(Level& L: levels){
    L.electrons = std::min(2 * L.degeneracy, remaining);
    remaining -= L.electrons;
  }
  return levels;
}

int homo_index(std::span<const Level> levels)
{
  for(int i = int(levels.size()) - 1; i >= 0; i--)
    if(levels[i].electrons > 0) return i;
  throw std::invalid_argument("hueckel::homo_index: no level is occupied");
}

// The LUMO of a filled HOMO -- the @pre has_lumo both readers share, checked
// because reading past the levels would return a plausible number instead of
// failing.
static const Level& lumo_of(std::span<const Level> levels, int homo, const char* who)
{
  if(homo + 1 >= int(levels.size()))
    throw std::logic_error(std::string("hueckel::") + who + ": the highest occupied level is the last "
                           "level, so there is no LUMO (@pre has_lumo violated)");
  return levels[homo + 1];
}

Shell classify_shell(std::span<const Level> levels, int homo)
{
  if(!filled(levels[homo])) return Shell::open;
  return bonding(lumo_of(levels, homo, "classify_shell"))? Shell::pseudo_closed : Shell::properly_closed;
}

double homo_lumo_gap(std::span<const Level> levels, int homo)
{
  return levels[homo].x_rep - lumo_of(levels, homo, "homo_lumo_gap").x_rep;
}

double pi_energy(std::span<const Level> levels)
{
  double E_pi = 0;
  for(const Level& L: levels) E_pi += L.x_rep * L.electrons;
  return E_pi;
}

Analysis analyze(std::vector<double> x, int n_electrons)
{
  for(size_t i = 1; i < x.size(); i++)
    if(x[i-1] < x[i])
      throw std::invalid_argument("hueckel::analyze: the spectrum is not descending -- x[" +
                                  std::to_string(i-1) + "] = " + std::to_string(x[i-1]) + " < x[" +
                                  std::to_string(i) + "] = " + std::to_string(x[i]));
  if(n_electrons <= 0 || int64_t(n_electrons) >= 2 * int64_t(x.size()))
    throw std::invalid_argument("hueckel::analyze: " + std::to_string(n_electrons) + " electrons in " +
                                std::to_string(x.size()) + " orbitals must be positive and leave the top "
                                "level unfilled (there would be no LUMO)");

  Analysis A;
  A.levels      = aufbau_fill(degenerate_levels(x, degeneracy_tol), n_electrons);
  A.homo_level  = homo_index(A.levels);
  A.shell       = classify_shell(A.levels, A.homo_level);
  // The database stores 0 for an open shell (spiral.f:292-296); the spectral
  // gap itself is defined only between a filled level and the one above it.
  A.gap           = A.shell == Shell::open? 0.0 : homo_lumo_gap(A.levels, A.homo_level);
  A.excitation_eV = excitation_in_eV(A.gap);
  A.E_pi          = pi_energy(A.levels);
  A.E_resonance   = resonance_energy(A.E_pi, n_electrons);
  A.TRE           = babic_TRE(A.E_pi, int(x.size()));   // per ATOM = per orbital of the pi system
  A.dTRE_C60      = A.TRE - TRE_C60;
  A.dTRE_C60_kcal = in_kcal(A.dTRE_C60);
  A.n_electrons   = n_electrons;
  A.x             = std::move(x);
  return A;
}

}  // namespace hueckel
