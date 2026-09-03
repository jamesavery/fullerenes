#include "fullerenes/graphview.hh"
#include "fullerenes/dense_linalg.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// The adjacency spectrum and the invariants that are functions of it alone
// (contracts on the declarations in graphview.hh).

vector<double> GraphView::adjacency_spectrum() const
{
  // The precondition is checked, not inferred: jacobi_eig's false means its
  // sweep budget tripped and nothing more, and an asymmetric matrix does not
  // trip it -- it converges quietly to numbers that are not the spectrum.
  if(!adjacency_is_symmetric())
    throw std::logic_error("GraphView::adjacency_spectrum: the adjacency is not symmetric "
                           "(@pre symmetric violated -- a graph bug, not an input error)");

  // A[u][v] = the NUMBER of u-v arcs, so a multigraph gets its own spectrum
  // and the Perron root stays the (multi-)degree.
  vector<double> A(size_t(N) * N, 0.0);
  for(node_t u = 0; u < N; u++)
    for(node_t v: nbrs(u)) A[size_t(u)*N + v] += 1.0;

  vector<double> lambda;
  if(!LinAlg::jacobi_eig(std::move(A), N, lambda, nullptr))
    throw std::runtime_error("GraphView::adjacency_spectrum: cyclic Jacobi did not converge within "
                             "its sweep budget on a symmetric " + to_string(N) + "x" + to_string(N) + " matrix");
  std::reverse(lambda.begin(), lambda.end());   // jacobi_eig sorts ascending
  return lambda;
}

std::array<double,16> spectral_moments(std::span<const double> spectrum)
{
  std::array<double,16> M{};
  for(double x: spectrum){
    double p = 1;
    for(size_t k = 0; k < M.size(); k++, p *= x) M[k] += p;
  }
  return M;
}

double estrada_index(std::span<const double> spectrum)
{
  double ee = 0;
  for(double x: spectrum) ee += std::exp(x);
  return ee;
}

// EE_even/EE: the even-walk share of the Estrada index.
double bipartivity(std::span<const double> spectrum)
{
  double EE_even = 0;
  for(double x: spectrum) EE_even += std::cosh(x);
  return EE_even / estrada_index(spectrum);
}
