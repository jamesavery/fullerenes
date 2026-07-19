#include "fullerenes/eisenstein.hh"

#include <cmath>
#include <string>

// (-1,1)  x___x (0,1)
//        / \ / \
// (-1,0)x---x---x (1,0)
//        \ / \ /
//  (0,-1) x---x (1,-1)
Eisenstein Eisenstein::unit[7] = {{1,0},{0,1},{-1,1},{-1,0},{0,-1},{1,-1},{1,0}};

// =====================================================================
// Number theory: norm-form representatives.
// =====================================================================

std::vector<Eisenstein> sector0_reps_of_norm(int N) {
  std::vector<Eisenstein> out;
  if (N == 0) { out.push_back(Eisenstein(0, 0)); return out; }
  if (N < 0) return out;   // not a norm: no reps (callers treat empty as "non-Loeschian")
  // a >= 0, b >= 0:  a^2 + a*b + b^2 == N
  for (int b = 0; 3 * b * b <= 4 * N; ++b) {
    const long disc = 4L * N - 3L * b * b;
    const long s    = (long)std::lround(std::sqrt((double)disc));
    for (long ds = -1; ds <= 1; ++ds) {
      const long st = s + ds;
      if (st < 0) continue;
      if (st * st != disc) continue;
      const long two_a = st - b;
      if (two_a < 0) continue;
      if (two_a % 2 != 0) continue;
      const int a = (int)(two_a / 2);
      out.push_back(Eisenstein(a, b));
      break;
    }
  }
  return out;
}

Eisenstein eisenstein_of_norm(int N) {
  auto all = sector0_reps_of_norm(N);
  if (all.empty())
    throw std::logic_error("eisenstein_of_norm: no solution for N=" + std::to_string(N));
  return all[0];
}

// =====================================================================
// Lattice isometries.
// =====================================================================

LatticeIsometry isometry_from_segments(Eisenstein a_from, Eisenstein b_from,
                                       Eisenstein a_to,   Eisenstein b_to) {
  const Eisenstein vf = b_from - a_from, vt = b_to - a_to;
  if (vf.norm2() == 0 || vf.norm2() != vt.norm2())
    throw std::logic_error("isometry_from_segments: segment vectors have unequal or zero norms");
  const Eisenstein u = vt.div(vf);
  if (u.norm2() != 1 || !(u * vf == vt))
    throw std::logic_error("isometry_from_segments: segment alignment is not an Eisenstein unit");
  return { u, a_to - u * a_from, /*reflect=*/false };
}

LatticeIsometry align(Eisenstein z_from, Eisenstein z_to) {
  const long N = z_from.norm2();
  if (N != z_to.norm2())
    throw std::logic_error(
        "align: norm mismatch z_from=(" + std::to_string(z_from.first) + ","
        + std::to_string(z_from.second) + ") |" + std::to_string(N)
        + "|, z_to=(" + std::to_string(z_to.first) + ","
        + std::to_string(z_to.second) + ") |" + std::to_string(z_to.norm2()) + "|");
  if (N == 0) {
    // Degenerate: both are zero; any unit works.
    return LatticeIsometry{ Eisenstein(1, 0), Eisenstein(0, 0), /*reflect=*/false };
  }

  // Rotation branch: U_rot = z_to * z_from.complex_conj() / N
  const Eisenstein num_rot = z_to * z_from.complex_conj();
  if (num_rot.first % N == 0 && num_rot.second % N == 0) {
    const Eisenstein U((int)(num_rot.first / N), (int)(num_rot.second / N));
    return LatticeIsometry{ U, Eisenstein(0, 0), /*reflect=*/false };
  }

  // Reflection branch: U_ref = z_to * z_from / N
  const Eisenstein num_ref = z_to * z_from;
  if (num_ref.first % N == 0 && num_ref.second % N == 0) {
    const Eisenstein U((int)(num_ref.first / N), (int)(num_ref.second / N));
    return LatticeIsometry{ U, Eisenstein(0, 0), /*reflect=*/true };
  }

  throw std::logic_error(
      "align: neither branch divisible by N=" + std::to_string(N)
      + " (z_from=(" + std::to_string(z_from.first) + "," + std::to_string(z_from.second)
      + "), z_to=(" + std::to_string(z_to.first) + "," + std::to_string(z_to.second)
      + "), num_rot=(" + std::to_string(num_rot.first) + "," + std::to_string(num_rot.second)
      + "), num_ref=(" + std::to_string(num_ref.first) + "," + std::to_string(num_ref.second)
      + "))");
}
