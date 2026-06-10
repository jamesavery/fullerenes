#include "fullerenes/eisenstein.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

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
  if (all.empty()) {
    std::fprintf(stderr, "eisenstein_of_norm: no solution for N=%d\n", N);
    std::abort();
  }
  return all[0];
}

// =====================================================================
// D6 alignment.
// =====================================================================

D6Affine align(Eisenstein z_from, Eisenstein z_to) {
  const long N = z_from.norm2();
  if (N != z_to.norm2()) {
    std::fprintf(stderr,
        "align: norm mismatch z_from=(%d,%d) |%ld|, z_to=(%d,%d) |%d|\n",
        z_from.first, z_from.second, N,
        z_to.first, z_to.second, z_to.norm2());
    std::abort();
  }
  if (N == 0) {
    // Degenerate: both are zero; any unit works.
    return D6Affine{ Eisenstein(1, 0), /*reflect=*/false };
  }

  // Rotation branch: U_rot = z_to * z_from.complex_conj() / N
  const Eisenstein num_rot = z_to * z_from.complex_conj();
  if (num_rot.first % N == 0 && num_rot.second % N == 0) {
    const Eisenstein U((int)(num_rot.first / N), (int)(num_rot.second / N));
    return D6Affine{ U, /*reflect=*/false };
  }

  // Reflection branch: U_ref = z_to * z_from / N
  const Eisenstein num_ref = z_to * z_from;
  if (num_ref.first % N == 0 && num_ref.second % N == 0) {
    const Eisenstein U((int)(num_ref.first / N), (int)(num_ref.second / N));
    return D6Affine{ U, /*reflect=*/true };
  }

  std::fprintf(stderr,
      "align: neither branch divisible by N=%ld "
      "(z_from=(%d,%d), z_to=(%d,%d), num_rot=(%d,%d), num_ref=(%d,%d))\n",
      N, z_from.first, z_from.second, z_to.first, z_to.second,
      num_rot.first, num_rot.second, num_ref.first, num_ref.second);
  std::abort();
}
