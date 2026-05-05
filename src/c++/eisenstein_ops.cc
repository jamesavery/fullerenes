#include "fullerenes/eisenstein_ops.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

Eisenstein from_cartesian(double x, double y) {
    // b = 2y/sqrt(3); a = x - b/2
    const double bf = 2.0 * y / std::sqrt(3.0);
    const long   b  = std::lround(bf);
    const double af = x - 0.5 * (double)b;
    const long   a  = std::lround(af);
    return Eisenstein((int)a, (int)b);
}

std::vector<Eisenstein> sector0_reps_of_norm(int N) {
    std::vector<Eisenstein> out;
    if (N == 0) { out.push_back(Eisenstein(0, 0)); return out; }
    if (N < 0) {
        std::fprintf(stderr, "sector0_reps_of_norm: negative norm %d\n", N);
        std::abort();
    }
    // a >= 0, b >= 0:  a^2 + a*b + b^2 = N
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

    // Rotation branch: U_rot = z_to * complex_conj(z_from) / N
    const Eisenstein num_rot = z_to * complex_conj(z_from);
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

IntBary integer_barycentric(Eisenstein p,
                            Eisenstein P0,
                            Eisenstein P1,
                            Eisenstein P2)
{
    IntBary r;
    const Eisenstein d01 = P1 - P0;
    const Eisenstein d02 = P2 - P0;
    r.denom = wedge(d01, d02);
    // n0 = wedge(P1-p, P2-p)
    // n1 = wedge(p-P0, P2-P0) = wedge(p-P0, d02)
    // n2 = wedge(P1-P0, p-P0) = wedge(d01, p-P0)
    r.n0 = wedge(P1 - p, P2 - p);
    r.n1 = wedge(p - P0, d02);
    r.n2 = wedge(d01, p - P0);
    return r;
}

static long gcd_long(long a, long b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { long t = a % b; a = b; b = t; }
    return a;
}

ReducedBary reduce_to_lowest_terms(IntBary bw) {
    long m0 = bw.n0, m1 = bw.n1, m2 = bw.n2, d = bw.denom;
    const int n_zero = (m0 == 0) + (m1 == 0) + (m2 == 0);
    if (n_zero == 1) {
        const long red = (m0 == 0) ? gcd_long(m1, m2)
                       : (m1 == 0) ? gcd_long(m0, m2)
                                   : gcd_long(m0, m1);
        if (red > 1) { m0 /= red; m1 /= red; m2 /= red; d /= red; }
    }
    return { m0, m1, m2, d };
}

coord3d barycentric_combine(ReducedBary b,
                            const coord3d& C0,
                            const coord3d& C1,
                            const coord3d& C2)
{
    coord3d v;
    for (int i = 0; i < 3; ++i)
        v[i] = ((double)b.m0 * C0[i]
              + (double)b.m1 * C1[i]
              + (double)b.m2 * C2[i]) / (double)b.denom;
    return v;
}
