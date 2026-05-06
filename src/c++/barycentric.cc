#include "fullerenes/barycentric.hh"

namespace {

long gcd_long(long a, long b) {
  if (a < 0) a = -a;
  if (b < 0) b = -b;
  while (b) { long t = a % b; a = b; b = t; }
  return a;
}

}  // namespace

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
