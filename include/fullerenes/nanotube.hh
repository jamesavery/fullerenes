#pragma once
// The (5,0) carbon nanotube capped by half-dodecahedra, as a fullerene
// dual: the thinnest capped nanotube, and the extremal isomer for the
// intrinsic diameter of the fullerene surface at its size (a closed curve
// encircling a fullerene surface is at least five hexagons long, so by the
// coarea inequality no isomer of the same area has a longer geodesic than
// this tube, up to the cap depth).  The series C20, C30, C40, ... : each
// ring adds ten carbons (five dual vertices).
//
// Dual layout, n_rings >= 0: vertex 0 the north pole, 1..5 the north cap
// (degree 5), then n_rings rings of five degree-6 vertices, then the south
// cap 5 vertices (degree 5) and the south pole.  n_rings == 0 is the
// icosahedron, the dual of C20.  Neighbour rings are built CCW as seen
// from outside -- consecutive neighbours u, w of v bound the face (v, u,
// w) -- and the construction is checked by is_consistently_oriented, which
// throws rather than returning a mis-oriented graph.  (The benchmark this
// was promoted from listed both pole rings in the opposite sense, behind
// an assertion compiled out in release builds.)

#include "fullerenes/graph.hh"

#include <stdexcept>

// The dual of the (5,0) nanotube with n_rings hexagon rings: 12 + 5 n_rings
// vertices, the dual of C(20 + 10 n_rings).
inline Graph nanotube_50_dual(int n_rings) {
  if (n_rings < 0) throw std::invalid_argument("nanotube_50_dual: n_rings < 0");
  const int N = 12 + 5 * n_rings;
  auto mod5 = [](int i) { return ((i % 5) + 5) % 5; };
  auto cn  = [&](int i) { return 1 + mod5(i); };                    // north cap
  auto rng = [&](int j, int i) { return 6 + 5 * j + mod5(i); };     // ring j
  auto cs  = [&](int i) { return 6 + 5 * n_rings + mod5(i); };      // south cap
  const int pole_N = 0, pole_S = 11 + 5 * n_rings, last = n_rings - 1;
  // The vertex below cap slot i of the north cap, and above cap slot i of
  // the south cap: the first / last ring, or the other cap when there is
  // no ring.
  auto below_north = [&](int i) { return n_rings ? rng(0, i) : cs(i); };
  auto above_south = [&](int i) { return n_rings ? rng(last, i) : cn(i); };

  Graph adj(N, GRAPH_DMAX);
  // The pole rings run against the cap index: cap slot i's ring lists
  // cn(i+1) before the vertices below it, so the face (pole_N, cn(i+1),
  // cn(i)) is the one both agree on.
  for (int i = 4; i >= 0; i--) adj.push_back(pole_N, cn(i));
  for (int i = 0; i < 5; i++)
    adj.assign_row(cn(i), {pole_N, cn(i + 1), below_north(i + 1), below_north(i), cn(i - 1)});
  for (int j = 0; j < n_rings; j++)
    for (int i = 0; i < 5; i++) {
      const int up_l  = j ? rng(j - 1, i - 1) : cn(i - 1);
      const int up_r  = j ? rng(j - 1, i) : cn(i);
      const int dn_r  = (j < last) ? rng(j + 1, i + 1) : cs(i + 1);
      const int dn_l  = (j < last) ? rng(j + 1, i) : cs(i);
      adj.assign_row(rng(j, i), {up_l, up_r, rng(j, i + 1), dn_r, dn_l, rng(j, i - 1)});
    }
  for (int i = 0; i < 5; i++)
    adj.assign_row(cs(i), {pole_S, cs(i - 1), above_south(i - 1), above_south(i), cs(i + 1)});
  for (int i = 0; i < 5; i++) adj.push_back(pole_S, cs(i));

  Graph G(adj);
  if (!G.is_consistently_oriented())
    throw std::logic_error("nanotube_50_dual: construction is not consistently oriented");
  return G;
}
