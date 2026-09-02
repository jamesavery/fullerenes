#pragma once
// ============================================================================
// delaunay_cyclotomic.hh -- CyclotomicMetric: the exact metric policy of
// the FLATTENED KIS SURFACE, the third regime of the iDT machinery
// (BandedFloatMetric / ExactIntegerMetric, delaunay_view.hh; layer 2 of
// the cyclotomic-algebra-and-idt-extension debt entry,
// claude-projects/delaunay/refactor-debt.md; mathematics in
// claude-projects/delaunay/cyclotomic-idt.tex, whose sec. 4-6 carry the
// predicate forms, the flip transport and this policy's design).
//
// THE CARRY.  Three exact arrays beside the DCEL's float shadows:
//   - lsq[h]     per half-edge (twin-paired): the squared length in
//                Z[gamma], gamma = 2 cos(pi/15), at the x kLsqScale
//                convention;
//   - f_wedge[f] per face slot: the face's delta-normalized wedge w, with
//                16 Area^2 = (2 - gamma_2) w^2 -- carried, never
//                re-extracted (no integer square root exists off the
//                lattice);
//   - curv_k[v]  per vertex: the curvature index, kappa = k * pi/15 (the
//                pentagon-incidence count of a cubic corner; 0 at every
//                face centre).  Purely combinatorial and flip-invariant:
//                the flatness word costs no arithmetic at all.
// Every predicate consults the constructor-validated cyclotomic::Diamond,
// so a corrupt carry REFUSES (InvariantViolated through the view's latch)
// instead of mis-classifying; flips transport (lsq, both face wedges)
// through Diamond::flipped's exact divisions.
//
// THE Q-FRAME DEVELOPMENT (prepare_star / ear / commit_star /
// first_tie_side).  Flat-star removal develops the star's rim into the
// plane.  The module's points live in (1/5) Z[zeta_30], but no module
// anchor must be FOUND (the recovery problem the wedge carry exists to
// avoid): develop in the frame scaled by the conjugate of the unknown
// anchor -- q_i := kLsqScale * P_i * conj(P_0) -- and every q_i is an
// INTEGRAL Zeta30 value computable by ring arithmetic alone:
//     q_0 = Ls(0)   (real!),
//     q_{i+1} = q_i * (d_i + delta w_i) / (2 Ls(i)),
// with d_i = Ls(i) + Ls(i+1) - Lr(i) the apex dot form, w_i the fan
// face's carried wedge, delta = zeta - conj(zeta) = 2 zeta - gamma, and
// the division exact by the module geometry (a failure falsifies the
// carry and trips).  All downstream reads are invariant under the frame's
// common rotation-scale |P_0|^2 > 0: ear CCW and sector signs directly,
// the ear diagonal's lsq and the new faces' wedges after one exact
// division by the frame scale Ls(0).  A flat apex forces closure --
// q_k == q_0 -- checked, like the lattice development's Unclosed trip;
// unlike the lattice development there is NO anchor-orbit scan.
//
// ENTRY BOUNDARY.  derive_cyclotomic_kis_carry builds the carry from a
// FRESH kis DCEL (fullerene kis complex: face centres first, then cubic
// vertices; delaunay_alexandrov.hh's kis_metric shape) and VERIFIES it
// loudly -- edge classes against the float shadows at the integrality
// band, exactly one centre per face, curvature indices against the float
// cone angles and Gauss-Bonnet (sum k = 60) -- mirroring
// derive_exact_lsq_carry's discipline: entering the exact regime on a
// metric it does not describe would be a silent wrong answer.
//
// Host tier (the sign oracle's exact rung is host-tier by design).
// ============================================================================

#include "cyclotomic.hh"
#include "delaunay_view.hh"

#include <cmath>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace cyclotomic {

// The regime's curvature quantum: kappa = k * pi/15 at every vertex (a
// cubic corner's pentagon count k; 0 at face centres).
inline constexpr double kCurvatureQuantum = std::numbers::pi_v<double> / 15;

// Does the CCW sector from direction `from` to direction `to` subtend at
// most pi?  The exact mirror of delaunay_detail::lattice_sector_at_most_pi
// (its banner carries the argument), over Zeta30 directions: wedge > 0 is
// (0, pi); wedge == 0 with negative real inner product is exactly pi (the
// accepted tie); wedge == 0 with positive dot would be a 0 / 2pi sector,
// excluded by the fan premise (defensive reject).  `ok` false = the sign
// oracle refused (caller trips).
inline bool sector_at_most_pi(const Zeta30& from, const Zeta30& to,
                              bool& ok) {
  const SignOr sw = sign_real(wedge(from, to));
  if (!sw) { ok = false; return false; }
  if (*sw == Sign::Positive) return true;
  if (*sw == Sign::Negative) return false;
  // wedge == 0: conj(from)*to is real, its x-part IS the inner product.
  const SignOr sd = sign_real((from.conj() * to).x);
  if (!sd) { ok = false; return false; }
  return *sd == Sign::Negative;
}

// ---------------------------------------------------------------------------
// The policy.  Spans are caller-owned (CyclotomicKisCarry below); the
// pending-transport members make the policy STATEFUL -- one instance must
// live across a whole reduction run (the drivers pass the metric by
// reference throughout).
// ---------------------------------------------------------------------------
// The carry as the policy sees it: the three exact arrays (file banner) and
// the two development scratch arrays, as caller-owned views.
struct CarryViews {
  std::span<Real30> lsq;                 // [nh_cap] per half-edge, twin-paired
  std::span<Real30> f_wedge;             // [nf_cap] per face slot
  std::span<const signed char> curv_k;   // [nv0] curvature index
  std::span<Zeta30> dev;                 // [k_max+1] Q-frame development
  std::span<Real30> diag_pend;           // [k_max] accepted-ear diagonal lsq
};

struct CyclotomicMetric : CarryViews {
  // Flip transport, armed by flipped() and applied by the set_edge_length
  // the SAME flip issues (flip_edge calls them back to back on one h; a
  // plan-refused flip leaves a stale entry that the next flipped() simply
  // re-arms).
  struct PendingFlip {
    int h = -1;
    Real30 f2, w_origin, w_far;
  } pend{};
  // Ear-diagonal FIFO: ear() pushes each ACCEPTED diagonal's exact lsq in
  // acceptance order; splice_fan's set_edge_length calls consume them in
  // the same order (ear_clip_fan records diagonals 1:1 with acceptances).
  int n_diag = 0, i_diag = 0;
  Real30 dev_scale{};                    // Ls(0) of the current development
  int dev_k = 0;

  // The float shadows of an exact squared length: the x kLsqScale
  // convention cleared, and its root.
  static double shadow_lsq(const Real30& q) { return q.value() / (double)kLsqScale; }
  static double shadow_len(const Real30& q) { return std::sqrt(shadow_lsq(q)); }

  // Does h's float shadow agree with its exact lsq, at the integrality
  // band?  The one spelling of the exact<->float length bridge: the entry
  // boundary refuses on it, the audits re-check it.
  bool shadow_agrees(const DelaunayView& V, int h) const {
    const double sq = V.he_length[h] * V.he_length[h];
    return std::abs(sq - shadow_lsq(lsq[h])) <=
           delaunay_detail::lsq_integrality_band * std::max(1.0, sq);
  }
  // Does v's float curvature agree with its curvature index, kappa =
  // k * pi/15?  The one spelling of the exact<->float curvature bridge.
  bool cone_agrees(const DelaunayView& V, int v) const {
    return std::abs(V.curvature(v) - curv_k[v] * kCurvatureQuantum) <=
           delaunay_detail::curvature_agreement_band;
  }

  // A sign the oracle refused is a corrupt carry: trip by name and hand
  // the refusal on, so each predicate below is one composition over a
  // decided sign (the benign value falls out of the nullopt; the latch is
  // terminal, so nothing downstream reads it).
  static SignOr decided(DelaunayView& V, SignOr s, const char* what, int h) {
    if (!s) V.trip(DelaunayView::Status::InvariantViolated, what, h);
    return s;
  }

  // The wedge-carrying diamond of h; a carry the validated constructor
  // refuses is corrupt by definition -- trip, never guess.
  std::optional<Diamond> diamond_of(DelaunayView& V, int h) const {
    const auto A = V.diamond_arcs(h);
    Refusal why = Refusal::None;
    auto D = Diamond::make(lsq[A.e], lsq[A.a], lsq[A.b], lsq[A.c], lsq[A.d],
                           f_wedge[V.he_face[h]],
                           f_wedge[V.he_face[V.twin(h)]], &why);
    if (!D)
      V.trip(DelaunayView::Status::InvariantViolated,
             "cyclotomic: diamond carry refused (corrupt lsq/wedge)", h);
    return D;
  }

  bool is_flat(const DelaunayView&, int v) const { return curv_k[v] == 0; }

  // The sign of h's Delaunay form, decided -- the one form behind delaunay
  // and cocircular.
  SignOr delaunay_sign(DelaunayView& V, int h) const {
    const auto D = diamond_of(V, h);
    return D ? decided(V, D->delaunay_form_sign(),
                       "cyclotomic delaunay form: sign refused", h)
             : SignOr{};
  }
  bool delaunay(DelaunayView& V, int h) const {
    const SignOr s = delaunay_sign(V, h);
    return s && *s != Sign::Negative;
  }
  bool cocircular(DelaunayView& V, int h) const {
    const SignOr s = delaunay_sign(V, h);
    return s && *s == Sign::Zero;
  }
  bool convex(DelaunayView& V, int h) const {
    const auto D = diamond_of(V, h);
    if (!D) return false;
    const SignOr u = decided(V, D->convex_at_origin_sign(),
                             "cyclotomic convex: sign refused", h);
    const SignOr w = decided(V, D->reversed().convex_at_origin_sign(),
                             "cyclotomic convex: sign refused", h);
    return u && w && *u == Sign::Positive && *w == Sign::Positive;
  }

  // The exact order of two squared lengths, three-way -- the exact_metric
  // word (the canonical completion's corner keys are ordered by it): the
  // ring's order, cyclotomic::compare, decided.  A refusal (poisoned or
  // overflowed element) is a corrupt carry -- trip, never guess; 0 is
  // returned so the caller's status check ends the walk.
  // @post on Ok: result == sign(lsq[a] - lsq[b])
  int compare_lsq(DelaunayView& V, int a, int b) const {
    const SignOr s = decided(V, compare(lsq[a], lsq[b]),
                             "cyclotomic compare_lsq: sign refused", a);
    return s ? (int)*s : 0;
  }

  std::optional<Length> flipped(DelaunayView& V, int h) {
    const auto D = diamond_of(V, h);
    if (!D) return std::nullopt;
    DivTrace dt;
    const auto g = D->flipped(&dt);
    if (!g) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "cyclotomic flip: exact division refused (module premise)", h);
      return std::nullopt;
    }
    pend = {h, g->f2, g->w_origin, g->w_far};
    return Length{shadow_len(g->f2), 0};
  }

  void set_edge_length(DelaunayView& V, int h, Length l) {
    const int t = V.twin(h);
    V.he_length[h] = V.he_length[t] = l.len;
    if (pend.h == h) {
      // flip_edge writes before rewiring faces: he_face[h] is still the
      // old upper slot, which the flip turns into the far-side face
      // (B,D,v); the twin's slot becomes the origin-side face (D,B,u).
      lsq[h] = lsq[t] = pend.f2;
      f_wedge[V.he_face[h]] = pend.w_far;
      f_wedge[V.he_face[t]] = pend.w_origin;
      pend.h = -1;
      return;
    }
    if (i_diag < n_diag) {   // splice_fan's diagonal writes, in ear order
      lsq[h] = lsq[t] = diag_pend[i_diag++];
      return;
    }
    V.trip(DelaunayView::Status::InvariantViolated,
           "cyclotomic set_edge_length: unpaired length write", h);
  }

  // The Q-frame development of the flat star (file banner).
  void prepare_star(DelaunayView& V, DelaunayWorkspace& ws, int v) {
    pend.h = -1;
    n_diag = i_diag = 0;
    const FanPolygon& fan = ws.fan;
    const int k = fan.k;
    if (k + 1 > (int)dev.size()) {
      V.trip(DelaunayView::Status::CapacityExceeded,
             "cyclotomic prepare_star: development scratch", k);
      return;
    }
    dev_k = k;
    dev_scale = lsq[fan.spoke_he[0]];
    dev[0] = Zeta30{dev_scale, Real30{}};
    for (int i = 0; i < k; i++) {
      const Real30& Ls_i = lsq[fan.spoke_he[i]];
      const Real30& Ls_n = lsq[fan.spoke_he[(i + 1) % k]];
      const Real30& Lr_i = lsq[fan.inner_rim[i]];
      const Real30& w_i = f_wedge[V.he_face[fan.spoke_he[i]]];
      const Real30 d_i = Ls_i + Ls_n - Lr_i;
      const Zeta30 rot{d_i - Real30::gamma() * w_i, 2 * w_i};
      const Zeta30 num = dev[i] * rot;
      const Real30 den = 2 * Ls_i;
      const auto qx = exact_div(num.x, den);
      const auto qy = exact_div(num.y, den);
      if (!qx || !qy) {
        V.trip(DelaunayView::Status::InvariantViolated,
               "cyclotomic prepare_star: development division refused", v);
        return;
      }
      const Zeta30 q{*qx, *qy};
      if (i + 1 < k) {
        dev[i + 1] = q;
      } else if (!(q == dev[0])) {
        V.trip(DelaunayView::Status::InvariantViolated,
               "cyclotomic prepare_star: flat star fails to close", v);
        return;
      }
    }
  }

  // Exact ear acceptance over the Q-frame development: strict CCW ear +
  // apex sector at most pi (the two conditions of
  // delaunay_detail::ear_diag_sq_if_acceptable); on acceptance the exact
  // diagonal lsq is queued for splice_fan's paired write.
  Length ear(DelaunayView& V, const FanPolygon&, int pp, int pi, int pn) {
    bool ok = true;
    const SignOr s = sign_real(wedge(dev[pi] - dev[pp], dev[pn] - dev[pp]));
    if (!s) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "cyclotomic ear: CCW sign refused", pi);
      return {0, 0};
    }
    if (*s != Sign::Positive) return {0, 0};
    if (!sector_at_most_pi(dev[pp], dev[pn], ok)) {
      if (!ok)
        V.trip(DelaunayView::Status::InvariantViolated,
               "cyclotomic ear: sector sign refused", pi);
      return {0, 0};
    }
    const auto dsq = exact_div((dev[pn] - dev[pp]).lsq(), dev_scale);
    if (!dsq) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "cyclotomic ear: diagonal descale refused", pi);
      return {0, 0};
    }
    if (n_diag >= (int)diag_pend.size()) {
      V.trip(DelaunayView::Status::CapacityExceeded,
             "cyclotomic ear: diagonal queue", n_diag);
      return {0, 0};
    }
    diag_pend[n_diag++] = *dsq;
    return Length{shadow_len(*dsq), 0};
  }

  // The removal commit: the ear faces exist now (splice_fan pushed them in
  // triangle order), so their wedges land -- from the development, one
  // exact descale each, strictly positive by the ear acceptance.
  void commit_star(DelaunayView& V, DelaunayWorkspace& ws, int v) {
    if (i_diag != n_diag) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "cyclotomic commit_star: unconsumed ear diagonals",
             n_diag - i_diag);
      return;
    }
    const auto faces = ws.new_faces.live();
    if ((int)faces.size() != ws.tri.n_triangles) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "cyclotomic commit_star: face/triangle count mismatch",
             (int)faces.size());
      return;
    }
    for (int ti = 0; ti < ws.tri.n_triangles; ti++) {
      const auto& t = ws.tri.triangles[ti];
      const auto w = exact_div(
          wedge(dev[t.v1] - dev[t.v0], dev[t.v2] - dev[t.v0]), dev_scale);
      if (!w) {
        V.trip(DelaunayView::Status::InvariantViolated,
               "cyclotomic commit_star: wedge descale refused", v);
        return;
      }
      const SignOr s = sign_real(*w);
      if (!s || *s != Sign::Positive) {
        V.trip(DelaunayView::Status::InvariantViolated,
               "cyclotomic commit_star: non-positive ear-face wedge", v);
        return;
      }
      f_wedge[faces[ti]] = *w;
    }
  }

  // The tie-break's side order (the exact mirror of
  // lattice_first_tie_side, same walk, same verdict): develop side 0's
  // sector in the Q-frame and test whether it subtends at most pi.
  int first_tie_side(DelaunayView& V, DelaunayWorkspace& ws, int h_loop,
                     const double* /*theta*/) {
    int n = 0;
    for (int g = h_loop;; g = V.ccw(g)) {
      if (n >= (int)ws.poly.size()) {
        V.trip(DelaunayView::Status::CapacityExceeded,
               "cyclotomic first_tie_side: sector-arc scratch", n);
        return 0;
      }
      ws.poly[n++] = g;
      if (g == V.twin(h_loop)) break;
    }
    const Zeta30 q0{lsq[ws.poly[0]], Real30{}};
    Zeta30 qt = q0;
    for (int t = 0; t + 1 < n; t++) {
      const Real30& Ls_t = lsq[ws.poly[t]];
      const Real30& Ls_n = lsq[ws.poly[t + 1]];
      const Real30& Lr_t = lsq[V.he_next[ws.poly[t]]];
      const Real30& w_t = f_wedge[V.he_face[ws.poly[t]]];
      const Real30 d_t = Ls_t + Ls_n - Lr_t;
      const Zeta30 rot{d_t - Real30::gamma() * w_t, 2 * w_t};
      const Zeta30 num = qt * rot;
      const Real30 den = 2 * Ls_t;
      const auto qx = exact_div(num.x, den);
      const auto qy = exact_div(num.y, den);
      if (!qx || !qy) {
        V.trip(DelaunayView::Status::InvariantViolated,
               "cyclotomic first_tie_side: development division refused",
               h_loop);
        return 0;
      }
      qt = Zeta30{*qx, *qy};
    }
    bool ok = true;
    const bool le = sector_at_most_pi(q0, qt, ok);
    if (!ok) {
      V.trip(DelaunayView::Status::InvariantViolated,
             "cyclotomic first_tie_side: sector sign refused", h_loop);
      return 0;
    }
    return le ? 0 : 1;
  }
};

// ---------------------------------------------------------------------------
// The owned carry + the verified entry boundary.
// ---------------------------------------------------------------------------
struct CyclotomicKisCarry {
  std::vector<Real30> lsq, f_wedge, diag_pend;
  std::vector<signed char> curv_k;
  std::vector<Zeta30> dev;

  // A fresh policy over this carry: the views, and the transport state at
  // its defaults (pend disarmed, the ear FIFO empty).
  CyclotomicMetric metric() {
    return {CarryViews{lsq, f_wedge, curv_k, dev, diag_pend}};
  }
};

// Derive + VERIFY the cyclotomic carry from a FRESH kis DCEL (see the
// file banner's entry boundary).  n_centres = the face-centre vertex
// count (kis ids < n_centres); centre_size[c] = that face's size, 5 or 6.
// Every mismatch throws: this boundary is loud by design.
inline CyclotomicKisCarry derive_cyclotomic_kis_carry(
    const DelaunayView& V, int n_centres, std::span<const int> centre_size,
    const char* op = "derive_cyclotomic_kis_carry") {
  auto fail = [&](const std::string& what, long id) {
    throw std::runtime_error(std::string(op) + ": " + what + " (id " +
                             std::to_string(id) + ")");
  };
  if (n_centres <= 0 || n_centres >= V.nv ||
      (int)centre_size.size() < n_centres)
    fail("centre bookkeeping does not match the DCEL", n_centres);

  CyclotomicKisCarry c;
  c.lsq.assign(V.he_length.size(), Real30{});
  c.f_wedge.assign(V.f_he.size(), Real30{});
  c.curv_k.assign((std::size_t)V.nv, 0);
  c.dev.assign(V.he_length.size() + 1, Zeta30{});
  c.diag_pend.assign(V.he_length.size(), Real30{});
  // The verifying reads go through the policy's own bridge words, on views
  // of the arrays just sized (nothing below resizes them).
  const CyclotomicMetric m = c.metric();

  const Real30 U = Real30::lsq_cubic_edge();
  const Real30 S = Real30::lsq_pentagon_spoke();
  const Real30 WH = Real30::wedge_hexagon_kis();
  const Real30 WP = Real30::wedge_pentagon_kis();

  // Edge classes, each verified against its float shadow.
  for (int h = 0; h < V.nh; h++) {
    if (!V.alive(h)) continue;
    const int u = V.he_origin[h], w = V.dest(h);
    const bool cu = u < n_centres, cw = w < n_centres;
    if (cu && cw) fail("centre-centre edge (not a kis complex)", h);
    Real30 L = U;
    if (cu || cw) {
      const int centre = cu ? u : w;
      const int sz = centre_size[centre];
      if (sz != 5 && sz != 6) fail("centre face size not 5 or 6", centre);
      L = (sz == 5) ? S : U;
    }
    c.lsq[h] = L;
    if (!m.shadow_agrees(V, h))
      fail("he_length disagrees with the kis edge class", h);
  }

  // Face classes: exactly one centre corner each; the wedge by its size.
  for (int f = 0; f < V.nf; f++) {
    if (V.f_he[f] < 0) continue;
    const auto hs = V.face_halfedges(f);
    int centre = -1, n_c = 0;
    for (int j = 0; j < 3; j++) {
      const int vtx = V.he_origin[hs[j]];
      if (vtx < n_centres) { centre = vtx; n_c++; }
    }
    if (n_c != 1) fail("kis face without exactly one centre corner", f);
    c.f_wedge[f] = (centre_size[centre] == 5) ? WP : WH;
  }

  // Curvature indices: centres flat; cubic corners count their deg-5
  // centres.  Gauss-Bonnet (sum k = 60 on a fullerene kis) and the float
  // cone angles close the loop.
  long k_total = 0;
  for (int v = 0; v < V.nv; v++) {
    if (V.v_out[v] < 0) continue;
    int k = 0;
    if (v >= n_centres) {
      for (int h : V.incident(v)) {
        const int d = V.dest(h);
        if (d < n_centres && centre_size[d] == 5) k++;
      }
      if (k > 3) fail("cubic corner with more than 3 pentagons", v);
    }
    c.curv_k[v] = (signed char)k;
    k_total += k;
    if (!m.cone_agrees(V, v))
      fail("cone angle disagrees with the curvature index", v);
  }
  if (k_total != 60) fail("total curvature is not 4 pi (sum k != 60)", k_total);
  return c;
}

}  // namespace cyclotomic
