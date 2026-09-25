#pragma once
// ============================================================================
// delaunay_cyclotomic.hh -- CyclotomicMetricT: the exact metric policy of
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
// Every predicate consults the constructor-validated cyclotomic::DiamondT,
// so a corrupt carry REFUSES (InvariantViolated through the view's latch)
// instead of mis-classifying; flips transport (lsq, both face wedges)
// through DiamondT::flipped's exact divisions.
//
// THE TWO WIDTHS.  The policy is templated on the STORAGE type S of a
// carried number and the RING R of the arithmetic (cyclotomic.hh's
// WIDTHS paragraph); every read widens S to R, every write narrows R to
// S, and a value the storage cannot hold trips CapacityExceeded by name.
// Two instantiations are used:
//   CyclotomicMetric32  = <Stored30<int32_t>, Real30>      the batch tier
//                         (device and host sweeps; through ~C15,000);
//   CyclotomicMetric64  = <Real30, Real30Wide>              single huge
//                         isomers (delaunay_cyclotomic_wide.hh, host).
// The unqualified names (CyclotomicMetric, CarryStore, CyclotomicKisCarry)
// are the 64-bit ring over its own storage, kept for the parallel port's
// arena until it moves to the batch tier; they are not a supported tier.
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
// ONE tier of code: the policy is the same code on the host and in a GPU
// work-item (cyclotomic.hh's banner) -- the cubic chain's device
// reduction runs it over a per-isomer carry arena (claude-projects/
// parallel-primitives, IdtBatch::launch_kis_reduce_exact).
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

// The storage of a development point at storage type S: a pair of
// coordinates at that width.  For a ring used as its own storage the pair
// IS the point ring (the parallel port's arena holds Zeta30 values).
template <class S>
struct StoredZeta {
  S x, y;
};
template <class S>
struct zeta_store {
  using type = StoredZeta<S>;
};
template <class Coef>
struct zeta_store<Real30T<Coef>> {
  using type = Zeta30T<Real30T<Coef>>;
};
template <class S>
using zeta_store_t = typename zeta_store<S>::type;

template <class R, class S>
inline Zeta30T<R> widen_point(const StoredZeta<S>& z) {
  return {widen<R>(z.x), widen<R>(z.y)};
}
template <class R, class Coef>
inline Zeta30T<R> widen_point(const Zeta30T<Real30T<Coef>>& z) {
  return {widen<R>(z.x), widen<R>(z.y)};
}
template <class Z, class R>
inline std::optional<Z> narrow_point(const Zeta30T<R>& z) {
  const auto x = narrow<decltype(Z::x)>(z.x), y = narrow<decltype(Z::y)>(z.y);
  if (!x || !y) return std::nullopt;
  return Z{*x, *y};
}

// Does the CCW sector from direction `from` to direction `to` subtend at
// most pi?  The exact mirror of delaunay_detail::lattice_sector_at_most_pi
// (its banner carries the argument), over Zeta30 directions: wedge > 0 is
// (0, pi); wedge == 0 with negative real inner product is exactly pi (the
// accepted tie); wedge == 0 with positive dot would be a 0 / 2pi sector,
// excluded by the fan premise (defensive reject).  A refused sign leaves
// its name in `tr` (caller trips).
template <class R>
inline bool sector_at_most_pi(const Zeta30T<R>& from, const Zeta30T<R>& to,
                              SignTrace& tr) {
  const SignOr sw = sign_real(wedge(from, to), &tr);
  if (!sw) return false;
  if (*sw == Sign::Positive) return true;
  if (*sw == Sign::Negative) return false;
  // wedge == 0: conj(from)*to is real, its x-part IS the inner product.
  const SignOr sd = sign_real((from.conj() * to).x, &tr);
  if (!sd) return false;
  return *sd == Sign::Negative;
}

// ---------------------------------------------------------------------------
// The policy.  Spans are caller-owned (CyclotomicKisCarryT below); the
// pending-transport members make the policy STATEFUL -- one instance must
// live across a whole reduction run (the drivers pass the metric by
// reference throughout).
// ---------------------------------------------------------------------------
// The carry as the policy sees it: the three exact arrays (file banner) and
// the two development scratch arrays, as caller-owned views at storage
// type S.
template <class S>
struct CarryViewsT {
  std::span<S> lsq;                      // [nh_cap] per half-edge, twin-paired
  std::span<S> f_wedge;                  // [nf_cap] per face slot
  std::span<const signed char> curv_k;   // [nv0] curvature index
  std::span<zeta_store_t<S>> dev;        // [k_max+1] Q-frame development
  std::span<S> diag_pend;                // [k_max] accepted-ear diagonal lsq
};

template <class S, class R>
struct CyclotomicMetricT : CarryViewsT<S> {
  using Views = CarryViewsT<S>;
  using Views::lsq;
  using Views::f_wedge;
  using Views::curv_k;
  using Views::dev;
  using Views::diag_pend;
  using ring = R;
  using storage = S;
  using Zs = zeta_store_t<S>;
  using Zr = Zeta30T<R>;

#ifdef FULLERENES_SIGN_FILTER_CHECK
  // The floating-point filter differential's counters for THIS isomer, or
  // null.  Bound by the caller into every trace the predicates below make,
  // so the check observes the same predicates production runs.  Instrument
  // only; absent from a production or device build.
  SignFilterCounters* sign_ctr = nullptr;
#endif

  // Flip transport, armed by flipped() and applied by the set_edge_length
  // the SAME flip issues (flip_edge calls them back to back on one h; a
  // plan-refused flip leaves a stale entry that the next flipped() simply
  // re-arms).
  struct PendingFlip {
    int h = -1;
    R f2, w_origin, w_far;
  } pend{};
  // Ear-diagonal FIFO: ear() pushes each ACCEPTED diagonal's exact lsq in
  // acceptance order; splice_fan's set_edge_length calls consume them in
  // the same order (ear_clip_fan records diagonals 1:1 with acceptances).
  int n_diag = 0, i_diag = 0;
  R dev_scale{};                         // Ls(0) of the current development
  int dev_k = 0;

  // The carry read at arithmetic width: squared length of h, wedge of face
  // f, development point i.
  R L(int h) const { return widen<R>(lsq[h]); }
  R W(int f) const { return widen<R>(f_wedge[f]); }
  Zr Qd(int i) const { return widen_point<R>(dev[i]); }

  // The carry written at storage width: a value the width cannot hold
  // trips by name (the storage-width refusal) and is not written.
  static bool put(DelaunayView& V, std::span<S> arr, int i, const R& v,
                  const char* what) {
    const auto n = narrow<S>(v);
    if (!n) {
      V.trip(DelaunayView::Status::CapacityExceeded, what, i);
      return false;
    }
    arr[i] = *n;
    return true;
  }
  static bool put_point(DelaunayView& V, std::span<Zs> arr, int i, const Zr& z,
                        const char* what) {
    const auto n = narrow_point<Zs>(z);
    if (!n) {
      V.trip(DelaunayView::Status::CapacityExceeded, what, i);
      return false;
    }
    arr[i] = *n;
    return true;
  }

  // The float shadows of an exact squared length: the x kLsqScale
  // convention cleared, and its root.
  static double shadow_lsq(const R& q) { return q.value() / (double)kLsqScale; }
  static double shadow_len(const R& q) { return std::sqrt(shadow_lsq(q)); }

  // Does h's float shadow agree with its exact lsq, at the integrality
  // band?  The one spelling of the exact<->float length bridge: the entry
  // boundary refuses on it, the audits re-check it.
  bool shadow_agrees(const DelaunayView& V, int h) const {
    const double sq = V.he_length[h] * V.he_length[h];
    return std::abs(sq - shadow_lsq(L(h))) <=
           delaunay_detail::lsq_integrality_band * std::max(1.0, sq);
  }
  // Does v's float curvature agree with its curvature index, kappa =
  // k * pi/15?  The one spelling of the exact<->float curvature bridge.
  bool cone_agrees(const DelaunayView& V, int v) const {
    return std::abs(V.curvature(v) - curv_k[v] * kCurvatureQuantum) <=
           delaunay_detail::curvature_agreement_band;
  }

  // A refusal of the ring, tripped by name.  A POISONED value is a
  // coefficient overflow: the carry has left the ring's guaranteed
  // envelope (carry_coeff_max), a CAPACITY refusal, the same kind as a
  // storage width exceeded.  Every other refusal (a non-divisible
  // quotient, a non-positive wedge, an inconsistent carry, an undecided
  // sign) falsifies the carry: InvariantViolated.
  static void refuse(DelaunayView& V, Refusal why, const char* what, int id) {
    V.trip(why == Refusal::Poisoned ? DelaunayView::Status::CapacityExceeded
                                    : DelaunayView::Status::InvariantViolated,
           what, id);
  }
  // A sign, decided or refused by name, so each predicate below is one
  // composition over a decided sign (the benign value falls out of the
  // nullopt; the latch is terminal, so nothing downstream reads it).
  static SignOr decided(DelaunayView& V, SignOr s, const SignTrace& tr,
                        const char* what, int h) {
    if (!s) refuse(V, tr.refusal, what, h);
    return s;
  }

  // The wedge-carrying diamond of h; a carry the validated constructor
  // refuses is corrupt (or, poisoned, outside the envelope) -- trip, never
  // guess.
  std::optional<DiamondT<R>> diamond_of(DelaunayView& V, int h) const {
    const auto A = V.diamond_arcs(h);
    Refusal why = Refusal::None;
    auto D = DiamondT<R>::make(L(A.e), L(A.a), L(A.b), L(A.c), L(A.d),
                               W(V.he_face[h]), W(V.he_face[V.twin(h)]), &why);
    if (!D) refuse(V, why, "cyclotomic: diamond carry refused (lsq/wedge)", h);
    return D;
  }

  bool is_flat(const DelaunayView&, int v) const { return curv_k[v] == 0; }

  // The sign of h's Delaunay form, decided -- the one form behind delaunay
  // and cocircular.
  SignOr delaunay_sign(DelaunayView& V, int h) const {
    const auto D = diamond_of(V, h);
    if (!D) return SignOr{};
    SignTrace tr;
    FULLERENES_BIND_SIGN_CTR(tr, sign_ctr);
    return decided(V, D->delaunay_form_sign(&tr), tr,
                   "cyclotomic delaunay form: sign refused", h);
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
    SignTrace tu, tw;
    FULLERENES_BIND_SIGN_CTR(tu, sign_ctr);
    FULLERENES_BIND_SIGN_CTR(tw, sign_ctr);
    const SignOr u = decided(V, D->convex_at_origin_sign(&tu), tu,
                             "cyclotomic convex: sign refused", h);
    const SignOr w = decided(V, D->reversed().convex_at_origin_sign(&tw), tw,
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
    SignTrace tr;
    FULLERENES_BIND_SIGN_CTR(tr, sign_ctr);
    const SignOr s = decided(V, compare(L(a), L(b), &tr), tr,
                             "cyclotomic compare_lsq: sign refused", a);
    return s ? (int)*s : 0;
  }

  std::optional<Length> flipped(DelaunayView& V, int h) {
    const auto D = diamond_of(V, h);
    if (!D) return std::nullopt;
    DivTrace dt;
    const auto g = D->flipped(&dt);
    if (!g) {
      refuse(V, dt.refusal, "cyclotomic flip: exact division refused", h);
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
      if (!put(V, lsq, h, pend.f2,
               "cyclotomic set_edge_length: squared length exceeds the storage width"))
        return;
      lsq[t] = lsq[h];
      if (!put(V, f_wedge, V.he_face[h], pend.w_far,
               "cyclotomic set_edge_length: wedge exceeds the storage width") ||
          !put(V, f_wedge, V.he_face[t], pend.w_origin,
               "cyclotomic set_edge_length: wedge exceeds the storage width"))
        return;
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
    dev_scale = L(fan.spoke_he[0]);
    if (!put_point(V, dev, 0, Zr{dev_scale, R{}},
                   "cyclotomic prepare_star: development point exceeds the storage width"))
      return;
    Zr q_i = Qd(0);
    for (int i = 0; i < k; i++) {
      const R Ls_i = L(fan.spoke_he[i]);
      const R Ls_n = L(fan.spoke_he[(i + 1) % k]);
      const R Lr_i = L(fan.inner_rim[i]);
      const R w_i = W(V.he_face[fan.spoke_he[i]]);
      const R d_i = Ls_i + Ls_n - Lr_i;
      const Zr rot{d_i - R::gamma() * w_i, 2 * w_i};
      const Zr num = q_i * rot;
      const R den = 2 * Ls_i;
      DivTrace dt;
      const auto qx = exact_div(num.x, den, &dt);
      const auto qy = qx ? exact_div(num.y, den, &dt) : std::optional<R>{};
      if (!qx || !qy) {
        refuse(V, dt.refusal, "cyclotomic prepare_star: development division refused", v);
        return;
      }
      const Zr q{*qx, *qy};
      if (i + 1 < k) {
        if (!put_point(V, dev, i + 1, q,
                       "cyclotomic prepare_star: development point exceeds the storage width"))
          return;
        q_i = q;
      } else if (!(q == Qd(0))) {
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
    const Zr qp = Qd(pp), qi = Qd(pi), qn = Qd(pn);
    SignTrace tr;
    FULLERENES_BIND_SIGN_CTR(tr, sign_ctr);
    const SignOr s = decided(V, sign_real(wedge(qi - qp, qn - qp), &tr), tr,
                             "cyclotomic ear: CCW sign refused", pi);
    if (!s) return {0, 0};
    if (*s != Sign::Positive) return {0, 0};
    SignTrace ts;
    FULLERENES_BIND_SIGN_CTR(ts, sign_ctr);
    if (!sector_at_most_pi(qp, qn, ts)) {
      if (ts.refusal != Refusal::None)
        refuse(V, ts.refusal, "cyclotomic ear: sector sign refused", pi);
      return {0, 0};
    }
    DivTrace dt;
    const auto dsq = exact_div((qn - qp).lsq(), dev_scale, &dt);
    if (!dsq) {
      refuse(V, dt.refusal, "cyclotomic ear: diagonal descale refused", pi);
      return {0, 0};
    }
    if (n_diag >= (int)diag_pend.size()) {
      V.trip(DelaunayView::Status::CapacityExceeded,
             "cyclotomic ear: diagonal queue", n_diag);
      return {0, 0};
    }
    if (!put(V, diag_pend, n_diag, *dsq,
             "cyclotomic ear: diagonal exceeds the storage width"))
      return {0, 0};
    n_diag++;
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
      const Zr q0 = Qd(t.v0);
      DivTrace dt;
      const auto w = exact_div(wedge(Qd(t.v1) - q0, Qd(t.v2) - q0), dev_scale, &dt);
      if (!w) {
        refuse(V, dt.refusal, "cyclotomic commit_star: wedge descale refused", v);
        return;
      }
      SignTrace tr;
      FULLERENES_BIND_SIGN_CTR(tr, sign_ctr);
      const SignOr s = decided(V, sign_real(*w, &tr), tr,
                               "cyclotomic commit_star: ear-face wedge sign refused", v);
      if (!s) return;
      if (*s != Sign::Positive) {
        V.trip(DelaunayView::Status::InvariantViolated,
               "cyclotomic commit_star: non-positive ear-face wedge", v);
        return;
      }
      if (!put(V, f_wedge, faces[ti], *w,
               "cyclotomic commit_star: wedge exceeds the storage width"))
        return;
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
    const Zr q0{L(ws.poly[0]), R{}};
    Zr qt = q0;
    for (int t = 0; t + 1 < n; t++) {
      const R Ls_t = L(ws.poly[t]);
      const R Ls_n = L(ws.poly[t + 1]);
      const R Lr_t = L(V.he_next[ws.poly[t]]);
      const R w_t = W(V.he_face[ws.poly[t]]);
      const R d_t = Ls_t + Ls_n - Lr_t;
      const Zr rot{d_t - R::gamma() * w_t, 2 * w_t};
      const Zr num = qt * rot;
      const R den = 2 * Ls_t;
      DivTrace dt;
      const auto qx = exact_div(num.x, den, &dt);
      const auto qy = qx ? exact_div(num.y, den, &dt) : std::optional<R>{};
      if (!qx || !qy) {
        refuse(V, dt.refusal, "cyclotomic first_tie_side: development division refused",
               h_loop);
        return 0;
      }
      qt = Zr{*qx, *qy};
    }
    SignTrace ts;
    FULLERENES_BIND_SIGN_CTR(ts, sign_ctr);
    const bool le = sector_at_most_pi(q0, qt, ts);
    if (ts.refusal != Refusal::None) {
      refuse(V, ts.refusal, "cyclotomic first_tie_side: sector sign refused", h_loop);
      return 0;
    }
    return le ? 0 : 1;
  }
};

// ---------------------------------------------------------------------------
// The carry's storage as caller-owned WRITABLE views (the derivation's
// output), the owned carry, and the verified entry boundary in two layers:
// a view body that fills a store and reports a refusal by name (no
// allocation, no exceptions -- callable from a GPU work-item over a
// per-isomer arena), and the owner that allocates, calls it, and converts
// a refusal to the documented throw.  Sizes: carry_capacities below.
// ---------------------------------------------------------------------------
template <class S, class R>
struct CarryStoreT {
  std::span<S> lsq, f_wedge, diag_pend;
  std::span<signed char> curv_k;
  std::span<zeta_store_t<S>> dev;

#ifdef FULLERENES_SIGN_FILTER_CHECK
  // Carried so a caller can arm the filter differential without changing
  // any signature between here and the predicates (instrument only).
  SignFilterCounters* sign_ctr = nullptr;
#endif

  CarryViewsT<S> views() const { return {lsq, f_wedge, curv_k, dev, diag_pend}; }
  // A fresh policy over this store: the views, and the transport state at
  // its defaults (pend disarmed, the ear FIFO empty).
  CyclotomicMetricT<S, R> metric() const {
    CyclotomicMetricT<S, R> m{views()};
    FULLERENES_BIND_SIGN_CTR_STORE(m, sign_ctr);
    return m;
  }
};

// A refused derivation: what failed and the id it failed on; ok() when the
// carry was derived and verified.
struct CarryRefusal {
  const char* what = nullptr;
  long id = 0;
  bool ok() const { return what == nullptr; }
};

// The carry's sizes over a DCEL with nv vertices, nh half-edges and nf
// faces: lsq and diag_pend one entry per half-edge, f_wedge one per face,
// curv_k one per vertex, dev one per half-edge plus one.
struct CarryCapacities {
  std::size_t lsq, f_wedge, diag_pend, curv_k, dev;
};
constexpr CarryCapacities carry_capacities(long nv, long nh, long nf) {
  return {std::size_t(nh), std::size_t(nf), std::size_t(nh), std::size_t(nv), std::size_t(nh + 1)};
}
// The carry's sizes at V's capacities: the storage a carry needs for the
// whole lifetime of V (the view never grows past its capacities).
inline CarryCapacities carry_capacities(const DelaunayView& V) {
  return carry_capacities(V.nv_cap, V.nh_cap, V.nf_cap);
}

// Derive + VERIFY the cyclotomic carry from a FRESH kis DCEL into a store
// (see the file banner's entry boundary).  n_centres = the face-centre
// vertex count (kis ids < n_centres); centre_size[c] = that face's size,
// 5 or 6 (any integer element type).  The first mismatch is returned by
// name; the store is then partially written and must not be used.
// @pre  every store span holds at least carry_capacities(V.nv, V.nh, V.nf);
//       V is a fresh kis DCEL
template <class Size, class S, class R>
inline CarryRefusal derive_cyclotomic_kis_carry_into(
    CarryStoreT<S, R> c, const DelaunayView& V, int n_centres,
    std::span<const Size> centre_size) {
  if (n_centres <= 0 || n_centres >= V.nv ||
      (long)centre_size.size() < n_centres)
    return {"centre bookkeeping does not match the DCEL", n_centres};
  const CarryCapacities need = carry_capacities(V.nv, V.nh, V.nf);
  if (c.lsq.size() < need.lsq || c.f_wedge.size() < need.f_wedge ||
      c.curv_k.size() < need.curv_k || c.dev.size() < need.dev ||
      c.diag_pend.size() < need.diag_pend)
    return {"carry store smaller than the DCEL", V.nh};
  for (auto& q : c.lsq) q = S{};
  for (auto& q : c.f_wedge) q = S{};
  for (auto& k : c.curv_k) k = 0;
  for (auto& z : c.dev) z = zeta_store_t<S>{};
  for (auto& q : c.diag_pend) q = S{};
  // The verifying reads go through the policy's own bridge words, on the
  // views of the store (nothing below resizes them).
  const CyclotomicMetricT<S, R> m = c.metric();

  // The constants table at storage width (|c| <= 100: every width holds it).
  const auto U = narrow<S>(R::lsq_cubic_edge());
  const auto Sp = narrow<S>(R::lsq_pentagon_spoke());
  const auto WH = narrow<S>(R::wedge_hexagon_kis());
  const auto WP = narrow<S>(R::wedge_pentagon_kis());
  if (!U || !Sp || !WH || !WP)
    return {"the constants table exceeds the storage width", 0};

  // Edge classes, each verified against its float shadow.
  for (int h = 0; h < V.nh; h++) {
    if (!V.alive(h)) continue;
    const int u = V.he_origin[h], w = V.dest(h);
    const bool cu = u < n_centres, cw = w < n_centres;
    if (cu && cw) return {"centre-centre edge (not a kis complex)", h};
    S L = *U;
    if (cu || cw) {
      const int centre = cu ? u : w;
      const int sz = (int)centre_size[centre];
      if (sz != 5 && sz != 6) return {"centre face size not 5 or 6", centre};
      L = (sz == 5) ? *Sp : *U;
    }
    c.lsq[h] = L;
    if (!m.shadow_agrees(V, h))
      return {"he_length disagrees with the kis edge class", h};
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
    if (n_c != 1) return {"kis face without exactly one centre corner", f};
    c.f_wedge[f] = ((int)centre_size[centre] == 5) ? *WP : *WH;
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
        if (d < n_centres && (int)centre_size[d] == 5) k++;
      }
      if (k > 3) return {"cubic corner with more than 3 pentagons", v};
    }
    c.curv_k[v] = (signed char)k;
    k_total += k;
    if (!m.cone_agrees(V, v))
      return {"cone angle disagrees with the curvature index", v};
  }
  if (k_total != 60) return {"total curvature is not 4 pi (sum k != 60)", k_total};
  return {};
}

template <class S, class R>
struct CyclotomicKisCarryT {
  std::vector<S> lsq, f_wedge, diag_pend;
  std::vector<signed char> curv_k;
  std::vector<zeta_store_t<S>> dev;

  CyclotomicKisCarryT() = default;
  // Storage at the given sizes, every entry at its zero value.
  explicit CyclotomicKisCarryT(const CarryCapacities& cap)
      : lsq(cap.lsq, S{}), f_wedge(cap.f_wedge, S{}), diag_pend(cap.diag_pend, S{}),
        curv_k(cap.curv_k, 0), dev(cap.dev, zeta_store_t<S>{}) {}

  CarryStoreT<S, R> store() { return {lsq, f_wedge, diag_pend, curv_k, dev}; }
  // A fresh policy over this carry: the views, and the transport state at
  // its defaults (pend disarmed, the ear FIFO empty).
  CyclotomicMetricT<S, R> metric() { return store().metric(); }
};

// The owning entry boundary: allocate the carry at the DCEL's capacities,
// derive and verify it (derive_cyclotomic_kis_carry_into), and convert a
// refusal to the documented throw.  Every mismatch throws: this boundary
// is loud by design.
template <class S = Real30, class R = Real30>
inline CyclotomicKisCarryT<S, R> derive_cyclotomic_kis_carry(
    const DelaunayView& V, int n_centres, std::span<const int> centre_size,
    const char* op = "derive_cyclotomic_kis_carry") {
  CyclotomicKisCarryT<S, R> c(carry_capacities(V));
  const CarryRefusal r =
      derive_cyclotomic_kis_carry_into<int, S, R>(c.store(), V, n_centres, centre_size);
  if (!r.ok())
    throw std::runtime_error(std::string(op) + ": " + r.what + " (id " +
                             std::to_string(r.id) + ")");
  return c;
}

// ---- The batch tier, 32/64 (file banner) ----
using Stored30Narrow = Stored30<int32_t>;
using CyclotomicMetric32 = CyclotomicMetricT<Stored30Narrow, Real30>;
using CarryStore32 = CarryStoreT<Stored30Narrow, Real30>;
using CyclotomicKisCarry32 = CyclotomicKisCarryT<Stored30Narrow, Real30>;

// ---- The 64-bit ring over its own storage: the parallel port's arena
// until it moves to the batch tier (not a supported tier, see banner) ----
using CarryViews = CarryViewsT<Real30>;
using CyclotomicMetric = CyclotomicMetricT<Real30, Real30>;
using CarryStore = CarryStoreT<Real30, Real30>;
using CyclotomicKisCarry = CyclotomicKisCarryT<Real30, Real30>;

}  // namespace cyclotomic
