#pragma once

// =====================================================================
// optim::LBFGS -- the limited-memory BFGS direction generator, a
// faithful transcription of minimize::lbfgs's direction computation
// (fullerenes/minimize.hh): ring buffer of (s, y) pairs, two-loop
// recursion with gamma = 1/(rho ytY) initial scaling, restart-on-
// non-descent, and the damped curvature skip
// s.y > 1e-10 sqrt(|s|^2 |y|^2).
//
// The arithmetic is kept operation-for-operation identical to the
// parent (same loop order, same expressions) because the wu
// re-expression gate requires bit-identical trajectories.  Deltahedron
// carries a slightly different L-BFGS arithmetic (gamma as ys/yy, deque
// storage); its migration gate is quality-distribution parity, not
// trajectory identity, and uses THIS policy.
// =====================================================================

#include "fullerenes/optim/core.hh"

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

namespace optim {

struct LBFGS {
  using paradigm = line_search_tag;

  int  m = 10;                    // history pairs
  bool unit_first_trial = false;  // true: first trial t0 = 1 always (Deltahedron);
                                  // false: minimize.hh's 1/gmax on an empty history (wu)

  struct State {
    std::vector<std::vector<double>> S, Y;   // ring buffers of (s, y)
    std::vector<double> rho, alpha;
    std::vector<double> q;                   // two-loop scratch
    int stored = 0, head = 0;

    State(std::size_t n, const LBFGS& p)
        : S(std::max(1, p.m), std::vector<double>(n)),
          Y(std::max(1, p.m), std::vector<double>(n)),
          rho(std::max(1, p.m)), alpha(std::max(1, p.m)), q(n) {}
  };

  void reset(State& st) const { st.stored = 0; }
  void invalidate(State& st) const { st.stored = 0; }   // curvature pairs straddle the move

  bool has_history(const State& st) const { return st.stored > 0; }

  // d = -H~ g by the two-loop recursion (minimize.hh lines: copy,
  // newest->oldest alpha loop, gamma scaling, oldest->newest beta loop).
  void direction(SeqCtx, State& st, std::span<const double> g,
                 std::span<double> d) const {
    const int M = (int)st.S.size();
    std::copy(g.begin(), g.end(), st.q.begin());
    for (int i = 0; i < st.stored; ++i) {                 // newest -> oldest
      const int j = ((st.head - 1 - i) % M + M) % M;
      st.alpha[j] = st.rho[j] * la::dot(st.S[j], st.q);
      la::axpy(st.q, -st.alpha[j], st.Y[j]);
    }
    if (st.stored > 0) {
      const int jn = ((st.head - 1) % M + M) % M;
      const double gamma = 1.0 / (st.rho[jn] * la::dot(st.Y[jn], st.Y[jn]));
      for (double& v : st.q) v *= gamma;
    }
    for (int i = st.stored - 1; i >= 0; --i) {            // oldest -> newest
      const int j = ((st.head - 1 - i) % M + M) % M;
      const double beta = st.rho[j] * la::dot(st.Y[j], st.q);
      la::axpy(st.q, st.alpha[j] - beta, st.S[j]);
    }
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = -st.q[i];
  }

  // Store the (s, y) pair when it carries curvature information.
  void accepted(SeqCtx, State& st, std::span<const double> s,
                std::span<const double> y, std::span<const double>) const {
    const int M = (int)st.S.size();
    auto& s_slot = st.S[st.head];
    auto& y_slot = st.Y[st.head];
    double sn2 = 0, yn2 = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
      s_slot[i] = s[i];
      y_slot[i] = y[i];
      sn2 += s[i] * s[i];
      yn2 += y[i] * y[i];
    }
    const double sy = la::dot(s_slot, y_slot);
    if (sy > 1e-10 * std::sqrt(sn2 * yn2)) {
      st.rho[st.head] = 1.0 / sy;
      st.head = (st.head + 1) % M;
      st.stored = std::min(st.stored + 1, M);
    }
  }

  // First-iteration trial step scaled by the gradient (minimize.hh t0).
  double first_trial(const State& st, double gmax) const {
    if (unit_first_trial) return 1.0;
    return (st.stored == 0) ? std::min(1.0, 1.0 / std::max(gmax, 1e-12))
                            : 1.0;
  }
};

}  // namespace optim
