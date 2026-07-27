#pragma once

// Dense linear algebra on small matrices and vectors — deliberately
// BLAS/LAPACK-free.
//
// At least one deployed OpenBLAS (0.3.20 with DYNAMIC_ARCH selecting
// Cooperlake AVX-512 kernels) silently returns WRONG dgesv solutions
// (info = 0, residual ~1e-2) for n ≳ 60, and wrong dgemm products at
// n ≥ 128.  The systems solved through this module are tiny (n = tens:
// cone counts, small Jacobians), so an in-house cyclic Jacobi and a
// partial-pivot LU are fast enough — and their correctness does not
// depend on which BLAS the host has installed.  Do NOT "optimize" these
// back to LAPACK calls.
//
// Vector arithmetic (v+w, v−w, v*s) is NOT defined here: the global
// template operators in auxiliary.hh already provide it.

#include "fullerenes/dense_linalg_view.hh"   // the allocation-free bodies
#include "fullerenes/matrix.hh"

#include <expected>
#include <vector>

// This header is the OWNER API: allocating conveniences over
// matrix<double>/vector.  The algorithm bodies live at view level
// (dense_linalg_view.hh: spans + row stride, allocation-free, device-legal,
// written once); each function here is a thin wrapper around its view body,
// and the failure-channel table in the view header maps every entry point's
// singular behaviour.  The vector reductions (dot, norm, sum_sq, max_abs,
// energy, is_usable_step, negate) have NO owner overloads: std::vector
// converts implicitly to std::span, so the view functions serve owner
// callers directly.  The symmetric-eigen family (jacobi_eig, SymEigen) is
// owner-only.

namespace LinAlg {

using V = std::vector<double>;

// The live block of an owner matrix as a view (matrix<T> stores row-major
// flat with stride n; empty matrices yield an empty view).
inline MatConstView view_of(const matrix<double>& A) {
  const std::size_t sz = (std::size_t)A.m * A.n;
  return {std::span<const double>(sz ? A.data() : nullptr, sz), A.m, A.n, A.n};
}

// (Vector reductions: see the view header -- no owner overloads, vectors
// convert to spans.)

// --- Symmetric eigendecomposition (cyclic Jacobi) ---

// Cyclic Jacobi eigendecomposition of a symmetric matrix (flat row-major
// A, consumed).  Eigenvalues ascending in lam; if V_out is non-null, row m
// of *V_out is the unit eigenvector paired with lam[m].  The sweep budget
// is a guard, not a knob: cyclic Jacobi converges quadratically once the
// off-diagonal is small, so a trip means non-symmetric input (caller bug).
// The bool IS the two-outcome encoding (true = converged, false = guard
// trip); callers symmetrize their input, so false cannot occur through
// sym_eigvals/SymEigen::decompose — they translate it to their empty-
// result outcome rather than propagate an exception.
bool jacobi_eig(std::vector<double> A, int n, std::vector<double>& lam,
                std::vector<double>* V_out);

// Eigenvalues of the symmetrized A, sorted ascending; empty on guard trip.
std::vector<double> sym_eigvals(const matrix<double>& A);

// --- LU with partial pivoting ---

// LU-with-sign solve: (x, det_sign) on success.
// det_sign ∈ {−1, +1} is sign(det A), the product of
//   (−1)^{# row swaps}  ×  product of signs of diag(U).
// Singular = an exact zero pivot.
struct LuSolved { V x; int det_sign; };
enum class LuFail { Singular };

std::expected<LuSolved, LuFail> solve_with_sign(const matrix<double>& A,
                                                const V& b);

// Signed determinant of A via the same partial-pivot LU as solve_with_sign:
// magnitude ∏|diag(U)|, sign (−1)^{# row swaps} × ∏ sign(diag(U)).  A singular
// matrix has det = 0 — a legitimate VALUE, not a failure — so an exact zero
// pivot returns 0.0 rather than reporting LuFail::Singular the way
// solve_with_sign does.
double det(const matrix<double>& A);

// Solve A·x = b.  Zero vector on failure (singular A).
// @pre b.size() == A.m.  The result has A.m entries; historically the
// pre-view implementation returned b.size() entries when b was longer --
// no caller relied on it, and the size contract is now explicit.
V solve(const matrix<double>& A, const V& b);

// Solve (A + λI)·x = b.  @pre b.size() == A.m.
V solve_shifted(const matrix<double>& A, const V& b, double lambda);

// Matrix-vector product A·v.
V matvec(const matrix<double>& A, const V& v);

// --- Symmetric truncated pseudoinverse ---
//
// For symmetric A: the spectral decomposition A = Q·diag(λ)·Qᵀ yields the
// pseudoinverse A⁺ = Q·diag(λ⁺)·Qᵀ with (λ⁺)ᵢ = 1/λᵢ when |λᵢ| > cutoff,
// else 0 — the unique minimum-norm solve that is automatically orthogonal
// to the (soft) kernel.
namespace SymEigen {

// Q has eigenvectors as columns; λ ascending.  Empty on eigensolver failure.
struct Decomp {
  matrix<double> Q;
  V              lambda;
};
Decomp decompose(const matrix<double>& A);

// Signature of a symmetric operator at a numerical "zero" tolerance:
// (#λ > tol, #λ < −tol, #|λ| ≤ tol).
struct Signature { int pos, neg, zero; };
Signature signature(const V& lambda, double tol);

// x = Σ_{|λᵢ|>cutoff} (qᵢᵀ b / λᵢ) qᵢ.  rank = #non-truncated modes;
// lambda_min_kept = smallest |λ| among them (0 if all truncated).
struct ApplyResult { V x; int rank; double lambda_min_kept; };
ApplyResult apply(const Decomp& d, const V& b, double cutoff);

// Standalone pseudoinverse solve x = A⁺·b with cutoff = rcond·max|λ|.
struct Solution {
  V         x;
  int       rank           = 0;
  double    lambda_abs_min = 0;
  double    lambda_abs_max = 0;
  Signature sig            = {0, 0, 0};
};
Solution solve(const matrix<double>& A, const V& b, double rcond);

}  // namespace SymEigen

}  // namespace LinAlg
