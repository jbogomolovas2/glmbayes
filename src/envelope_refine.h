// envelope_refine.h -- optional refinement of envelope tangency points.
//
// Nygren & Nygren (2006) state that the optimal tangency for a restricted
// likelihood-subgradient density satisfies a self-consistency condition: the
// tangency equals the expectation under that restricted density.  glmbayes
// currently uses the closed-form normal-calibrated positions z* +/- omega as
// an approximation to it, which is excellent for near-Gaussian posteriors and
// can be catastrophically misplaced otherwise.
//
// This header exposes a short damped fixed-point iteration toward that
// condition, using the existing positions as starting values.  The partition,
// the standardization, the mixture construction and the accept/reject sampler
// are all untouched.
//
// VALIDITY IS NOT AT STAKE.  Any tangent to a globally concave log-likelihood
// is a valid upper envelope; optimality only affects the envelope's mass, i.e.
// the rejection constant.  Refinement can change speed but never correctness,
// and a cell whose refinement fails simply keeps the tangency glmbayes would
// have used anyway.

#ifndef GLMBAYES_ENVELOPE_REFINE_H
#define GLMBAYES_ENVELOPE_REFINE_H

#include "RcppArmadillo.h"
#include <string>

namespace glmbayes {
namespace env {

// Mills ratio M(x) = (1 - Phi(x)) / phi(x), x >= 0.  Continued fraction in the
// tail, direct evaluation near the origin.
double mills_ratio(double x);

// Mean of a standard normal truncated to (a, b).  Stable into the 1e5 tails,
// where forming exp(log phi - log tail) loses ~9 significant digits.
double trunc_norm_mean(double a, double b);

// log( Phi(b) - Phi(a) ), branch-selected, matching Set_Grid.cpp.
double log_box_prob(double a, double b);

// log(exp(x) + exp(y)), overflow-safe.  -Inf is the identity, +Inf absorbs.
double logaddexp_safe(double x, double y);

// Draw a mixture component from NORMALIZED LOG weights, using one uniform.
// Equivalent to the cumulative scan over exp(logw), but without ever
// exponentiating an individual weight, so a cell whose probability is below
// the smallest representable double keeps it instead of underflowing to zero.
// Bounded at n-1: the last cell absorbs any residual, so cumulative roundoff
// cannot run past the end of the array.
//
// Exposed so it can be tested directly on supplied log weights, with no
// envelope and no family involved.
int select_cell_log(const double* logw, int n, double u);

// Safeguarded Newton refinement of the fixed-point equation r = m(r).
// Finite-difference Jacobians accelerate modest grids; larger grids retain
// damped fixed-point updates. Every accepted cell update must have finite,
// nonincreasing envelope mass (up to roundoff for root corrections). accelerate=false retains the legacy iteration.
// Optional diagnostics report final-state residuals, per-cell damping, masses,
// and update counts without requiring verbose output.
//
//   G4        l1 x l2, parameters x grid points.  MODIFIED IN PLACE.
//   GIndex    l2 x l1, per-cell per-axis interval code (1,2,3 for a split axis,
//             4 for an axis EnvelopeSize left unsplit), as used by Set_Grid.
//   Lint      2 x l1, per-axis interval cutpoints.
//   maxit     ceiling on evaluated update passes. Convergence depends on the
//             model and grid; there is no universal contraction factor.
//   rho0      initial step fraction.  1.0 takes the full fixed-point step and
//             lets per-cell backtracking damp only the cells that need it;
//             smaller values merely guarantee slow geometric approach.
//   tol       stop when max|m(r) - r| < tol.
//   min_gain  stop when a pass improves the log total envelope mass by less
//             than this, so well-placed families do not pay for the
//             pathological iteration budget.
//
// dlogW_out receives the improvement in the LOG TOTAL ENVELOPE MASS,
// logsumexp_j(log W_j) before minus after (positive means the envelope
// shrank).  Note this is a log-sum-exp over cells, not a sum of logs.
// resid_out receives the final max|m(r) - r|.  Both are diagnostics, so we can
// tell whether the iteration budget is adequate and whether ordinary families
// can exit early.  Returns the number of passes performed.
int refine_tangency(Rcpp::NumericMatrix        G4,
                    const Rcpp::NumericMatrix& GIndex,
                    const arma::mat&           Lint,
                    const Rcpp::NumericVector& y,
                    const Rcpp::NumericMatrix& x,
                    const Rcpp::NumericMatrix& mu,
                    const Rcpp::NumericMatrix& P,
                    const Rcpp::NumericVector& alpha,
                    const Rcpp::NumericVector& wt,
                    const std::string&         family,
                    const std::string&         link,
                    bool                       use_opencl,
                    bool                       verbose,
                    int                        maxit    = 60,
                    double                     rho0     = 1.0,
                    double                     tol      = 1e-6,
                    double                     min_gain = 1e-3,
                    double*                    dlogW_out = NULL,
                    double*                    resid_out = NULL,
                    bool                       accelerate = true,
                    Rcpp::List*                diagnostics = NULL);

} // namespace env
} // namespace glmbayes

#endif
