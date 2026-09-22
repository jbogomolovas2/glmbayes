// envelope_refine.cpp -- see envelope_refine.h for the rationale.
//
// THE FIXED POINT.  On a cell A with tangency r and c = grad f2(r) (the
// gradient of the negative log posterior, i.e. what EnvelopeEval returns as
// cbars), the restricted envelope proposal is
//
//     q_r(z) prop exp(-c'z) phi(z) 1_A(z)   =>   Z ~ N(-c, I) restricted to A,
//
// which is exactly what rnorm_ct(logrt, loglt, -cbars, 1, logU) draws.  The cell
// envelope mass is
//
//     log W_A(r) = -NegLL(r) + r'c + 0.5 c'c + sum_i log[Phi(hi_i) - Phi(lo_i)]
//
// matching Set_LogP.cpp's logP(i,1) = logP(i,0) - NegLL + 0.5 c'c + r'c, with
// lo_i, hi_i the shifted bounds Set_Grid.cpp builds from Lint and cbars.
//
// Differentiating, with H = grad^2 L(r) negative semidefinite by concavity,
//
//     grad_r log W_A = H (m(r) - r),      m(r) = E_q[Z | A]
//
// so the stationary condition is r = m(r), and the step direction
// d = m(r) - r satisfies d' grad_r log W_A = d' H d <= 0.  It is a descent
// direction for the cell envelope mass by concavity alone, whatever the
// starting tangency.  That is what makes the damped update with backtracking
// safe rather than merely hopeful.
//
// Componentwise, since Z + c is a standard normal truncated to the fixed
// interval [Lint(0,i) + c_i, Lint(1,i) + c_i],
//
//     m_i(r) = -c_i + R(Lint(0,i) + c_i, Lint(1,i) + c_i)
//
// with R the truncated standard normal mean.  On an axis EnvelopeSize left
// unsplit (GIndex code 4, cell = the whole line) this is m_i = -c_i exactly.
//
// CONVERGENCE. Plain fixed-point updates can require very small damping
// on unsplit axes. For modest grids we solve a finite-difference Newton
// system for r-m(r), project into the cell, and accept only finite updates
// with nonincreasing cell mass. Larger grids use the original iteration.
// No universal convergence rate is assumed.
//
// NUMERICS.  R(a,b) has to survive arguments of order 1e5, which is exactly
// where the pathological cells live.  Computing it as
// exp(log phi(a) - log(1 - Phi(a))) fails there: both logs are about -5e9 and
// differ by ~11.5, so the subtraction discards ~9 significant digits (returns
// 100000.06367 where the truth is 100000.00001).  The stable route factors
// phi(a) out of numerator and denominator:
//
//     R(a,b) = (1 - e) / (M(a) - e M(b)),  e = exp(-(b^2-a^2)/2) = phi(b)/phi(a)
//     R(a,inf) = 1 / M(a)
//
// M(x) uses a continued fraction in the tail; the asymptotic series is only
// good to ~1e-4 at x = 5, far too coarse at the branch point.
//
// NOTE ON ||c||.  A large tangency gradient is NOT by itself a failure
// diagnostic.  After refinement the pathological case still carries
// ||c|| ~ 7e3 while its envelope is excellent.  What matters is whether the
// tangency and its cell are matched, i.e. the integrated cell mass, not the
// size of the tilt.

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include "envelope_refine.h"
#include "Envelopefuncs.h"

#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

using namespace Rcpp;

namespace glmbayes {
namespace env {

// ---------------------------------------------------------------- primitives

double mills_ratio(double x) {
  if (x < 3.0) {
    // logs are O(1)-O(10) here, no cancellation
    double lt = R::pnorm(x, 0.0, 1.0, /*lower*/0, /*log*/1);
    double ld = R::dnorm(x, 0.0, 1.0, /*log*/1);
    return std::exp(lt - ld);
  }
  long double cf = 0.0L;
  for (int k = 120; k >= 1; --k)
    cf = (long double)k / ((long double)x + cf);
  return (double)(1.0L / ((long double)x + cf));
}

namespace {

double R_pos(double a, double b) {
  if (!R_finite(b)) return 1.0 / mills_ratio(a);
  // (b-a)(b+a) rather than b*b - a*a: better relative accuracy when the
  // endpoints are large and close together
  double d = 0.5 * (b - a) * (b + a);
  double e = std::exp(-d);
  double num = -std::expm1(-d);
  double den = mills_ratio(a) - e * mills_ratio(b);
  // With the factored formula this should not happen.  If it does, fail safe:
  // return NaN so the caller skips this coordinate and keeps the tangency
  // glmbayes would have used.  Returning a midpoint would be inventing an
  // answer, and we have already shown the midpoint is badly wrong for a
  // narrow deep-tail interval once a*h = O(1).
  if (!(den > 0.0) || !R_finite(den) || !R_finite(num)) return R_NaReal;
  return num / den;
}

// log(exp(x) + exp(y)) without overflow
inline double logaddexp(double x, double y) {
  // -Inf is the identity; +Inf absorbs.  The naive !R_finite test returns the
  // OTHER argument for +Inf, which is exactly wrong in a routine whose whole
  // purpose is catastrophic envelope weights.
  if (ISNAN(x) || ISNAN(y)) return R_NaReal;
  if (x == R_PosInf || y == R_PosInf) return R_PosInf;
  if (x == R_NegInf) return y;
  if (y == R_NegInf) return x;
  double M = std::max(x, y), m = std::min(x, y);
  return M + std::log1p(std::exp(m - M));
}

} // anonymous namespace

double trunc_norm_mean(double a, double b) {
  if (!(a < b)) return a;
  if (!R_finite(a) && !R_finite(b)) return 0.0;      // unsplit axis: m = -c
  if (b <= 0.0) return -trunc_norm_mean(-b, -a);     // reflect the left tail
  if (a >= 0.0) return R_pos(a, b);                  // stable Mills branch
  // straddles zero: phi and the box probability are both O(1)
  double pa = R_finite(a) ? R::dnorm(a, 0.0, 1.0, 0) : 0.0;
  double pb = R_finite(b) ? R::dnorm(b, 0.0, 1.0, 0) : 0.0;
  double lo = R_finite(a) ? R::pnorm(a, 0.0, 1.0, 1, 0) : 0.0;
  double hi = R_finite(b) ? R::pnorm(b, 0.0, 1.0, 1, 0) : 1.0;
  double den = hi - lo;
  if (!(den > 0.0)) return R_NaReal;   // fail safe, do not invent a value
  return (pa - pb) / den;
}

double log_box_prob(double a, double b) {
  if (!(a < b)) return R_NegInf;
  if (!R_finite(a) && !R_finite(b)) return 0.0;
  if (!R_finite(a)) return R::pnorm(b, 0.0, 1.0, 1, 1);
  if (!R_finite(b)) return R::pnorm(a, 0.0, 1.0, 0, 1);
  double g1 = -a, g2 = b;
  if (g1 >= g2) {                                    // upper-anchored
    double lup = R::pnorm(b, 0.0, 1.0, 1, 1);
    double llo = R::pnorm(a, 0.0, 1.0, 1, 1);
    return lup + std::log(-std::expm1(llo - lup));
  } else {                                           // lower-anchored
    double lgt_a = R::pnorm(a, 0.0, 1.0, 0, 1);
    double lgt_b = R::pnorm(b, 0.0, 1.0, 0, 1);
    return lgt_a + std::log(-std::expm1(lgt_b - lgt_a));
  }
}

// -------------------------------------------------- log-space mixture draw

double logaddexp_safe(double x, double y) {
  if (ISNAN(x) || ISNAN(y)) return R_NaReal;
  if (x == R_PosInf || y == R_PosInf) return R_PosInf;
  if (x == R_NegInf) return y;
  if (y == R_NegInf) return x;
  double M = std::max(x, y), m = std::min(x, y);
  return M + std::log1p(std::exp(m - M));
}

int select_cell_log(const double* logw, int n, double u) {
  if (n <= 1) return 0;
  if (!(u > 0.0)) return 0;          // u == 0 selects the first cell
  const double lu = std::log(u);
  double acc = R_NegInf;
  int J = 0;
  while (J < n - 1) {                // last cell is the fallback
    acc = logaddexp_safe(acc, logw[J]);
    if (lu <= acc) return J;
    ++J;
  }
  return n - 1;
}

// ------------------------------------------------------------- cell bounds
// Exactly the codes Set_Grid.cpp uses to build Down/Up, so the package has one
// convention rather than two.  Code 4 is an axis EnvelopeSize left unsplit.
namespace {
inline void cell_bounds(int code, double lint0, double lint1, double cb,
                        double* lo, double* hi) {
  switch (code) {
    case 1:  *lo = R_NegInf;     *hi = lint0 + cb;  break;
    case 2:  *lo = lint0 + cb;   *hi = lint1 + cb;  break;
    case 3:  *lo = lint1 + cb;   *hi = R_PosInf;    break;
    default: *lo = R_NegInf;     *hi = R_PosInf;    break;   // code 4
  }
}
} // anonymous namespace

// ------------------------------------------------------------- refinement

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
                    int                        maxit,
                    double                     rho0,
                    double                     tol,
                    double                     min_gain,
                    double*                    dlogW_out,
                    double*                    resid_out,
                    bool                       accelerate,
                    Rcpp::List*                diagnostics)
{
  const int l1 = G4.nrow();      // parameters
  const int l2 = G4.ncol();      // grid cells

  if (dlogW_out) *dlogW_out = 0.0;
  if (resid_out) *resid_out = R_NaReal;
  maxit = std::max(0, maxit);

  std::vector<double> logW_cur(l2), rho(l2, rho0);
  std::vector<double> cb_cur((size_t)l2 * l1);
  Rcpp::NumericMatrix G4_try(l1, l2);

  // total envelope mass is a LOG-SUM-EXP over cells, not a sum of logs:
  // sum_j W_j, reported in logs.  Summing logW_j would be the log of the
  // product, which is not a meaningful quantity here.
  double lse_0 = R_NegInf, lse_prev = R_NegInf;

  // --- evaluate at the incoming tangencies -------------------------------
  {
    Rcpp::List ev = EnvelopeEval(G4, y, x, mu, P, alpha, wt,
                                 family, link, use_opencl, false);
    Rcpp::NumericVector NegLL = ev["NegLL"];
    arma::mat cbars = Rcpp::as<arma::mat>(ev["cbars"]);   // l2 x l1

    for (int j = 0; j < l2; ++j) {
      double quad = 0.0, lin = 0.0, box = 0.0;
      for (int i = 0; i < l1; ++i) {
        double cb = cbars(j, i);
        cb_cur[(size_t)j*l1 + i] = cb;
        quad += 0.5 * cb * cb;
        lin  += G4(i, j) * cb;
        double lo, hi;
        cell_bounds((int) GIndex(j, i), Lint(0, i), Lint(1, i), cb, &lo, &hi);
        box += log_box_prob(lo, hi);
      }
      logW_cur[j] = -NegLL[j] + lin + quad + box;
      lse_0 = logaddexp(lse_0, logW_cur[j]);
    }
    lse_prev = lse_0;

    // The proposal evaluation has finite checks; the initial one did not.
    // If the incoming envelope cannot even be scored, leave G4 untouched --
    // the promise is that a failed refinement retains glmbayes's tangencies.
    if (ISNAN(lse_0) || lse_0 == R_PosInf) {
      if (verbose)
        Rcpp::Rcout << "[EnvelopeBuild:refine_tangency] initial envelope mass "
                       "not usable; tangencies left unchanged\n";
      return 0;
    }
  }

  // m(r) is cheap once the likelihood gradient is known. Differentiate the
  // fixed-point equation F(r)=r-m(r), rather than the potentially enormous
  // cell mass: J_F delta = m(r)-r. Numerical derivatives change only proposed
  // tangencies; the original envelope mass check remains authoritative.
  auto targets = [&](const NumericMatrix& cbars) {
    NumericMatrix target(l1, l2);
    for (int j = 0; j < l2; ++j) for (int i = 0; i < l1; ++i) {
      double cb = cbars(j, i), lo, hi;
      cell_bounds((int) GIndex(j, i), Lint(0, i), Lint(1, i), cb, &lo, &hi);
      target(i, j) = -cb + trunc_norm_mean(lo, hi);
    }
    return target;
  };
  auto current_targets = [&]() {
    NumericMatrix cbars(l2, l1);
    for (int j = 0; j < l2; ++j) for (int i = 0; i < l1; ++i)
      cbars(j, i) = cb_cur[(size_t)j*l1 + i];
    return targets(cbars);
  };
  auto residual = [&](const NumericMatrix& target) {
    double result = 0.0;
    for (int j = 0; j < l2; ++j) for (int i = 0; i < l1; ++i) {
      const double delta = target(i,j) - G4(i,j);
      if (!R_finite(delta)) return R_PosInf;
      result = std::max(result, std::fabs(delta));
    }
    return result;
  };
  // Bound Jacobian storage (~16 MB) and the number of extra likelihood
  // evaluations per pass. Large grids retain the original low-memory path.
  const bool use_newton = accelerate && l1 <= 16 &&
      (double) l1 * l1 * l2 <= 2000000.0;
  NumericVector accepted(l2), rejected(l2), newton_steps(l2);
  int iter = 0;
  double resid = R_PosInf;
  NumericMatrix target = current_targets();

  while (iter < maxit) {
    Rcpp::checkUserInterrupt();
    resid = residual(target);
    if (resid < tol) break;
    NumericMatrix step(l1, l2);
    for (int j = 0; j < l2; ++j) for (int i = 0; i < l1; ++i)
      step(i,j) = R_finite(target(i,j)) ? target(i,j)-G4(i,j) : 0.0;
    std::vector<bool> used_newton(l2, false);
    if (use_newton) {
      std::vector<arma::mat> jac(l2, arma::eye<arma::mat>(l1,l1));
      for (int axis = 0; axis < l1; ++axis) {
        NumericMatrix perturbed = clone(G4);
        NumericVector h(l2);
        for (int j = 0; j < l2; ++j) {
          h[j] = 1e-5 * std::max(1.0, std::fabs(G4(axis,j)));
          perturbed(axis,j) += h[j];
        }
        List ev = EnvelopeEval(perturbed,y,x,mu,P,alpha,wt,
                               family,link,use_opencl,false);
        NumericMatrix pt = targets(as<NumericMatrix>(ev["cbars"]));
        for (int j = 0; j < l2; ++j) for (int i = 0; i < l1; ++i)
          jac[j](i,axis) -= (pt(i,j)-target(i,j))/h[j];
      }
      for (int j = 0; j < l2; ++j) {
        arma::vec rhs(l1), delta;
        for (int i = 0; i < l1; ++i) rhs[i] = step(i,j);
        if (jac[j].is_finite() && rhs.is_finite() &&
            arma::solve(delta,jac[j],rhs,arma::solve_opts::no_approx) &&
            delta.is_finite()) {
          for (int i = 0; i < l1; ++i) step(i,j) = delta[i];
          used_newton[j] = true;
        }
      }
    }
    bool any_move = false;
    for (int j = 0; j < l2; ++j) for (int i = 0; i < l1; ++i) {
      double cand = G4(i,j) + rho[j]*step(i,j);
      if (!R_finite(cand)) cand = G4(i,j);
      const int code = (int) GIndex(j,i);
      // Project Newton candidates into their fixed cell. The legacy convex
      // combination already stays there. Neither path changes the partition.
      if (code == 1 || code == 2)
        cand = std::min(cand, Lint(code == 1 ? 0 : 1,i));
      if (code == 2 || code == 3)
        cand = std::max(cand, Lint(code == 2 ? 0 : 1,i));
      G4_try(i,j) = cand;
      if (cand != G4(i,j)) any_move = true;
    }
    if (!any_move) break;
    ++iter;
    if (verbose)
      Rcpp::Rcout << "    [refine] pass " << iter << ": resid = " << resid
                  << ", newton = " << use_newton << "\n";

    // --- evaluate the proposal ------------------------------------------
    Rcpp::List ev2 = EnvelopeEval(G4_try, y, x, mu, P, alpha, wt,
                                  family, link, use_opencl, false);
    Rcpp::NumericVector NegLL2 = ev2["NegLL"];
    arma::mat cbars2 = Rcpp::as<arma::mat>(ev2["cbars"]);

    // --- accept per cell only if its envelope mass did not increase ------
    for (int j = 0; j < l2; ++j) {
      double quad = 0.0, lin = 0.0, box = 0.0;
      bool ok = R_finite(NegLL2[j]);
      for (int i = 0; ok && i < l1; ++i) {
        double cb = cbars2(j, i);
        if (!R_finite(cb)) { ok = false; break; }
        quad += 0.5 * cb * cb;
        lin  += G4_try(i, j) * cb;
        double lo, hi;
        cell_bounds((int) GIndex(j, i), Lint(0, i), Lint(1, i), cb, &lo, &hi);
        double lb = log_box_prob(lo, hi);
        if (!R_finite(lb)) { ok = false; break; }
        box += lb;
      }
      double logW_new = ok ? (-NegLL2[j] + lin + quad + box) : R_PosInf;

      // Close to the root, the remaining mass improvement is quadratic
      // in the residual and can fall below floating-point resolution. Permit
      // a Newton correction only within a summation-roundoff allowance and
      // only when it at least halves an already small cell residual. This
      // changes the tangent, never the likelihood or rejection acceptance rule.
      double old_resid = 0.0, new_resid = 0.0;
      for (int i = 0; i < l1; ++i) {
        old_resid = std::max(old_resid, std::fabs(target(i,j)-G4(i,j)));
        double lo, hi, cb = cbars2(j,i);
        cell_bounds((int) GIndex(j,i),Lint(0,i),Lint(1,i),cb,&lo,&hi);
        double delta = -cb + trunc_norm_mean(lo,hi) - G4_try(i,j);
        new_resid = R_finite(delta) ? std::max(new_resid,std::fabs(delta)) : R_PosInf;
      }
      const double roundoff = 16.0 * std::numeric_limits<double>::epsilon() *
          (1.0 + std::fabs(NegLL2[j]) + std::fabs(lin) + std::fabs(quad) + std::fabs(box));
      const bool root_correction = used_newton[j] && old_resid < 1e-3 &&
          new_resid < 0.5*old_resid && logW_new <= logW_cur[j] + roundoff;
      if (ok && R_finite(logW_new) && (logW_new <= logW_cur[j] || root_correction)) {
        for (int i = 0; i < l1; ++i) {
          G4(i, j) = G4_try(i, j);
          cb_cur[(size_t)j*l1 + i] = cbars2(j, i);
        }
        logW_cur[j] = logW_new;
        accepted[j] += 1;
        if (used_newton[j]) {
          newton_steps[j] += 1;
          rho[j] = std::min(rho0, 2.0*rho[j]);
        }
      } else {
        rejected[j] += 1;
        rho[j] *= 0.5;              // shorter step for this cell next pass
      }
    }

    target = current_targets();
    resid = residual(target);

    // --- second stopping rule: total envelope mass has stopped improving --
    double lse = R_NegInf;
    for (int j = 0; j < l2; ++j) lse = logaddexp(lse, logW_cur[j]);
    double gain = lse_prev - lse;                 // positive = envelope shrank
    lse_prev = lse;
    if (verbose)
      Rcpp::Rcout << "    [refine] pass " << iter << ": gain = " << gain
                  << ", lse = " << lse << "\n";
    // the gain rule must never fire while the fixed point is still far
    // away: logsumexp is dominated by the worst cell, so the total can
    // look flat for several passes while individual cells improve hugely
    if (R_finite(gain) && gain >= 0.0 && gain < min_gain && resid < 1e-3) break;
  }

  resid = residual(current_targets());
  double lse_fin = R_NegInf;
  for (int j = 0; j < l2; ++j) lse_fin = logaddexp(lse_fin, logW_cur[j]);

  if (diagnostics) {
    *diagnostics = List::create(_["iterations"]=iter, _["residual"]=resid,
      _["converged"]=(resid < tol), _["accelerated"]=use_newton,
      _["log_mass_initial"]=lse_0, _["log_mass_final"]=lse_fin,
      _["cell_log_mass"]=wrap(logW_cur), _["step_fraction"]=wrap(rho),
      _["accepted_updates"]=accepted, _["rejected_updates"]=rejected,
      _["newton_updates"]=newton_steps);
  }
  if (resid_out) *resid_out = resid;
  if (dlogW_out) *dlogW_out = lse_0 - lse_fin;     // positive = improvement

  if (verbose) {
    // per-cell logW at exit, so the internal quantity can be checked against
    // the package's own PLSD after Set_Grid/Set_LogP have run
    double mn = R_PosInf, mx = R_NegInf;
    for (int j = 0; j < l2; ++j)
      if (R_finite(logW_cur[j])) { mn = std::min(mn, logW_cur[j]);
                                   mx = std::max(mx, logW_cur[j]); }
    Rcpp::Rcout << "[EnvelopeBuild:refine_tangency] logW range at exit: ["
                << mn << ", " << mx << "]\n";
    Rcpp::Rcout << "[EnvelopeBuild:refine_tangency] iterations = " << iter
                << ", log total envelope mass: " << lse_0 << " -> " << lse_fin
                << " (improvement " << (lse_0 - lse_fin) << ")"
                << ", final max|m(r)-r| = " << resid << "\n";
  }

  return iter;
}

} // namespace env
} // namespace glmbayes
