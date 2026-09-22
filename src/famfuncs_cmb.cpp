// famfuncs_cmb.cpp -- Conway-Maxwell-binomial family for glmbayes.
//
// CMB in its natural parametrisation:
//   log f(y | theta, nu) = theta*y + nu*t(y) - kappa(theta, nu)
//   t(y)   = log C(m, y)
//   kappa  = log sum_{j=0}^{m} exp(theta*j + nu*t(j))
// Globally concave in (beta, gamma) because kappa is a log-sum-exp of affine
// functions.  Likelihood, gradient and Hessian all come from the same
// (m+1)-term sum: no inner solve, no implicit derivatives.
//
// TWO LINEAR PREDICTORS.  Both theta and nu carry covariates:
//   theta_i = x_i' beta,   nu_i = z_i' gamma,
// with identity links on both (theta is already the natural parameter, and nu
// must be free to go negative for the boundary regime).  They are supplied as
// one STACKED design of 2N rows:
//
//     x = [ X  0 ]   rows 1..N     -> theta_1..theta_N
//         [ 0  Z ]   rows N+1..2N  -> nu_1..nu_N
//
// so that eta = alpha + x*b splits by ROW into (theta, nu).  A scalar nu is
// the special case Z = 1_N.
//
// WHY ROWS AND NOT COORDINATES.  glmb_Standardize_Model maps b -> L*b and
// x -> x*LInv, and EnvelopeBuild/rNormalGLM_std receive the rotated x with
// the raw y, alpha, wt.  A GLM is unaffected because its likelihood sees b
// only through x*b.  CMB sees b through two functionals; carrying nu as a
// coordinate b[nu_index] is not invariant, since after rotation
// nu = e'LInv b is a dense combination.  The envelope would then evaluate the
// wrong nu while still satisfying max log h <= 0 -- silent and undetectable
// downstream.  Expressed as rows, both functionals survive any
// right-multiplication of x exactly.
//
// CALLING CONTRACT (mirrors f2_f3_binomial_logit):
//   b      l2 x m1  parameters x grid points   (b = G4 = trans(G3))
//   y      l1 = 2N  PROPORTIONS in rows 1..N; rows N+1..2N are placeholders
//   x      l1 x l2  stacked design as above
//   alpha  l1       offset added to the linear predictor
//   wt     l1 = 2N  trial counts m in rows 1..N; rows N+1..2N placeholders
//   grad returned m1 x l2.

#include <cmath>
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <limits>
#include <RcppParallel.h>
#include <Rmath.h>
#include "famfuncs.h"
#include "progress_utils.h"

#include <map>
#include <vector>
#include <numeric>

using namespace Rcpp;
using namespace RcppParallel;
using namespace glmbayes::progress;

namespace glmbayes {
namespace fam {

namespace {

struct CMBData {
  int N;
  std::vector<int>    m;
  std::vector<double> yc;
  std::vector<double> ty;
  std::map<int, std::vector<double> > tab;
};

// Templated on the container so the same code serves NumericVector and
// RcppParallel::RVector<double>.  Uses only R::lchoose, pure Rmath with no
// allocation, so it is safe from a worker thread.
template <typename VY, typename VW>
CMBData cmb_prep_t(const VY& y, const VW& wt, int N) {
  CMBData d;
  d.N = N;
  d.m.resize(N);
  d.yc.resize(N);
  d.ty.resize(N);

  bool all_one = true;

  for (int j = 0; j < N; ++j) {
    int mj = (int) std::lround(wt[j]);
    if (mj < 1) Rcpp::stop("CMB: trial count m < 1");
    if (mj != 1) all_one = false;

    double ycj = (double) std::lround(wt[j] * y[j]);
    if (ycj < 0.0 || ycj > (double) mj)
      Rcpp::stop("CMB: implied count outside [0, m]");

    d.m[j]  = mj;
    d.yc[j] = ycj;

    if (d.tab.find(mj) == d.tab.end()) {
      std::vector<double> t(mj + 1);
      for (int k = 0; k <= mj; ++k)
        t[k] = R::lchoose((double) mj, (double) k);
      d.tab[mj] = t;
    }
    d.ty[j] = d.tab[mj][(int) ycj];
  }

  if (all_one)
    Rcpp::stop("CMB: all m = 1; t(y) is identically 0 and nu is unidentified");

  return d;
}

inline void cmb_moments(double theta, double nu,
                        const std::vector<double>& t,
                        double* lZ, double* EY, double* Et) {
  const int mj = (int) t.size() - 1;

  double M = R_NegInf;
  for (int k = 0; k <= mj; ++k) {
    double lw = theta * (double) k + nu * t[k];
    if (lw > M) M = lw;
  }

  double s = 0.0, sy = 0.0, st = 0.0;
  for (int k = 0; k <= mj; ++k) {
    double e = std::exp(theta * (double) k + nu * t[k] - M);
    s += e;
    if (EY != NULL) { sy += e * (double) k; st += e * t[k]; }
  }

  *lZ = M + std::log(s);
  if (EY != NULL) { *EY = sy / s; *Et = st / s; }
}

inline int cmb_N(int l1) {
  if (l1 % 2 != 0)
    Rcpp::stop("CMB: nrow(x) must be even (2N: N theta rows then N nu rows). "
               "Build the design with cmb_augment().");
  return l1 / 2;
}

} // anonymous namespace

// ------------------------------------------------------------------- f2/f3

Rcpp::List f2_f3_cmb(
    Rcpp::NumericMatrix  b,
    Rcpp::NumericVector  y,
    Rcpp::NumericMatrix  x,
    Rcpp::NumericMatrix  mu,
    Rcpp::NumericMatrix  P,
    Rcpp::NumericVector  alpha,
    Rcpp::NumericVector  wt,
    int                  progbar
) {
  const int l1 = x.nrow();
  const int l2 = x.ncol();
  const int m1 = b.ncol();
  const int N  = cmb_N(l1);

  if (N < 1) Rcpp::stop("CMB: empty design");
  if (b.nrow() != l2) Rcpp::stop("CMB: nrow(b) must equal ncol(x)");
  if (y.size() != l1 || wt.size() != l1 || alpha.size() != l1)
    Rcpp::stop("CMB: y, wt and alpha must have length nrow(x)");

  CMBData d = cmb_prep_t(y, wt, N);

  arma::mat x2    (x.begin(),     l1, l2, false);
  arma::mat alpha2(alpha.begin(), l1, 1,  false);
  arma::mat mu2   (mu.begin(),    l2, 1,  false);
  arma::mat P2    (P.begin(),     l2, l2, false);

  Rcpp::NumericMatrix b2temp(l2, 1);
  Rcpp::NumericMatrix bmu(l2, 1);
  arma::mat bmu2(bmu.begin(), l2, 1, false);

  Rcpp::NumericVector qf(m1);

  Rcpp::NumericMatrix out(l2, m1);
  arma::mat out2(out.begin(), l2, m1, false);

  arma::colvec eta(l1);
  arma::colvec sc(l1);

  for (int i = 0; i < m1; ++i) {
    Rcpp::checkUserInterrupt();
    if (progbar == 1) {
      progress_bar(i, m1 - 1);
      if (i == m1 - 1) Rcpp::Rcout << "" << std::endl;
    }

    b2temp = b(Range(0, l2 - 1), Range(i, i));
    arma::mat b2(b2temp.begin(), l2, 1, false);

    bmu2 = b2 - mu2;
    double prior = 0.5 * arma::as_scalar(bmu2.t() * P2 * bmu2);

    eta = alpha2 + x2 * b2;

    double ll = 0.0;
    for (int j = 0; j < N; ++j) {
      const double th = eta(j);
      const double nu = eta(N + j);
      double lZ, EY, Et;
      cmb_moments(th, nu, d.tab[d.m[j]], &lZ, &EY, &Et);

      ll += th * d.yc[j] + nu * d.ty[j] - lZ;
      sc(j)     = d.yc[j] - EY;
      sc(N + j) = d.ty[j] - Et;
    }

    qf(i) = -ll + prior;

    Rcpp::NumericMatrix::Column outcol = out(_, i);
    arma::mat outtemp2(outcol.begin(), l2, 1, false);
    outtemp2 = P2 * bmu2 - x2.t() * sc;
  }

  arma::mat grad = trans(out2);

  return Rcpp::List::create(
    Rcpp::Named("qf")   = qf,
    Rcpp::Named("grad") = grad
  );
}

// ----------------------------------------------------------------------- f1

NumericVector f1_cmb(NumericMatrix b, NumericVector y, NumericMatrix x,
                     NumericVector alpha, NumericVector wt) {
  const int l1 = x.nrow();
  const int l2 = x.ncol();
  const int m1 = b.ncol();
  const int N  = cmb_N(l1);

  CMBData d = cmb_prep_t(y, wt, N);

  arma::mat x2    (x.begin(),     l1, l2, false);
  arma::mat alpha2(alpha.begin(), l1, 1,  false);

  Rcpp::NumericMatrix b2temp(l2, 1);
  NumericVector res(m1);
  arma::colvec eta(l1);

  for (int i = 0; i < m1; ++i) {
    b2temp = b(Range(0, l2 - 1), Range(i, i));
    arma::mat b2(b2temp.begin(), l2, 1, false);

    eta = alpha2 + x2 * b2;

    double ll = 0.0;
    for (int j = 0; j < N; ++j) {
      double lZ;
      cmb_moments(eta(j), eta(N + j), d.tab[d.m[j]], &lZ, NULL, NULL);
      ll += eta(j) * d.yc[j] + eta(N + j) * d.ty[j] - lZ;
    }
    res(i) = -ll;
  }
  return res;
}

// ------------------------------------------------------------------- f2, f3

NumericVector f2_cmb(NumericMatrix b, NumericVector y, NumericMatrix x,
                     NumericMatrix mu, NumericMatrix P, NumericVector alpha,
                     NumericVector wt, int progbar) {
  Rcpp::List r = f2_f3_cmb(b, y, x, mu, P, alpha, wt, progbar);
  return Rcpp::as<NumericVector>(r["qf"]);
}

arma::mat f3_cmb(NumericMatrix b, NumericVector y, NumericMatrix x,
                 NumericMatrix mu, NumericMatrix P, NumericVector alpha,
                 NumericVector wt, int progbar) {
  Rcpp::List r = f2_f3_cmb(b, y, x, mu, P, alpha, wt, progbar);
  return Rcpp::as<arma::mat>(r["grad"]);
}

// ------------------------------------------------- thread-safe rmat variant
// Used by the parallel sampler's accept step.  No R allocation: only Rmath
// and Armadillo views over RcppParallel memory.

arma::vec f2_cmb_rmat(
    const RMatrix<double>& b,
    const RVector<double>& y,
    const RMatrix<double>& x,
    const RMatrix<double>& mu,
    const RMatrix<double>& P,
    const RVector<double>& alpha,
    const RVector<double>& wt,
    const int progbar
) {
  std::size_t l1 = x.nrow();
  std::size_t l2 = x.ncol();
  std::size_t m1 = b.ncol();
  int N = cmb_N((int) l1);

  CMBData d = cmb_prep_t(y, wt, N);

  arma::mat b2full(const_cast<double*>(&*b.begin()),     l2, m1, false);
  arma::mat x2    (const_cast<double*>(&*x.begin()),     l1, l2, false);
  arma::mat mu2   (const_cast<double*>(&*mu.begin()),    l2, 1,  false);
  arma::mat P2    (const_cast<double*>(&*P.begin()),     l2, l2, false);
  arma::mat alpha2(const_cast<double*>(&*alpha.begin()), l1, 1,  false);

  arma::vec res(m1, arma::fill::none);
  arma::mat bmu(l2, 1, arma::fill::none);

  for (std::size_t i = 0; i < m1; ++i) {
    arma::mat b_i(b2full.colptr(i), l2, 1, false);

    bmu = b_i - mu2;
    double mahal = 0.5 * arma::as_scalar(bmu.t() * P2 * bmu);

    double ll = 0.0;
    for (int j = 0; j < N; ++j) {
      double th = alpha2(j, 0)     + arma::dot(x2.row(j),     b_i);
      double nu = alpha2(N + j, 0) + arma::dot(x2.row(N + j), b_i);
      double lZ;
      cmb_moments(th, nu, d.tab[d.m[j]], &lZ, NULL, NULL);
      ll += th * d.yc[j] + nu * d.ty[j] - lZ;
    }

    res(i) = -ll + mahal;
  }

  return res;
}

} // namespace fam
} // namespace glmbayes
