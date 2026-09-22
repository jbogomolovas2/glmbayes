// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>
#include <algorithm>
#include "rng_utils.h"


#include <Rmath.h>       // For R::qnorm


// Don't add mutex if EMSCRIPTEN

#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
#include <tbb/mutex.h>   // For thread locking
tbb::mutex qnorm_mutex;  // Local mutex for this file
#endif


using namespace Rcpp;
using namespace glmbayes::rng;


double safe_qnorm_logp(double logp, double mu, double sigma, bool lower_tail) {
#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
    tbb::mutex::scoped_lock lock(qnorm_mutex);
#endif  
  return R::qnorm(logp, mu, sigma, lower_tail, true);  // log.p = TRUE
}

namespace glmbayes {
namespace rng {

// Compatibility for callers with only the two tail probabilities. These
// cannot encode both endpoints in extreme tails; new envelopes supply logU.
double rnorm_ct(double lgrt, double lglt, double mu, double sigma) {
  const bool lower = lgrt >= lglt;
  const double anchor = lower ? lglt : lgrt;
  const double opposite = lower ? lgrt : lglt;
  const double other_endpoint = std::log(-std::expm1(opposite));
  const double log_mass = anchor + std::log(-std::expm1(other_endpoint - anchor));
  if (std::isnan(log_mass)) return R_NaN;
  return rnorm_ct(lgrt, lglt, mu, sigma, log_mass);
}

double rnorm_ct(double lgrt, double lglt, double mu, double sigma, double log_mass) {
  if (std::isnan(log_mass)) return rnorm_ct(lgrt, lglt, mu, sigma);
  const bool lower = lgrt >= lglt;
  const double anchor = lower ? lglt : lgrt;
  // Subtract U times the interval mass from the larger endpoint probability,
  // using the stable tail. Never recover an endpoint by complementing a tail
  // that may already have rounded to zero (e.g. pnorm(100, log.p=TRUE)).
  double U;
  do { U = runif_safe(); } while (U == 0.0);
  const double fraction = std::exp(std::min(0.0, log_mass - anchor));
  const double logp = anchor + std::log1p(-U * fraction);
  return safe_qnorm_logp(logp, mu, sigma, lower);
}

} // namespace rng
} // namespace glmbayes
