// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <random>

#include <Rmath.h>   // libR Mathlib: Rf_pgamma, Rf_qgamma
#include "rng_utils.h"

namespace {
// Storage is thread-local; stream identity is assigned per output row.
thread_local std::mt19937_64 safe_rng_engine;
}

using namespace glmbayes::rng;


// log P(X <= d) for X ~ InvGamma(shape, rate)
// X = 1/Y, Y ~ Gamma(shape, rate)
// P(X <= d) = P(Y >= 1/d) = upper-tail Gamma at 1/d
double log_p_inv_gamma_safe(double dispersion,
                            double shape,
                            double rate) {
  double y = 1.0 / dispersion;
  
  return ::Rf_pgamma(y,
                      shape,
                      1.0 / rate,
                      /*lower_tail=*/0,
                      /*log_p=*/1);  // log upper tail
}


// Safe inverse-gamma CDF using libR Mathlib (Rf_pgamma)
double p_inv_gamma_safe(double dispersion,
                        double shape,
                        double rate) {
  // For X ~ InvGamma(shape, rate), Y = 1/X ~ Gamma(shape, rate)
  // So P(X <= d) = P(Y >= 1/d) = 1 - F_Y(1/d)
  double y = 1.0 / dispersion;
  
  double Fy = ::Rf_pgamma(y, shape, 1.0 / rate, /*lower_tail=*/1, /*log_p=*/0);
  
  return 1.0 - Fy;
}


double q_inv_gamma_safe(double u,
                        double shape,
                        double rate,
                        double disp_upper,
                        double disp_lower) {
  // log CDF at bounds
  double lg_low = log_p_inv_gamma_safe(disp_lower, shape, rate);
  double lg_upp = log_p_inv_gamma_safe(disp_upper, shape, rate);
  
  // ensure lg_upp >= lg_low (it should be, but be defensive)
  if (lg_upp <= lg_low) {
    // interval numerically collapsed: fall back to a boundary
    return disp_upper;
  }
  
  // We want p1 = p_low + u * (p_upp - p_low)
  // with p_low = exp(lg_low), p_upp = exp(lg_upp).
  //
  // Write:
  //   p1 = exp(lg_low) + u * (exp(lg_upp) - exp(lg_low))
  //       = exp(lg_low) * [1 + u * (exp(lg_upp - lg_low) - 1)]
  //
  // So:
  //   log p1 = lg_low + log(1 + u * (exp(lg_upp - lg_low) - 1))
  
  double diff   = lg_upp - lg_low;          // > 0
  double expdif = std::exp(diff);           // exp(lg_upp - lg_low)
  double inner  = 1.0 + u * (expdif - 1.0); // in (1, expdif]
  
  double lg_p1  = lg_low + std::log(inner);
  
  // Now p2 = 1 - p1, but we stay in log-space:
  // log p2 = log(1 - exp(lg_p1))
  double lg_p2 = std::log1p(-std::exp(lg_p1));
  
  double y = ::Rf_qgamma(lg_p2,
                         shape,
                         1.0 / rate,
                         /*lower_tail=*/1,
                         /*log_p=*/1);

  if (!std::isfinite(y) || y <= 0.0) {
    return disp_upper;  // safe fallback, never 0
  }

  return 1.0 / y;
}





namespace glmbayes {

namespace rng {

std::uint64_t seed_from_R() {
  // Only this function touches R's RNG. Exactly two draws per sampler job,
  // regardless of thread count, rejection work, or pilot/calibration size.
  const auto lo = static_cast<std::uint32_t>(::unif_rand() * 4294967296.0);
  const auto hi = static_cast<std::uint32_t>(::unif_rand() * 4294967296.0);
  return (static_cast<std::uint64_t>(hi) << 32) | lo;
}

void seed_draw(std::uint64_t job_seed, std::uint64_t draw_index) {
  std::seed_seq seed {
    static_cast<std::uint32_t>(job_seed),
    static_cast<std::uint32_t>(job_seed >> 32),
    static_cast<std::uint32_t>(draw_index),
    static_cast<std::uint32_t>(draw_index >> 32),
    std::uint32_t(0x676c6d62)
  };
  safe_rng_engine.seed(seed);
}

double runif_safe() {
  // Explicit 53-bit conversion avoids implementation-dependent behavior of
  // uniform_real_distribution. Exclude both endpoints for inverse CDFs.
  std::uint64_t bits;
  do { bits = safe_rng_engine() >> 11; } while (bits == 0);
  return static_cast<double>(bits) * (1.0 / 9007199254740992.0);
}


double rinvgamma_ct_safe(double shape,
                       double rate,
                       double disp_upper,
                       double disp_lower) {
  // draw uniform(0,1) from thread‑local RNG
  double p = runif_safe();
  
  // invert CDF at p to get inverse‑gamma draw
  // q_inv_gamma must be pure C++ math, no R calls
  return q_inv_gamma_safe(p, shape, rate, disp_upper, disp_lower);
  
  
}

// log P(disp_lower <= X <= disp_upper) for X ~ InvGamma(shape, rate)
// Consistent with: log_p_inv_gamma_safe, p_inv_gamma_safe,
// q_inv_gamma_safe, rinvgamma_ct_safe.
double log_p_inv_gamma_ct_safe(double disp_lower,
                               double disp_upper,
                               double shape,
                               double rate
                               )
{
  // log P(X <= d)
  double lg_low = log_p_inv_gamma_safe(disp_lower, shape, rate);
  double lg_upp = log_p_inv_gamma_safe(disp_upper, shape, rate);
  
  // Valid interval requires lg_upp > lg_low
  if (!(lg_upp > lg_low)) {
    return -std::numeric_limits<double>::infinity();
  }
  
  // Stable log-difference:
  // log( exp(lg_upp) - exp(lg_low) )
  // = lg_upp + log(1 - exp(lg_low - lg_upp))
  double diff = lg_low - lg_upp;
  return lg_upp + std::log1p(-std::exp(diff));
}

}
}
