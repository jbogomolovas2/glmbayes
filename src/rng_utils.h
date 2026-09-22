// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/**
 * @file rng_utils.h
 * @brief Thread‑safe uniform sampling and truncated‑distribution utilities
 *        used throughout glmbayes.
 *
 * @namespace glmbayes::rng
 * @brief Random‑number generators for uniform draws and truncated Normal
 *        and truncated inverse‑Gamma distributions, including numerically
 *        stable “central truncation” variants.
 *
 * @section ImplementedIn
 *   These declarations are implemented in:
 *     - rng_utils.cpp
 *     - (optionally) supporting truncated‑distribution source files
 *
 * @section UsedBy
 *   These functions are consumed by:
 *     - RcppParallel‑based samplers (Normal, Normal–Gamma, and GLM samplers)
 *     - Envelope construction routines requiring truncated inverse‑Gamma draws
 *     - Diagnostic and simulation utilities across the glmbayes backend
 *
 * @section Responsibilities
 *   Provides:
 *     - a thread‑safe uniform RNG (`runif_safe`) suitable for parallel contexts,
 *     - truncated inverse‑Gamma samplers (`rinvgamma_ct_safe`, `rinvgamma_ct`)
 *       for dispersion parameters restricted to \([d_{\min}, d_{\max}]\),
 *     - truncated Normal sampler (`rnorm_ct`) for draws restricted to
 *       \([l_{\mathrm{lower}}, l_{\mathrm{upper}}]\),
 *     - log‑density evaluation for truncated inverse‑Gamma distributions
 *       (`log_p_inv_gamma_ct_safe`).
 *
 *   The suffix `_ct` denotes “central truncation,” meaning all functions
 *   operate on a parameter restricted to a finite interval and are designed
 *   to remain numerically stable even when the bounds are tight or the shape
 *   and rate parameters induce heavy tails.
 */

#ifndef RNG_UTILS_H
#define RNG_UTILS_H

#include <cstdint>

namespace glmbayes{

namespace rng {

// Draw a job key from R's RNG on the main thread, inside an RNGScope.
std::uint64_t seed_from_R();
// Reset the current thread's stream for one output row. Workers share the
// immutable job key; stream identity depends on the row, never the thread.
void seed_draw(std::uint64_t job_seed, std::uint64_t draw_index);
// Thread-safe uniform RNG (0, 1), used only after seed_draw().
double runif_safe();

double rinvgamma_ct_safe(double shape,
                       double rate,
                       double disp_upper,
                       double disp_lower);



double  rnorm_ct(double lgrt,double lglt,double mu,double sigma);
// The interval mass preserves both boundaries when a complementary tail
// rounds to zero. Envelope samplers should use this overload.
double rnorm_ct(double lgrt, double lglt, double mu, double sigma, double log_mass);

double rinvgamma_ct(double shape,double rate,double disp_upper,double disp_lower);

double log_p_inv_gamma_ct_safe(double disp_lower,
                               double disp_upper,
                               double shape,
                               double rate);

}
}
#endif
