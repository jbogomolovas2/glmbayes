# Maintainer review: CMB contribution and tangency refinement

### Review at a glance

This review covers contribution `b971596511c90f33de290497b2c7e8feb1a4e580` against base `7959a42d32f822a4b65ddfaf94927190ad5db536`. The four sections below retain the complete inventory and supporting evidence.

1. **[Changed signatures](#1-changed-signatures).** Public `EnvelopeSort` adds `logPLSD=NULL`. The internal/native `EnvelopeBuild` chain adds `refine=TRUE` and `refine_maxit=60`; worker constructors and a truncated-normal overload also change. Most public fitting signatures are unchanged. **Compatibility defect:** the manual envelope wrapper still passes 16 arguments to an 18-argument native entry and crashes on the public non-Gaussian path.

2. **[Additional returned items](#2-additional-returned-items).** Regular envelopes add `logPLSD` and nested `refinement` diagnostics; fit/summary objects add `diagnostics`. Independent Normal-Gamma fits now retain their envelope; direct Gaussian fits return `Envelope=NULL`. Sorted envelope arrays and refinement vectors use different row orders, linked by `cell_grid_index`. Types, availability and consumers are listed below.

3. **[New functions and files](#3-new-functions-and-files).** CMB adds `cmb`, `cmb_augment`, `cmb_fitted`, `cmb_nu`, and native likelihood/derivative kernels. Separate shared proposals add tangency refinement (`envelope_refine.h/.cpp`), tail/log-weight utilities, reproducible RNG streams and diagnostics. The full 51-path ledger labels generated files and validation code and distinguishes CMB requirements from these broader changes.

4. **[CMB implementation](#4-cmb-implementation).** The finite-support likelihood uses two natural predictors, `theta=X beta` and `nu=Z gamma`, represented by a stacked design. R/C++ derivatives feed the existing normal-prior sampler; fitted helpers calculate finite-sum means. Limitations include incomplete formula/GLM-family support, omitted offsets in fitted helpers, fractional-count validation and an undeclared validation dependency.

**Data and model specifications:** [six synthetic datasets and fitted models](cases/README.md), including the exact zero-heavy input and paired settings, address the additional dataset request. These comparisons do not isolate non-normality versus correlation as a cause of refinement failure.

**Strongest CMB evidence:** the generated zero-heavy case completed all three measured 20,000-draw runs with refinement (total median **0.594 s**, range 0.593–0.628); all three disabled runs hit the 200,000-proposal cap before the first accepted draw. The same revision, inputs and numerical safeguards were used in both settings. The other five CMB cases completed either way, with modest proposal improvements and largely overlapping timings. Poisson slowed from a median 0.330 s to 0.344 s.

**Artificial stress check:** the deliberately displaced stiff-axis example also needs refinement to complete under these limits. It demonstrates sensitivity to a poor starting tangency, not typical CMB or Gaussian performance. The separate Gaussian envelope comparison forces a path that ordinary Gaussian fits do not use.

All 33 completed measured CMB fits pass the existing integration tolerances; capped runs are failures. The public suite reports 302 passing assertions and 11 OpenCL skips. The paper mapping separates envelope validity, finite-sample efficiency and asymptotics; preservation of the asymptotic guarantees remains unestablished. **This review leaves the published snapshot unchanged; that snapshot enables refinement by default in the shared GLM sampler.**

### Provenance and reading order

Reviewed contribution: **b971596511c90f33de290497b2c7e8feb1a4e580**. Base: **7959a42d32f822a4b65ddfaf94927190ad5db536**. Read-only remote verification on 2026-09-28 PDT found `jbogomolovas2/glmbayes:refs/heads/cmb-public-prep` exactly at the contribution; its merge base is exactly the recorded base. The diff is 51 paths, 2,583 insertions and 234 deletions. All contribution source links below use the full reviewed commit, never a moving branch.

The clean managed review checkout had its HEAD and empty `git status --porcelain` verified; its machine-specific location is omitted from this public copy. A `git archive` of that commit was installed into a separate review library. No untracked exploratory code, private observations, pre-existing object files or pre-installed fork was used. Compilation products are in the separate build copy, not the clean checkout. Only the six generated CMB inputs, synthetic stiff-axis input, and public upstream example data enter the comparisons.

Reading order: start with the overview above; use sections 1–4 as the detailed reference for each request. Section 4 includes the paper mapping, full timing table and reproduction instructions. [Evidence summary](EVIDENCE_SUMMARY.md), [full diff](reviewed.diff), [machine-readable function inventory](r-function-inventory.csv), [results](results/summary.csv), and [publication notes](PUBLICATION_NOTES.md) accompany this note.

## 1. Changed signatures

### Public R functions

Only `EnvelopeSort` changes the formal argument list of an existing exported R function. It appends `logPLSD = NULL` to `(l1, l2, GIndex, G3, cbars, logU, logrt, loglt, logP, LLconst, PLSD, a1, E_draws, lg_prob_factor = NULL, UB2min = NULL)`. Existing positional/named calls remain valid; the new optional vector is reordered with `PLSD` and returned only when supplied. Both regular and independent Normal-Gamma envelope builders call it; the latter need not supply the new field. [Source](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R)

The following existing public functions have **unchanged arguments and defaults**, but changed behavior. The exact old/new formals for each are recorded in the signature appendix below and CSV.

| Functions | Callers and compatibility implications |
|---|---|
| `dNormal`, its internal `plinks` closure | Prior constructors now admit CMB/identity; other prior-family dispatch is unchanged. |
| `glmbfamfunc`, family-local `f1`, `f2`, `f3`, `f4`, `f7` | Used by mode optimization, deviance and information calculations; adds CMB implementations of the existing callback contracts. The branch-local definitions are new implementations, not changes to existing families' callback arguments. |
| `rNormal_reg` | Called through `rglmb`/prior simulation dispatch; admits CMB, skips `glmb.wfit` for it, removes `fit` by assigning NULL. Existing post-fit consumers that require an IRLS fit cannot assume it exists. |
| `rindepNormalGamma_reg` | Called by `rglmb`/`rlmb`; now returns the actual `Envelope`, previously NULL. More retained memory; callers testing `is.null(Envelope)` observe a change. |
| `rglmb`, `rlmb` | Called by users and formula interfaces; add diagnostics and structured proposal-cap errors via `.glmb_run_simulation`. Public controls still do not expose refinement. Existing parallel/serial controls and defaults are unchanged. |
| `glmb`, `lmb` | Formula interfaces copy simulation diagnostics. No complete CMB formula interface is added: the supported CMB entry is augmented matrices with `rglmb`. |
| `residuals.rglmb` | User/S3 residual calls explicitly error for CMB; per-observation deviance residuals remain unsupported. |
| `summary.rglmb` | User/S3 summary calls use CMB finite-sum fitted means; suppress GLM MLE/SE values with NA and return NULL model for CMB. Carries diagnostics. |
| `summary.glmb`, `print.summary.glmb`, `print.summary.rglmb` | Carry/print structured diagnostics when present; otherwise retain the old candidate-count print. Printed output changes without argument changes. |

The internal `.uni_lmb` also copies diagnostics without changing its formals. No other existing R function definition changes in the parsed diff. The parser records the final same-named family closure in each version; the CMB `f1`–`f7` rows therefore denote branch additions, not replacements of earlier branches. The full diff is authoritative for these branches.

### R wrapper → registration → C++ trace, including the defect

The unchanged public formals are `EnvelopeBuild(bStar, A, y, x, mu, P, alpha, wt, family="binomial", link="logit", Gridtype=2L, n=1L, n_envopt=NULL, sortgrid=FALSE, use_opencl=FALSE, verbose=FALSE)` in both revisions. The unchanged manual `.EnvelopeBuild_cpp(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)` has no defaults in either revision. Its unchanged signature is nevertheless affected by the native arity change.

- `rglmb` → `rNormal_reg` → `.rNormalGLM_cpp` → generated `_glmbayes_rNormalGLM_cpp_export` (18 arguments, unchanged) → `rNormalGLM_cpp_export` → `glmbayes::sim::rNormalGLM` → `glmbayes::env::EnvelopeBuild(..., true, 60)`. This path bypasses the defective R envelope wrapper. [R dispatcher](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R), [native wrapper](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/export_wrappers.cpp), [sampler](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp).
- Internal generated `EnvelopeBuild_cpp_export` changes from 16 arguments to the same 16 plus `refine = TRUE, refine_maxit = 60L`. Its `.Call`, generated C trampoline, native registration (16 → 18), export wrapper and C++ implementation all agree. A call to this R function using the old 16 arguments receives its R defaults and is source-compatible. [Generated R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/RcppExports.R), [generated C++](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/RcppExports.cpp).
- **Defect D1:** public `EnvelopeBuild` has unchanged formals; for non-Gaussian families it calls `.EnvelopeBuild_cpp`, whose unchanged manual wrapper passes only 16 arguments directly to the now-18-argument native symbol. C++ default arguments do not apply to an R `.Call` pointer invocation. Native registration reports 18; `length(formals(.EnvelopeBuild_cpp))` is 16. On this build the isolated reproduction **segfaulted** at this call, rather than producing a catchable R error. Do not claim compatibility for this path. Gaussian public `EnvelopeBuild` instead uses the separate independent Normal-Gamma builder. [Manual wrapper](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rcpp_wrappers.R), [public wrapper](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R), [crash log](regression-export.log), [safe arity inspection](defects.log).
- The comparison harness uses the consistent generated export explicitly, not a repaired namespace binding. The reviewed package remains unchanged.

### Existing C++ interfaces

| Function / layer | Old → new arguments/defaults | Callers / implications |
|---|---|---|
| `glmbayes::env::EnvelopeBuild` | Append `bool refine=true, int refine_maxit=60` after `verbose=false`; all earlier parameters/defaults unchanged. | Native GLM sampler explicitly passes true/60; export wrapper forwards settings. Source callers including the header can omit them. C++ ABI changes; recompile all dependents. |
| `EnvelopeBuild_cpp_export` | Append `bool refine=true, int refine_maxit=60` after mandatory `verbose`; preceding 16 arguments unchanged. | Rcpp-generated R/C++ wrappers; generated native entry must receive 18. |
| `_glmbayes_EnvelopeBuild_cpp_export` | 16 SEXP parameters → 18, appended `refineSEXP`, `refine_maxitSEXP`; no C defaults. | Registration count changes accordingly. Manual wrapper defect above. |
| `rNormalGLM_worker` constructor | Add mandatory `uint64_t job_seed_` after `n_`; add `RVector<double> logPLSD_r_`, `bool use_log_sel_` after `PLSD_r_`; add `RMatrix<double> logU_r_` after `logrt_r_`. | Internal parallel/pilot construction sites updated. No R API change. `operator()` arguments stay `(size_t begin,size_t end)`; per-row streams, cap failure and CMB scoring change its body. |
| `rIndepNormalGammaReg_worker` constructor | Add mandatory `uint64_t job_seed_` after `n_`, `RMatrix<double> logU_r_` after `logrt_r_`. | Internal parallel caller updated. `operator()` signature unchanged. |
| `glmbayes::rng::rnorm_ct` | Preserve `(double lgrt,double lglt,double mu,double sigma)`; add overload with mandatory final `double log_mass`. No defaults. | Normal-GLM and independent Normal-Gamma serial/parallel paths supply `logU`; legacy four-argument callers remain possible but cannot encode both extreme-tail endpoints reliably. |
| `runif_safe()` | Signature unchanged; engine changes from OS-seeded `mt19937`/standard distribution to per-draw-seeded `mt19937_64`, explicit nonzero 53-bit uniform. | Shared envelope samplers and truncated draws; results differ from the base's RNG streams. |
| `f2_f3_non_opencl`, `rNormalGLM`, `rNormalGLM_std`, `rNormalGLM_std_parallel`, `run_rcppparallel_pilot` | Signatures/defaults unchanged. | CMB dispatch; log-weight/tail safeguards; proposal cap now aborts without force acceptance; reproducible CPU streams. Pilot checks report cap failure after worker return. |
| `rIndepNormalGammaReg`, `rIndepNormalGammaReg_std`, `rIndepNormalGammaReg_std_parallel` | Signatures/defaults unchanged. | Return envelope; use seeded row streams and five-argument truncated normal. |
| `rNormalReg` | Signature/defaults unchanged. | Adds explicit `Envelope=NULL` to direct Gaussian output; diagnostics can identify direct sampling. |

[Envelope declaration](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/Envelopefuncs.h), [Normal GLM implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp), [independent Normal-Gamma implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rIndepNormalGammaReg.cpp), [RNG declarations](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rng_utils.h), [tail sampler](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rnorm_ct.cpp). There are no new CMB `.Call` registrations: new family kernels are reached by existing dispatch.

Exact old/new C++ formal lists (including worker constructors and unchanged native interfaces) are in [cpp-function-inventory.csv](cpp-function-inventory.csv). Declarations record defaults; definitions do not repeat header defaults.

### Exact R signature appendix

Rows below include all changed existing R definitions, including body-only changes. An argument without a default is required. New definitions are catalogued in section 3 and the CSV.

| Function | Old formals | New formals |
|---|---|---|
| [`EnvelopeBuild_cpp_export`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/RcppExports.R) | `(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)` | `(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose, refine = TRUE, refine_maxit = 60L)` |
| [`glmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/glmb.R) | `(formula, family = binomial, pfamily = dNormal(mu, Sigma, dispersion = 1), n = 1000, data, weights, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, subset, offset, na.action, Gridtype = 2, n_envopt = NULL, start = NULL, etastart, mustart, control = list(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...)` | unchanged (same arguments/defaults) |
| [`lmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/lmb.R) | `(formula, pfamily, n = 1000, data, subset, weights, na.action, method = "qr", model = TRUE, x = TRUE, y = TRUE, qr = TRUE, singular.ok = TRUE, contrasts = NULL, offset, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, ...)` | unchanged (same arguments/defaults) |
| [`.uni_lmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/lmb.R) | `(formula, pfamily, n = 1000, data, subset, weights, na.action, method = "qr", model = TRUE, x = TRUE, y = TRUE, qr = TRUE, singular.ok = TRUE, contrasts = NULL, offset, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, ...)` | unchanged (same arguments/defaults) |
| [`dNormal`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R) | `(mu, Sigma, dispersion = NULL)` | unchanged (same arguments/defaults) |
| [`dNormal::plinks`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R) | `(family)` | unchanged (same arguments/defaults) |
| [`residuals.rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/residuals.glmb.R) | `(object, ysim = NULL, ...)` | unchanged (same arguments/defaults) |
| [`rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rglmb.R) | `(n = 1, y, x, family = gaussian(), pfamily, offset = NULL, weights = 1, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)` | unchanged (same arguments/defaults) |
| [`rlmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rlmb.R) | `(n = 1, y, x, pfamily, offset = rep(0, nobs), weights = NULL, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, progbar = FALSE)` | unchanged (same arguments/defaults) |
| [`rindepNormalGamma_reg`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R) | `(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(), Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, progbar = TRUE)` | unchanged (same arguments/defaults) |
| [`rNormal_reg`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R) | `(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(), Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, progbar = FALSE)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(family, lik_shape = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f1`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f2`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, mu, P, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f3`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, mu, P, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f4`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, alpha = 0, wt = 1, dispersion = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f7`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, mu, P, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`EnvelopeSort`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(l1, l2, GIndex, G3, cbars, logU, logrt, loglt, logP, LLconst, PLSD, a1, E_draws, lg_prob_factor = NULL, UB2min = NULL)` | `(l1, l2, GIndex, G3, cbars, logU, logrt, loglt, logP, LLconst, PLSD, a1, E_draws, lg_prob_factor = NULL, UB2min = NULL, logPLSD = NULL)` |
| [`summary.glmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.glmb.R) | `(object, ...)` | unchanged (same arguments/defaults) |
| [`print.summary.glmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.glmb.R) | `(x, digits = max(3, getOption("digits") - 3), ...)` | unchanged (same arguments/defaults) |
| [`summary.rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.rglmb.R) | `(object, ...)` | unchanged (same arguments/defaults) |
| [`print.summary.rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.rglmb.R) | `(x, digits = max(3, getOption("digits") - 3), ...)` | unchanged (same arguments/defaults) |

## 2. Additional returned items

### Envelopes and sampler returns

| Field | Type / availability | Meaning and consumers |
|---|---|---|
| `Envelope$logPLSD` | Numeric vector, one entry per cell, regular `EnvelopeBuild` sorted, unsorted and sort-memory-fallback paths | Normalized log mixture weights. CPU serial/parallel Normal-GLM selection uses these when a finite log weight exponentiates to zero. Ordinary weights still use the bounded ordinary scan. `EnvelopeSort` returns it only if supplied; no blanket promise for other builders. |
| `Envelope$refinement` | List, regular builder on all return paths; empty when disabled (also potentially incomplete after unusable initial mass) | Detailed refinement state below. Fit diagnostics summarize it; public tests and this harness inspect it. |
| `Envelope` in native `rIndepNormalGammaReg` and R `rindepNormalGamma_reg` | Actual list `Env3`, previously absent/native or NULL/R | Shared diagnostics can now distinguish rejection sampling from direct sampling. It is a different envelope builder: do not assume regular-builder refinement fields exist. |
| `Envelope` in native `rNormalReg` | Explicit NULL | Conjugate Gaussian path, no envelope/refinement. |
| Existing `thetabars`, `cbars`, `LLconst`, `logP`, `PLSD` | Same names; values now reflect refined tangencies | `G3` is synchronized with refined `G4` before `setlogP_C2`, avoiding stale tangencies paired with new gradients. `logP` returned here is the unnormalized box log probability (first internal column), not the full cell envelope mass. |

The full cell log mass used in this review is `logP - LLconst + rowSums(cbars^2)/2`; total mass is its log-sum-exp. These are masses against the standardized normal base, with the likelihood's retained constants. They are comparable between on/off for a given case, not as evidence across different datasets. Acceptance prediction uses the matching integration normalization.

`refinement` fields (from [implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.cpp)):

| Field(s) | Type / meaning |
|---|---|
| `iterations` | Integer evaluated update passes, at most configured maxit (60 here). |
| `residual` | Numeric scalar maximum absolute `m(r)-r` over axes/cells, final state. |
| `converged` | Logical `residual < 1e-6`; efficiency diagnostic, not an MCMC convergence test or accuracy certificate. |
| `accelerated` | Logical eligibility for Newton path; does not mean every cell accepted a Newton update. |
| `log_mass_initial`, `log_mass_final` | Numeric scalars, log-sum-exp of cell envelope masses before/after. |
| `cell_log_mass` | Numeric vector of final per-cell log masses. |
| `step_fraction` | Numeric vector of final per-cell damping factors. |
| `accepted_updates`, `rejected_updates`, `newton_updates` | Numeric vectors of per-cell update counts; accepted Newton steps are a subset of accepted updates. |
| `cell_grid_index` | Numeric matrix encoding integer region codes; builder attaches original `GIndex` to identify the diagnostic cells. |

**Sorted-path caveat:** `EnvelopeSort` reorders main arrays and `logPLSD` by the same permutation. The builder attaches `refinement` afterward; its vectors remain in **original grid order**. Join them to sorted envelope rows using `cell_grid_index` versus `GridIndex`; do not zip them by position. Unsorted and memory-fallback returns carry both fields as well. [Harness check](harness-check.log) verifies normal sorted/unsorted alignment; the allocation-failure fallback is inspected in source, not forced at runtime. The fallback rebuilds an unsorted list and does not retain a `sort_ok` field.

### Fit, summary and error diagnostics

`rglmb`, `rlmb`, `glmb`, `lmb`/`.uni_lmb`, `summary.glmb`, and `summary.rglmb` add/copy `diagnostics`, a list of class `glmb_diagnostics`. The low-level `rNormal_reg` called directly does not itself invoke the new diagnostic wrapper. The two summary printers consume this object; `print.glmb_diagnostics` is a new S3 method. [Definitions](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/fit_diagnostics.R)

| Field(s) | Type / availability / meaning |
|---|---|
| `status`, `method`, `sampling_mode` | Character scalars: completed/refinement_incomplete; direct/envelope_rejection; serial/parallel/not_applicable. Parallel describes the requested CPU mode, not an independently measured device path. |
| `requested_parallel`, `requested_opencl` | Logical requested flags; no proof of actual OpenCL execution. |
| `draws`, `total_proposals`, `proposals_per_draw`, `acceptance_rate`, `max_proposals` | Numeric/integer counts or rates; proposal fields NA if no envelope or `iters` is not a positive integral count per returned draw. |
| `envelope$status`, `$cells`, `$iterations`, `$residual` | Character status and integer/numeric summaries. Status can be converged, incomplete, not_reported, not_applicable; absent refinement values become NA. |
| `build$version`, `$library`, `$R_version`, `$platform` | Character build identity; library is local installation path. |
| `build$fingerprint` | Named character MD5 vector `compiled_code`, `r_code`; NA for unavailable files. Cached once per namespace. |
| `rng_kind` | Character vector from `RNGkind()`. |

Cap errors passing through `.glmb_run_simulation` have classes `glmbayes_sampling_error`, `error`, `condition`, and fields `message` (character), `call`, `parent` (original error), `diagnostics`. Error diagnostics include `status='failed'`, `reason='proposal_limit'`, `draws_returned=0L`, numeric `proposal_limit`, `total_proposals=NA`, character `sampling_mode`, `action`, and `build`. Native low-level harness calls receive the original error, not this R wrapper class. No partial samples are returned. Other errors are rethrown unchanged. Failed comparisons are not passes.

### New CMB return contracts

`cmb()` returns a class-`family` list with `family='cmb'`, `link='identity'`, `valideta` and `validmu` functions (both return TRUE). Despite that class it is not a complete `stats::family` object. `cmb_augment` returns numeric `x` (2N × (p+q)), `y`, `weights`, `offset` (each length 2N); details in section 4. `cmb_fitted` and `cmb_nu` each return draw × N numeric matrices. `glmbfamfunc(cmb())` returns the existing callback names, implemented for CMB. Native `f2_f3_cmb` returns `qf` (one numeric value per parameter column) and `grad` (grid-points × parameters matrix), consumed by `EnvelopeEval` as `NegLL` and `cbars`.

For CMB `summary.rglmb` additionally changes existing values: `fitted.values` is draw × N, `linear.predictors` still contains 2N augmented entries; `family` is the minimal CMB family; `model` is NULL; MLE coefficient/SE comparison columns are NA. These are deliberate limitations, not missing posterior draws.

## 3. New functions and files

### New functions and dependencies

| Group | New names, purpose and visibility |
|---|---|
| CMB public R | `cmb(link='identity')`, `cmb_augment(X,yc,m,Z=NULL)`, `cmb_fitted(object,scale=c('proportion','count'))`, `cmb_nu(object)` are exported. Base R finite sums and matrix operations; sampler requires existing normal prior and standardization. |
| CMB internal R | `.is_cmb`, `.cmb_nobs`, `.cmb_m`; `glmbfamfunc` closures `.cmb_N`, `.cmb_tab`, `.cmb_mom`, `.cmb_dat`, `.cmb_eta`, `.cmb_sat` and its `EYf`; CMB-specific `f1/f2/f3/f4/f7`. Layout checks, sufficient-statistic moments, likelihood/prior derivatives, saturated-likelihood solve. |
| CMB C++ | `f1_cmb`, `f2_cmb`, `f3_cmb`, `f2_f3_cmb`, `f2_cmb_rmat`; private `CMBData`, `cmb_prep_t`, `cmb_moments`, `cmb_N`. Internal namespace kernels, declared in `famfuncs.h`; Rcpp/Armadillo/Rmath/RcppParallel container dependencies. No new exported Rcpp entry points. |
| Shared refinement | `refine_tangency`; private bounds/target/residual helpers and lambdas. Internal C++ API; depends on `EnvelopeEval`, Rmath tails and Armadillo linear solves. |
| Numerical stability | `mills_ratio`, `trunc_norm_mean`, `log_box_prob`, `logaddexp_safe`, `select_cell_log`; private `R_pos`, `logaddexp`, `cell_bounds`; new `rnorm_ct` overload. Shared C++ utilities, not R exports. |
| Shared RNG/cap | `seed_from_R`, `seed_draw`, private `check_proposal_limit` and cap constant. Seed extraction only on R thread; per-row C++ streams. |
| Diagnostics | `.glmb_build_info`, `.glmb_fit_diagnostics`, `.glmb_run_simulation`, cache environment; exported S3 `print.glmb_diagnostics`. Base R build hashes and count checking; no new package dependency. |
| Validation | `quadratic_envelope` test helper; matrix script `lse`, `quad`, `importance` and local callbacks; six new test files and Python subprocess runner. Tests use testthat; importance fallback also requires **mvtnorm**, used via `::` but not declared in DESCRIPTION. |

CMB's finite-support likelihood, two-predictor augmentation, dispatch and mean conversion are its requirements. Refinement, log mixture weights, tail stability, row-seeded RNG, shared proposal failure handling and diagnostics are **broader proposed improvements**. The controlled comparison isolates refinement from those other changes; it does not estimate each safeguard's independent effect or justify a combined upstream merge.

### Complete revision file ledger

Every path in the immutable diff appears below, including documentation, generated interfaces, validation and the deleted machine-specific build file. `A/M/D` means added/modified/deleted; generated files are explicitly labeled. This table and [diff-inventory.tsv](diff-inventory.tsv) cover the full revision rather than only CMB source.

| Status | File | Area / visibility and dependencies |
|---|---|---|
| M | [.Rbuildignore](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/.Rbuildignore) | **validation / packaging** — Exclude local analyses, backups and generated machine Makevars from builds; R build tooling. |
| M | [.gitignore](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/.gitignore) | **validation / packaging** — Ignore local analyses, editor/build files; Git only. |
| A | [CMB_CONTRIBUTION.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/CMB_CONTRIBUTION.md) | **validation** — Contribution scope and historical public results; documentation only. |
| M | [DESCRIPTION](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/DESCRIPTION) | **CMB / packaging** — Development version 0.9.76.9001 and X-CMB-Release v1; no added dependencies. |
| M | [NAMESPACE](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/NAMESPACE) | **CMB / diagnostics** — Generated-style namespace file with CMB exports and appended diagnostic S3 registration. |
| M | [R/RcppExports.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/RcppExports.R) | **shared sampling** — GENERATED Rcpp R wrapper; refinement arguments/defaults; Rcpp native symbol. |
| A | [R/fit_diagnostics.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/fit_diagnostics.R) | **diagnostics** — New fit/build/error diagnostics and S3 printer; base R hashes, RNGkind. |
| M | [R/glmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/glmb.R) | **diagnostics** — Formula GLM diagnostic propagation; rglmb. |
| M | [R/lmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/lmb.R) | **diagnostics** — Formula linear/univariate helper diagnostic propagation; rlmb. |
| M | [R/pfamily.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R) | **CMB** — Normal-prior dispatch, public augmentation/fitted helpers and family; finite sums. |
| M | [R/residuals.glmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/residuals.glmb.R) | **CMB** — Explicit unsupported CMB residual error; S3 methods. |
| M | [R/rglmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rglmb.R) | **diagnostics / shared sampling** — Diagnostic wrapper and reproducibility docs; simfunction. |
| M | [R/rlmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rlmb.R) | **diagnostics / shared sampling** — Diagnostic wrapper; simfunction. |
| M | [R/simfunction.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R) | **CMB / diagnostics** — CMB normal-prior dispatch, IRLS bypass, ING envelope retention. |
| M | [R/simulationpipeline.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | **CMB / numerical stability** — CMB callback closures, sorted log weights; base R, native builders. |
| M | [R/summary.glmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.glmb.R) | **diagnostics** — Summary carries diagnostics; printer uses new S3 method. |
| M | [R/summary.rglmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.rglmb.R) | **CMB / diagnostics** — CMB fitted means/MLE limitations and diagnostics; cmb_fitted. |
| M | [README.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/README.md) | **CMB / validation** — Public usage/scope and validation pointers; documentation. |
| A | [inst/examples/Ex_CMB.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/examples/Ex_CMB.R) | **CMB / validation** — Generated-data public example; augmented matrices and dNormal. |
| A | [inst/validation/CMB_V1.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/CMB_V1.md) | **validation** — Release scope; documentation. |
| A | [inst/validation/CMB_VALIDATION.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/CMB_VALIDATION.md) | **validation** — Public suite and integration protocol; testthat/mvtnorm/Python. |
| A | [inst/validation/RNG_REPRODUCIBILITY.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/RNG_REPRODUCIBILITY.md) | **shared sampling / validation** — CPU RNG contract and limitations; documentation. |
| A | [inst/validation/cmb-matrix.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/cmb-matrix.R) | **CMB / validation** — Six generated cases, quadrature/importance references and fits; mvtnorm fallback. |
| A | [inst/validation/run-cmb-matrix.py](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/run-cmb-matrix.py) | **validation** — Sequential case timeout/CSV runner; Python standard library and Rscript. |
| A | [man/cmb.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/cmb.Rd) | **CMB** — GENERATED roxygen help for four CMB public helpers. |
| A | [man/print.glmb_diagnostics.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/print.glmb_diagnostics.Rd) | **diagnostics** — GENERATED roxygen S3 printer help. |
| M | [man/rglmb.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/rglmb.Rd) | **diagnostics / shared sampling** — GENERATED roxygen return/RNG help. |
| M | [man/rlmb.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/rlmb.Rd) | **diagnostics** — GENERATED roxygen diagnostic return help. |
| M | [man/simfuncs.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/simfuncs.Rd) | **diagnostics** — GENERATED roxygen ING envelope return help. |
| M | [src/EnvelopeBuild.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeBuild.cpp) | **shared sampling / numerical stability** — Refine fixed tangencies, synchronize G3, add log weights/diagnostics; EnvelopeEval/refine/EnvelopeSort. |
| M | [src/EnvelopeEval.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeEval.cpp) | **CMB** — CPU CMB joint likelihood/gradient dispatch; famfuncs_cmb. |
| M | [src/Envelopefuncs.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/Envelopefuncs.h) | **shared sampling** — Append defaulted refinement arguments to builder declaration. |
| D | [src/Makevars](https://github.com/jbogomolovas2/glmbayes/commit/b971596511c90f33de290497b2c7e8feb1a4e580#diff-57dba0e6774ac3caeb9123b3600bddbcfe0e095d1f5b2f2d7115e5e739b0ad2e) | **validation / packaging** — DELETED machine-specific generated build configuration; configure creates environment-appropriate file. |
| M | [src/RcppExports.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/RcppExports.cpp) | **shared sampling** — GENERATED Rcpp trampoline/registration arity 18; Rcpp/export_wrappers. |
| A | [src/envelope_refine.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.cpp) | **shared sampling / numerical stability** — New refinement, stable means/mass and log categorical selection; Rmath, Armadillo, EnvelopeEval. |
| A | [src/envelope_refine.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.h) | **shared sampling / numerical stability** — New internal declarations, defaults and rationale. |
| M | [src/export_wrappers.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/export_wrappers.cpp) | **shared sampling** — Handwritten Rcpp export wrapper forwards refinement arguments. |
| M | [src/famfuncs.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/famfuncs.h) | **CMB** — New internal family kernel declarations; Rcpp/Armadillo/RcppParallel. |
| A | [src/famfuncs_cmb.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/famfuncs_cmb.cpp) | **CMB** — New finite-support normalizer and likelihood/gradient kernels; Rmath/Armadillo. |
| M | [src/rIndepNormalGammaReg.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rIndepNormalGammaReg.cpp) | **shared sampling / diagnostics** — Row streams, interval-mass-aware normal draws, returned envelope. |
| M | [src/rNormalGLM.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp) | **CMB / shared sampling** — CMB accept score, refinement, cap failure, log selection, row streams and logU. |
| M | [src/rNormalReg.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalReg.cpp) | **diagnostics** — Direct Gaussian sampler returns Envelope=NULL. |
| M | [src/rng_utils.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rng_utils.cpp) | **shared sampling** — Job seed from R, per-row stream assignment, explicit uniform conversion. |
| M | [src/rng_utils.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rng_utils.h) | **shared sampling / numerical stability** — Seed function and log_mass overload declarations. |
| M | [src/rnorm_ct.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rnorm_ct.cpp) | **numerical stability** — Endpoint-preserving tail inverse-CDF draw; shared RNG/Rmath. |
| A | [tests/testthat/test-cmb-public.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-cmb-public.R) | **validation** — Finite-sum likelihood/gradient/binomial identity and fitted means. Depends on testthat and built package. |
| A | [tests/testthat/test-envelope-refinement.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-envelope-refinement.R) | **validation** — Stiff unsplit, already-fixed and large-p fallback examples. Depends on testthat and built package. |
| A | [tests/testthat/test-fit-diagnostics.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-fit-diagnostics.R) | **validation** — Fit/build/summary diagnostics and failure conditions. Depends on testthat and built package. |
| A | [tests/testthat/test-rejection-cap.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-rejection-cap.R) | **validation** — No force acceptance or partial samples on cap. Depends on testthat and built package. |
| A | [tests/testthat/test-rng-reproducibility.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-rng-reproducibility.R) | **validation** — CPU row-stream repeatability and serial/parallel equivalence. Depends on testthat and built package. |
| A | [tests/testthat/test-truncated-normal-tails.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-truncated-normal-tails.R) | **validation** — Extreme-tail endpoint and interval behavior. Depends on testthat and built package. |

## 4. CMB implementation

### Likelihood, parameterization, augmentation and dispatch

For count `yc_i` in `{0,…,m_i}`, the implementation uses

`log p(yc_i | theta_i, nu_i) = theta_i*yc_i + nu_i*log choose(m_i,yc_i) - log sum_{j=0}^{m_i} exp(theta_i*j + nu_i*log choose(m_i,j))`.

Here `theta = X beta + offset_theta`, `nu = Z gamma + offset_nu`, both with identity links. `nu` is unrestricted (including negative values) on finite support. At `nu=1`, `theta` is the binomial logit; for other `nu`, applying `plogis(theta)` is not the mean. The normalizer is a stabilized finite log-sum-exp. The log likelihood is concave in the natural parameters because the normalizer is convex; linear predictors preserve this. A proper positive-definite normal prior supplies posterior curvature, including when the likelihood is weakly identified. All `m=1` inputs are rejected because the second sufficient statistic vanishes. [R family/helpers](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R), [native kernels](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/famfuncs_cmb.cpp).

`cmb_augment` constructs `x = rbind(cbind(X,0),cbind(0,Z))`, defaults `Z` to an intercept, sets `y=c(yc/m,0_N)`, `weights=c(m,1_N)`, and zero offsets of length 2N. The final N response/weight entries are placeholders, not observations. Both predictor functionals are rows of the design; they survive the sampler's right-multiplication/rotation of `x`. Treating `nu` as a single fixed coordinate after rotation would be wrong. Explicit offsets can be supplied to `rglmb` after augmentation; the fitted helper defect below affects their reporting.

For `T=(Y,log choose(m,Y))`, the score is the observed statistic minus its finite-sum expectation. The negative posterior gradient is `-x'*(observed-expected)+P*(b-mu)`. R `f7` constructs information from `Var(Y)`, `Var(log choose(m,Y))` and their covariance, with cross-block terms; it is likelihood information without the prior. `f2/f3` drive BFGS mode finding and `optim`'s numerical Hessian drives standardization. Native `f2_f3_cmb` jointly computes likelihood and gradient; native acceptance scoring uses `f2_cmb` or `f2_cmb_rmat`. There is no inner normalizing solve; work scales with the sum of the finite support sizes, number of evaluated parameter points and repeated evaluations. Tables are built per call; a cheap closed-form binomial kernel is not substituted.

`f4` computes total deviance using a saturated likelihood at the draw's own `nu`; interior observations use bracketed Newton/bisection with tolerance 1e-10 and at most 200 iterations, endpoints use limiting saturated log likelihood zero. That stopping ceiling is not accompanied by a general convergence certificate. CMB dispatch is CPU-only here: `EnvelopeEval` adds CMB to `f2_f3_non_opencl`, and serial/parallel Normal-GLM acceptance adds CMB scoring. `rNormalGLM` forces dispersion to 1 for CMB so trial-count weights are not rescaled. [Callbacks](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R), [CPU dispatch](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeEval.cpp), [sampler](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp).

Limitations: no complete GLM family apparatus, CMB formula interface, IRLS fit, ordinary inverse link, generic influence/residual method, or CMB OpenCL certification. Use `cmb_fitted`/`cmb_nu` for zero-offset fitted summaries; generic inherited methods need individual assessment. Validation covers small 2–4 coefficient CMB examples with moderate trial sizes, not high-dimensional or very large-support performance.

### Snapshot defects and qualification of claims

1. **D1 — native wrapper mismatch/crash:** described in section 1; public non-Gaussian `EnvelopeBuild` is broken. Regression tests miss it because refinement tests use the generated export. No fix was applied.
2. **D2 — offsets omitted by fitted helpers:** `cmb_fitted` and `cmb_nu` use only `x %*% coefficients`, ignoring fitted `offset`/`offset2`. The supplied synthetic check with theta offset 1 and nu offset .5 returns count mean 2 instead of 2.710088, and nu 1 instead of 1.5. Affects CMB summaries using those helpers. [Reproducer](check-defects.R), [result](defects.log).
3. **D3 — insufficient count validation:** `cmb_augment` accepts fractional counts/trial sizes (e.g. 1.2 of 4.4); downstream R/C++ code rounds implied counts and m, silently changing the stated data. It also does not comprehensively reject zero/nonfinite trial sizes at the public boundary. The finite-support model requires integer m≥1 and integer counts in range. Generated validation inputs satisfy this.
4. **D4 — validation dependency metadata:** `cmb-matrix.R` calls `mvtnorm` for difficult references, but DESCRIPTION does not declare it. It was present here (version recorded in reference metadata); a reproduction needs it explicitly.
5. **Documentation qualifications:** comments saying failed refinement retains the original tangency are too broad: a rejected update retains the last accepted point; the algorithm does not roll the whole cell back to initialization. Nested diagnostic order differs from sorted envelope order. Header claims that refinement can “never” affect correctness assume exact valid evaluation; floating-point/inverse-CDF limitations remain. The unchanged A08 vignette transcription omits square roots in the omega formula and presents Theorem 3 as an inequality; the paper's displayed Theorem 3 is a limit. This is inherited documentation, not a new contribution defect. The C++ initialization matches the paper's square-root formula.

### Original initialization and refinement algorithm

For standardized axis i with curvature `a_i=diag(A)_i`, initialization is

`omega_i = (sqrt(2)-exp(-1.20491-.7321*sqrt(.5+a_i)))/sqrt(1+a_i)`.

Tangencies start at `bStar_i + {-omega_i,0,omega_i}`, and cutpoints at `bStar_i ± omega_i/2`. Gridtype 1 selects one/three points via a curvature threshold, 2 uses `EnvelopeOpt`'s cost criterion with `n_envopt`, 3 always splits, 4 leaves axes unsplit. Products form cells; unsplit axes have region code 4. Refinement leaves these cutpoints, grid size and standardization fixed. [Builder](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeBuild.cpp)

Writing `F(r)` for the negative residual log likelihood including the shifted prior quadratic and `c=grad F(r)`, the restricted proposal is `N(-c,I)` within the fixed cell. Its unnormalized log envelope mass is

`log W_A(r) = -F(r) + r'c + .5*c'c + log Pr[N(-c,I) in A]`.

The target equation is `r=m(r)`, the restricted proposal mean. The code's differentiation gives `grad log W_A = Hessian(log likelihood)*(m-r)`. Under its concavity assumptions, the raw `m-r` direction is a descent direction. Singular curvature means stationarity need not imply this particular fixed point; neither this calculation nor the implementation establishes global uniqueness/convergence.

Each pass computes target means using continued-fraction Mills ratios in positive tails and reflection on the left. For at most 16 parameters and `p*p*cells <= 2,000,000`, it forms forward finite-difference Jacobians of `r-m(r)` (step `1e-5*max(1,abs(r_axis))`) and solves per-cell Newton systems. Failed/nonfinite solves retain the fixed-point proposal. Larger problems use damped fixed-point steps without Jacobian storage. Candidates are projected into their original cell; nonfinite coordinates stay at their current value.

A candidate must have finite evaluated likelihood, gradient and box mass. Its cell log mass must not increase, except a narrowly defined Newton root correction: old residual <1e-3, new residual less than half, and mass increase within `16*epsilon*(1+abs(F)+abs(linear)+abs(quadratic)+abs(box))`. Rejection halves that cell's damping; accepted Newton steps can double it back toward 1. Last accepted values are retained. An unusable initial total mass (NaN/+Inf) returns immediately without moves or full diagnostics; it does not guarantee a usable envelope afterward.

Stopping: maximum residual <1e-6, no representable move, or 60 evaluated passes in this builder. The helper has a `min_gain` default 1e-3, but **the builder passes 0**, disabling the nonnegative-gain early-exit rule (`gain < min_gain`). There is no universal iteration bound for convergence. Cost includes initial evaluation, up to p extra full-grid gradient evaluations per Newton pass, one trial evaluation per pass, per-cell dense solves and bounded Jacobian storage; backtracking can consume further passes. Construction and total time must therefore accompany acceptance rates. [Implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.cpp)

### Paper assumptions mapped to code

The supplied full text of Nygren & Nygren (2006), *Likelihood Subgradient Densities*, JASA 101:1144–1156, DOI [10.1198/016214506000000357](https://doi.org/10.1198/016214506000000357), was read directly, including sections 2, 3.1–3.3, Remarks 7–19 and the relevant appendix. PDF pages 3–8 correspond to journal pages 1145–1150 (the file includes a cover page). The formulas on journal pages 1147–1150 were checked visually.

| Paper result / scope | Implementation mapping and limits |
|---|---|
| Definition 2, Theorem 1, section 2 (pp.1145–1146): valid global subgradient bound, finite tilted-prior MGF | Proper normal prior and globally concave natural CMB likelihood supply the mathematical setting; `EnvelopeEval` computes tangents. Changing a tangency does not by itself invalidate that exact-arithmetic bound. Runtime numerical evaluation is a separate concern. |
| Claim 2, Remarks 5–6, section 3.1 (p.1147): combine restricted envelopes over a partition | `Set_Grid`, `setlogP_C2`, `PLSD/logPLSD`, restricted-normal draws and acceptance scoring implement the mixture. Refinement keeps cells fixed and changes their tangent masses. The text immediately after Remark 6 states the tangency/restricted-mean condition used here. It does not give this safeguarded Newton algorithm or its convergence rate. |
| Example 2, Remarks 7–8 (p.1148): diagonal normal restriction, inverse transform, tail accuracy | Rectangular standardized cells yield product probabilities and means. Mills-ratio and `logU` safeguards address numerical risk, but finite-precision tail/categorical correctness is not proved for all inputs. A 53-bit uniform cannot make arbitrarily tiny cell probabilities practically selectable merely because logs remain finite. |
| Section 3.2, Theorem 2 (p.1148) | The normal-data three-interval construction has limit `2/sqrt(pi)` as data precision grows. The calibrated omega is the starting construction, not a general finite-sample CMB bound. |
| Definition 3, Remarks 11–15, section 3.3 (p.1149) | Existing Cholesky/rotation and prior-quadratic redistribution achieve the intended standard form when positive-definiteness and smoothness assumptions hold. CMB's two row-block predictors preserve its likelihood under that transform. |
| Theorem 3 (p.1149), Remarks 10/18 (pp.1148/1150) | The displayed normal-data theorem gives limit `(2/sqrt(pi))^k` for fixed dimension and specified diagonal precision scaling. General log-concave statements rely on applicable Bayesian normal approximation. We have not proved the needed regularity/uniform envelope control for growing dimension, weak identification, changing priors, extreme CMB regimes, adaptive grid selection or numerical refinement. |
| Remarks 9/17 and 19 (pp.1148–1150) | Off-mode placement can preserve validity while harming efficiency; setup costs can favor smaller blocks. The deliberately stiff test illustrates placement sensitivity, not a counterexample to the paper's correctly centered normal construction. |

**Three separate questions:** (i) envelope validity under the paper's assumptions; (ii) finite-sample acceptance and total cost on specific inputs; (iii) asymptotic behavior along a specified sequence of models. The measurements below address (ii), and reference checks support sampled moments/rates on those inputs. They do not prove (i) numerically for all inputs or establish (iii). Reducing valid cell masses at a fixed partition is an efficiency argument conditional on valid computation, not a new proof or a universal guarantee.

### Controlled comparison protocol

[compare.R](compare.R) calls the pinned generated builder with explicit `refine=TRUE/FALSE`, then the same pinned CPU `rNormalGLM_std_cpp_export`. All numerical safeguards remain active. No source patch or default change is used. One intercept case was cross-checked against public `rglmb`: identical proposal counts and coefficient differences at most 4.89e-12 (relative `all.equal` tolerance 1e-12). Matrix operation ordering accounts for floating-point differences; see [check](harness-check.log).

Ten cases × two settings × (one warm-up + three measured repetitions) = **80 fits**, all targeting **20,000 accepted draws**. Seed 719 generates the unchanged six public synthetic datasets; sampling seeds are 9100 for warm-up, 9101–9103 for measured runs, paired between settings. Order alternates on/off by repetition. Each fit runs alone in a fresh R process with serial native sampling, OpenCL false, RcppParallel threads 1 and BLAS/OpenMP thread environment limits 1. Gridtype 2 and `n_envopt=1000` stay fixed except the existing stiff test's Gridtype 4. Parent subprocess timeout is 300 seconds per entire fit process; native cap stays 200,000 proposals per accepted draw.

The six CMB cases and priors come verbatim from [cmb-matrix.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/cmb-matrix.R). The stiff test uses `curvature=c(10,1000), center=1`, zero predictor and a single unsplit cell from [the existing regression example](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-envelope-refinement.R), extended to sampling. This intentional off-mode construction is not a representative correctly standardized Gaussian fit.

Public upstream examples: binomial uses MASS Menarche's centered age design and informative prior exactly as [Ex_rglmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/examples/Ex_rglmb.R); Poisson uses its Dobson counts/outcome/treatment and `Prior_Setup`. Gaussian uses the Dobson plant-weight data and design from [Ex_03_Dobsonlinreg.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/demo/Ex_03_Dobsonlinreg.R), with `Prior_Setup` normal prior and a fixed plug-in dispersion `var(resid(lm(weight~group)))`. These public examples predate the contribution. Gaussian's ordinary dNormal route is direct: the paired **envelope-only Gaussian diagnostic** deliberately forces the shared builder/sampler, and is not a speed claim about ordinary Gaussian fitting or the example's unknown-dispersion prior.

Construction time covers BFGS/standardization and envelope building (stiff test supplies standard quantities). Sampling time covers the native sampling call. Total time includes construction, sampling, transformation and checkpoint/result serialization inside the harness; it excludes package/data setup and R process launch. Consequently total is not exactly the sum of the two phase timings. No reference integration or off-setting residual evaluation is timed. Disabled residuals are evaluated afterward with `refine=TRUE, refine_maxit=0`; assertions verify unchanged tangencies/weights. Raw runs retain unavailable residuals; `runs-with-residuals.csv` adds these measured initial residuals.

### Results and tradeoffs

Each table value below is **median [minimum, maximum] over the three measured runs**; warm-ups are preserved in raw results but excluded. Every completed run returned 20,000 draws. No run timed out. Failed/capped times are time-to-failure, not equivalent completed workloads. Native failures expose no partial draws/count vector; their messages identify cap at draw 1, establishing zero accepts and 200,000 proposals in these particular failures. Raw inaccessible counts are retained as NA; [failure-counts.csv](results/failure-counts.csv) separately records the exact zero accepts / 200,000 proposals established by these first-draw cap messages. [Machine-readable counts, acceptance, passes, cells and all timing ranges](results/summary.csv)

| Case | Refinement | Complete/cap/timeout | Construction s | Sampling s | Total s | Proposals | log mass | Residual |
|---|---|---|---|---|---|---|---|---|
| intercepts | TRUE | 3/0/0 | 0.031 [0.03, 0.032] | 0.306 [0.306, 0.309] | 0.353 [0.351, 0.356] | 25513 [25395, 25540] | -128.55 [-128.55, -128.55] | 4.63989e-07 [4.63989e-07, 4.63989e-07] |
| intercepts | FALSE | 3/0/0 | 0.029 [0.029, 0.029] | 0.308 [0.307, 0.308] | 0.352 [0.351, 0.353] | 25768 [25710, 25835] | -128.541 [-128.541, -128.541] | 0.00614947 [0.00614947, 0.00614947] |
| theta_slope | TRUE | 3/0/0 | 0.035 [0.035, 0.036] | 0.349 [0.34, 0.354] | 0.406 [0.396, 0.414] | 28838 [28702, 28848] | -128.253 [-128.253, -128.253] | 2.31869e-10 [2.31869e-10, 2.31869e-10] |
| theta_slope | FALSE | 3/0/0 | 0.033 [0.032, 0.034] | 0.35 [0.343, 0.351] | 0.405 [0.397, 0.407] | 29252 [29230, 29439] | -128.238 [-128.238, -128.238] | 0.011722 [0.011722, 0.011722] |
| nu_slope | TRUE | 3/0/0 | 0.034 [0.033, 0.035] | 0.338 [0.338, 0.343] | 0.393 [0.391, 0.398] | 28877 [28737, 29037] | -133.582 [-133.582, -133.582] | 3.46712e-10 [3.46712e-10, 3.46712e-10] |
| nu_slope | FALSE | 3/0/0 | 0.033 [0.031, 0.033] | 0.342 [0.34, 0.342] | 0.396 [0.392, 0.396] | 29270 [29250, 29582] | -133.566 [-133.566, -133.566] | 0.00992152 [0.00992152, 0.00992152] |
| strong_prior | TRUE | 3/0/0 | 0.04 [0.04, 0.042] | 0.369 [0.369, 0.374] | 0.436 [0.436, 0.444] | 31951 [31789, 32129] | -124.475 [-124.475, -124.475] | 1.2255e-09 [1.2255e-09, 1.2255e-09] |
| strong_prior | FALSE | 3/0/0 | 0.034 [0.033, 0.034] | 0.374 [0.374, 0.383] | 0.436 [0.436, 0.443] | 32539 [32266, 32659] | -124.457 [-124.457, -124.457] | 0.118226 [0.118226, 0.118226] |
| weak_prior | TRUE | 3/0/0 | 0.042 [0.041, 0.045] | 0.387 [0.38, 0.454] | 0.459 [0.451, 0.524] | 33208 [33086, 33319] | -136.903 [-136.903, -136.903] | 7.71018e-08 [7.71018e-08, 7.71018e-08] |
| weak_prior | FALSE | 3/0/0 | 0.036 [0.035, 0.037] | 0.39 [0.386, 0.392] | 0.455 [0.449, 0.457] | 34130 [34070, 34163] | -136.871 [-136.871, -136.871] | 0.00641066 [0.00641066, 0.00641066] |
| zero_heavy | TRUE | 3/0/0 | 0.129 [0.129, 0.13] | 0.437 [0.436, 0.468] | 0.594 [0.593, 0.628] | 45569 [45388, 45644] | -49.7923 [-49.7923, -49.7923] | 4.45723e-08 [4.45723e-08, 4.45723e-08] |
| zero_heavy | FALSE | 0/3/0 | 0.036 [0.036, 0.037] | 1.217 [1.215, 1.222] | 1.275 [1.275, 1.28] | NA | 751131 [751131, 751131] | 929.771 [929.771, 929.771] |
| stiff_axis | TRUE | 3/0/0 | 0.003 [0.002, 0.003] | 4.642 [4.626, 4.819] | 4.671 [4.651, 4.844] | 2.10052e+06 [2.09023e+06, 2.10521e+06] | -0.918939 [-0.918939, -0.918939] | 1.04592e-08 [1.04592e-08, 1.04592e-08] |
| stiff_axis | FALSE | 0/3/0 | 0.003 [0.002, 0.003] | 0.431 [0.43, 0.438] | 0.437 [0.435, 0.443] | NA | 500554 [500554, 500554] | 1001 [1001, 1001] |
| gaussian | TRUE | 3/0/0 | 0.003 [0.003, 0.003] | 0.213 [0.211, 0.214] | 0.234 [0.233, 0.235] | 25441 [25366, 25513] | -25.2209 [-25.2209, -25.2209] | 3.25469e-08 [3.25469e-08, 3.25469e-08] |
| gaussian | FALSE | 3/0/0 | 0.003 [0.002, 0.003] | 0.214 [0.212, 0.214] | 0.234 [0.234, 0.235] | 25441 [25367, 25514] | -25.2209 [-25.2209, -25.2209] | 0.000113381 [0.000113381, 0.000113381] |
| binomial | TRUE | 3/0/0 | 0.003 [0.003, 0.003] | 0.226 [0.224, 0.229] | 0.246 [0.246, 0.251] | 25473 [25416, 25539] | -62.6083 [-62.6083, -62.6083] | 5.44887e-09 [5.44887e-09, 5.44887e-09] |
| binomial | FALSE | 3/0/0 | 0.003 [0.002, 0.003] | 0.23 [0.229, 0.234] | 0.251 [0.251, 0.256] | 25491 [25408, 25539] | -62.6077 [-62.6077, -62.6077] | 0.00268295 [0.00268295, 0.00268295] |
| poisson | TRUE | 3/0/0 | 0.017 [0.017, 0.019] | 0.29 [0.289, 0.296] | 0.344 [0.343, 0.352] | 36843 [36396, 36855] | -36.0385 [-36.0385, -36.0385] | 8.8669e-11 [8.8669e-11, 8.8669e-11] |
| poisson | FALSE | 3/0/0 | 0.003 [0.003, 0.004] | 0.29 [0.287, 0.291] | 0.33 [0.327, 0.331] | 36941 [36770, 37087] | -36.0326 [-36.0326, -36.0326] | 0.0186753 [0.0186753, 0.0186753] |


Refinement is **critical under these limits for zero_heavy**: all three enabled runs complete in median 0.594 s, with median 45,569 proposals and acceptance about 0.439. All three disabled runs cap at the first draw. Log mass falls from 751130.940 to -49.7923; residual falls from 929.771 to 4.46e-8, with 31 passes on 81 cells. This is practical necessity for this fixed setup, not for all CMB data or all possible grids.

**Artificial stress check, separate from the CMB evidence:** the deliberately displaced stiff-axis test completes with refinement in median 4.671 s; disabled runs cap at draw 1. Refinement is necessary for completion in this constructed setup under the chosen limits. Even with refinement, acceptance is only about 0.00952 (~105 proposals/draw): a single unsplit cell remains inefficient. One Newton pass addresses the deliberately bad tangency, not the coarse partition. This is not evidence that typical Gaussian fits require refinement.

For the other five CMB cases, both settings finish and reference-check successfully. Cell mass improvements imply only about 0.9–3.2% reductions in expected proposals; total times mostly overlap. Weak-prior enabled time is 0.459 [0.451,0.524] s versus 0.455 [0.449,0.457] s disabled, so improved acceptance does not establish a speed benefit. Gaussian is effectively unchanged. Menarche has a small measured timing improvement but almost unchanged proposal counts; three subsecond repetitions are too little evidence for a general speed assertion. Poisson is slower with refinement: 0.344 [0.343,0.352] s versus 0.330 [0.327,0.331] s, despite slightly fewer proposals. Construction grows from median .003 to .017 s; sampling medians are both .290 s.

### Regression, references, environment and reproduction

The existing public suite was run **once**, unchanged, on the pinned installed package: **302 passed assertions, 0 failed, 0 error tests, 0 test warnings, 11 skipped tests, 68 tests total**. All eleven skips were OpenCL-unavailable cases (see exact reasons in the retained log/CSV). [Raw log](regression.log), [CSV](regression.csv), [R result object](regression.rds). The initial CSV export failed on its list column after testing was complete; a separate export from the saved RDS removed that column, without rerunning the suite. A later isolated wrapper reproduction segfault is retained separately from the suite result; a passing suite does not erase D1.

[references.R](references.R) reuses only the original integration portion of `cmb-matrix.R`. It does not rerun its serial/parallel fitting loop. For the first five cases, independently scaled tensor Gauss-Hermite rules agree within 1e-5 in log mass, mean and covariance. Zero-heavy does not meet that threshold even at 52 nodes; the unchanged fallback uses two seeded Student-t integrations, 32 × 32,768 evaluations each. They agree within 1.147 combined SE, with ESS 579,786 and 579,058 (>100,000). The reference status is resolved for all six. Reference generation had a separate 30-minute ceiling per case, outside all fit timings.

All **33 completed measured CMB fits** pass the original six-SE mean, covariance and predicted proposal-rate tests, including reference uncertainty. Maximum scores across them are 3.3614 (mean), 2.6049 (covariance), and 2.0312 (proposal rate), below 6. Refinement convergence is reported separately: disabled valid sampling is not required to satisfy the refinement fixed-point tolerance. Enabled CMB runs all converge. The three capped zero-heavy measured runs have no moment/rate pass. Integration shares the compiled likelihood, so it checks sampling/integration rather than independently validating that kernel; the public finite-sum/derivative tests address complementary likelihood checks. [Per-run checks](results/reference-checks.csv), [reference logs/status and RDS files](references/status.csv).

Environment: macOS Ventura 13.0.1, arm64; R 4.5.2, glmbayes 0.9.76.9001; Apple clang 14.0.3; Accelerate BLAS and R LAPACK 3.12.1; Rcpp 1.1.2, RcppParallel 6.2.1, opencltools/nmathopencl 0.8.3. Full session/build fingerprints and dependency versions are in [environment.txt](environment.txt) and install/reference logs. Locale startup warnings were retained (C locale). CPU results do not certify OpenCL. Subsecond timing differences should be read in light of clock granularity, process/environment overhead and three repetitions; there is no new speed guarantee.

Follow the [reproduction instructions](README.md#inspect-or-reproduce). The publication runner checks the pinned revision and dependencies, then archives the immutable source into a new run directory, builds a private library, runs the public suite, 80 serial fits, six reference integrations, reference checks and summaries. Saved evidence is not overwritten. The measured fitting and integration algorithms are unchanged; publication-only harness changes are listed in [publication notes](PUBLICATION_NOTES.md). `reproduce-wrapper-defect.R` is an optional separate-process **crash reproducer**, deliberately excluded from the normal runner. [Checksums](SHA256SUMS) identify delivery files; full raw draws/envelopes are retained for independent reanalysis.

**Subsequent upstream is separate:** read-only remote HEAD was `dc2d905d2244350c76068609abe031bb01421132`, with eight commits listed by the GitHub comparison against the base. They cover OpenCL-core delegation, an independent Normal-Gamma prior guard, output gating/method return values, progress utilities, DIC/residual/summary alignment with glmbayesCore, importing/reexporting Prior_Setup, a merge, and further migration updates. Overlap exists in DESCRIPTION/NAMESPACE, `glmb`, `simfunction`, residuals/summary, envelope and independent Normal-Gamma code. They are not in the compared binary and need separate integration review. [Exact commit/file metadata](upstream-compare.json) and [upstream ledger](UPSTREAM_CHANGES.md) preserve that distinction.

This review leaves the published contribution and its defaults unchanged. The reviewed contribution itself enables refinement in the shared GLM sampler. The measurement/review stage did not post a GitHub comment, email, push, PR or publication, or restructure the contribution. The subsequent documentation publication is described in [publication notes](PUBLICATION_NOTES.md).
