# Public CMB validation

Run from the repository root with the installed fork selected by `R_LIBS`:

```sh
Rscript -e 'library(glmbayes); testthat::test_dir("tests/testthat", reporter="summary")'
python3 inst/validation/run-cmb-matrix.py /tmp/cmb-validation
```

Six generated-data cases cover intercepts, theta slopes, nu slopes, strong and
weak priors, and zero-heavy responses. Fixed seeds, 2--4 coefficients and normal
prior SDs 0.5--20 are used. Each case fits 20,000 draws in serial and parallel
mode: 12 fits and 240,000 posterior draws in total. The wrapper applies a
per-case timeout and exits nonzero on failure.

The reference integrates the posterior using tensor Gauss-Hermite rules at
two scales, increasing resolution and recentering if needed. References must
agree within 1e-5 for log mass, means and covariance. If unresolved, two seeded
Student-t importance integrations are used (32 batches of 32,768 evaluations
each), requiring effective sample sizes over 100,000 and agreement within six
combined SE. Reference uncertainty is included in the six-SE posterior and
proposal-rate comparisons. A failed reference is not a passing sampler test.

These checks share the compiled likelihood with the sampler; they independently
check integration and sampling, not the likelihood formula. The separate
`test-cmb-public.R` checks direct finite-support PMF calculations, derivatives,
the binomial special case and fitted-mean conversion on synthetic inputs.

Shared tests cover bounded failures, stiff refinement, extreme-tail boundaries
and exact CPU serial/parallel RNG reproducibility. See `RNG_REPRODUCIBILITY.md`.
The 200,000-proposal cap bounds proposals per Normal-GLM draw, not elapsed time
or construction. Refinement convergence concerns efficiency, not proof of
posterior accuracy. OpenCL is not certified by CPU-only tests.

Private emu observations and their data-specific tests/references are excluded
from this public branch. Historical local test totals therefore differ from
the public suite. Public results are recorded in `CMB_CONTRIBUTION.md`.
