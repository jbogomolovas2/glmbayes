# CPU envelope RNG reproducibility

Call `set.seed()` before a fit. For the same model, seed, package build, R RNG
kind, and numerical environment, CPU envelope draws repeat. Normal GLM
(including CMB) and independent Normal-Gamma samplers use the same stream
assignment in serial and parallel execution. Worker count and scheduling do
not determine which stream a draw receives.

The implementation obtains one 64-bit job key from two R uniform draws on the
main thread, under the calling Rcpp RNG scope. It uses `std::seed_seq` to seed
an `mt19937_64` engine from that key and the 64-bit output-row index. Engine
storage is thread-local, but its identity is reset for each output row.
Mixture selection, truncated proposals, dispersion proposals, and rejection
uniforms all use that row's stream. No worker accesses R's RNG.

Pilot and calibration passes replay their assigned rows. They do not consume
additional R uniforms or perturb the final draws. Rejection work for one row
does not move another row's stream. Uniforms use an explicit 53-bit conversion
and exclude both endpoints, rather than relying on a library-specific
`uniform_real_distribution` mapping.

Exactly two uniforms advance R's RNG per low-level sampler invocation. Other
R operations in a larger analysis may consume their own random numbers.
Consecutive calls without resetting the seed therefore receive different keys.
The main thread's R RNG state is the same after serial or parallel execution
of a given invocation, even when worker timing changes.

Increasing the requested number of draws preserves the existing prefix only
when the envelope and other numerical inputs remain the same. Pin `n_envopt`
when making that comparison. The new streams intentionally differ from the
previous OS-seeded implementation. No promise is made of bitwise equality
across R/library/compiler/platform changes, or of OpenCL reproducibility.
Independent stream initialization is not a mathematical proof of non-overlap
among all possible pseudorandom sequences.

`test-rng-reproducibility.R` checks repeated calls, different seeds, consecutive
calls, draw prefixes, subsequent R state, Mersenne-Twister and L'Ecuyer-CMRG
keys, and 1/2/4-thread runs. It tests exact serial/parallel equality for CMB and
independent Normal-Gamma coefficients, dispersion, and proposal counts.
Posterior-accuracy and rejection-cap tests are retained separately so exact
repeatability cannot conceal a consistently wrong sampler.
