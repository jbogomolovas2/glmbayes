# What the comparisons show

**Refinement made the difference in the zero-heavy CMB example.** All three measured runs with refinement produced 20,000 draws, taking a median 0.594 seconds. Without it, every run reached 200,000 proposals before accepting a draw.

The other five CMB examples worked with either setting. Refinement reduced the expected number of proposals by about 1–3%, but gave little overall speed improvement. Poisson was slightly slower: 0.344 versus 0.330 seconds. The separate stiff-axis example also needed refinement, but it deliberately starts from a poor tangency and is a stress test.

Whether refinement preserves the paper's asymptotic guarantees remains an open question. These comparisons also do not establish whether non-normality or correlation causes the difficult CMB case.

## How the comparison was run

Ten cases were tested with refinement enabled and disabled, using matching inputs, seeds and numerical safeguards. Each setting had one warm-up and three measured runs, requesting 20,000 draws on a single CPU thread. Each draw had a 200,000-proposal limit; each fit had a five-minute timeout.

Of 80 runs, 72 finished and eight reached the proposal limit, including warm-ups. None timed out. **All 33 completed measured CMB runs passed the numerical reference checks.** Runs that reached the cap were counted as failures.

The package regression suite passed 302 assertions. Eleven OpenCL tests were skipped, so these results cover CPU sampling only. See the [full comparison details](MAINTAINER_REVIEW.md#controlled-comparison-protocol) and [timing table](results/summary.md).

## What needs fixing

- The public envelope wrapper passes 16 arguments where the native function expects 18; the reproduction crashed R.
- CMB fitted-value helpers omit offsets.
- Fractional counts and trial sizes are accepted and later rounded.
- The validation scripts use `mvtnorm`, which is missing from the package dependencies.

These defects remain in the reviewed code. Refinement is currently enabled by default in the shared sampler.

## Data and supporting material

[Datasets and model specifications](cases/README.md) · [Four-part review](MAINTAINER_REVIEW.md) · [Individual runs](results/runs.csv) · [Reference checks](results/reference-checks.csv) · [Reproduction instructions](README.md#inspect-or-reproduce)

The results refer to contribution `b971596…` against base `7959a42…`. [Environment details](environment.txt), [later upstream changes](UPSTREAM_CHANGES.md), and the full revision identifiers are recorded in the review.
