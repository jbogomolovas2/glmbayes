# Publication notes

This documentation-only branch publishes the review prepared on September 28, 2026, for contribution `b971596511c90f33de290497b2c7e8feb1a4e580` against base `7959a42d32f822a4b65ddfaf94927190ad5db536`. The contribution remains the measured code revision. Later upstream changes are recorded separately and were not integrated.

## Public contents

The original delivery checksum manifest supplied the explicit file allowlist. The public copy retains the review, inventories, exact source diff, reproduction scripts, environment metadata, tables, logs, numerical integration references, posterior draws and envelopes. The unsent GitHub reply draft is excluded. A landing page, these publication notes, and six CSV inputs with explicit model specifications are added. The integer responses and trial counts are exact extracts; covariates are exported with 17 significant digits (base R CSV parsing differs by at most one machine epsilon). Full binary-precision inputs remain in the unchanged RDS files. The generating inputs were verified against the pinned generator and all 48 CMB comparison inputs and six references, addressing the additional dataset/model request in the current issue comment. The duplicate ZIP, build/source copies, installed library, supplied paper, and private/exploratory analysis directories are excluded.

The saved comparisons remain the original 80 runs. No model, sampler setting, timing, count, posterior sample, envelope or integration reference was changed during publication. No new performance claim is based on publication checks. The source diff is byte-for-byte unchanged, including the already-public upstream build paths in the deleted `src/Makevars` file.

## Redactions and documentation edits

- The local checkout path in the review and revision metadata is replaced by a descriptive placeholder.
- The host name and local review-library prefix in the environment and installation log are redacted.
- In `regression.rds`, 19 shared test source environments have their working-directory prefix replaced with `[review-root]/`. Regression outcome and timing tables were checked identical before and after redaction. Test messages, source text and line references are retained.
- Public links replace references to the excluded reply draft. Statements about the review stage are distinguished from this later publication stage. Checksums were regenerated for this public copy.

All saved numeric CSV tables, posterior/envelope/reference RDS files, original source diff and function inventories remain byte-identical to the original delivery. The regression result object is the sole RDS file with a metadata-only redaction.

## Reproduction harness adjustments

The measured fitting and integration scripts (`compare.R`, `references.R`, their Python runners and `validate-results.R`) are unchanged. Publication changes are limited to:

- `reproduction.sh`: resolves the repository from its own location, checks the pinned commit and R dependencies (including `mvtnorm`), archives from the repository root, and creates a new output directory. An existing destination is refused. It checks comparison/reference outcomes after the run; documented capped comparisons remain recorded failures.
- `regression.R`: retains its saved results and also returns an error when the public suite has failed assertions or error tests.
- `inventory.R`: writes its R signature inventory beside the script rather than to the former local review-directory path.
- `verify-package.py`: checks the complete public manifest and relative links, the original result totals, regression totals, and recorded summary statistics. It reads saved evidence without rerunning fits.

The README documents the fresh-clone commands. Syntax, dependency preflight, refusal to overwrite existing output, and integrity/reanalysis checks were performed for publication. The full 80-fit/reference-generation pipeline was not rerun for this documentation publication. The original build/test/comparison logs remain the evidence for the reported experiments.

## Outstanding code work

The wrapper argument-count crash, omitted CMB offsets, count validation and undeclared validation dependency remain documented defects of the reviewed snapshot. Public controls/defaults, upstream integration and contribution structure remain subjects for the maintainer discussion. This publication does not itself open a PR or post an issue comment.
