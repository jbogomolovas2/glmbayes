# Publication notes

The review concerns contribution `b971596511c90f33de290497b2c7e8feb1a4e580`, based on `7959a42d32f822a4b65ddfaf94927190ad5db536`. The saved measurements were made on September 28, 2026. This branch adds documentation and evidence only.

## Saved evidence

The results, draws, envelopes, references, source diff and function inventories are preserved from the original review. The six CSV datasets were extracted from the saved inputs and checked against the generator, all 48 CMB comparison inputs and six references. Counts and trial sizes are exact; covariates have 17 significant digits, with full precision retained in the RDS files.

Private study data, the supplied paper, installed libraries, build copies and the unsent reply draft are excluded. Local machine paths and the hostname were removed from public metadata. In `regression.rds`, 19 test source-directory references were redacted; outcomes and timings are unchanged. Already-public upstream paths in `reviewed.diff` are retained.

## Reproduction scripts

The original fitting and integration scripts are unchanged. Publication adjustments make the runner check dependencies and use a fresh output directory, stop on regression failures, and write the signature inventory beside its script. The verifier checks file hashes, links and saved result totals. [Reproduction instructions](README.md#inspect-or-reproduce)

Publication checks reproduced the reference-check summaries from saved draws; the full sampling experiment was not rerun. Subsequent wording edits shorten the presentation without changing the measurements or package code. [SHA256SUMS](SHA256SUMS) covers the current public files.
