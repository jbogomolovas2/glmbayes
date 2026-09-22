#!/usr/bin/env python3
"""Run from repository root, with R_LIBS selecting the build to validate."""
import argparse
import csv
from pathlib import Path
import subprocess

CASES = ['intercepts', 'theta_slope', 'nu_slope', 'strong_prior', 'weak_prior',
         'zero_heavy']
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=Path)
parser.add_argument('--timeout', type=float, default=240)
parser.add_argument('--cases', nargs='+', choices=CASES, default=CASES)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
failed = []
for case in args.cases:
    # Stream logs so reference progress is inspectable while a case runs.
    with (args.output / f'{case}.log').open('w') as log:
        try:
            result = subprocess.run(
                ['Rscript', '--vanilla', 'inst/validation/cmb-matrix.R',
                 case, str(args.output.resolve())],
                stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
            if result.returncode:
                failed.append(case)
            print(case, 'PASS' if result.returncode == 0 else 'FAIL', flush=True)
        except subprocess.TimeoutExpired:
            failed.append(case)
            log.write(f'\nTIMEOUT after {args.timeout} seconds\n')
            print(case, 'TIMEOUT', flush=True)
rows = []
for case in args.cases:
    path = args.output / f'{case}.csv'
    if path.exists() and case not in failed:
        with path.open() as src:
            rows.extend(csv.DictReader(src))
if rows:
    with (args.output / 'summary.csv').open('w', newline='') as dst:
        writer = csv.DictWriter(dst, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
raise SystemExit(1 if failed else 0)
