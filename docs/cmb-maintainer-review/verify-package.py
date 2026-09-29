"""Verify the published saved evidence without rebuilding or rerunning fits."""
from pathlib import Path
from urllib.parse import unquote, urlsplit
import csv
import hashlib
import math
import re
import statistics

root = Path(__file__).resolve().parent
pin = 'b971596511c90f33de290497b2c7e8feb1a4e580'
def table(path):
    with (root / path).open(newline='') as stream:
        return list(csv.DictReader(stream))

manifest = {}
for row in (root / 'SHA256SUMS').read_text().splitlines():
    sha, path = row.split('  ', 1)
    target = (root / path).resolve()
    assert root in target.parents and path not in manifest, path
    assert target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == sha, path
    manifest[path] = sha
actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name != 'SHA256SUMS'}
assert actual == set(manifest), (actual - set(manifest), set(manifest) - actual)
assert 'GITHUB_REPLY_DRAFT.md' not in actual
assert not any(Path(p).suffix in {'.zip', '.o', '.so', '.dylib', '.pdf'} for p in actual)

runs = table('results/runs.csv')
assert len(runs) == 80
assert len({(x['case'], x['refine'], x['repetition']) for x in runs}) == 80
cases = {'intercepts','theta_slope','nu_slope','strong_prior','weak_prior','zero_heavy','stiff_axis','gaussian','binomial','poisson'}
assert {x['case'] for x in runs} == cases
for x in runs:
    assert int(x['n_target']) == 20000
    capped = x['case'] in {'zero_heavy', 'stiff_axis'} and x['refine'] == 'FALSE'
    assert x['status'] == ('proposal_cap' if capped else 'completed'), x
    if not capped:
        assert int(float(x['accepted'])) == 20000
assert sum(x['status'] == 'completed' for x in runs) == 72
assert sum(x['status'] == 'completed' and x['repetition'] != '0' for x in runs) == 54

checks = table('results/reference-checks.csv')
refs = table('references/status.csv')
assert len(refs) == 6 and all(x['status'] == 'resolved' for x in refs)
cmb = {x['case'] for x in refs}
expected = {(x['case'],x['refine'],x['repetition']) for x in runs
            if x['case'] in cmb and x['status'] == 'completed' and x['repetition'] != '0'}
assert len(checks) == 33 and {(x['case'],x['refine'],x['repetition']) for x in checks} == expected
assert all(x['passed'] == 'TRUE' and all(math.isfinite(float(x[m])) and float(x[m]) < 6
           for m in ['mean_z','cov_z','proposal_z']) for x in checks)
for case in sorted(cmb):
    data = table('cases/' + case.replace('_','-') + '-data.csv')
    assert len(data) == 60 and [int(x['observation']) for x in data] == list(range(1,61))
    assert all(int(x['trials']) == 10 and 0 <= int(x['successes']) <= 10 for x in data)
    assert all(math.isclose(float(x['v']), -1 + 2*i/59, abs_tol=1e-14) for i,x in enumerate(data))
    if case == 'zero_heavy':
        assert sum(int(x['successes']) == 0 for x in data) == 48
        assert sum(int(x['successes']) for x in data) == 18
regression = table('regression.csv')
assert len(regression) == 68
assert sum(int(x['passed']) for x in regression) == 302
assert sum(x['skipped'] == 'TRUE' for x in regression) == 11
assert all(int(x['failed']) == 0 and x['error'] == 'FALSE' and int(x['warning']) == 0 for x in regression)
summary = table('results/summary.csv')
raw = table('results/runs-with-residuals.csv')
assert len(summary) == 20
for row in summary:
    selected = [x for x in raw if x['case'] == row['case'] and x['refine'] == row['refine'] and x['repetition'] != '0']
    assert len(selected) == int(row['measured']) == 3
    for status in ['completed','proposal_cap','timeout']:
        assert sum(x['status'] == status for x in selected) == int(row[status])
    for metric in ['construction_seconds','sampling_seconds','total_seconds','accepted','proposals','acceptance','log_mass','residual','iterations','cells']:
        values = [float(x[metric]) for x in selected if x.get(metric) not in (None, '', 'NA')]
        for suffix, operation in [('median',statistics.median),('min',min),('max',max)]:
            value = row[metric + '_' + suffix]
            assert (value == 'NA' if not values else math.isclose(float(value), operation(values), rel_tol=1e-12, abs_tol=1e-12)), (row['case'],metric,suffix)
paths = [x.split('\t')[1] for x in (root / 'diff-inventory.tsv').read_text().splitlines()]
review = (root / 'MAINTAINER_REVIEW.md').read_text()
assert len(paths) == 51 and all(path in review for path in paths)
assert re.findall(r'^## (\d)\. (.+)$', review, re.M) == [
    ('1','Changed signatures'),('2','Additional returned items'),
    ('3','New functions and files'),('4','CMB implementation')]

def anchors(text):
    found = set()
    counts = {}
    for heading in re.findall(r'^#{1,6}\s+(.+)$', text, re.M):
        heading = re.sub(r'\[([^\]]+)\]\([^)]+\)', r'\1', heading)
        key = re.sub(r'[^\w\- ]', '', heading.lower()).replace(' ', '-')
        count = counts.get(key, 0)
        found.add(key + ('-' + str(count) if count else ''))
        counts[key] = count + 1
    return found

for file in root.rglob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)', file.read_text()):
        uri = urlsplit(link)
        if uri.scheme:
            if '/blob/' in uri.path:
                assert '/blob/' + pin + '/' in uri.path, (file,link)
            continue
        target = (file.parent / unquote(uri.path)).resolve() if uri.path else file
        assert target.exists(), (file,link)
        if uri.fragment and target.suffix == '.md':
            assert unquote(uri.fragment) in anchors(target.read_text()), (file,link)
print(f'Verified {len(manifest)} files/checksums and all relative links; 51-path inventory, 80 runs, 33 passing completed CMB checks, six references, 302 regression assertions and 11 skips; summary statistics match saved runs.')
