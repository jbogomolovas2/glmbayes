import csv,statistics,pathlib
root=pathlib.Path(__file__).resolve().parent
r=list(csv.DictReader((root/'results/runs-with-residuals.csv').open()))
metrics=['construction_seconds','sampling_seconds','total_seconds','accepted','proposals','acceptance','log_mass','residual','iterations','cells']
rows=[]
for case in dict.fromkeys(x['case'] for x in r):
 for ref in ('TRUE','FALSE'):
  a=[x for x in r if x['case']==case and x['refine']==ref and x['repetition']!='0']
  d=dict(case=case,refine=ref,measured=len(a),completed=sum(x['status']=='completed' for x in a),proposal_cap=sum(x['status']=='proposal_cap' for x in a),timeout=sum(x['status']=='timeout' for x in a))
  for m in metrics:
   vals=[float(x[m]) for x in a if x.get(m) not in (None,'','NA')]
   for suffix,f in [('median',statistics.median),('min',min),('max',max)]:d[m+'_'+suffix]=f(vals) if vals else 'NA'
  rows.append(d)
with (root/'results/summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
with (root/'results/summary.md').open('w') as f:
 f.write('| Case | Refinement | Complete/cap/timeout | Construction s | Sampling s | Total s | Proposals | log mass | Residual |\n|---|---|---|---|---|---|---|---|---|\n')
 def fmt(d,m):
  v=d[m+'_median'];return 'NA' if v=='NA' else f'{v:.6g} [{d[m+"_min"]:.6g}, {d[m+"_max"]:.6g}]'
 for d in rows:f.write('| '+ ' | '.join([d['case'],d['refine'],f"{d['completed']}/{d['proposal_cap']}/{d['timeout']}"]+[fmt(d,m) for m in ['construction_seconds','sampling_seconds','total_seconds','proposals','log_mass','residual']])+' |\n')

raw=list(csv.DictReader((root/'results/runs.csv').open()))
with (root/'results/failure-counts.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['case','refine','repetition','status','accepted_before_failure','proposals_before_failure','draws_returned','evidence']);w.writeheader()
 for a in raw:
  if a['status']=='proposal_cap' and 'zero acceptances for draw 1 of 20000' in a['message']:
   w.writerow(dict(case=a['case'],refine=a['refine'],repetition=a['repetition'],status=a['status'],accepted_before_failure=0,proposals_before_failure=200000,draws_returned=0,evidence=a['message']))
