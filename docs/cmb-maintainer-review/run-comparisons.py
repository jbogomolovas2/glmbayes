#!/usr/bin/env python3
"""Run from this directory after installing the pinned source into library/."""
import os, subprocess, csv, time, pathlib, json
root=pathlib.Path(__file__).resolve().parent
os.chdir(root)
env=os.environ.copy();env.update(R_LIBS=str(root/'library'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',RCPP_PARALLEL_NUM_THREADS='1')
(root/'results').mkdir(exist_ok=True)
rows=[]
for case in ['intercepts','theta_slope','nu_slope','strong_prior','weak_prior','zero_heavy','stiff_axis','gaussian','binomial','poisson']:
 for rep in range(4): # one complete warm-up per setting, then three measured fits
  for refine in (True,False) if rep%2==0 else (False,True):
   key=f'{case}-{int(refine)}-{rep}';dest=root/'results'/f'{key}.csv'
   start=time.monotonic()
   with (root/'results'/f'{key}.log').open('w') as log:
    try:
     proc=subprocess.run(['Rscript','--vanilla','compare.R',case,str(refine).upper(),str(rep),str(dest)],env=env,stdout=log,stderr=subprocess.STDOUT,timeout=300)
     status='process_error' if proc.returncode else None
    except subprocess.TimeoutExpired:status='timeout'
   row=next(csv.DictReader(dest.open())) if dest.exists() else dict(case=case,refine=refine,repetition=rep)
   if status:row.update(status=status,total_seconds=time.monotonic()-start,message='External 300-second per-fit timeout' if status=='timeout' else 'See process log')
   rows.append(row)
   keys=list(dict.fromkeys(k for r in rows for k in r))
   with (root/'results'/'runs.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
   print(key,row.get('status'),row.get('total_seconds'),flush=True)
