import subprocess,os,pathlib,csv
root=pathlib.Path(__file__).resolve().parent;os.chdir(root)
env=os.environ.copy();env.update(R_LIBS=str(root/'library'),OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
(root/'references').mkdir(exist_ok=True)
rows=[]
for case in ['intercepts','theta_slope','nu_slope','strong_prior','weak_prior','zero_heavy']:
 with (root/'references'/f'{case}.log').open('w') as f:
  try:r=subprocess.run(['Rscript','--vanilla','references.R',case],env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800);status='resolved' if r.returncode==0 else 'error'
  except subprocess.TimeoutExpired:status='unresolved_timeout'
 rows.append(dict(case=case,status=status));print(case,status,flush=True)
 with (root/'references'/'status.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=['case','status']);w.writeheader();w.writerows(rows)
