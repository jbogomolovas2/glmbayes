import re,subprocess,csv,pathlib
root=pathlib.Path(__file__).resolve().parent
base='7959a42d32f822a4b65ddfaf94927190ad5db536';pin='b971596511c90f33de290497b2c7e8feb1a4e580'
spec={'src/Envelopefuncs.h':['EnvelopeBuild'],'src/export_wrappers.cpp':['EnvelopeBuild_cpp_export'],'src/RcppExports.cpp':['_glmbayes_EnvelopeBuild_cpp_export'],'src/rNormalGLM.cpp':['rNormalGLM_worker','run_rcppparallel_pilot','rNormalGLM_std_parallel','rNormalGLM_std','rNormalGLM'],'src/rIndepNormalGammaReg.cpp':['rIndepNormalGammaReg_worker','rIndepNormalGammaReg_std','rIndepNormalGammaReg_std_parallel','rIndepNormalGammaReg'],'src/EnvelopeEval.cpp':['f2_f3_non_opencl'],'src/rNormalReg.cpp':['rNormalReg'],'src/rng_utils.cpp':['runif_safe'],'src/rng_utils.h':['rnorm_ct']}
def extract(rev,path,name):
 s=subprocess.check_output(['git','show',rev+':'+path],text=True);s=re.sub(r'/\*.*?\*/|//[^\n]*','',s,flags=re.S)
 vals=[]
 for m in re.finditer(r'\b'+re.escape(name)+r'\s*\(',s):
  i=m.end();depth=1
  while i<len(s) and depth:
   if s[i]=='(':depth+=1
   if s[i]==')':depth-=1
   i+=1
  tail=s[i:].lstrip()
  if not tail.startswith(('{',':',';','BEGIN_RCPP')):continue
  # Header declarations only; source forward prototypes are excluded except trampoline.
  if tail.startswith(';') and not path.endswith('.h'):continue
  before=s[s.rfind('\n',0,m.start())+1:m.start()].strip()
  if '=' in before or before in ('return',''): # constructors are bare names
   if name not in ('rNormalGLM_worker','rIndepNormalGammaReg_worker'):continue
  args=' '.join(s[m.end():i-1].split());value=name+'('+args+')'
  if value not in vals:vals.append(value)
 return '\n'.join(vals)
rows=[]
for path,names in spec.items():
 for name in names:
  old=extract(base,path,name);new=extract(pin,path,name)
  assert old and new,(path,name,old,new)
  rows.append(dict(file=path,symbol=name,old=old,new=new,kind='unchanged signature' if old==new else 'changed signature'))
with (root/'cpp-function-inventory.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
