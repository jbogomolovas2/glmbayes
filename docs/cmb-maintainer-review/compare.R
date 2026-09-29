# Validation-only harness: pinned generated exports; no package mutation.
suppressPackageStartupMessages(library(glmbayes))
args <- commandArgs(TRUE); case<-args[1]; refine<-as.logical(args[2]); repid<-as.integer(args[3]); out<-args[4]
RcppParallel::setThreadOptions(numThreads=1)
n<-20000L; seed<-9100L+repid; set.seed(719)
sourcefile<-file.path('source','inst','validation','cmb-matrix.R')
cmbcases<-c('intercepts','theta_slope','nu_slope','strong_prior','weak_prior','zero_heavy')
if(case %in% cmbcases) {
 name<-case
 lines<-readLines(sourcefile)
 eval(parse(text=lines[grep('^set.seed',lines)[1]:(grep('^ff<-',lines)[1]-1)]))
 family<-cmb()
} else if(case=='stiff_axis') {
 p<-2; a<-list(y=0,x=matrix(0,1,p),offset=0,weights=1);mu<-matrix(0,p,1);P<-diag(c(10,1000));family<-gaussian()
} else if(case=='binomial') {
 data(menarche,package='MASS'); a<-list(y=menarche$Menarche/menarche$Total,x=cbind(1,menarche$Age-13),offset=rep(0,nrow(menarche)),weights=menarche$Total)
 mu<-matrix(c(0,log(9)/3),2,1);Sigma<-diag(c((log(9))^2,(3*mu[2]/2)^2));P<-solve(Sigma);family<-binomial()
} else if(case=='poisson') {
 counts<-c(18,17,15,20,10,20,25,13,12);outcome<-gl(3,1,9);treatment<-gl(3,3)
 ps<-Prior_Setup(counts~outcome+treatment,family=poisson())
 a<-list(y=ps$y,x=as.matrix(ps$x),offset=rep(0,9),weights=rep(1,9));mu<-matrix(ps$mu);P<-solve(ps$Sigma);family<-poisson()
} else if(case=='gaussian') {
 ctl<-c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14);trt<-c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
 weight<-c(ctl,trt);group<-gl(2,10,20);ps<-Prior_Setup(weight~group,family=gaussian())
 a<-list(y=weight,x=model.matrix(~group),offset=rep(0,20),weights=rep(1/var(resid(lm(weight~group))),20));mu<-matrix(ps$mu);P<-solve(ps$Sigma);family<-gaussian()
} else stop('unknown case')
p<-ncol(a$x);ff<-glmbayes:::glmbfamfunc(family)
row<-data.frame(case=case,refine=refine,repetition=repid,seed=seed,n_target=n,status='started',construction_seconds=NA_real_,sampling_seconds=NA_real_,total_seconds=NA_real_,accepted=NA_real_,proposals=NA_real_,acceptance=NA_real_,log_mass=NA_real_,residual=NA_real_,iterations=NA_real_,converged=NA,cells=NA_integer_,message='')
write_row<-function() write.csv(row,out,row.names=FALSE)
write_row();t0<-proc.time()[3];sampling_start<-NA_real_
tryCatch({
 if(case=='stiff_axis') {
  s<-list(bstar2=rep(1,p),A=P,x2=a$x,mu2=mu,P2=P,L2Inv=diag(p),L3Inv=diag(p));alpha<-a$offset
 } else {
  alpha<-as.vector(a$offset+a$x%*%mu);zero<-matrix(0,p,1)
  o<-optim(rep(0,p),ff$f2,ff$f3,y=a$y,x=a$x,mu=zero,P=P,alpha=alpha,wt=a$weights,method='BFGS',hessian=TRUE)
  stopifnot(o$convergence==0)
  s<-glmbayes:::glmb_Standardize_Model_cpp_export(a$y,a$x,P,matrix(o$par,p,1),o$hessian)
 }
 E<-glmbayes:::EnvelopeBuild_cpp_export(s$bstar2,s$A,a$y,s$x2,s$mu2,s$P2,alpha,a$weights,family$family,family$link,if(case=='stiff_axis')4L else 2L,n,1000L,TRUE,FALSE,FALSE,refine,60L)
 row$construction_seconds<-proc.time()[3]-t0
 lw<-as.vector(E$logP)-as.vector(E$LLconst)+.5*rowSums(E$cbars^2);lse<-function(v){b<-max(v);b+log(sum(exp(v-b)))}
 row$log_mass<-lse(lw);row$cells<-length(E$PLSD)
 if(length(E$refinement)) {row$residual<-E$refinement$residual;row$iterations<-E$refinement$iterations;row$converged<-E$refinement$converged}
 # Persist envelope before sampling, including failed/timed-out fits.
 saveRDS(list(envelope=E,standard=s,input=a,mu=mu,P=P,alpha=alpha),sub('.csv$','-envelope.rds',out));write_row()
 set.seed(seed);sampling_start<-proc.time()[3]
 fit<-glmbayes:::rNormalGLM_std_cpp_export(n,a$y,s$x2,s$mu2,s$P2,alpha,a$weights,ff$f2,E,family$family,family$link,0L,FALSE)
 row$sampling_seconds<-proc.time()[3]-sampling_start
 b<-sweep(as.matrix(fit$out)%*%t(s$L2Inv%*%s$L3Inv),2,as.vector(mu),'+')
 row$status<-'completed';row$accepted<-nrow(b);row$proposals<-sum(fit$draws);row$acceptance<-nrow(b)/row$proposals
 saveRDS(list(b=b,draws=fit$draws,envelope=E),sub('.csv$','.rds',out))
},error=function(e){row$status<<-if(grepl('reached 200000 proposals',conditionMessage(e)))'proposal_cap' else 'error';row$message<<-conditionMessage(e);if(is.finite(sampling_start))row$sampling_seconds<<-proc.time()[3]-sampling_start})
row$total_seconds<-proc.time()[3]-t0;write_row();print(row)
