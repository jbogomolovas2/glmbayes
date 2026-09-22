# Run from repository root: Rscript inst/validation/cmb-matrix.R CASE OUTPUT_DIR
# Select the installed build with R_LIBS. Each case is run in a separate process.
suppressPackageStartupMessages(library(glmbayes))
args <- commandArgs(TRUE); stopifnot(length(args)==2)
name <- args[1]; out <- args[2]
stopifnot(name %in% c('intercepts','theta_slope','nu_slope','strong_prior',
                     'weak_prior','zero_heavy'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(719)
v <- seq(-1,1,length.out=60); one <- matrix(1,60,1)
X <- cbind(1,v); Z <- X; sdprior <- 5
m <- rep(10,60); theta <- -.3+.5*v; nu <- 1+.3*v
if(name=='intercepts') { X<-Z<-one;theta<-rep(-.3,60);nu<-rep(1,60) }
if(name=='theta_slope') Z<-one
if(name=='nu_slope') X<-one
if(name=='strong_prior') sdprior<-.5
if(name=='weak_prior') sdprior<-20
if(name=='zero_heavy') {theta<--3+.3*v;nu<-rep(.5,60)}
y <- vapply(seq_along(m),function(i){j<-0:m[i];w<-theta[i]*j+nu[i]*lchoose(m[i],j);sample(j,1,prob=exp(w-max(w)))},numeric(1))
a <- cmb_augment(X,y,m,Z=Z);p<-ncol(a$x);mu<-matrix(0,p,1);P<-diag(p)/sdprior^2
ff<-glmbayes:::glmbfamfunc(cmb())
o<-optim(rep(0,p),fn=ff$f2,gr=ff$f3,y=a$y,x=a$x,mu=mu,P=P,alpha=a$offset,wt=a$weights,method='BFGS',hessian=TRUE)
stopifnot(o$convergence==0,min(eigen(o$hessian,symmetric=TRUE)$values)>0)
s<-glmbayes:::glmb_Standardize_Model_cpp_export(a$y,a$x,P,matrix(o$par,p,1),o$hessian)
inv<-s$L2Inv%*%s$L3Inv;L0<-t(chol(solve(diag(p)+s$A)));center<-s$bstar2
lse<-function(v){b<-max(v);b+log(sum(exp(v-b)))}
quad<-function(nodes,scale){
 J<-matrix(0,nodes,nodes);J[cbind(1:(nodes-1),2:nodes)]<-sqrt(1:(nodes-1));e<-eigen(J+t(J),symmetric=TRUE)
 logw<-2*log(abs(e$vectors[1,]));L<-scale*L0;logZ<--Inf;mn<-rep(0,p);cross<-matrix(0,p,p)
 for(start in seq(0,nodes^p-1,by=5000)){
  id<-seq.int(start,min(nodes^p-1,start+4999));idx<-sapply(0:(p-1),function(k)(id%/%nodes^k)%%nodes+1L)
  u<-matrix(e$values[idx],nrow(idx),p);z<-sweep(u%*%t(L),2,center,'+')
  ev<-glmbayes:::EnvelopeEval_cpp_export(t(z),a$y,s$x2,s$mu2,s$P2,a$offset,a$weights,'cmb','identity',FALSE,FALSE)
  terms<-rowSums(matrix(logw[idx],nrow(idx),p))-ev$NegLL-.5*rowSums(z^2)+sum(log(diag(L)))+.5*rowSums(u^2)
  nxt<-lse(c(logZ,lse(terms)));old<-exp(logZ-nxt);w<-exp(terms-nxt);b<-z%*%t(inv)
  mn<-old*mn+colSums(b*w);cross<-old*cross+crossprod(b,b*w);logZ<-nxt
 }
 list(logZ=logZ,mean=mn,cov=cross-tcrossprod(mn),nodes=nodes,scale=scale)
}
for(nodes in c(26L,36L,52L)){
 q<-quad(nodes,1);q2<-quad(nodes,1.5)
 err<-max(abs(q$logZ-q2$logZ),max(abs(q$mean-q2$mean)),max(abs(q$cov-q2$cov)))
 cat(name,'quadrature nodes',nodes,'agreement',err,'\n');flush.console()
 if(err<1e-5) break
 if(nodes==26L) {
  # Pilot moments only select a better integration proposal; they are not
  # accepted as a reference until two independently scaled rules agree.
  center<-as.vector(solve(inv,q$mean))
  V<-solve(inv)%*%q$cov%*%t(solve(inv));L0<-t(chol(V))
 }
}
# For strongly non-Gaussian posteriors, check an independently sampled,
# heavy-tailed integration proposal and retain its uncertainty explicitly.
importance <- function(seed,pilot) {
 set.seed(seed); B<-32L; N<-32768L; df<-8
 centerz<-as.vector(solve(inv,pilot$mean))
 V<-solve(inv)%*%pilot$cov%*%t(solve(inv)); V<-V*(df-2)/df
 batches<-lapply(seq_len(B),function(k){
  z<-mvtnorm::rmvt(N,sigma=V,df=df,delta=centerz,type='shifted')
  ev<-glmbayes:::EnvelopeEval_cpp_export(t(z),a$y,s$x2,s$mu2,s$P2,a$offset,a$weights,'cmb','identity',FALSE,FALSE)
  lw<--ev$NegLL-.5*rowSums(z*z)-p/2*log(2*pi)-mvtnorm::dmvt(z,delta=centerz,sigma=V,df=df,log=TRUE)
  zsum<-lse(lw);w<-exp(lw-zsum);b<-z%*%t(inv);mn<-colSums(b*w);cross<-crossprod(b,b*w)
  list(logZ=zsum-log(N),mean=mn,cross=cross,cov=cross-tcrossprod(mn),logsum2=lse(2*lw))
 })
 logs<-vapply(batches,function(b)b$logZ,numeric(1));logZ<-lse(logs)-log(B)
 w<-exp(logs-logZ)/B;mn<-Reduce('+',Map(function(b,w)b$mean*w,batches,w));cross<-Reduce('+',Map(function(b,w)b$cross*w,batches,w))
 list(logZ=logZ,mean=mn,cov=cross-tcrossprod(mn),method='Student-t importance',
  logZ_se=sd(exp(logs-logZ))/sqrt(B),
  mean_se=apply(sapply(batches,function(b)b$mean),1,sd)/sqrt(B),
  cov_se=matrix(apply(sapply(batches,function(b)as.vector(b$cov)),1,sd)/sqrt(B),p,p),
  ess=exp(2*(logZ+log(B*N))-lse(vapply(batches,function(b)b$logsum2,numeric(1)))))
}
if(err>=1e-5) {
 pilot<-q;q<-importance(381,pilot);q2<-importance(992,pilot)
 agreement<-max(abs(q$logZ-q2$logZ)/sqrt(q$logZ_se^2+q2$logZ_se^2),
   abs(q$mean-q2$mean)/sqrt(q$mean_se^2+q2$mean_se^2),
   abs(q$cov-q2$cov)/sqrt(q$cov_se^2+q2$cov_se^2))
 cat(name,'independent importance agreement (SE)',agreement,'ESS',q$ess,q2$ess,'\n')
 stopifnot(agreement<6,min(q$ess,q2$ess)>1e5)
 err<-NA_real_
} else {q$logZ_se<-0;q$mean_se<-rep(0,p);q$cov_se<-matrix(0,p,p)}
saveRDS(list(case=name,a=a,sdprior=sdprior,reference=q,reference_check=q2,reference_error=err),file.path(out,paste0(name,'-reference.rds')))

rows<-list();n<-20000L
for(parallel in c(FALSE,TRUE)){
 t0<-proc.time()[3]
 f<-rglmb(n=n,y=a$y,x=a$x,offset=a$offset,weights=a$weights,family=cmb(),pfamily=dNormal(mu=mu,Sigma=diag(p)*sdprior^2),n_envopt=1000L,use_parallel=parallel,verbose=FALSE)
 b<-as.matrix(f$coefficients);E<-f$Envelope;S<-cov(b);centered<-sweep(b,2,q$mean,'-');se<-matrix(0,p,p)
 for(i in 1:p)for(j in 1:p)se[i,j]<-sd(centered[,i]*centered[,j])/sqrt(n)
 mz<-max(abs(colMeans(b)-q$mean)/sqrt(diag(S)/n+q$mean_se^2));cz<-max(abs(S-q$cov)/sqrt(se^2+q$cov_se^2))
 lw<-as.vector(E$logP)-as.vector(E$LLconst)+.5*rowSums(E$cbars^2);alpha<-exp(lse(lw)-q$logZ)
 pz<-abs(mean(f$iters)-alpha)/sqrt(alpha*(alpha-1)/n+(alpha*q$logZ_se)^2)
 passed<-all(is.finite(b))&&mz<6&&cz<6&&pz<6&&isTRUE(E$refinement$converged)
 row<-data.frame(case=name,parallel=parallel,p=p,zeros=sum(y==0),prior_sd=sdprior,n=n,cells=length(E$PLSD),passes=E$refinement$iterations,residual=E$refinement$residual,seconds=proc.time()[3]-t0,proposals=mean(f$iters),predicted=alpha,mean_z=mz,cov_z=cz,proposal_z=pz,passed=passed)
 print(row);flush.console();rows[[length(rows)+1]]<-row
 write.csv(do.call(rbind,rows),file.path(out,paste0(name,'.csv')),row.names=FALSE)
 saveRDS(list(row=row,mean=colMeans(b),cov=S,reference=q,envelope=E),file.path(out,paste0(name,'-',parallel,'.rds')))
}
stopifnot(all(vapply(rows,function(r)r$passed,logical(1))))
