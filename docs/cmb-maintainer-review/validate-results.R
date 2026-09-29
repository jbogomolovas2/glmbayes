library(glmbayes)
r<-read.csv('results/runs.csv');rows<-list()
for(i in seq_len(nrow(r))) {
 z<-r[i,];key<-paste(z$case,as.integer(z$refine),z$repetition,sep='-');ep<-paste0('results/',key,'-envelope.rds')
 if(!file.exists(ep))next
 e<-readRDS(ep);E<-e$envelope
 if(!z$refine) {
  s<-e$standard;a<-e$input
  EE<-glmbayes:::EnvelopeBuild_cpp_export(s$bstar2,s$A,a$y,s$x2,s$mu2,s$P2,e$alpha,a$weights,if(z$case %in% c('stiff_axis','gaussian'))'gaussian' else if(z$case %in% c('binomial','poisson'))z$case else 'cmb',if(z$case=='binomial')'logit' else if(z$case=='poisson')'log' else 'identity',if(z$case=='stiff_axis')4L else 2L,20000L,1000L,TRUE,FALSE,FALSE,TRUE,0L)
  stopifnot(isTRUE(all.equal(E$thetabars,EE$thetabars)),isTRUE(all.equal(E$PLSD,EE$PLSD)))
  r$residual[i]<-EE$refinement$residual;r$iterations[i]<-0L;r$converged[i]<-NA
 }
 if(z$status!='completed' || z$repetition==0 || !z$case %in% c('intercepts','theta_slope','nu_slope','strong_prior','weak_prior','zero_heavy'))next
 refpath<-paste0('references/',z$case,'-reference.rds')
 if(!file.exists(refpath)){rows[[length(rows)+1]]<-data.frame(case=z$case,refine=z$refine,repetition=z$repetition,status='unresolved_reference',mean_z=NA,cov_z=NA,proposal_z=NA,passed=FALSE);next}
 q<-readRDS(refpath)$reference;f<-readRDS(paste0('results/',key,'.rds'));b<-f$b;n<-nrow(b);p<-ncol(b);S<-cov(b)
 centered<-sweep(b,2,q$mean,'-');se<-matrix(0,p,p)
 for(k in 1:p)for(j in 1:p)se[k,j]<-sd(centered[,k]*centered[,j])/sqrt(n)
 mz<-max(abs(colMeans(b)-q$mean)/sqrt(diag(S)/n+q$mean_se^2));cz<-max(abs(S-q$cov)/sqrt(se^2+q$cov_se^2))
 alpha<-exp(z$log_mass-q$logZ);pz<-abs(mean(f$draws)-alpha)/sqrt(alpha*(alpha-1)/n+(alpha*q$logZ_se)^2)
 rows[[length(rows)+1]]<-data.frame(case=z$case,refine=z$refine,repetition=z$repetition,status='checked',mean_z=mz,cov_z=cz,proposal_z=pz,passed=all(is.finite(b))&&all(is.finite(c(mz,cz,pz)))&&mz<6&&cz<6&&pz<6)
}
write.csv(r,'results/runs-with-residuals.csv',row.names=FALSE)
write.csv(do.call(rbind,rows),'results/reference-checks.csv',row.names=FALSE)
