library(glmbayes)
e<-readRDS('results/intercepts-1-1-envelope.rds');h<-readRDS('results/intercepts-1-1.rds');a<-e$input
set.seed(9101)
f<-rglmb(n=20000,y=a$y,x=a$x,offset=a$offset,weights=a$weights,family=cmb(),pfamily=dNormal(mu=e$mu,Sigma=solve(e$P)),n_envopt=1000,use_parallel=FALSE,use_opencl=FALSE,verbose=FALSE)
cat('Maximum absolute coefficient difference:',max(abs(f$coefficients-h$b)), '\n');print(all.equal(unname(f$coefficients),unname(h$b)));print(all.equal(as.vector(f$iters),as.vector(h$draws)));stopifnot(isTRUE(all.equal(unname(f$coefficients),unname(h$b),tolerance=1e-12)),identical(as.vector(f$iters),as.vector(h$draws)))
cat('Pinned public rglmb and validation harness: intercept CMB coefficients agree to 1e-12 and proposal counts are identical, 20,000 draws.\n')
# Sorting and unsorted diagnostics retain cell identity; match by GridIndex.
s<-e$standard
build<-function(sort)glmbayes:::EnvelopeBuild_cpp_export(s$bstar2,s$A,a$y,s$x2,s$mu2,s$P2,e$alpha,a$weights,'cmb','identity',2L,20000L,1000L,sort,FALSE,FALSE,TRUE,60L)
x<-build(FALSE);y<-build(TRUE);key<-function(m)apply(m,1,paste,collapse=',');idx<-match(key(y$GridIndex),key(x$GridIndex))
stopifnot(isTRUE(all.equal(y$logPLSD,x$logPLSD[idx])),identical(y$refinement$cell_grid_index,x$refinement$cell_grid_index))
cat('Sorted logPLSD aligns with GridIndex; nested refinement vectors retain original order with cell_grid_index.\n')
