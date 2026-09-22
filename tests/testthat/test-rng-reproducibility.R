# Per-output-row streams must not depend on worker scheduling or calibration.
# Compare sample values and R RNG state, not timing-dependent fit metadata.
with_rng_threads <- function(code) {
  old <- Sys.getenv('RCPP_PARALLEL_NUM_THREADS',unset=NA_character_)
  on.exit(if(is.na(old)) Sys.unsetenv('RCPP_PARALLEL_NUM_THREADS') else
            Sys.setenv(RCPP_PARALLEL_NUM_THREADS=old))
  force(code)
}
rng_cmb_fit <- function(seed=312,parallel=TRUE,threads=2L,n=400L) {
  RcppParallel::setThreadOptions(numThreads=threads)
  set.seed(seed)
  x <- cbind(1,seq(-1,1,length.out=30))
  a <- cmb_augment(x,rep(c(0,1,2,4,6),6),6,Z=matrix(1,30,1))
  f <- rglmb(n=n,y=a$y,x=a$x,weights=a$weights,offset=a$offset,
             family=cmb(),pfamily=dNormal(mu=matrix(0,3,1),Sigma=25*diag(3)),
             n_envopt=1000L,use_parallel=parallel,verbose=FALSE)
  list(beta=f$coefficients,iters=f$iters,state=.Random.seed,next_random=runif(5))
}
test_that('CMB output and R RNG state repeat across scheduling choices', {
  with_rng_threads({
    a <- rng_cmb_fit(parallel=FALSE)
    expect_identical(rng_cmb_fit(parallel=FALSE),a)
    for(threads in c(1L,2L,4L)) {
      b <- rng_cmb_fit(threads=threads)
      expect_identical(b,a)
    }
    expect_false(identical(rng_cmb_fit(seed=313)$beta,a$beta))
    expect_equal(nrow(unique(a$beta)),nrow(a$beta))
  })
})
test_that('draw prefixes and R state are independent of requested draw count', {
  with_rng_threads({
    a <- rng_cmb_fit(n=400L)
    b <- rng_cmb_fit(n=700L,threads=4L)
    expect_identical(b$beta[seq_len(400),,drop=FALSE],a$beta)
    expect_identical(b$iters[seq_len(400)],as.vector(a$iters))
    expect_identical(b$state,a$state)
    expect_identical(b$next_random,a$next_random)
  })
})
test_that('independent Normal-Gamma draws repeat across thread counts', {
  with_rng_threads({
    fit <- function(parallel,threads=2L) {
      RcppParallel::setThreadOptions(numThreads=threads);set.seed(97)
      d <- data.frame(y=c(4.17,5.58,5.18,6.11,4.5,4.61,5.17,4.53,5.33,5.14),x=seq(-1,1,length.out=10))
      f <- lmb(y~x,data=d,n=300L,
               pfamily=dIndependent_Normal_Gamma(mu=c(5,0),Sigma=diag(2),shape=3,rate=2),
               verbose=FALSE,use_parallel=parallel)
      stopifnot(length(f$dispersion)==300L,length(f$iters)==300L)
      list(beta=f$coefficients,dispersion=f$dispersion,iters=f$iters,state=.Random.seed)
    }
    a <- fit(FALSE)
    expect_identical(fit(FALSE),a)
    b <- fit(TRUE,1L)
    expect_identical(fit(TRUE,4L),b)
    # Both algorithms use the same per-row stream. Check their public outputs.
    expect_identical(b,a)
  })
})

test_that('R RNG advancement is fixed and distinct calls get new keys', {
  draw <- function() glmbayes:::rNormalGLM_std_cpp_export(
    n=10L,y=0,x=matrix(0,1,1),mu=matrix(0,1,1),P=matrix(0,1,1),
    alpha=0,wt=1,f2=identity,Envelope=list(PLSD=1,logPLSD=0,
      cbars=matrix(0,1,1),loglt=matrix(0,1,1),logrt=matrix(0,1,1),
      logU=matrix(0,1,1),LLconst=.5*log(2*pi)),
    family='gaussian',link='identity',progbar=0L)
  set.seed(204);runif(2);expected_state <- .Random.seed
  set.seed(204);a <- draw()
  expect_identical(.Random.seed,expected_state)
  b <- draw()
  expect_false(identical(a$out,b$out))
  set.seed(204)
  expect_identical(draw(),a)
})

test_that('the L Ecuyer R generator also supplies reproducible stream keys', {
  old <- RNGkind();on.exit(do.call(RNGkind,as.list(old)))
  RNGkind("L'Ecuyer-CMRG")
  with_rng_threads({
    a <- rng_cmb_fit(parallel=FALSE)
    expect_identical(rng_cmb_fit(threads=4L),a)
  })
})
