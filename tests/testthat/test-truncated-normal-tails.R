# Force acceptance to inspect the production proposal independently of the
# posterior and rejection rule. Shifts put a narrow cell 100--1000 SD into
# either tail, where one of the two legacy log-tail probabilities is zero.
test_that('normal proposals preserve narrow cells in extreme tails', {
  cb <- c(100, -100, 1000, -1000, 0)
  lo <- c(0, -.01, 0, -.001, -1)
  hi <- c(.01, 0, .001, 0, 1)
  down <- lo+cb; up <- hi+cb
  lt <- pnorm(up,log.p=TRUE)
  rt <- pnorm(down,lower.tail=FALSE,log.p=TRUE)
  lm <- vapply(seq_along(cb),function(j) {
    if (down[j]>=0) {
      a <- pnorm(down[j],lower.tail=FALSE,log.p=TRUE)
      b <- pnorm(up[j],lower.tail=FALSE,log.p=TRUE)
    } else {
      a <- pnorm(up[j],log.p=TRUE)
      b <- pnorm(down[j],log.p=TRUE)
    }
    a+log(-expm1(b-a))
  },numeric(1))
  p <- length(cb); n <- 20000L
  E <- list(PLSD=1,logPLSD=0,cbars=matrix(cb,1),
            loglt=matrix(lt,1),logrt=matrix(rt,1),logU=matrix(lm,1),
            LLconst=1e100)
  ans <- glmbayes:::rNormalGLM_std_cpp_export(
    n=n,y=0,x=matrix(0,1,p),mu=matrix(0,p,1),P=matrix(0,p,p),
    alpha=0,wt=1,f2=identity,Envelope=E,family='gaussian',
    link='identity',progbar=0L)
  expect_equal(ans$draws,rep(1,n))
  expect_true(all(is.finite(ans$out)))
  for(j in seq_len(p)) {
    v <- ans$out[,j]
    expect_true(all(v>=lo[j]-1e-10 & v<=hi[j]+1e-10))
    # Independent integration on the actual cell, without normal CDF tails.
    density <- function(x) exp(-cb[j]*x-x*x/2)
    z <- integrate(density,lo[j],hi[j])$value
    m <- integrate(function(x) x*density(x),lo[j],hi[j])$value/z
    vref <- integrate(function(x) (x-m)^2*density(x),lo[j],hi[j])$value/z
    expect_lt(abs(mean(v)-m),6*sqrt(vref/n))
    expect_lt(abs(mean((v-m)^2)-vref),6*sd((v-m)^2)/sqrt(n))
  }
})
