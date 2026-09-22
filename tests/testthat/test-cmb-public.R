# Generated/toy inputs only; independent finite-support probability calculation.
test_that('CMB likelihood and derivatives agree with direct finite sums', {
  m <- c(4, 7, 5, 8); y <- c(0, 3, 5, 2)
  X <- cbind(1, c(-1, -.2, .4, 1))
  a <- cmb_augment(X, y, m)
  ff <- glmbayes:::glmbfamfunc(cmb())
  mu <- matrix(0, 3, 1); P <- diag(3)/5
  direct <- function(b) {
    theta <- as.vector(X %*% b[1:2]); nu <- b[3]
    ll <- vapply(seq_along(y), function(i) {
      k <- 0:m[i]; lw <- theta[i]*k + nu*lchoose(m[i], k)
      theta[i]*y[i] + nu*lchoose(m[i], y[i]) - max(lw) - log(sum(exp(lw-max(lw))))
    }, numeric(1))
    -sum(ll) + .5*as.numeric(crossprod(b, P %*% b))
  }
  for (nu in c(-.4, .7, 1, 2)) {
    b <- c(-.3, .5, nu)
    actual <- ff$f2(b,a$y,a$x,mu,P,a$offset,a$weights)
    expect_equal(as.numeric(actual),direct(b),tolerance=1e-10)
    numerical <- vapply(seq_along(b),function(j) {
      delta <- rep(0,3); delta[j] <- 1e-6
      (direct(b+delta)-direct(b-delta))/2e-6
    },numeric(1))
    expect_equal(as.numeric(ff$f3(b,a$y,a$x,mu,P,a$offset,a$weights)),numerical,tolerance=1e-6)
  }
  b <- c(-.3,.5,1)
  binomial <- -sum(dbinom(y,m,plogis(X %*% b[1:2]),log=TRUE)) + .5*sum(b*(P%*%b))
  expect_equal(direct(b),as.numeric(binomial),tolerance=1e-10)
})

test_that('fitted CMB means integrate the distribution for each draw', {
  X <- cbind(1,c(-1,0,1)); m <- c(4,6,8)
  a <- cmb_augment(X,c(0,3,7),m)
  b <- rbind(c(-.2,.4,1),c(-.2,.4,.6))
  object <- list(family=cmb(),x=a$x,prior.weights=a$weights,coefficients=b)
  expected <- t(vapply(seq_len(nrow(b)),function(j) vapply(seq_along(m),function(i) {
    k <- 0:m[i]; lw <- sum(X[i,]*b[j,1:2])*k+b[j,3]*lchoose(m[i],k)
    w <- exp(lw-max(lw)); sum(k*w)/sum(w)
  },numeric(1)),numeric(3)))
  expect_equal(cmb_fitted(object,'count'),expected)
  expect_equal(cmb_fitted(object,'proportion'),sweep(expected,2,m,'/'))
  expect_equal(as.numeric(cmb_fitted(object)[1,]),as.numeric(plogis(X%*%b[1,1:2])))
  expect_equal(unname(cmb_nu(object)),matrix(c(1,.6),2,3))
})
