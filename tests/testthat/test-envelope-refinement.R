# A Gaussian residual target yields an affine fixed-point equation. Starting
# away from its root exposes the stiff unsplit-axis failure without sampling.
quadratic_envelope <- function(curvature, center = 1, grid = 4L) {
  p <- length(curvature)
  glmbayes:::EnvelopeBuild_cpp_export(
    bStar = matrix(rep(center,p),p,1), A = diag(curvature,p),
    y = 0, x = matrix(0,1,p), mu = matrix(0,p,1), P = diag(curvature,p),
    alpha = 0, wt = 1, family = 'gaussian', link = 'identity',
    Gridtype = grid, n = 1L, n_envopt = 1L,
    sortgrid = FALSE, use_opencl = FALSE, verbose = FALSE)
}

test_that('Newton refinement resolves stiff unsplit axes', {
  E <- quadratic_envelope(c(10,1000))
  d <- E$refinement
  expect_true(d$accelerated)
  expect_true(d$converged)
  expect_lte(d$iterations,5)
  expect_lt(d$residual,1e-6)
  expect_equal(as.vector(E$thetabars),c(0,0),tolerance=1e-6)
  expect_lte(d$log_mass_final,d$log_mass_initial)
  expect_equal(E$PLSD,1)
})

test_that('a fixed point requires no update passes', {
  E <- quadratic_envelope(c(2,4),center=0)
  expect_equal(E$refinement$iterations,0L)
  expect_equal(E$refinement$residual,0)
  expect_true(E$refinement$converged)
})

test_that('large parameter counts retain the bounded-memory fallback', {
  E <- quadratic_envelope(rep(.5,17))
  expect_false(E$refinement$accelerated)
  expect_true(E$refinement$converged)
  expect_lt(E$refinement$residual,1e-6)
})
