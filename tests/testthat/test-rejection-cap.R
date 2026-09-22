# Exercise the compiled serial sampler without optimization or refinement.
# The residual target is constant; the proposal is an untruncated normal.
constant_target_draws <- function(intercept, n = 4L) {
  glmbayes:::rNormalGLM_std_cpp_export(
    n = n, y = 0, x = matrix(0, 1, 1), mu = matrix(0, 1, 1),
    P = matrix(0, 1, 1), alpha = 0, wt = 1, f2 = identity,
    Envelope = list(PLSD = 1, logPLSD = 0, cbars = matrix(0, 1, 1),
                    loglt = matrix(0, 1, 1), logrt = matrix(0, 1, 1),
                    LLconst = intercept),
    family = 'gaussian', link = 'identity', progbar = 0L)
}

test_that('a healthy envelope returns only accepted draws', {
  set.seed(103)
  ans <- constant_target_draws(0.5 * log(2 * pi), 100L)
  expect_equal(ans$draws, rep(1, 100))
  expect_true(all(is.finite(ans$out)))
  expect_equal(nrow(ans$out), 100L)
})

test_that('a nonaccepting envelope fails instead of returning partial draws', {
  expect_error(constant_target_draws(-1e300),
               '200000 proposals with zero acceptances.*no draws were returned')
})

test_that('NaN acceptance calculations cannot bypass the proposal counter', {
  expect_error(constant_target_draws(NaN),
               '200000 proposals with zero acceptances.*no draws were returned')
})
