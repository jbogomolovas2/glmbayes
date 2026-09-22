# Synthetic example; no study data or local paths required.
library(glmbayes)
set.seed(2026)
clutch <- rep(8L, 30)
size <- seq(-1, 1, length.out = length(clutch))
theta <- -0.3 + 0.4 * size
nu <- 0.8 + 0.2 * size
hatch <- vapply(seq_along(clutch), function(i) {
  y <- 0:clutch[i]
  logw <- theta[i] * y + nu[i] * lchoose(clutch[i], y)
  sample(y, 1, prob = exp(logw - max(logw)))
}, numeric(1))
X <- Z <- cbind(intercept = 1, size = size)
a <- cmb_augment(X, hatch, clutch, Z = Z)
p <- ncol(a$x)
fit <- rglmb(n = 2000L, y = a$y, x = a$x, offset = a$offset,
             weights = a$weights, family = cmb(),
             pfamily = dNormal(mu = matrix(0, p, 1), Sigma = 4 * diag(p)),
             n_envopt = 1000L, use_parallel = TRUE, verbose = FALSE)
print(fit$diagnostics)
expected_counts <- cmb_fitted(fit, scale = 'count')
expected_proportions <- cmb_fitted(fit, scale = 'proportion')
print(head(data.frame(size = size, observed = hatch,
                      expected = colMeans(expected_counts),
                      proportion = colMeans(expected_proportions))))
