# Synthetic CMB data and fitted model specifications

The six CSVs below contain the data used in the CMB comparisons. Each has 60 observations, 10 trials per observation, and `v` equally spaced from -1 to 1; `successes` is the integer response. The generator uses seed 719. Inputs were checked against the saved fits and references; the `.rds` files retain full precision.

## Fitted models

The natural-parameter likelihood is

$$
P(Y_i=y\mid\theta_i,\nu_i,m_i)
=\frac{\exp\{\theta_i y+\nu_i\log {m_i\choose y}\}}
{\sum_{j=0}^{m_i}\exp\{\theta_i j+\nu_i\log {m_i\choose j}\}},\qquad y=0,\ldots,m_i.
$$

Both predictors use identity links. Coefficients in the theta and nu predictors jointly have independent zero-mean normal priors with the common SD shown below. All offsets are zero; there are no random effects.

| Case and data | Fitted theta | Fitted nu | Prior SD | Coefficients |
|---|---|---|---|---|
| [Intercepts](intercepts-data.csv) | beta0 | gamma0 | 5 | 2 |
| [Theta slope](theta-slope-data.csv) | beta0 + beta1 v | gamma0 | 5 | 3 |
| [Nu slope](nu-slope-data.csv) | beta0 | gamma0 + gamma1 v | 5 | 3 |
| [Strong prior](strong-prior-data.csv) | beta0 + beta1 v | gamma0 + gamma1 v | 0.5 | 4 |
| [Weak prior](weak-prior-data.csv) | beta0 + beta1 v | gamma0 + gamma1 v | 20 | 4 |
| [Zero-heavy](zero-heavy-data.csv) | beta0 + beta1 v | gamma0 + gamma1 v | 5 | 4 |

The intercepts data were generated with theta = -0.3 and nu = 1. The zero-heavy data were generated with theta = -3 + 0.3 v and nu = 0.5. The other four datasets were generated with theta = -0.3 + 0.5 v and nu = 1 + 0.3 v. Thus, the theta-slope and nu-slope labels describe the **fitted** predictor restrictions; their generating distribution still varies both predictors. The strong- and weak-prior cases share that same realized dataset. See the [pinned generator](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/cmb-matrix.R) for the complete generation and reference-integration code.

## Zero-heavy example: the practical refinement failure

This dataset has **48 zeros among 60 observations**, 18 total successes, and observed counts ranging from 0 to 4. Both fitted predictors include an intercept and `v`; all four coefficients have independent Normal(0, 25) priors (SD 5). The nu slope is estimated even though its generating value is zero.

The paired runs use the same inputs and prior, the same mode optimization and standardization, `Gridtype=2`, `n_envopt=1000`, sorted envelopes, serial CPU execution, and no OpenCL. They request 20,000 accepted draws, with one warm-up (seed 9100) and three measured runs (seeds 9101, 9102, 9103) per setting. Only the generated builder's refinement flag changes; `refine_maxit=60`, the other numerical safeguards, the 200,000-proposal-per-draw cap, and the 300-second process timeout are retained. Both settings produce an 81-cell grid.

The [comparison script](../compare.R) uses the generated builder to toggle refinement, because the public fitting interface does not expose that control.

With refinement, all three measured runs finish (median total 0.594 seconds). Without it, all three hit the proposal cap before accepting their first draw. [Paired results](../results/summary.md) and the saved [enabled envelope/input](../results/zero_heavy-1-1-envelope.rds), [disabled envelope/input](../results/zero_heavy-0-1-envelope.rds), and [enabled draws](../results/zero_heavy-1-1.rds) permit inspection of the original inputs, standardized matrices, tangencies, and posterior draws using `readRDS()`.

The comparison shows where refinement is needed under these settings. It does not yet establish whether non-normality or correlation causes the difficulty, or how to detect such cases automatically.

## Other comparison inputs

The [same harness](../compare.R) specifies the public menarche binomial data, the nine-count Poisson example, and the two-group Gaussian example. The separate stiff-axis example is synthetic: a deliberately displaced starting tangency, two coordinates, precision diagonal (10, 1000), and a single-cell grid (`Gridtype=4`). It is an artificial placement stress test. The Gaussian comparison forces an envelope path that ordinary direct Gaussian fits do not use.
