diagnostics_case <- function(parallel=FALSE) {
  X <- matrix(1,20,1);a <- cmb_augment(X,rep(c(0,1,3,5),5),6,Z=X)
  rglmb(n=200L,y=a$y,x=a$x,offset=a$offset,weights=a$weights,family=cmb(),
        pfamily=dNormal(mu=matrix(0,2,1),Sigma=25*diag(2)),
        n_envopt=1000L,use_parallel=parallel,verbose=FALSE)
}
test_that('successful envelope diagnostics agree with actual sampling', {
  set.seed(614);f <- diagnostics_case(TRUE);d <- f$diagnostics
  expect_s3_class(d,'glmb_diagnostics')
  expect_identical(d$status,'completed')
  expect_identical(d$sampling_mode,'parallel')
  expect_equal(d$draws,nrow(f$coefficients))
  expect_equal(d$total_proposals,sum(f$iters))
  expect_equal(d$proposals_per_draw,mean(f$iters))
  expect_equal(d$acceptance_rate,nrow(f$coefficients)/sum(f$iters))
  expect_equal(d$max_proposals,max(f$iters))
  expect_equal(d$envelope$residual,f$Envelope$refinement$residual)
  expect_equal(d$envelope$cells,length(f$Envelope$PLSD))
  expect_match(d$build$fingerprint[['compiled_code']],'^[a-f0-9]{32}$')
  expect_identical(summary(f)$diagnostics,d)
  expect_output(print(d),'Envelope: converged')
})
test_that('incomplete or absent refinement is never reported as convergence', {
  f <- diagnostics_case();f$Envelope$refinement$converged <- FALSE
  d <- glmbayes:::.glmb_fit_diagnostics(f,list(use_parallel=FALSE))
  expect_identical(d$status,'refinement_incomplete')
  expect_output(print(d),'Envelope: incomplete')
  f$Envelope$refinement <- NULL
  d <- glmbayes:::.glmb_fit_diagnostics(f,list())
  expect_identical(d$envelope$status,'not_reported')
  expect_true(is.na(d$envelope$residual))
  f$Envelope$refinement <- list(converged=NA)
  f$Envelope$PLSD <- NULL
  d <- glmbayes:::.glmb_fit_diagnostics(f,list())
  expect_identical(d$envelope$status,'not_reported')
  expect_true(is.na(d$envelope$cells))
  f$iters <- NULL
  expect_true(is.na(glmbayes:::.glmb_fit_diagnostics(f,list())$total_proposals))
})
test_that('formula fits preserve diagnostics for direct and envelope sampling', {
  d <- data.frame(y=c(1,2,1,4,3,2,5,4),x=seq(-1,1,length.out=8))
  f <- glmb(y~x,data=d,n=100L,family=poisson(),
            pfamily=dNormal(mu=matrix(0,2,1),Sigma=diag(2)),verbose=FALSE)
  expect_s3_class(f$diagnostics,'glmb_diagnostics')
  expect_identical(summary(f)$diagnostics,f$diagnostics)
  g <- lmb(y~x,data=d,n=100L,
           pfamily=dIndependent_Normal_Gamma(mu=c(2,0),Sigma=diag(2),shape=3,rate=2),verbose=FALSE)
  expect_identical(g$diagnostics$method,'envelope_rejection')
  expect_gt(g$diagnostics$envelope$cells,0)
  expect_equal(g$diagnostics$total_proposals,sum(g$iters))
  h <- lmb(y~x,data=d,n=100L,
           pfamily=dNormal(mu=matrix(0,2,1),Sigma=diag(2),dispersion=1),verbose=FALSE)
  expect_identical(h$diagnostics$envelope$status,'not_applicable')
  expect_true(is.na(h$diagnostics$acceptance_rate))
})
test_that('real sampler exhaustion supplies a structured condition without draws', {
  sample_bad <- function() glmbayes:::rNormalGLM_std_cpp_export(
    n=2L,y=0,x=matrix(0,1,1),mu=matrix(0,1,1),P=matrix(0,1,1),
    alpha=0,wt=1,f2=identity,Envelope=list(PLSD=1,logPLSD=0,
      cbars=matrix(0,1,1),loglt=matrix(0,1,1),logrt=matrix(0,1,1),
      logU=matrix(0,1,1),LLconst=-1e100),
    family='gaussian',link='identity',progbar=0L)
  e <- tryCatch(glmbayes:::.glmb_run_simulation(sample_bad(),list(use_parallel=FALSE)),
                glmbayes_sampling_error=function(e)e)
  expect_s3_class(e,'glmbayes_sampling_error')
  expect_identical(e$diagnostics$status,'failed')
  expect_equal(e$diagnostics$proposal_limit,200000)
  expect_identical(e$diagnostics$draws_returned,0L)
  expect_null(e$coefficients)
  expect_match(conditionMessage(e),'predictor scaling')
  # Unrelated errors retain their original class and message.
  original <- structure(list(message='unrelated failure',call=NULL),class=c('test_error','error','condition'))
  expect_error(glmbayes:::.glmb_run_simulation(stop(original),list()),class='test_error')
})
