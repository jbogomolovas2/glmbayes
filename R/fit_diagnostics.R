# Shared diagnostics for rglmb/rlmb and their formula interfaces.
.glmb_diagnostics_cache <- new.env(parent=emptyenv())
.glmb_build_info <- function() {
  if (!is.null(.glmb_diagnostics_cache$build)) return(.glmb_diagnostics_cache$build)
  path <- getNamespaceInfo(asNamespace('glmbayes'),'path')
  dll <- getLoadedDLLs()[['glmbayes']]
  files <- c(compiled_code=if(is.null(dll))NA_character_ else dll[['path']],
             r_code=file.path(path,'R','glmbayes.rdb'))
  hashes <- setNames(rep(NA_character_,length(files)),names(files))
  ok <- !is.na(files) & file.exists(files)
  hashes[ok] <- unname(tools::md5sum(files[ok]))
  info <- list(version=as.character(getNamespaceVersion('glmbayes')),
               library=path,fingerprint=hashes,R_version=R.version.string,
               platform=R.version$platform)
  .glmb_diagnostics_cache$build <- info
  info
}
.glmb_fit_diagnostics <- function(fit,args) {
  E <- fit$Envelope; envelope <- is.list(E)
  ref <- if(envelope) E$refinement else NULL
  convergence <- if(!envelope)'not_applicable' else if(isTRUE(ref$converged))
    'converged' else if(identical(ref$converged,FALSE))'incomplete' else 'not_reported'
  n <- if(!is.null(fit$coefficients))nrow(as.matrix(fit$coefficients)) else length(fit$dispersion)
  counts <- as.numeric(fit$iters)
  valid_counts <- envelope && length(counts)==n && n>0 &&
    all(is.finite(counts) & counts>=1 & counts==floor(counts))
  total <- if(valid_counts)sum(counts) else NA_real_
  d <- list(status=if(convergence=='incomplete')'refinement_incomplete' else 'completed',
    method=if(envelope)'envelope_rejection' else 'direct',
    sampling_mode=if(!envelope)'not_applicable' else if(isTRUE(args$use_parallel))'parallel' else 'serial',
    requested_parallel=isTRUE(args$use_parallel),requested_opencl=isTRUE(args$use_opencl),
    draws=n,total_proposals=total,
    proposals_per_draw=if(valid_counts)total/n else NA_real_,
    acceptance_rate=if(valid_counts)n/total else NA_real_,
    max_proposals=if(valid_counts)max(counts) else NA_real_,
    envelope=list(status=convergence,cells=if(envelope && !is.null(E$PLSD))length(E$PLSD) else NA_integer_,
      iterations=if(is.null(ref$iterations))NA_integer_ else ref$iterations,
      residual=if(is.null(ref$residual))NA_real_ else ref$residual),
    build=.glmb_build_info(),rng_kind=RNGkind())
  class(d) <- 'glmb_diagnostics'
  d
}
.glmb_run_simulation <- function(expr,args) {
  fit <- tryCatch(force(expr),error=function(e) {
    message <- conditionMessage(e)
    if(!grepl('rejection sampler reached [0-9]+ proposals',message)) stop(e)
    cap <- as.numeric(sub('.*rejection sampler reached ([0-9]+) proposals.*','\\1',message))
    advice <- 'Inspect the envelope and predictor scaling before retrying; requesting more draws will not improve acceptance.'
    d <- structure(list(status='failed',reason='proposal_limit',draws_returned=0L,
      proposal_limit=cap,total_proposals=NA_real_,
      sampling_mode=if(isTRUE(args$use_parallel))'parallel' else 'serial',
      action=advice,build=.glmb_build_info()),class='glmb_diagnostics')
    stop(structure(list(message=paste(message,advice),call=conditionCall(e),
                        parent=e,diagnostics=d),
                   class=c('glmbayes_sampling_error','error','condition')))
  })
  fit$diagnostics <- .glmb_fit_diagnostics(fit,args)
  fit
}

#' Print sampling diagnostics
#' @param x A fit's diagnostics object.
#' @param ... Unused additional arguments.
#' @return The diagnostics object, invisibly.
#' @export
print.glmb_diagnostics <- function(x,...) {
  if(identical(x$status,'failed')) {
    cat('Sampling failed: proposal limit reached; no draws returned.\n',x$action,'\n',sep='')
    return(invisible(x))
  }
  if(identical(x$method,'direct')) {
    cat('Sampling: direct; envelope diagnostics not applicable.\n')
  } else {
    cat(sprintf('Sampling: %s; %s draws; %s proposals (%.3f per draw).\n',
                x$sampling_mode,format(x$draws,trim=TRUE),
                format(x$total_proposals,trim=TRUE),x$proposals_per_draw))
    cat(sprintf('Envelope: %s; %s cells; %s passes; residual %s.\n',
                x$envelope$status,x$envelope$cells,x$envelope$iterations,
                format(x$envelope$residual,digits=3)))
    if(identical(x$envelope$status,'incomplete'))
      cat('Refinement stopped before convergence; this concerns efficiency, not an MCMC convergence test.\n')
  }
  invisible(x)
}
