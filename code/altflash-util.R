

flash_init_to_altflash <- function(init,
                                   L_ebnm_fn,
                                   F_ebnm_fn){

  obj <- list(X = init$flash_fit$Y,
              EL = init$L_pm,
              EL2 = init$L_pm^2 + init$L_psd^2,
              ELtL = crossprod(init$L_pm),
              EF = init$F_pm,
              EF2 = init$F_pm^2 + init$F_psd^2,
              EFtF = crossprod(init$F_pm),
              tau = init$flash_fit$tau,
              L_ghat = init$L_ghat,
              F_ghat = init$F_ghat,
              L_ebnm_fn = L_ebnm_fn,
              F_ebnm_fn = F_ebnm_fn,
              tol = init$flash_fit$conv.tol,
              elbo = init$elbo)

  diag(obj$ELtL) <- diag(obj$ELtL) + colSums(init$L_psd^2)
  diag(obj$EFtF) <- diag(obj$EFtF) + colSums(init$F_psd^2)
  obj$trXtX <- sum(obj$X^2)

  return(obj)
}


# modify 'flash_backfit' function to record elbo and time
flash_backfit_modified <- function (flash, kset = NULL, extrapolate = TRUE, warmstart = TRUE,
          maxiter = 500, tol = NULL, verbose = NULL) {
  flash <- flashier:::get.fit(flash)
  tol <- flashier:::handle.tol.param(tol, flash)
  verbose.lvl <- flashier:::handle.verbose.param(verbose, flash)
  if (flashier:::is.timed.out(flash)) {
    report.timeout.no.backfit(verbose.lvl)
    verbose.lvl <- 0
  }
  if (is.null(kset)) {
    if (flashier:::get.n.factors(flash) > 0) {
      kset <- 1:flashier:::get.n.factors(flash)
    }
    else {
      flashier:::announce.no.backfit(verbose.lvl)
      verbose.lvl <- 0
    }
  }
  else {
    must.be.valid.kset(flash, kset)
  }
  flashier:::must.be.integer(maxiter, lower = 1, allow.null = FALSE)
  flashier:::must.be.integer(verbose.lvl, lower = -1, upper = 3, allow.null = FALSE)
  if (extrapolate) {
    method <- "extrapolate"
  }
  else {
    method <- "sequential"
  }
  if (method == "parallel") {
    flashier:::check.parallel.ok(flash, kset)
    if (missing(conv.crit.fn)) {
      conv.crit.fn <- function(new, old, k) {
        return(abs(flashier:::calc.obj.diff(new, old, k)))
      }
    }
  }
  flash <- flashier:::set.warmstart(flash, warmstart)
  verbose.fns <- flashier:::get.verbose.fns(flash)
  verbose.colnames <- flashier:::get.verbose.colnames(flash)
  verbose.colwidths <- flashier:::get.verbose.colwidths(flash)
  conv.crit.fn <- flashier:::get.conv.crit.fn(flash)
  conv.crit <- rep(Inf, flashier:::get.n.factors(flash))
  conv.crit[setdiff(1:flashier:::get.n.factors(flash), kset)] <- 0
  flashier:::announce.backfit(verbose.lvl, n.factors = length(kset), tol)
  flashier:::print_table.header(verbose.lvl, verbose.colnames, verbose.colwidths,
                     backfit = TRUE)
  if (method == "parallel") {
    kset <- setdiff(kset, which(is.zero(flash)))
    kset <- setdiff(kset, which.k.fixed(flash))
    cl <- parallel::makeCluster(getOption("cl.cores", 2L),
                                type = getOption("cl.type", "PSOCK"), useXDR = FALSE)
  }
  else if (method == "extrapolate") {
    extrapolate.control <- getOption("extrapolate.control",
                                     list())
    extrapolate.param <- flashier:::set.extrapolate.param(extrapolate.control)
  }
  iter <- 0
  old.obj <- flashier:::get.obj(flash)
  next.tol.target <- NULL
  if (method == "extrapolate") {
    extrapolate.param <- flashier:::init.beta(extrapolate.param)
    old.f <- flash
  }

  # record time and elbo
  init.time <- Sys.time()
  vec.time <- c(0)
  vec.elbo <- old.obj


  while (iter < maxiter && max(conv.crit) > tol && !flashier:::is.timed.out(flash)) {
    iter <- iter + 1
    kset <- flashier:::get.next.kset(method, kset, conv.crit, tol)
    if (!(method %in% c("parallel", "extrapolate"))) {
      for (k in kset) {
        old.f <- flash
        flash <- flashier:::update.one.factor(flash, k, iter, verbose.lvl)
        info <- flashier:::calc.update.info(flash, old.f, conv.crit.fn,
                                 verbose.fns, k)
        conv.crit[k] <- flashier:::get.conv.crit(info)
        flashier:::print_table.entry(verbose.lvl, verbose.colwidths,
                          iter, info, k = k, backfit = TRUE)
      }
    }
    else {
      if (method == "parallel") {
        old.f <- flash
        flash <- flashier:::update.factors.parallel(flash, kset,
                                         cl)
      }
      else if (method == "extrapolate") {
        proposed.f <- flashier:::extrapolate.f(flash, old.f, extrapolate.param)
        proposed.f <- flashier:::update.factors.in.kset(proposed.f,
                                             kset)
        old.f <- flash
        if (flashier:::get.obj(proposed.f) - flashier:::get.obj(flash) < tol) {
          flash <- flashier:::update.factors.in.kset(flash, kset)
          extrapolate.param <- flashier:::decelerate(extrapolate.param)
        }
        else {
          flash <- proposed.f
          extrapolate.param <- flashier:::accelerate(extrapolate.param)
        }
      }
      info <- flashier:::calc.update.info(flash, old.f, conv.crit.fn,
                               verbose.fns)
      conv.crit <- flashier:::get.conv.crit(info)
      flashier:::print_table.entry(verbose.lvl, verbose.colwidths,
                        iter, info, k = "all", backfit = TRUE)
    }
    if (is.null(next.tol.target) && max(conv.crit) > 0 &&
        max(conv.crit) < Inf) {
      next.tol.target <- 10^floor(log10(max(conv.crit)))
    }
    else if (!is.null(next.tol.target) && max(conv.crit) <
             next.tol.target) {
      flashier:::report.backfit.progress(verbose.lvl, next.tol.target)
      next.tol.target <- next.tol.target/10
    }

    # record time and elbo ###################
    vec.time <- c(vec.time, as.numeric(Sys.time() - init.time, units = "secs"))
    vec.elbo <- c(vec.elbo, flashier:::get.obj(flash))
    ##########################################

  }
  if (iter > 0 && flashier:::is.timed.out(flash)) {
    t.diff <- Sys.time() - flashier:::get.timeout.set.time(flash)
    flashier:::report.timeout.reached(verbose.lvl, t.diff)
    flash <- flashier:::set.timeout.reached.flag(flash)
  }
  if (iter == maxiter) {
    flashier:::report.maxiter.reached(verbose.lvl)
    flash <- flashier:::set.max.backfit.iter.reached.flag(flash)
  }
  else {
    flash <- flashier:::clear.max.backfit.iter.reached.flag(flash)
  }
  if (method == "parallel") {
    parallel::stopCluster(cl)
  }
  if (flashier:::get.obj(flash) > old.obj) {
    flashier:::report.backfit.complete(verbose.lvl, get.obj(flash))
  }
  flashier:::announce.wrapup(verbose.lvl)
  flash <- flashier:::wrapup.flash(flash, output.lvl = 3L)
  flashier:::report.completion(verbose.lvl)

  flash <- c(flash, list(vec.elbo = vec.elbo, vec.time = vec.time))
  return(flash)
}
