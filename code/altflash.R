
altflash <- function(obj,
                     solQb_nupdates = 1,
                     maxiter = 500,
                     tol = NULL
                     ){
  if(!is.null(tol)){
    obj$tol <- tol
  }

  # record time and elbo
  init.time <- Sys.time()
  vec.time <- c(0)
  vec.elbo <- obj$elbo

  for(i in 1:maxiter){

    obj <- solQb(obj, L_side = TRUE, nupdates = solQb_nupdates)
    obj <- solQb(obj, L_side = FALSE, nupdates = solQb_nupdates)
    obj <- solTAU(obj)

    # record time and elbo ###################
    vec.time <- c(vec.time, as.numeric(Sys.time() - init.time, units = "secs"))
    vec.elbo <- c(vec.elbo, obj$elbo)
    ##########################################

    # check convergence
    if(abs(vec.elbo[i+1] - vec.elbo[i]) < obj$tol){
      break
    }
  }

  obj <- c(obj, list(vec.elbo = vec.elbo, vec.time = vec.time))

  return(obj)
}



altflash_normal <- function(obj,
                            maxiter = 500,
                            tol = NULL
                            ){
  if(!is.null(tol)){
    obj$tol <- tol
  }

  obj$L_ghat_tau <- sapply(obj$L_ghat, function(x) 1/x$sd^2)
  obj$F_ghat_tau <- sapply(obj$F_ghat, function(x) 1/x$sd^2)

  # record time and elbo
  init.time <- Sys.time()
  vec.time <- c(0)
  vec.elbo <- obj$elbo

  for(i in 1:maxiter){

    obj <- solG_normal(obj)
    obj <- solQ_normal(obj, L_side=TRUE)
    obj <- solQ_normal(obj, L_side=FALSE)
    obj <- solTAU(obj)

    # record time and elbo ###################
    vec.time <- c(vec.time, as.numeric(Sys.time() - init.time, units = "secs"))
    vec.elbo <- c(vec.elbo, obj$elbo)
    ##########################################

    # check convergence
    if(abs(vec.elbo[i+1] - vec.elbo[i]) < obj$tol){
      break
    }
  }

  obj <- c(obj, list(vec.elbo = vec.elbo, vec.time = vec.time))

  return(obj)
}













