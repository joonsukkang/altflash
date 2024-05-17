solG <- function(obj){

  # L side
  for (k in 1:length(obj$L_ghat)){
    x <- c((as.vector(obj$X%*%obj$EF[,k]) - obj$EL[,-k, drop=FALSE]%*%obj$EFtF[-k,k])/obj$EFtF[k,k])
    s <- 1/sqrt(obj$tau*obj$EFtF[k,k])
    obj$L_ghat[[k]] <- obj$L_ebnm_fn(x=x, s=s,
                                     g_init = obj$L_ghat[[k]],
                                     output='fitted_g')$fitted_g
  }

  # F side
  for (k in 1:length(obj$F_ghat)){
    x <- c((as.vector(t(obj$X)%*%obj$EL[,k]) - obj$EF[,-k, drop=FALSE]%*%obj$ELtL[-k,k])/obj$ELtL[k,k])
    s <- 1/sqrt(obj$tau*obj$ELtL[k,k])
    obj$F_ghat[[k]] <- obj$F_ebnm_fn(x=x, s=s,
                                     g_init = obj$F_ghat[[k]],
                                     output='fitted_g')$fitted_g
  }

  return(obj)
}


solG_normal <- function(obj){

  obj$L_ghat_tau <- nrow(obj$EL)/diag(obj$ELtL)
  obj$F_ghat_tau <- nrow(obj$EF)/diag(obj$EFtF)

  return(obj)
}
