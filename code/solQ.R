

solQb <- function(obj, L_side = TRUE, nupdates=1){

  if(L_side){

    obj$XEF <- obj$X %*% obj$EF
    for (i in 1:nupdates){

      obj$L_KL <- 0
      for (k in 1:length(obj$L_ghat)){
        x <- c((as.vector(obj$XEF[,k]) - obj$EL[,-k, drop=FALSE]%*%obj$EFtF[-k,k])/obj$EFtF[k,k])
        s <- rep(1/sqrt(obj$tau*obj$EFtF[k,k]), times = length(x))
        out <- obj$L_ebnm_fn(x=x, s=s,
                             g_init = obj$L_ghat[[k]],
                             output=c('posterior_mean',
                                      "posterior_second_moment",
                                      "log_likelihood",
                                      'fitted_g'))
        obj$EL[,k] <- out$posterior$mean
        obj$EL2[,k] <- out$posterior$second_moment
        obj$L_KL <- obj$L_KL -(c(out$log_likelihood) +
                                 sum(log(2*pi*s^2)/2) +
                                 sum((x^2 - 2*x*obj$EL[,k] + obj$EL2[,k])/(2*s^2)))

        obj$L_ghat[[k]] <- out$fitted_g
      }
    }
    obj$ELtL <- crossprod(obj$EL);
    diag(obj$ELtL) <- colSums(obj$EL2)

  }else{

    obj$XtEL <- t(obj$X) %*% obj$EL
    for (i in 1:nupdates){

      obj$F_KL <- 0
      for (k in 1:length(obj$F_ghat)){
        x <- c((as.vector(obj$XtEL[,k]) - obj$EF[,-k, drop=FALSE]%*%obj$ELtL[-k,k])/obj$ELtL[k,k])
        s <- rep(1/sqrt(obj$tau*obj$ELtL[k,k]), times = length(x))
        out <- obj$F_ebnm_fn(x=x, s=s,
                             g_init = obj$F_ghat[[k]],
                             output=c('posterior_mean',
                                      "posterior_second_moment",
                                      "log_likelihood",
                                      'fitted_g'))
        obj$EF[,k] <- out$posterior$mean
        obj$EF2[,k] <- out$posterior$second_moment
        obj$F_KL <- obj$F_KL -(c(out$log_likelihood) +
                                 sum(log(2*pi*s^2)/2) +
                                 sum((x^2 - 2*x*obj$EF[,k] + obj$EF2[,k])/(2*s^2)))

        obj$F_ghat[[k]] <- out$fitted_g

      }
    }
    obj$EFtF <- crossprod(obj$EF);
    diag(obj$EFtF) <- colSums(obj$EF2)
  }

  return(obj)
}


solQ_normal <- function(obj, L_side = TRUE){

  if(L_side){
    EL <- t(solve(obj$EFtF + 1/obj$tau * diag(obj$L_ghat_tau), t(obj$EF) %*% t(obj$X)))
    EL2 <- EL^2 + matrix(1, nrow=nrow(EL), ncol=1) %*% t(1/(obj$L_ghat_tau + obj$tau * diag(obj$EFtF)))

    ELtL <- crossprod(EL);
    diag(ELtL) <- colSums(EL2)

    obj$EL <- EL
    obj$EL2 <- EL2
    obj$ELtL <- ELtL

    obj$L_KL <- nrow(obj$EL) * sum(log(1/sqrt(obj$L_ghat_tau))-1/2) +
      -sum(log(sqrt(obj$EL2-obj$EL^2))) + sum(diag(obj$ELtL)/(2/obj$L_ghat_tau))
  }else{
    EF <- t(solve(obj$ELtL + 1/obj$tau * diag(obj$F_ghat_tau), t(obj$EL) %*% obj$X))
    EF2 <- EF^2 + matrix(1, nrow=nrow(EF), ncol=1) %*% t(1/(obj$F_ghat_tau + obj$tau * diag(obj$ELtL)))

    EFtF <- crossprod(EF);
    diag(EFtF) <- colSums(EF2)

    obj$EF <- EF
    obj$EF2 <- EF2
    obj$EFtF <- EFtF

    obj$F_KL <- nrow(obj$EF) * sum(log(1/sqrt(obj$F_ghat_tau))-1/2) +
      -sum(log(sqrt(obj$EF2-obj$EF^2))) + sum(diag(obj$EFtF)/(2/obj$F_ghat_tau))
  }

  return(obj)
}
