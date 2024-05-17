solTAU <- function(obj){

  obj$EResidualSq <- obj$trXtX +
    -2 * sum((obj$X %*% obj$EF)*obj$EL) + sum(obj$ELtL * obj$EFtF)

  obj$tau <- prod(dim(obj$X))/obj$EResidualSq

  obj$elbo <- prod(dim(obj$X))/2 * (log(obj$tau) - log(2*pi))+
    - obj$tau/2 * obj$EResidualSq +
    - obj$L_KL - obj$F_KL


  return(obj)
}


# computeELBO <- function(obj){
#
#   # obj$EResidualSq <- obj$trXtX +
#   #   -2 * sum((obj$X %*% obj$EF)*obj$EL) + sum(obj$ELtL * obj$EFtF)
#
#   # obj$L_KL <- nrow(obj$EL) * sum(log(1/sqrt(obj$L_ghat_tau))-1/2) +
#   #   -sum(log(sqrt(obj$EL2-obj$EL^2))) + sum(diag(obj$ELtL)/(2/obj$L_ghat_tau))
#
#   # obj$F_KL <- nrow(obj$EF) * sum(log(1/sqrt(obj$F_ghat_tau))-1/2) +
#   #   -sum(log(sqrt(obj$EF2-obj$EF^2))) + sum(diag(obj$EFtF)/(2/obj$F_ghat_tau))
#
#   # obj$elbo <- prod(dim(obj$X))/2 * (log(obj$tau) - log(2*pi))+
#   #   - obj$tau/2 * obj$EResidualSq +
#   #   - obj$L_KL - obj$F_KL
#
#   return(obj)
# }
