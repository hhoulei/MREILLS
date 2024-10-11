#' @title BOOT_MREILLS
#'
#' @description Bootstrap for MREILLS.
#'
#' @param fdata_all A list of GWAS summary statistics for E populations,
#' each element is a list consists of betaGX, sebetaGX, betaGY and sebetaGY.
#' betaGX is nSNP*nexposure matrix of G-X association (beta coefficient),
#' sebetaGX is the corresponding standard error matrix,
#' betaGY is vector of G-Y association (beta coefficient),
#' sebetaGY is the corresponding standard error vector.
#'
#' @param r1 hyper-parameter, gamma in the MR-EILLS model.
#' @param lambda hyper-parameter, lambda in the MR-EILLS model.
#'
#' @param meth Optimization method to use, including "Nelder-Mead", "BFGS",
#' "CG", "L-BFGS-B", "SANN" and "Brent".
#'
#' @param numBoot the number of bootstrap.
#'
#'
#' @return A vector of standard error for causal effect estimation of
#' each exposure on the outcome.
#'
#' @export
#'
BOOT_MREILLS <- function(fdata_all,r1,meth,lambda,numBoot){


  estboot <- NULL
  for(ojo in 1:numBoot){

    locboot <- sample(1:nrow(fdata_all[[1]]$betaGX),replace = T)

    fdata_boot <- list()
    for(obo in 1:length(fdata_all)){
      fdata_once <- fdata_all[[obo]]
      fdata_boot[[obo]] <- list(betaGX=fdata_once$betaGX[locboot,],
                                sebetaGX=fdata_once$sebetaGX[locboot,],
                                betaGY=fdata_once$betaGY[locboot],
                                sebetaGY= fdata_once$sebetaGY[locboot])
    }

    x.inv <- try(MREILLS(fdata_boot,r1,meth,lambda),
                 silent = T)
    if('try-error' %in% class(x.inv)){
      next
    }else{
      estboot <- cbind(estboot,x.inv$par)
    }

  }

  se <- apply(estboot,1,function(x) sd(x,na.rm = T))
  return(se)
}
