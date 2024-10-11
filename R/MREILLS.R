#'
#' @title MREILLS
#'
#' @description Univariable and Multivariable Mendelian Randomization Method
#' Integrating Multiple Heterogeneous GWAS Summary Datasets (Lei Hou)
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
#' @importFrom MASS mvrnorm
#' @importFrom ggplot2 ggplot
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggridges theme_ridges
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#'
#' @return A list object, the results of causal effect estimation of
#' each exposure on the outcome.
#'
#' @examples
#' # generate data
#' #library(ggplot2)
#' #library(ggridges)
#'
#' Nvar <- 8
#' E <- 3
#' g <- 100
#' fdata_all <- list()
#' for(oi in 1:E){
#'
#'   Sigma_episx <- diag(8)
#'   diag(Sigma_episx) <- 0.01
#'   episx <- MASS::mvrnorm(g,rep(0,Nvar),Sigma = Sigma_episx)
#'
#'   alpha <- matrix(runif(g*(Nvar+1),0.05,0.2),ncol=Nvar+1)
#'   alpha[,Nvar+1] <- 0
#'
#'   beta <- matrix(runif((Nvar+1)*(Nvar+1),-1,1),ncol=Nvar+1)
#'   diag(beta) <- 0
#'
#'   betaGX <- (alpha[,-(Nvar+1)] + episx) %*% solve(diag(Nvar)-beta[-(Nvar+1),-(Nvar+1)])
#'   sebetaGX <- matrix(runif((Nvar)*g,0.01,0.05),nrow=g)
#'
#'   episy <- rnorm(g,0,0.01)
#'   betaGY <- betaGX %*% beta[-(Nvar+1),(Nvar+1)]+alpha[,Nvar+1]+episy
#'   betaGY <- c(betaGY)
#'   sebetaGY <- runif(g,0.01,0.05)
#'
#'   rownames(betaGX) <- paste0('rs',1:g)
#'   rownames(sebetaGX) <- paste0('rs',1:g)
#'   names(betaGY) <- paste0('rs',1:g)
#'   names(sebetaGY) <- paste0('rs',1:g)
#'
#'   data <- list(betaGX=betaGX,
#'                sebetaGX=sebetaGX,
#'                betaGY=betaGY,
#'                sebetaGY=sebetaGY)
#'
#'   fdata_all[[oi]] <- data
#' }
#' res_eills <- MREILLS(fdata_all,r1=0.1,meth='L-BFGS-B',lambda=100)
#' res_eills
#' boot_eills <- BOOT_MREILLS(fdata_all,r1=0.1,meth='L-BFGS-B',lambda=100,numBoot=3)
#' boot_eills
#' res_lambda <- CHO_lambda(fdata_all)
#' res_lambda$plot
#'
#' @export
#'
MREILLS <- function(fdata_all,r1,meth,lambda){


  Nv <- ncol(fdata_all[[1]]$betaGX)
  Ng <- nrow(fdata_all[[1]]$betaGX)
  Qloss <- function(bb) {

    # IV selection
    QSj <- NULL
    for(sj in 1:Ng){
      QSj1 <- NULL
      Xe <- NULL
      for(sje in 1:length(fdata_all)){
        oncee <- fdata_all[[sje]]
        QSj1 <- c(QSj1,
                  oncee$betaGY[sj]-matrix(oncee$betaGX[sj,],nrow=1)%*%bb)
        Xe <- cbind(Xe,oncee$betaGX[sj,])
      }

      Xe <- c(abs(Xe)%*%abs(QSj1))
      QSj <- c(QSj,sum(abs(QSj1))+sum(Xe))
    }

    Sj <- (QSj<lambda)


    #bb <- c(0,0,0,0,0,0,0)
    Re <- unlist(lapply(fdata_all,function(x) mean((x$sebetaGY[Sj]^(-2)/sum(x$sebetaGY[Sj]^(-2)))*((x$betaGY[Sj]-x$betaGX[Sj,]%*%bb)^2))))

    Jpe <- NULL
    weight <- NULL
    for(oe in 1:length(fdata_all)){
      once <- fdata_all[[oe]]
      epi <- once$betaGY[Sj]-once$betaGX[Sj,]%*%bb

      Jpe_once <- NULL
      for(oj in 1:Nv){
        Jpe_once <- c(Jpe_once,mean(once$betaGX[Sj,oj]*epi)^2)
      }
      Jpe <- cbind(Jpe,Jpe_once)

      weight <- c(weight,sum(once$sebetaGY[Sj]^(-2)))
    }

    weight <- weight/sum(weight)

    Jp <- Jpe%*%weight
    J <- sum(Jp)

    R <- sum(weight*Re)

    return(R+r1*J)
  }

  result <- optim(par = rep(0,Nv),
                  fn = Qloss,
                  method = meth
  )
  return(result)
}
