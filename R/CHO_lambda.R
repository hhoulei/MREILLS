#' @title CHO_lambda
#'
#' @description plot ridge plot for choosing lambda.
#'
#' @param fdata_all A list of GWAS summary statistics for E populations,
#' each element is a list consists of betaGX, sebetaGX, betaGY and sebetaGY.
#' betaGX is nSNP*nexposure matrix of G-X association (beta coefficient),
#' sebetaGX is the corresponding standard error matrix,
#' betaGY is vector of G-Y association (beta coefficient),
#' sebetaGY is the corresponding standard error vector.
#'
#'
#' @return the ridge plot.
#'
#' @export
#'

CHO_lambda <- function(fdata_all){

  Nv <- ncol(fdata_all[[1]]$betaGX)
  SNPall <- list()
  for(i in 1:length(fdata_all)){
    SNPall[[i]] <- rownames(fdata_all[[i]]$betaGX)
  }
  SNPall1 <- unique(unlist(SNPall))

  bb <- rep(0,Nv)
  QSj <- NULL

  for(sj in 1:length(SNPall1)){
    SNP_once <- SNPall1[sj]
    QSj1 <- NULL
    QSj2 <- NULL
    for(je in 1:length(fdata_all)){
      oncee <- fdata_all[[je]]
      loc <- which(rownames(oncee$betaGX) == SNP_once)
      if(length(loc)>0){
        epi1 <- oncee$betaGY[loc]-matrix(oncee$betaGX[loc,],nrow=1)%*%bb
        epi1 <- c(epi1)
        QSj1 <- c(QSj1,abs(epi1))
        QSj2 <- c(QSj2,sum(abs(epi1*oncee$betaGX[loc,])))
      }
    }
    QSj <- c(QSj, sum(QSj1)+sum(QSj2))
  }

  dt <- data.frame(QSj=QSj,
                   y=rep(1,length(QSj)))
  dt$y <- as.factor(dt$y)

  p1 <- ggplot(dt,aes(x = QSj, y= y,fill =  y)) +
    geom_density_ridges(scale=1)+
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0))+
    theme_ridges()+
    theme()+
    theme(panel.background = element_rect(fill = 'transparent'))+
    theme(legend.position = "none")+
    ylab('')+
    xlab('')
  #p1

  return(list(QSj=QSj,
              plot=p1))

}

