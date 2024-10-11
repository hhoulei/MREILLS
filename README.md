# MREILLS

***Installation***  
`devtools::install_github("hhoulei/MREILLS")`  

***Toy Example***  
`betaXG <- NULL`  
`sebetaXG <- NULL`  
`sebetaYG <- NULL`  
`for(i in 1:5){`  
`  betaXG <- cbind(betaXG,runif(100,0.05,0.2))`  
`  sebetaXG <- cbind(sebetaXG,runif(100,0.01,0.1))`  
`  sebetaYG <- cbind(sebetaYG,runif(100,0.01,0.1))`  
`}`  
`betaYG <- betaXG*0.2 + rnorm(100)`  
`race_name <- c('target1','target2','target3','target4','auxiliary')`  
`colnames(betaXG) <- race_name`  
`colnames(betaYG) <- race_name`  
`colnames(sebetaXG) <- race_name`  
`colnames(sebetaYG) <- race_name`  
`rho <- NULL`  
`meth <- 'CG'`  
`RESULT <- TEMR(betaXG,betaYG,sebetaXG,sebetaYG,rho,meth)`  
`RESULT`  
`# generate data`  
`#library(ggplot2)`  
`#library(ggridges)`  

`Nvar <- 8`  
`E <- 3`  
`g <- 100`  
`fdata_all <- list()`  
`for(oi in 1:E){`  

`  Sigma_episx <- diag(8)`  
`  diag(Sigma_episx) <- 0.01`  
`  episx <- MASS::mvrnorm(g,rep(0,Nvar),Sigma = Sigma_episx)`  

`  alpha <- matrix(runif(g*(Nvar+1),0.05,0.2),ncol=Nvar+1)`  
`  alpha[,Nvar+1] <- 0`  

`  beta <- matrix(runif((Nvar+1)*(Nvar+1),-1,1),ncol=Nvar+1)`  
`  diag(beta) <- 0`  

`  betaGX <- (alpha[,-(Nvar+1)] + episx) %*% solve(diag(Nvar)-beta[-(Nvar+1),-(Nvar+1)])`  
`  sebetaGX <- matrix(runif((Nvar)*g,0.01,0.05),nrow=g)`  

`  episy <- rnorm(g,0,0.01)`  
`  betaGY <- betaGX %*% beta[-(Nvar+1),(Nvar+1)]+alpha[,Nvar+1]+episy`  
`  betaGY <- c(betaGY)`  
`  sebetaGY <- runif(g,0.01,0.05)`  

`  rownames(betaGX) <- paste0('rs',1:g)`  
`  rownames(sebetaGX) <- paste0('rs',1:g)`  
`  names(betaGY) <- paste0('rs',1:g)`  
`  names(sebetaGY) <- paste0('rs',1:g)`  

`  data <- list(betaGX=betaGX,`  
`               sebetaGX=sebetaGX,`  
`               betaGY=betaGY,`  
`               sebetaGY=sebetaGY)`  

`  fdata_all[[oi]] <- data`  
`}`  
`res_eills <- MREILLS(fdata_all,r1=0.1,meth='L-BFGS-B',lambda=100)`  
`res_eills`  
`boot_eills <- BOOT_MREILLS(fdata_all,r1=0.1,meth='L-BFGS-B',lambda=100,numBoot=3)`  
`boot_eills`  
`res_lambda <- CHO_lambda(fdata_all)`  
`res_lambda$plot`  

***Citation***:  
Univariable and Multivariable Mendelian Randomization Method Integrating Multiple Heterogeneous GWAS Summary Datasets
Lei Hou, Xiaohua Zhou

Please contact houlei@bicmr.pku.edu.cn for any questions. We will continue to update this R package and reduce the problems that may be encountered during its installation.
