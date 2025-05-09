% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BOOT_MREILLS.R
\name{BOOT_MREILLS}
\alias{BOOT_MREILLS}
\title{BOOT_MREILLS}
\usage{
BOOT_MREILLS(fdata_all, r1, meth, lambda, numBoot)
}
\arguments{
\item{fdata_all}{A list of GWAS summary statistics for E populations,
each element is a list consists of betaGX, sebetaGX, betaGY and sebetaGY.
betaGX is nSNP*nexposure matrix of G-X association (beta coefficient),
sebetaGX is the corresponding standard error matrix,
betaGY is vector of G-Y association (beta coefficient),
sebetaGY is the corresponding standard error vector.}

\item{r1}{hyper-parameter, gamma in the MR-EILLS model.}

\item{meth}{Optimization method to use, including "Nelder-Mead", "BFGS",
"CG", "L-BFGS-B", "SANN" and "Brent".}

\item{lambda}{hyper-parameter, lambda in the MR-EILLS model.}

\item{numBoot}{the number of bootstrap.}
}
\value{
A vector of standard error for causal effect estimation of
each exposure on the outcome.
}
\description{
Bootstrap for MREILLS.
}
