\name{GeoVarestbootstrap}  
\alias{GeoVarestbootstrap}
\encoding{UTF-8}
\title{Computes stderr estimation using parametric bootstrap}
\description{
  The procedure return stderr estimation for compositelikelihood estimation}
\usage{GeoVarestbootstrap(fit,K=100,sparse=FALSE,GPU=NULL,local=c(1,1))}
\arguments{
  \item{fit}{A fitted object obtained from the
    \code{\link{GeoFit}}.}
     \item{K}{The number of simulations in the parametric bootstrap.}
       \item{sparse}{Logical; if \code{TRUE} then  cholesky decomposition is performed
  using sparse matrices algorithms (spam packake).}
       \item{GPU}{Numeric; if \code{NULL} (the default) 
      no OpenCL computation is performed. The user can choose the device to be used. Use \code{DeviceInfo()} function to see available devices, only double precision devices are allowed} 
        \item{local}{Numeric; number of local work-items of the OpenCL setup}
}

\value{  
  \item{stderr}{Stderr estimation with parametric bootstrap}
  \item{varcov}{The matrix of the variance-covariance of the estimates;}
  \item{claic}{The estimated composite information criterion}
   \item{clbic}{The estimated bayesian composite information criterion}
}


\seealso{\code{\link{GeoFit}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}
}




\keyword{Composite}
