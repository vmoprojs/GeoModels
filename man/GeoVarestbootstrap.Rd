\name{GeoVarestbootstrap}  
\alias{GeoVarestbootstrap}
\encoding{UTF-8}
\title{Update a \code{GeoFit} object   using parametric bootstrap for std error estimation}
\description{
  The procedure update a \code{GeoFit} object  estimating stderr estimation using parametric bootstrap.}
\usage{GeoVarestbootstrap(fit,K=100,sparse=FALSE,
GPU=NULL,local=c(1,1),optimizer="Nelder-Mead",
                     lower=NULL, upper=NULL, 
                     method="cholesky",alpha=0.95, M=30,L=500,seed=1)}
\arguments{
  \item{fit}{A fitted object obtained from the
    \code{\link{GeoFit}}.}
     \item{K}{The number of simulations in the parametric bootstrap.}
       \item{sparse}{Logical; if \code{TRUE} then  cholesky decomposition is performed
  using sparse matrices algorithms (spam packake).}
       \item{GPU}{Numeric; if \code{NULL} (the default) 
      no OpenCL computation is performed. The user can choose the device to be used. Use \code{DeviceInfo()} function to see available devices, only double precision devices are allowed} 
        \item{local}{Numeric; number of local work-items of the OpenCL setup}
        \item{optimizer}{The type of optimization algorithm. See  \code{\link{GeoFit}} for details.}
        \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize}.}
       \item{upper}{An optional named list giving the values  for the upper bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize}.}
       \item{method}{String; The method of simulation. Default is \code{cholesky}. For large data set two options are \code{Vecchia} or \code{TB}}
         \item{alpha}{ Numeric; The level of the confidence interval. }
         \item{M}{Numeric; the number of neighboords in the Vecchia method.} 
      \item{L}{Numeric; the number of lines in the  turning band method.} 
        \item{seed}{Numeric; The seed used in the  n-fold  kriging cross-validation. Default is 1.}
     }
\details{ The function update a \code{GeoFit} object  estimating stderr estimation 
and confidence interval using parametric bootstrap.}
\value{  
  Returns an (updated) object of class \code{GeoFit}.
}


\seealso{\code{\link{GeoFit}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}




\keyword{Composite}
