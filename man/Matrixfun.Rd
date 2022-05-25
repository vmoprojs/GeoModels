\name{MatSqrt, MatInv, MatLogDet}
\alias{MatSqrt}
\alias{MatInv}
\alias{MatLogDet}
\encoding{UTF-8}
\title{Square root, inverse  and  log determinant of a (semi)positive definite matrix, given a matrix decomposition.}
\description{
 Square root, inverse  and  log determinant of a (semi)positive definite matrix, given a matrix decomposition.
}
\usage{
MatSqrt(mat.decomp,method) 
MatInv(mat.decomp,method)
MatLogDet(mat.decomp,method)
}

\value{The function returna a  square root or inverse  or  log determinant of a (semi)positive definite matrix, given a matrix decomposition.}

\arguments{
  \item{mat.decomp}{numeric; a matrix decomposition.}
  \item{method}{string; the type of matrix decomposition. Two possible choices: \code{cholesky} and \code{svd}.}
}

\seealso{\code{\link{MatDecomp}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
 library(GeoModels)
 ################################################################
 ###
 ### Example 1. Inverse of Covariance matrix associated to
 ### a Matern correlation model
 ###
 ###############################################################
 # Define the spatial-coordinates of the points:
 x <- runif(15, 0, 1)
 y <- runif(15, 0, 1)
 coords <- cbind(x,y)
 # Matern Parameters
 param=list(smooth=0.5,sill=1,scale=0.2,nugget=0)
 a=matrix <- GeoCovmatrix(coordx=coords, corrmodel="Matern", param=param)

 ## decomposition with cholesky method  
 b=MatDecomp(a$covmat,method="cholesky")
 ## inverse of covariance matrix
 inverse=MatInv(b,method="cholesky")
}
\keyword{Composite}