\name{GeoCovDisplay}
\alias{GeoCovDisplay}
\encoding{UTF-8}
\title{Image plot displaying the pattern of the sparsness  of a  covariance matrix.}
\description{
  Image plot displaying the pattern of the sparsness  of a  covariance matrix.
}
\usage{
GeoCovDisplay(covmatrix,limits=FALSE,pch=2)
}
\arguments{
  \item{covmatrix}{An object of class matrix. See the Section \bold{Details}.}
  \item{limits}{Logical; If TRUE  and the covariance matrix is spatiotemporal or spatial bivariate
  then vertical and horizontal lines are added  to the image plot.}      
    \item{pch}{Type of symbols to use in the image plot.}
}

\details{ 
  For a given covariance matrix object (\code{\link{GeoCovmatrix}})
  the function diplays the  pattern of   the sparsness  of a  covariance matrix
  where the white color represents 0 entries and black color represents  non zero entries}
\value{Returns an image plot.}


\seealso{\code{\link{GeoCovmatrix}}}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\examples{

library(GeoModels)


  # Define the spatial-coordinates of the points:
x <- runif(100, 0, 2)
y <- runif(100, 0, 2)
coords=cbind(x,y)
matrix1 <- GeoCovmatrix(coordx=coords, corrmodel="GenWend", param=list(smooth=0,
                      power2=4,sill=1,scale=0.2,nugget=0))
 
GeoCovDisplay(matrix1)

}
\keyword{Sparsness pattern}
