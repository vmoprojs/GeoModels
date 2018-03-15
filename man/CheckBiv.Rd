\name{CheckBiv}
\alias{CheckBiv}
\encoding{UTF-8}
\title{Checking Bivariate covariance models}
\description{
  The procedure control if the correlation model is bivariate.
}
\usage{
CheckBiv(numbermodel)
}
\arguments{
  \item{numbermodel}{numeric; the number associated to a given correlation model.}
}

\value{Returns TRUE or FALSE depending if the correlation model is bivariate or not.} 




\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales-Oñate, \email{victor.morales@uv.cl}
}




\examples{

library(GeoModels)
CheckBiv(CkCorrModel("Bi_matern_sep"))

}

\keyword{Composite}