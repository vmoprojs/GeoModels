\name{CheckST}
\alias{CheckST}
\encoding{UTF-8}
\title{Checking SpaceTime covariance models}
\description{
  The procedure control if the correlation model is spacetime.
}
\usage{
CheckST(numbermodel)
}
\arguments{
  \item{numbermodel}{numeric; the number associated to a given correlation model.}
}

\value{Returns TRUE or FALSE depending if the correlation model is spatiotemporal or not.} 



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}



\examples{

library(GeoModels)
CheckST(CkCorrModel("gneiting"))

}

\keyword{Composite}