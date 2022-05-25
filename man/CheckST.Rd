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

\details{ The function check if  the correlation model is spacetime.} 

\value{Returns TRUE or FALSE depending if the correlation model is spacetime or not.} 



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}



\examples{

library(GeoModels)
CheckST(CkCorrModel("gneiting"))

}

\keyword{Composite}