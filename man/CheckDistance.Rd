\name{CheckDistance}
\alias{CheckDistance}
\encoding{UTF-8}
\title{Checking Distance}
\description{
  The procedure controls the type of distance.
}
\usage{
CheckDistance(distance)
}
\arguments{
  \item{distance}{String; the type of distance, for the
    description see \code{\link{GeoCovmatrix}}.
    Default is \code{Eucl}. Other possible values are \code{Geod} 
    and \code{Chor} that is euclidean, geodesic and chordal distance.}
}

\details{ The function check if  the type of distance is valid.} 


\value{Returns 0,1,2 for  euclidean,geodesic, chordal distances 
respectively. Otherwise returns NULL.} 





\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Composite}