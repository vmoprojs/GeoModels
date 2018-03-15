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

\value{Returns 0,1,2 for  euclidean,geodesic, chordal distances 
respectively. Otherwise returns NULL.} 

\seealso{\code{\link{GeoFit}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales-Oñate, \email{victor.morales@uv.cl}
}

\keyword{Composite}