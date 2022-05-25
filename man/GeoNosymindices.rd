\name{GeoNosymindices}
\alias{GeoNosymindices}
\encoding{UTF-8}
\title{GeoNosymindices.}
\description{
 Given a matrix of indices and associated distances the function return  a matrix of indices 
 and associated distances, deleting the symmetric indices.
}
\usage{
GeoNosymindices(X,Y)
}
\arguments{
  \item{X}{A matrix of indices}
  \item{Y}{Associated distances}

}

\details{
 The function return the  matrix of indices and associated distances, deleting the symmetric indices. 
}


\value{
  Returns a list containing the following components:

  \item{xy}{Matrix of indices}
  \item{d}{Associated distance}
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\keyword{Scatterplot}
