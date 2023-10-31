\name{CkModel}
\alias{CkModel}
\encoding{UTF-8}
\title{Checking Random Field type}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of random field inserted by the users.
}
\usage{
CkModel(model)
}
\arguments{
   \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
}

\value{The function returns a numeric positive integer, or NULL if the model is invalid.}
\details{The function controls the type of random field inserted by the users.}
\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Composite}
