\name{CkLikelihood}
\alias{CkLikelihood}
\encoding{UTF-8}
\title{Checking Composite-likelihood Type}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of the composite-likelihood inserted
  by the users.
}
\usage{
CkLikelihood(likelihood)
}
\arguments{
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default.}
}

\value{The function returns a numeric positive integer, or NULL if the likelihood is invalid.}

\details{  The function controls the type of the composite-likelihood inserted
  by the users.}
\seealso{\code{\link{GeoFit}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Composite}
