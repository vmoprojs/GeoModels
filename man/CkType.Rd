\name{CkType}
\alias{CkType}
\encoding{UTF-8}
\title{Checking Likelihood Objects}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of likelihood objects inserted by the users.
}
\usage{
CkType(type)
}
\arguments{
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
}

\seealso{\code{\link{GeoFit}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}
}

\keyword{Composite}
