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

\value{The function returns a numeric positive integer, or NULL if the type of likelihood is invalid.}

\details{The procedure checks the  likelihood Object}


\seealso{\code{\link{GeoFit}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Composite}
