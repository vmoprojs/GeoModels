\name{CkVarType}
\alias{CkVarType}
\encoding{UTF-8}
\title{Checking Variance Estimates Type}
\description{
  Subroutine called by InitParam.
  The procedure controls the method used to compute the estimates' variances.
}
\usage{
CkVarType(type)
}
\arguments{
  \item{type}{String; the method used to compute the estimates'
    variances. If \code{SubSamp}  the estimates' variances
    are computed by the sub-sampling method, see \code{\link{GeoFit}}.}
}

\details{The procedure controls the method used to compute the estimates' variances}
\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\keyword{Composite}
