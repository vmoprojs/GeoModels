\name{SimCE}
\alias{SimCE}
\encoding{UTF-8}
\title{Circulant embeeding simulation}
\description{
  Subroutine called by GeoSimapprox.
  The procedure return a simulation on a regular grid from a standard spatial Gaussian random field
  with a specified correlation model
}
\usage{
SimCE(M,N,x,y,z,corrmodel,param,mean.val, max.ext)
}
\arguments{
  \item{M}{Numeric; the dimension of x}
  \item{N}{Numeric; the dimension of y}
  \item{x}{A numeric \eqn{M}{M}-dimensional vector giving 1-dimension of spatial coordinates.}
  \item{y}{A numeric \eqn{N}{N}-dimensional vector giving 1-dimension of spatial coordinates.}
    \item{z}{A numeric \eqn{N}{N}-dimensional vector giving 1-dimension of spatial coordinates.}
  \item{corrmodel}{String; the name of a correlation model.}
  \item{param}{A list of parameter values required in the simulation procedure.}
  \item{mean.val}{The mean of the random field.}
  \item{max.ext}{The maximum extension of the simulation window.}
}

\value{Return a list from an \code{optim} call.}

\seealso{\code{\link{GeoSimapprox}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Simulation}
