\name{CkInput}
\alias{CkInput}
\encoding{UTF-8}
\title{Checking Input}
\description{
  Subroutine called by the fitting procedures.
  The procedure controls the the validity of the input inserted by the users.
}
\usage{
CkInput(coordx, coordy, coordt, coordx_dyn,corrmodel, data, distance, 
           fcall, fixed, grid,likelihood, maxdist, maxtime, 
            model, n,  optimizer, param, radius,
           start, taper, tapsep,  type, varest, vartype, 
           weighted,X)
}
\arguments{
   \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see \code{\link{GeoFit}}.}
    \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }  
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl},
               the euclidean distance. See the Section \bold{Details}.}
  \item{fcall}{String; \code{Fitting} to call the fitting procedure and
    \code{simulation} to call the simulation.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite-likelihood
    computation.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal lag separation in the composite-likelihood.}
   \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates.
   The default is 6371, the radius of the earth.}  
  \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
  \item{n}{Numeric; the number of trials in a binomial random fields.
         Default is \eqn{1}.}    
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{param}{A numeric vector of parameters, needed only in
    simulation. See \code{\link{GeoSim}}.}
  \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default.}
  \item{taper}{String; the name of the tapered correlation
    function.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the space time quasi taper (see \bold{Details}).}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
  the estimate variances, see \code{\link{GeoFit}}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted. If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
  \item{X}{Numeric; Matrix of space-time covariates in the linear mean specification.}
}

\seealso{\code{\link{GeoFit}}}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\keyword{Composite}
