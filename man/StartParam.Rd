\name{StartParam}
\alias{StartParam}
\encoding{UTF-8}
\title{Initializes the Parameters for Estimation Procedures}
\description{
  Subroutine called by the fitting procedures.
  The procedure initializes the parameters for the fitting procedure.
}
\usage{
StartParam(coordx, coordy, coordt,coordx_dyn, corrmodel, data, distance, fcall, 
          fixed, grid,likelihood,  maxdist, maxtime, model, n,  param, 
          parscale, paramrange, radius,  start, taper, tapsep,
          type,typereal,varest, vartype, weighted, winconst,
          winstp, winconst_t, winstp_t,X)
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
    \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }  
  \item{corrmodel}{String; the name of a correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl},
            the euclidean distance. See the Section \bold{Details}.}
  \item{fcall}{String; "fitting" to call the fitting procedure and
    "simulation" to call the simulation procedure.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{likelihood}{String; the configuration of the composite
    likelihood.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite-likelihood
    computation.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal lag considered in the composite-likelihood
    computation.} 
  \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates.
   The default is 6371, the radius of the earth.}  
  \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
   \item{n}{Numeric; number of trials for binomial random fields.}  
  \item{param}{A numeric vector of parameter values required in the
    simulation procedure of random fields.}
  \item{parscale}{A numeric vector of scaling factor to improve the
    maximizing procedure, see \code{\link{optim}}.}
  \item{paramrange}{A numeric vector of parameters ranges, see
    \code{\link{optim}}.}
   \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure.}
  \item{taper}{String; the name of the type of covariance matrix.
     It can be   \code{Standard} (the default value) or
    \code{Tapering} for taperd covariance matrix.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the space time 
   adaptive taper (see \bold{Details}).}
  \item{type}{String; the type of likelihood objects. Temporary value set
    to be "WLeastSquare" (weigthed least-square) in order to compute the starting values.}
  \item{typereal}{String; the real type of likelihood objects. See \code{\link{GeoFit}}.}
  \item{varest}{Logical; if \code{TRUE} the estimates' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
    the estimate variances, see the Section \bold{Details}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see \code{\link{GeoFit}}.}
       \item{winconst}{Numeric; a positive value for computing the spatial  sub-window
      in the sub-sampling procedure.}
 \item{winstp}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the spatial  sub-sampling procedure.}
   \item{winconst_t}{Numeric; a positive value for computing the temporal sub-window
         in the sub-sampling procedure.}
 \item{winstp_t}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the temporal  sub-sampling procedure.}
    \item{X}{Numeric; Matrix of space-time covariates.}   
}

\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}
}

\keyword{Composite}
