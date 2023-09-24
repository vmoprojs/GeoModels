\name{WlsStart}
\alias{WlsStart}
\encoding{UTF-8}
\title{Computes Starting Values based on Weighted Least Squares}
\description{
  Subroutine called by GeoFit.
  The function returns opportune starting values for the
  composite-likelihood fitting procedure based on weigthed least squares.
}
\usage{
WlsStart(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, fcall,
        fixed, grid,likelihood, maxdist, neighb,maxtime, model, n, param, 
        parscale, paramrange, radius,start, taper, tapsep, type, varest,
        vartype, weighted, winconst,winconst_t, winstp_t, winstp,copula,X,memdist,nosym)
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
  \item{corrmodel}{String; the name of a correlation model, for the
    description.}
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
    \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhood  in the composite likelihood
    computation. See the Section \bold{Details} for more information.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal separation considered in the composite-likelihood
    computation.} 
  \item{model}{String; the name of the model. Here the default is
    \code{NULL}.}
    \item{n}{Numeric; number of trials in a binomial random field.}  
  \item{param}{A numeric vector of parameter values required in the
    simulation procedure of random fields.}
  \item{parscale}{A numeric vector with scaling values for improving the
    maximisation routine.}
  \item{paramrange}{A numeric vector with the range of the parameter
    space.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere when using the great circle distance.
          Default value is the radius of the earth in Km (i.e. 6371)}   
  \item{start}{A numeric vector with starting values.}
   \item{taper}{String; the name of the type of covariance matrix.
     It can be \code{Standard} (the default value) or
    \code{Tapering} for taperd covariance matrix.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the space time quasi taper (see \bold{Details}).}
  \item{type}{String; the type of estimation method.}
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
         \item{copula}{The type of copula.}
    \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}  
        \item{memdist}{ Logical; if \code{TRUE} then  the distances in the  composite likelihood 
      are computed before  the optmization. }
      \item{nosym}{Logical; if TRUE two simmetric weights are not considered }
}

\value{A list with  a set of useful informations in the estimation procedure.}


\details{Internal function called by  \code{\link{GeoFit}}.}

\seealso{\code{\link{GeoFit}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\keyword{LeastSquare}
