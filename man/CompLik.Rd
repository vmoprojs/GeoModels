\name{CompLik}
\alias{CompLik}
\encoding{UTF-8}
\title{Optimizes the Composite log-likelihood}
\description{
  Subroutine called by GeoFit.
  The procedure estimates the model parameters by maximisation of the
  composite log-likelihood.
}
\usage{
CompLik(bivariate, coordx, coordy ,coordt, coordx_dyn, corrmodel, data, distance,  
               flagcorr, flagnuis, fixed, GPU, grid,likelihood,local,lower, 
                 model, n, namescorr, namesnuis, namesparam,
               numparam, numparamcorr, optimizer, onlyvar, param,
               spacetime, type,upper, varest, vartype,
               weigthed, winconst, winstp,winconst_t, winstp_t, ns, X)
}
\arguments{
 \item{bivariate}{Logical; if \code{TRUE} then the data come froma a bivariate random field. 
                Otherwise from a univariate random field.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates 
    or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL}
    then a spatial random field is expected.}  
    \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl},
            the euclidean distance. See the Section \bold{Details}.}
  \item{flagcorr}{A numeric vector of binary values denoting which
    paramerters of the correlation function will be estimated.}
  \item{flagnuis}{A numeric vector of binary values denoting which
    nuisance paramerters will be estimated.}
  \item{fixed}{A numeric vector of parameters that
    will be considered as known values.}
  \item{GPU}{Numeric; if \code{NULL} (the default) 
      no GPU computation is performed.} 

  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}

  \item{likelihood}{String; the configuration of the
    compositelikelihood, see \code{\link{GeoFit}}.}

 \item{local}{Numeric; number of local work-items of the GPU}
  \item{lower}{A numeric vector with the lower bounds of the parameters' ranges.}
  \item{model}{Numeric; the id value of the density associated to the likelihood
    objects.}
   \item{n}{Numeric; number of trials in a binomial random fields.}  
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{namesparam}{String; the names of the parameters to be maximised.}
  \item{numparam}{Numeric; the number of parameters to be maximised.}
  \item{numparamcorr}{Numeric; the number of correlation parameters.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{onlyvar}{Logical; if \code{TRUE}  (and varest is  \code{TRUE}) only
   the variance covariance matrix is computed
  without optimizing. \code{FALSE} is the default.}   
  \item{param}{A numeric vector of parameters' values.}
  \item{spacetime}{Logical; if \code{TRUE} the random field is
    spatial-temporal otherwise is a spatial field.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
  \item{upper}{A numeric vector with the upper bounds of the parameters'
    ranges.}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
  the estimate variances, see \code{\link{GeoFit}}.}
  \item{weigthed}{Logical; if \code{TRUE} then decreasing weigths coming 
  from a compactly supported correlation function
  with compact support \code{maxdist} (\code{maxtime})are used.}
  \item{winconst}{Numeric; a positive value for computing the spatial  sub-window
      in the sub-sampling procedure.}
 \item{winstp}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the spatial  sub-sampling procedure.}
   \item{winconst_t}{Numeric; a positive value for computing the temporal sub-window
         in the sub-sampling procedure.}
 \item{winstp_t}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the temporal  sub-sampling procedure.}
   \item{ns}{Numeric; Number of (dynamical) temporal instants. }   
      \item{X}{Numeric; Matrix of space-time covariates in the linear mean specification.}
}

\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}
}

\keyword{Composite}
