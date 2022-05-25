\name{CompIndLik2}
\alias{CompIndLik2}
\encoding{UTF-8}
\title{Optimizes the Composite indipendence log-likelihood}
\description{
  Subroutine called by GeoFit.
  The procedure estimates the model parameters by maximisation of the indipendence 
  composite log-likelihood.
}
\usage{
CompIndLik2(bivariate, coordx, coordy ,coordt,
coordx_dyn, data, flagcorr, flagnuis, fixed,grid,
              lower, model, n, namescorr, namesnuis, 
              namesparam,
              numparam, optimizer, onlyvar, parallel, 
              param, spacetime, type,
              upper, namesupper, varest, ns, X,
              sensitivity,copula,MM)
}
\arguments{
 \item{bivariate}{Logical; if \code{TRUE} then the data come froma a bivariate random field. 
                Otherwise from a univariate random field.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates 
    or a numeric vector assigning 1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL}
    then a spatial random field is expected.}
    \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL} }
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{flagcorr}{A numeric vector of binary values denoting which
    paramerters of the correlation function will be estimated.}
  \item{flagnuis}{A numeric vector of binary values denoting which
    nuisance paramerters will be estimated.}
  \item{fixed}{A numeric vector of parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}

   \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize}. The 
    names of the list must be the same of the names in the   \code{start} list.}
  \item{model}{Numeric; the id value of the density associated to the likelihood objects.}
  \item{n}{Numeric; number of trials in a binomial random fields.}  
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{namesparam}{String; the names of the parameters to be maximised.}
  \item{numparam}{Numeric; the number of parameters to be maximised.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). \code{Nelder-Mead} is the
    default. Other possible choices are \code{ucminf},\code{nlm}, \code{BFGS}  \code{L-BFGS-B} and \code{nlminb}. 
    In these last two cases upper and lower bounds can be passed by the user. In the case
    of one-dimensional optimization, the function \code{optimize} is used.}  
  \item{onlyvar}{Logical; if \code{TRUE}  (and varest is  \code{TRUE}) only
   the variance covariance matrix is computed
  without optimizing. \code{FALSE} is the default.}  
   \item{parallel}{Logical; if \code{TRUE} optmization is performed using optimParallel 
   using the maximum number of cores,  when optimizer is 
          \code{L-BFGS-B}.\code{FALSE} is the default.}  
  \item{param}{A numeric vector of parameters values.}
  \item{spacetime}{Logical; if \code{TRUE} the random field is spatial-temporal otherwise is a spatial field.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by pairwise marginal likelihoods.}
  \item{upper}{An optional named list giving the values  for the upper bound 
      of the space parameter when the optimizer is or \code{L-BFGS-B} or  \code{nlminb} or \code{optimize}. 
       The names of the list must be the same of the names in the   \code{start} list.}
   \item{namesupper}{String; the names of the upper limit of the parameters.}
  \item{varest}{Logical; if \code{TRUE} the estimate variances and standard errors are returned.
    \code{FALSE} is the default.}
   \item{ns}{Numeric; Number of (dynamical) temporal instants. }   
  \item{X}{Numeric; Matrix of space-time covariates in the linear mean specification.}
  \item{sensitivity}{Logical; if \code{TRUE} then  the sensitivy  matrix is computed}
   \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
          \item{MM}{Numeric;a non constant fixed mean }
}

\value{Return a list from an \code{optim} call.}

\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Composite}
