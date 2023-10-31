\name{Lik}
\alias{Lik}
\encoding{UTF-8}
\title{Optimizes the Log Likelihood}
\description{
  Subroutine called by GeoFit.
  The procedure estimates the model parameters by maximization of the
  log-likelihood.
}
\usage{
Lik(copula,bivariate,coordx,coordy,coordt,
coordx_dyn,corrmodel,data,fixed,flagcor,flagnuis,
           grid,lower,mdecomp,model,namescorr,
           namesnuis,namesparam,numcoord,
           numpairs,numparamcor,numtime,optimizer,
           onlyvar,parallel,param,radius,setup,
           spacetime,sparse,varest,taper,type,
           upper,ns,X,neighb,MM,aniso)
}
\arguments{
  \item{copula}{String; the type of copula. It can be "Beta" or "Gaussian"}
  \item{bivariate}{Logical; if \code{TRUE} then the data come froma a bivariate random field.
        Otherwise from a univariate random field.}
    \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) assigning 2-dimensions of spatial coordinates 
    or a numeric \eqn{d}{d}-dimensional vector assigning
    1-dimension of spatial coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, 
    the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL}
    then a spatial random field is expected.}
    \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }  
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{flagcor}{A numeric vector of flags denoting which correlation
    parameters have to be estimated.}
  \item{flagnuis}{A numeric verctor of flags denoting which nuisance
    parameters have to estimated.}
  \item{fixed}{A numeric vector of parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize}. The 
    names of the list must be the same of the names in the   \code{start} list.}
  \item{model}{Numeric; the id value of the density associated to the likelihood
    objects.}
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{namesparam}{String; the names of the parameters to be maximised.}
  \item{numcoord}{Numeric; the number of coordinates.}
  \item{numpairs}{Numeric; the number of pairs.}
  \item{numparamcor}{Numeric; the number of the correlation parameters.}
  \item{numtime}{Numeric; the number of temporal observations.}
   \item{mdecomp}{String; the type of matrix decomposition used in the simulation. Default is cholesky.
    The other possible choices is \code{svd} (Singular values decomposition).}     
         \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). \code{Nelder-Mead} is the
    default. Other possible choices are \code{nlm}, \code{BFGS}  \code{L-BFGS-B} and \code{nlminb}. 
    In these last two cases upper and lower bounds can be passed by the user. In the case
    of one-dimensional optimization, the function \code{optimize} is used.}  
   \item{onlyvar}{Logical; if \code{TRUE}  (and varest is  \code{TRUE}) only the variance covariance matrix is computed
  without optimizing. \code{FALSE} is the default.}  
  \item{parallel}{Logical; if \code{TRUE} optmization is performed using optimParallel 
   using the maximum number of cores,  when optimizer is  
          \code{L-BFGS-B}.\code{FALSE} is the default.} 
  \item{param}{A numeric vector of parameters.}
  \item{sparse}{Logical; if \code{TRUE} then  maximum likelihood is computed
  using sparse matrices algorithms.\code{FALSE} is the default.}
    \item{radius}{Numeric; the radius of the sphere when considering data on a sphere.}
    \item{ns}{Numeric: vector of number of location sites for each temporal instants}  
  \item{setup}{A List of useful components for the estimation based on the maximum tapered likelihood.}
  \item{spacetime}{Logical; if the random field is spatial
    (\code{FALSE}) or spatio-temporal (\code{TRUE}).}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{taper}{String; the name of the taper correlation function.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
\item{upper}{An optional named list giving the values  for the upper bound 
      of the space parameter when the optimizer is or \code{L-BFGS-B} or  \code{nlminb} or \code{optimize}. 
       The names of the list must be the same of the names in the   \code{start} list.}
    \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}
  \item{neighb}{Numeric;parameter for vecchia approximation using GPvecchia package}
    \item{MM}{Numeric;a non constant fixed mean }
           \item{aniso}{Logical; should anisotropy be considered?}
}

\value{Return a list from an \code{optim} call.}

\seealso{\code{\link{GeoFit}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\keyword{Composite}
