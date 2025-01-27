\name{GeoNeighbSelect}
\alias{GeoNeighbSelect}
\encoding{UTF-8}
\title{A brute force algorithm for spatial or spatiotemoral optimal neighboord selection
for pairwise composite likelihood estimation.}
\description{
  The procedure performs different pairwise composite likelihood estimation using 
  user's specified  spatial or spatiotemporal neighboords in the weight function. 
  The neighbor minimizing the sum of the squared differences
   between the estimated semivariogram  and the  empirical semivariogram is selected. 
   The procedure needs an object obtained using the GeoVariogram function.
}
\usage{
GeoNeighbSelect(data, coordx, coordy=NULL,coordz=NULL, coordt=NULL, coordx_dyn=NULL,
    copula=NULL,corrmodel=NULL, distance="Eucl",fixed=NULL,anisopars=NULL,
    est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal',lower=NULL,
    neighb=c(1,2,3,4,5),maxtime=Inf, memdist=TRUE,model='Gaussian',
    n=1, ncores=NULL,optimizer='Nelder-Mead', parallel=FALSE, 
    bivariate=FALSE,radius=6371, start=NULL,type='Pairwise', upper=NULL, 
    weighted=FALSE,X=NULL,nosym=FALSE,spobj=NULL,spdata=NULL,vario=NULL)
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation)  or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on 
 regular grid) 
 or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)   or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid).
   For the description see the Section \bold{Details}.}
      \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix  or  (\eqn{d \times 3}{d x 3})-matrix 
     Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; Optional argument, the default is \code{NULL}.}
      \item{coordz}{A numeric vector giving 1-dimension of
    spatial coordinates;  Optional argument, the default is \code{NULL}.}

     \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates.  Optional argument, the default is \code{NULL}
    then a spatial random fields is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}}
  \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl}, the euclidean distance.
    See the Section \bold{Details}.}
  \item{fixed}{An optional named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated.}
  \item{anisopars}{ A list of two elements: "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
  \item{est.aniso}{ A bivariate logical vector  providing which  anisotropic parameters must be estimated.}
  
 
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
     of non-equispaced spatial sites (irregular grid).}
  
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default (see Section \bold{Details} in \link{GeoFit}).}
   
   \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb} or \code{bobyqa}  or \code{optimize}. The 
    names of the list must be the same of the names in the   \code{start} list.}
  
    
      \item{neighb}{Numeric; a vector of positive integers  indicating the 
    order of neighborhood  in the weight function of  composite likelihood 
    (see Section \bold{Details} in \link{GeoFit}).}
 
    \item{maxtime}{Numeric; an optional positive integer indicating the 
    order of temporal neighborhood  in the composite likelihood
    computation.}


    \item{memdist}{ Logical; if \code{TRUE} then  all  the distances  useful in the  composite likelihood 
      estimation are computed before  the optimization. \code{FALSE} is deprecated. }  
  
  \item{model}{String; the type of random fields and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section \bold{Details} in \link{GeoFit}.}
   
     \item{n}{Numeric; number of trials in a binomial random fields; number of successes in a negative binomial random fields }
     \item{ncores}{Numeric;  the number of cores involved in the parallelization }  
  
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). \code{Nelder-Mead} is the
    default. Other possible choices are \code{nlm}, \code{BFGS}, \code{SANN}, \code{L-BFGS-B} and \code{nlminb} and
    \code{bobyqa}.
    In these last three cases upper and lower bounds can be passed by the user. In the case
    of one-dimensional optimization, the function \code{optimize} is used.}  
  
  \item{parallel}{Logical; if \code{TRUE} the procedure is parallelized using dofuture.} 
  \item{bivariate}{Logical; if \code{TRUE} the bivariate case is considered.}  
 
   \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates.
   The default is 6371, the radius of the earth.}
 
  \item{start}{An optional named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see Section \bold{Details} in \link{GeoFit}).}
 
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods (see Section \bold{Details} in \link{GeoFit}).}
 
  \item{upper}{An optional named list giving the values  for the upper bound 
      of the space parameter when the optimizer is or \code{L-BFGS-B} or  \code{bobyqa} or   \code{nlminb} or \code{optimize}. 
       The names of the list must be the same of the names in the   \code{start} list.}  
 
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted (see Section \bold{Details} in \link{GeoFit}). If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
 
  \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}
  \item{nosym}{Logical; if TRUE simmetric weights are not considered. This allows a faster but less efficient CL estimation. }
 
   \item{spobj}{An object of class sp or spacetime}
   \item{spdata}{Character:The name of data in the sp or spacetime object}
   \item{vario}{An object of the class GeoVariogram obtained using the \link{GeoVariogram} function }
}
\details{
 The procedure performs different pairwise composite likelihood estimation using 
  user's specified  spatial or spatiotemporal neighboords in the weight function. 
  The neighbor  minimizing the sum of the squared differences between the estimated 
  semivariogram  and the  empirical semivariogram is selected.
    The procedure needs an object obtained using the GeoVariogram function.
}
\value{
  Returns a list with information on the best selected neighbor.
}
\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)


######### spatial case
set.seed(32)
N=500  # number of location sites
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
mean <- 0.2
# Set the covariance model's parameters:
corrmodel <- "Matern"
sill <- 1;nugget <- 0
scale <- 0.2/3;smooth=0.5

model="Gaussian"
param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)
# Simulation 
data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param,model=model)$data
I=Inf
fixed<-list(nugget=nugget)
start<-list(mean=mean,scale=scale,smooth=smooth,sill=sill)
lower<-list(mean=-I,scale=0,sill=0,smooth=0)
upper<-list(mean=I,scale=I,sill=I,smooth=I)

vario = GeoVariogram(coordx=coords,data=data,maxdist=0.3,numbins=15)

neighb=c(1,2,3,4) ## trying  different neighbs
selK <- GeoNeighbSelect(vario=vario,data=data,coordx=coords,corrmodel=corrmodel, 
                        model=model,neighb=neighb,
                        likelihood="Conditional",type="Pairwise",parallel=FALSE,
                        optimizer="nlminb",lower=lower,upper=upper,
                        start=start,fixed=fixed)
print(selK$best_neighb) ## selected neighbor




}
\keyword{Composite}
