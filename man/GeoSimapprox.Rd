\name{GeoSimapprox}
\alias{GeoSimapprox}
\encoding{UTF-8}
\title{Fast simulation of Gaussian and non Gaussian Random Fields.}
\description{
  Simulation of Gaussian and some non Gaussian spatial, spatio-temporal and spatial bivariate 
 random fields using approximate methods of simulation (circulant embeeding and turning bands) (see Examples). 
}
\usage{
GeoSimapprox(coordx, coordy=NULL, coordt=NULL, 
coordx_dyn=NULL,corrmodel, distance="Eucl",GPU=NULL, 
grid=FALSE,local=c(1,1),max.ext=1,
method="TB", L=2000,model='Gaussian',
n=1,param,anisopars=NULL, radius=6371,X=NULL,spobj=NULL,nrep=1)
}
\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates.
     Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates.  Optional argument, the default is \code{NULL}
    then a spatial RF is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details}  of \code{\link{GeoFit}}.}
    \item{GPU}{Numeric; if \code{NULL} (the default) 
      no GPU computation is performed.} 

  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
          \item{local}{Numeric; number of local work-items of the GPU}
      \item{max.ext}{Numeric; The maximum extension of the simulation window (for the spatial CE method).}
     \item{method}{String; the type of approximation method. Default is \code{TB}  that is the 
     turning band method. The other possible choice is 
     and \code{CE} (circular embeeding).} 

      \item{L}{Numeric; the number of lines in the  turning band method.} 
  \item{model}{String; the type of RF and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{n}{Numeric; the number of trials for binomial RFs. The number of successes in the negative Binomial RFs. Default is \eqn{1}.} 
  \item{param}{A list of parameter values required in the
    simulation procedure of RFs, see \bold{Examples}.}
      \item{anisopars}{ A list of two elements "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
    \item{radius}{Numeric; a value indicating  the radius of the sphere when using the great 
    circle distance. Default value is the radius of the earth in Km (i.e. 6371)}  
  \item{X}{Numeric; Matrix of space-time covariates.}
    \item{spobj}{An object of class sp or spacetime}
    \item{nrep}{Numeric; Numbers of indipendent replicates.}
}


\value{
  Returns an object of class \code{GeoSim}.
  An object of class \code{GeoSim} is a list containing
  at most the following components:
  
  \item{bivariate}{Logical:\code{TRUE} if the Gaussian RF is bivariate, otherwise  \code{FALSE};}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{coordx_dyn}{A  list of dynamical (in time) spatial coordinates;}
  \item{corrmodel}{The correlation model; see \code{\link{GeoCovmatrix}}.}
  \item{data}{The vector or matrix or array of data, see
    \code{\link{GeoFit}};}
  \item{distance}{The type of spatial distance;}
  \item{method}{The method of simulation}
  \item{model}{The type of RF, see \code{\link{GeoFit}}.}
  \item{n}{The number of trial for Binomial RFs;the number of successes in a negative Binomial RFs;} 
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number the temporal realisations of the RF;}
  \item{param}{The vector of  parameters' estimates;}
  \item{radius}{The radius of the sphere if coordinates are passed in lon/lat format;}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if
  spatial RF;}
    \item{nrep}{The number of indipendent replicates;}
}

\references{

T. Gneiting, H. Sevcikova, D. B. Percival, M. Schlather and Y. Jiang (2006)
Fast and Exact Simulation of Large Gaussian Lattice Systems in R2: Exploring the Limits
\emph{Journal of Computational and Graphical Statistics} 15 (3)

D. Arroyo, X. Emery (2020) An R Implementation of a Continuous Spectral Algorithm for Simulating Vector Gaussian Random Fields in Euclidean Spaces \emph{ACM Transactions on Mathematical Software}  47(1)

J. Guiness (2018) Permutation and Grouping Methods for Sharpening Gaussian Process Approximations \emph{Technometrics} 60(4) 415-429.

}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)


################################################################
###
### Example 1. Simulation of a large spatial Gaussian RF 
###            with  Matern  covariance model
###            using circulant embeeding method
###            It works only for regular grid
###            for any typpe of correlation model
###############################################################
set.seed(68)
x = seq(0,1,0.006)
y = seq(0,1,0.006)
param=list(smooth=1.5,mean=0,sill=1,scale=0.2/3,nugget=0)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSimapprox(coordx=x,coordy=y, grid=TRUE,corrmodel="Matern", model="Gaussian",
                      method="CE",param=param)$data
fields::image.plot( matrix(data1, length(x), length(y), byrow = TRUE) )

################################################################
###
### Example 2. Simulation of a large spatial Gaussian RF 
###            with  Matern  covariance model
###            using Turning band method
###            It works for (ir)regular grid
###            for the Matern model
###############################################################
set.seed(68)
x = runif(8000)
y = runif(8000)
coords=cbind(x,y)
param=list(smooth=0.5,mean=0,sill=1,scale=0.1,nugget=0)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSimapprox(coords, corrmodel="Matern", model="Gaussian",
                      method="TB",param=param)$data
quilt.plot(coords,data1)



################################################################
###
### Example 3. Simulation of a large spacetime Gaussian RF 
###            with separable matern  covariance model
###            using  Circular embeeding method
###            It works  for (large) regular time grid
###############################################################
set.seed(68)
coordt <- (0:100)
coords <- cbind( runif(100, 0 ,1), runif(100, 0 ,1))
param <- list(mean  = 0, sill = 1, nugget = 0.25,
              scale_s = 0.05, scale_t = 2, 
              smooth_s = 0.5, smooth_t = 0.5)
# Simulation of a spatial Gaussian RF with Matern correlation function
param<-list(nugget=0,mean=0,scale_s=0.2/3,scale_t=2/3,sill=1,smooth_s=0.5,smooth_t=0.5)

data <- GeoSimapprox(coordx=coords, coordt=coordt, corrmodel="Matern_Matern",
                     model="Gaussian",method="CE",param=param)$data
dim(data)

################################################################
###
### Example 4. Simulation of a large spacetime Gaussian RF 
###            with separable GenWend covariance model
###            using  Circular embeeding method in time
###############################################################
set.seed(68)
# Simulation of a spatial Gaussian RF with Matern correlation function
param<-list(nugget=0,mean=0,scale_s=0.2,scale_t=3,sill=1,
             smooth_s=0,smooth_t=0, power2_s=4,power2_t=4)

data <- GeoSimapprox(coordx=coords, coordt=coordt, corrmodel="GenWend_GenWend",
                     model="Gaussian",method="CE",param=param)$data
dim(data)


################################################################
###
### Example 6. Simulation of a large bivariate Gaussian RF
### with  bivariate Matern correlation model
###
###############################################################

# Define the spatial-coordinates of the points:
#x <- runif(20000, 0, 2)
#y <- runif(20000, 0, 2)
#coords <- cbind(x,y)

# Simulation of a bivariate spatial Gaussian RF:
# with a  Bivariate Matern
#set.seed(12)
#param=list(mean_1=4,mean_2=2,smooth_1=0.5,smooth_2=0.5,smooth_12=0.5,
#           scale_1=0.12,scale_2=0.1,scale_12=0.15,
#           sill_1=1,sill_2=1,nugget_1=0,nugget_2=0,pcol=0.5)
#data <- GeoSimapprox(coordx=coords,corrmodel="Bi_matern",
#              param=param,method="TB",L=1000)$data
#opar=par(no.readonly = TRUE)
#par(mfrow=c(1,2))
#quilt.plot(coords,data[1,],col=terrain.colors(100),main="1",xlab="",ylab="")
#quilt.plot(coords,data[2,],col=terrain.colors(100),main="2",xlab="",ylab="")







}
\keyword{Simulation}
