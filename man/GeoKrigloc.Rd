\name{GeoKrigloc}
\alias{GeoKrigloc}
\encoding{UTF-8}
\title{Spatial (bivariate) and spatio temporal optimal linear  local  prediction 
          for Gaussian and non Gaussian RFs.}
\description{
For a  given set of  spatial location sites and temporal instants,
the function computes optimal linear prediction and associated mean square error 
for the Gaussian and non Gaussian case using a spatial (temporal) neighborhood
computed using the  function \code{\link{GeoNeighborhood}}
}
\usage{
GeoKrigloc(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, 
           loc, neighb=NULL, maxdist=NULL, maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE,  
           param, anisopars=NULL,radius=6371, sparse=FALSE, time=NULL, type="Standard",type_mse=NULL,
            type_krig="Simple",weigthed=TRUE, which=1, copula=NULL,X=NULL,Xloc=NULL)
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)  or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) giving the data used for prediction.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates used for prediction. qn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates.   Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates used for prediction; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates used for prediction; the default is \code{NULL}
    then a spatial random field is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
   \item{distance}{String; the name of the spatial distance. The default
     is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details} of \code{\link{GeoFit}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data used for prediction
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
  \item{loc}{A numeric (\eqn{n \times 2}{n x 2})-matrix (where
    \code{n} is the number of spatial sites) giving 2-dimensions of spatial coordinates to be predicted.}
   \item{neighb}{Numeric; an optional positive integer indicating the 
    order of the  neighborhood.}
  \item{maxdist}{Numeric; an optional positive value indicating the distance in the spatial neighborhood.}
  \item{maxtime}{Numeric; an optional positive integer value indicating the order of the temporal neighborhood.}
     \item{method}{String; the type of matrix decomposition used in the simulation. Default is \code{cholesky}. 
     The other possible choices is \code{svd}.}
    \item{n}{Numeric; the number of trials in a binomial random fields.
         Default is \eqn{1}.}   
     \item{nloc}{Numeric; the number of trials  of the locations sites to be predicted in the binomial random field.
         If missing then a rounded mean of n is considered.} 
    \item{mse}{Logical; if \code{TRUE} (the default) MSE of the kriging predictor  is computed}    
   \item{model}{String; the type of RF and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}      
   \item{param}{A list of parameter values required  for the correlation model.See the Section \bold{Details}.}
       \item{anisopars}{ A list of two elements: "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
    \item{radius}{Numeric: the radius of the sphere if coordinates are passed in lon/lat format;}
     \item{sparse}{Logical; if \code{TRUE} kriging is computed with sparse matrices algorithms 
          using spam package. Default is FALSE. It should be used with compactly supported covariances.}        
    \item{time}{A numeric (\eqn{m \times 1}{m x 1}) vector (where
    \code{m} is the number of temporal instants) giving the temporal instants to be predicted; the default is \code{NULL}
    then only spatial prediction is performed.}
   \item{type}{String; if \code{Standard}   then standard kriging is performed;if \code{Tapering}   
   then  kriging with covariance tapering is performed;if \code{Pairwise} then pairwise
     kriging is performed}
    \item{type_mse}{String; if  \code{Theoretical} then theoretical MSE pairwise kriging is computed.
       If \code{SubSamp} then  an estimation based on subsampling is computed.}  
   \item{type_krig}{String; the type of kriging. If \code{Simple} (the
    default) then simple kriging is performed. (See the Section \bold{Details}).}
   \item{weigthed}{Logical; if \code{TRUE} then decreasing weigths coming 
  from a compactly supported correlation function
  with compact support \code{maxdist} (\code{maxtime})are used in the pairwise kriging.}
      \item{which}{Numeric; In the case of bivariate (tapered) cokriging it indicates which variable to predict.
           It can be 1 or 2}
    \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
  \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}
    \item{Xloc}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification
     associated to predicted locations.}
}

\details{
This function use the \code{\link{GeoKrig}} with a
spatial or spatio-temporal neighborhood computed using  the function \code{\link{GeoNeighborhood}}.
The neighborhood is specified with the option \code{maxdist}  and \code{maxtime}.}

\value{
  Returns an object of class \code{Kg}.
  An object of class \code{Kg} is a list containing
  at most the following components:
    \item{bivariate}{\code{TRUE} if spatial bivariate cokriging is performed, otherwise  \code{FALSE};}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates used for prediction;}
  \item{corrmodel}{String: the correlation model;}
  \item{covmatrix}{The covariance matrix if \code{type} is \code{Standard}. An object of class spam if \code{type} is \code{Tapering}   }
  \item{data}{The vector or matrix or array of data used for prediction}
  \item{distance}{String: the type of spatial distance;}
  \item{grid}{\code{TRUE} if the spatial data used for prediction are observed in a regular grid, otherwise  \code{FALSE};}
  \item{loc}{A (\eqn{n \times 2}{n x 2})-matrix of spatial locations to be predicted.}
  \item{n}{The number of trial for Binomial RFs} 
  \item{nozero}{In the case of tapered simple kriging the percentage of non zero values in the covariance matrix.
  Otherwise is NULL.}
  \item{numcoord}{Numeric:he number \eqn{d}{d} of spatial coordinates used for prediction;}
  \item{numloc}{Numeric: the number \eqn{n}{n} of spatial coordinates to be predicted;}
  \item{numtime}{Numeric: the number \eqn{d}{d}  of the temporal instants used for prediction;}
  \item{numt}{Numeric: the number \eqn{m}{m}  of the temporal instants to be predicted;}
  \item{model}{The type of RF, see \code{\link{GeoFit}}.}
  \item{param}{Numeric: The covariance parameters;}
  \item{pred}{A (\eqn{m \times n}{n x m})-matrix of spatio or spatio temporal kriging prediction;}
    \item{radius}{Numeric: the radius of the sphere if coordinates are pssed in lon/lat format;}
  \item{spacetime}{\code{TRUE} if spatio-temporal kriging and \code{FALSE} if
  spatial kriging;}
   \item{tapmod}{String: the taper model if \code{type} is \code{Tapering}. Otherwise is NULL.}
  \item{time}{A \eqn{m}{m}-dimensional vector of temporal coordinates to be predicted;}
  \item{type}{String:  the type of kriging (Standard or Tapering).}
  \item{type_krig}{String: the type of kriging.}
  \item{mse}{A (\eqn{m \times n}{n x m})-matrix of spatio or spatio temporal mean square error kriging prediction;}

}


\references{
  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.
  Furrer R.,  Genton, M.G. and  Nychka D. (2006). \emph{Covariance
  Tapering for Interpolation of Large Spatial Datasets}.
  Journal of Computational and Graphical Statistics, \bold{15-3}, 502--523.


}

\seealso{\code{\link{GeoCovmatrix}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\examples{

################################################################
############### Examples of Spatial local kriging  #############
################################################################
require(GeoModels)
require(fields)
####
model="Gaussian"

# Define the spatial-coordinates of the points:
set.seed(79)
x = runif(500, 0, 1)
y = runif(500, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "Exponential"
mean=0
sill=1
nugget=0
scale=0.3/3
param=list(mean=mean,sill=sill,nugget=nugget,scale=scale)


# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,
              param=param)$data

# locations to predict
xx=seq(0,1,0.025)
loc_to_pred=as.matrix(expand.grid(xx,xx))

################################################################
###
### Example 1. Comparing spatial kriging with local kriging for
### a Gaussian random field with exponential correlation.
### 
###############################################################
pr=GeoKrig(loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,
      model=model,param= param, data=data,mse=TRUE)

pr_loc=GeoKrigloc(data=data,loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,
      model=model,neighb=50, param= param,mse=TRUE)

colour = rainbow(100)
par(mfrow=c(2,2))
zlim=c(-2.6,2.6)
# simple kriging map prediction
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=colour,
           zlim=zlim,xlab="",ylab="",
           main="Kriging with  exponential model ")

# simple kriging MSE map prediction variance
image.plot(xx, xx, matrix(pr$mse,ncol=length(xx)),col=colour,
           xlab="",ylab="",main="Std error")

# simple  kriging map prediction
image.plot(xx, xx, matrix(pr_loc$pred,ncol=length(xx)),col=colour,
           zlim=zlim,xlab="",ylab="",main="Local Kriging with exponential model")

# simple kriging MSE map prediction variance
image.plot(xx, xx, matrix(pr_loc$mse,ncol=length(xx)),col=colour,
           xlab="",ylab="",main="Std error")


############################################################
#### Example: spatio temporal  Gaussian local kriging ######
############################################################


require(GeoModels)
require(fields)
set.seed(78)
coords=cbind(runif(50),runif(50))
coordt=seq(0,5,0.25)
corrmodel="Exp_Exp"
param=list(nugget=0,mean=0,scale_s=0.2/3,scale_t=0.25/3,sill=2)

data = GeoSim(coordx=coords, coordt=coordt,corrmodel="Exp_Exp", 
                         param=param)$data
##  four location to predict
loc_to_pred=matrix(runif(8),4,2)
## three temporal instants to predict
time=c(0.5,1.2,3.7)


pr=GeoKrig(data=data,loc=loc_to_pred,time=time,coordx=coords,coordt=coordt,corrmodel=corrmodel,
  model="Gaussian", param= param, mse=TRUE)
pr_loc=GeoKrigloc(data=data,loc=loc_to_pred,time=time,coordx=coords,coordt=coordt,corrmodel=corrmodel,
  neigh=25,maxtime=1,model="Gaussian", param= param, mse=TRUE)

##  full and local prediction 
pr$pred
pr_loc$pred


############################################################
### Example 7. Bivariate local cokriging  ##################
############################################################
require(GeoModels)
require(fields)
# Define the spatial-coordinates of the points:
set.seed(128)
x = runif(200, 0, 1)
y = runif(200, 0, 1)
coords=cbind(x,y)

param=list(scale=0.3/3,mean_1=0,mean_2=0,sill_1=1,sill_2=1,
          nugget_1=0,nugget_2=0,pcol=0.7,smooth=0.5)
          
data = GeoSim(coordx=coords, corrmodel="Bi_matern_sep", param=param)$data

fixed=list(nugget_1=0,nugget_2=0,smooth=0.5,mean_1=0,mean_2=0)
start=list(sill_1=var(data[,1]),sill_2=var(data[,2]),scale=0.3/3,
      pcol=cor(data[1,],data[2,]))
# Maximum Composite likelihood fitting of the random field:
fitcl= GeoFit(data, coordx=coords, corrmodel="Bi_matern_sep",
        likelihood="Marginal",type="Pairwise",neighb=c(3,3,3),
        optimizer="BFGS", start=start,fixed=fixed)


# locations to predict

xx=seq(0,1,0.025)
loc_to_pred=as.matrix(expand.grid(xx,xx))

pr=GeoKrig(loc=loc_to_pred,coordx=coords,corrmodel="Bi_matern_sep",
       param= append(fitcl$param,fitcl$fixed), data=data,which=1,mse=TRUE)

pr_loc=GeoKrigloc(loc=loc_to_pred,coordx=coords,corrmodel="Bi_matern_sep",maxdist=0.15,param= append(fitcl$param,fitcl$fixed), 
  data=data,which=1,mse=TRUE)

par(mfrow=c(2,2))

colour = rainbow(100)
zlim=c(-2.6,2.6)
# simple kriging map prediction
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=colour,
           zlim=zlim,xlab="",ylab="",
           main="CoKriging with Bivariate Matern model")

# simple kriging MSE map prediction variance
image.plot(xx, xx, matrix(pr$mse,ncol=length(xx)),col=colour,
           xlab="",ylab="",main="Std error")

# simple  kriging map prediction
image.plot(xx, xx, matrix(pr_loc$pred,ncol=length(xx)),col=colour,
           zlim=zlim,xlab="",ylab="",main="Local CoKriging with Bivariate Matern model")

# simple kriging MSE map prediction variance
image.plot(xx, xx, matrix(pr_loc$mse,ncol=length(xx)),col=colour,
           xlab="",ylab="",main="Std error")


}
\keyword{Composite}
