\name{GeoKrigloc}
\alias{GeoKrigloc}
\encoding{UTF-8}
\title{Spatial (bivariate) and spatio temporal optimal linear  local  prediction 
          for Gaussian and non Gaussian RFs.}
\description{
For a  given set of  spatial location sites (and temporal instants),
the function computes optmal local  linear prediction and the associated mean squared error 
for the Gaussian and non Gaussian case using a spatial (temporal) neighborhood
computed using the  function \code{\link{GeoNeighborhood}}
}
\usage{
GeoKrigloc(estobj=NULL,data, coordx, coordy=NULL, coordz=NULL,coordt=NULL,
    coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, 
    loc, neighb=NULL, maxdist=NULL, 
    maxtime=NULL, method="cholesky",
    model="Gaussian", n=1,nloc=NULL, mse=FALSE,  
    param, anisopars=NULL,radius=6371,
    sparse=FALSE, time=NULL, type="Standard",type_mse=NULL,
    type_krig="Simple",weigthed=TRUE, 
    which=1, copula=NULL,X=NULL,Xloc=NULL,
    Mloc=NULL,spobj=NULL,spdata=NULL,parallel=FALSE,ncores=NULL)
}
\arguments{
   \item{estobj}{An object of class Geofit that includes information about data, model and estimates.}
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)  or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) giving the data used for prediction.}
   \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix  or  (\eqn{d \times 3}{d x 3})-matrix 
     Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; Optional argument, the default is \code{NULL}.}
      \item{coordz}{A numeric vector giving 1-dimension of
    spatial coordinates;  Optional argument, the default is \code{NULL}.}

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
    default) then simple kriging is performed. If \code{Optim} then optimal kriging is performed for some non-Gaussian RFs}
   \item{weigthed}{Logical; if \code{TRUE} then decreasing weigths coming 
  from a compactly supported correlation function
  with compact support \code{maxdist} (\code{maxtime})are used in the pairwise kriging.}
      \item{which}{Numeric; In the case of bivariate (tapered) cokriging it indicates which variable to predict.
           It can be 1 or 2}
    \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
  \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}
    \item{Xloc}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification
     associated to predicted locations.}
    \item{Mloc}{Numeric; Vector of spatio(temporal) estimated means 
     associated to predicted locations.}
                \item{spobj}{An object of class sp or spacetime}
   \item{spdata}{Character:The name of data in the sp or spacetime object}
       \item{parallel}{Logical; if \code{TRUE} then 
      the prediction computation is  parallelized 
      }
       \item{ncores}{Numeric; number of cores involved in parallelization.} 
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
   \item{coordz}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
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



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

################################################################
############### Examples of Spatial local kriging  #############
################################################################
require(GeoModels)
####
model="Gaussian"

# Define the spatial-coordinates of the points:
set.seed(759)
x = runif(1000, 0, 1)
y = runif(1000, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "GenWend"
mean=0; sill=1
nugget=0; scale=0.2
param=list(mean=mean,sill=sill,nugget=nugget,smooth=0,
scale=scale,power2=4)

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,
              param=param)$data

# Maximum pairwise likelihood fitting of the space time random field:

start=list(scale=scale,sill=sill,mean=mean)
fixed=list(power2=4,smooth=0,nugget=0)
fit = GeoFit(data, coordx=coords, corrmodel=corrmodel, 
              start=start,fixed=fixed,
             likelihood='Conditional', type='Pairwise',
              neighb=3)

# locations to predict
loc_to_pred=matrix(runif(8),4,2)
################################################################
###
### Example 1. Comparing spatial kriging with local kriging for
### a Gaussian random field with GenWend correlation.
### 
###############################################################
param=append(fit$param,fit$fixed)
pr=GeoKrig(fit,loc=loc_to_pred,mse=TRUE)

pr_loc=GeoKrigloc(fit,loc=loc_to_pred,neighb=100,mse=TRUE)

pr$pred;
pr_loc$pred


############################################################
#### Example: spatio temporal  Gaussian local kriging ######
############################################################


require(GeoModels)
require(fields)
set.seed(78)
coords=cbind(runif(100),runif(100))
coordt=seq(0,5,0.25)
corrmodel="Matern_Matern"
param=list(nugget=0,mean=0,scale_s=0.2/3,scale_t=0.25/3,sill=2,
  smooth_s=0.5,smooth_t=0.5)

data = GeoSim(coordx=coords, coordt=coordt,
     corrmodel=corrmodel, param=param)$data


# Maximum pairwise likelihood fitting of the space time random field:
start = list(scale_s=0.2/3,scale_t=0.25,sill=2,mean=0)
fixed = list(smooth_s=0.5,smooth_t=0.5,nugget=0)
I=Inf
lower=list(scale_s=0,scale_t=0,sill=0,mean=-I)
upper=list(scale_s=I,scale_t=I,sill=I,mean=I)
fit = GeoFit(data, coordx=coords, coordt=coordt, model=model, corrmodel=corrmodel, 
             likelihood='Conditional', type='Pairwise',start=start,fixed=fixed,
             optimizer="nlminb",lower=lower,upper=upper,
              neighb=3,maxtime=1)

##  four location to predict
loc_to_pred=matrix(runif(8),4,2)
## three temporal instants to predict
time=c(0.5,1.5,3.5)


pr=GeoKrig(fit,loc=loc_to_pred,time=time,mse=TRUE)
pr_loc=GeoKrigloc(fit,loc=loc_to_pred,time=time,
  neigh=25,maxtime=1, mse=TRUE)

##  full and local prediction 
pr$pred
pr_loc$pred



############################################################
#### Example: spatio bivariate  Gaussian local cokriging ######
############################################################
#set.seed(6)
#NN=1500 # number of spatial locations 
#x = runif(NN, 0, 1); 
#y = runif(NN, 0, 1) 
#coords=cbind(x,y) 

## setting parameters
#mean_1 = 2; mean_2= -1
#nugget_1 =0;nugget_2=0
#sill_1 =0.5; sill_2 =1; 

### correlation parameters
#CorrParam("Bi_Matern")
#scale_1=0.2/3; scale_2=0.15/3; scale_12=0.5*(scale_2+scale_1) 
#smooth_1=smooth_2=smooth_12=0.5
#pcol = -0.4
#param= list(nugget_1=nugget_1,nugget_2=nugget_2,
#          sill_1=sill_1,sill_2=sill_2,
#          mean_1=mean_1,mean_2=mean_2,
#          smooth_1=smooth_1, smooth_2=smooth_2,smooth_12=smooth_12,
#          scale_1=scale_1, scale_2=scale_2,scale_12=scale_12,
#          pcol=pcol)

## simulation
#data = GeoSim(coordx=coords, corrmodel="Bi_Matern",model=model,param=param)$data

#fixed=list(mean_1=mean_1,mean_2=mean_2, nugget_1=nugget_1,nugget_2=nugget_2, 
#       smooth_1=smooth_1, smooth_2=smooth_2,smooth_12=smooth_12)

#start=list( sill_1=sill_1,sill_2=sill_2,
#        scale_1=scale_1,scale_2=scale_2,scale_12=scale_12, pcol=pcol)

## estimation with maximum likelihood 
#fit = GeoFit(data=data,coordx=coords, corrmodel="Bi_Matern",
 # likelihood="Marginal",type="Pairwise",optimizer="BFGS",neighb=5,
  #start=start,fixed=fixed)

###### co-kriging for the fist component ##############
#xx=seq(0,1,0.022)
#loc_to_pred=as.matrix(expand.grid(xx,xx))
#pr1 = GeoKrigloc(fit,which=1,mse=TRUE,loc=loc_to_pred,neighb=100)

#opar=par(no.readonly = TRUE)
#par(mfrow=c(1,2))
#zlim=c(-2.5,2.5)
#colour = rainbow(100)
#quilt.plot(coords,data[1,] ,col=colour,main = paste(" Fist component"))   
#quilt.plot(loc_to_pred,pr1$pred,col=colour,
#          main = paste(" Kriging first component"),ylab="")
#par(opar)

}
\keyword{Composite}
