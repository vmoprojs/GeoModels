\name{GeoKrig}
\alias{GeoKrig}
\encoding{UTF-8}
\title{Spatial (bivariate) and spatio temporal optimal linear prediction 
          for Gaussian and non Gaussian random fields.}
\description{
For a  given set of  spatial location sites (and temporal instants),
the function computes optimal linear prediction and associated mean square error 
for the Gaussian and non Gaussian case.
}
\usage{
GeoKrig(data, coordx, coordy=NULL, coordt=NULL, 
coordx_dyn=NULL, corrmodel,distance="Eucl",
    grid=FALSE, loc, maxdist=NULL, maxtime=NULL,
    method="cholesky", model="Gaussian", n=1,nloc=NULL,mse=FALSE, 
    lin_opt=TRUE,  param, anisopars=NULL,radius=6371, sparse=FALSE,
    taper=NULL,tapsep=NULL, time=NULL, type="Standard",type_mse=NULL,
     type_krig="Simple",weigthed=TRUE,which=1,
     copula=NULL, X=NULL,Xloc=NULL,Mloc=NULL)
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)  or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) giving the data used for prediction.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving 1-dimension of spatial coordinates used for prediction. \eqn{d}{d}-dimensional vector giving 1-dimension of spatial coordinates.   Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates used for prediction; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
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
  \item{lin_opt}{Logical;If TRUE (default) then  optimal (pairwise) linear  kriging is computed. Otherwise  optimal (pairwise) kriging is computed
  in the mean square sense.}  
  \item{loc}{A numeric (\eqn{n \times 2}{n x 2})-matrix (where
    \code{n} is the number of spatial sites) giving 2-dimensions of spatial coordinates to be predicted.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial compact support in the case of covariance tapering kriging.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal compact support in the case of covasriance tapering kriging.}
     \item{method}{String; the type of matrix decomposition used in the simulation. Default is \code{cholesky}. 
     The other possible choices is \code{svd}.}
    \item{n}{Numeric; the number of trials in a binomial random fields.
         Default is \eqn{1}.}   
     \item{nloc}{Numeric; the number of trials  of the locations sites to be predicted in a binomial random fields type II.
         Default is \eqn{1}.} 
    \item{mse}{Logical; if \code{TRUE} (the default) MSE of the kriging predictor  is computed}    
   \item{model}{String; the type of RF and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}      
   \item{param}{A list of parameter values required  for the correlation model.See the Section \bold{Details}.}
      \item{anisopars}{ A list of two elements: "angle" and "ratio" i.e. the anisotropy angle and the anisotropy ratio, respectively.}
    \item{radius}{Numeric: the radius of the sphere if coordinates are passed in lon/lat format;}
     \item{sparse}{Logical; if \code{TRUE} kriging is computed with sparse matrices algorithms 
          using spam package. Default is FALSE. It should be used with compactly supported covariances.}        
   \item{taper}{String; the name of the taper correlation
    function, see the Section \bold{Details}.}
    \item{tapsep}{Numeric; an optional  value indicating the separabe
      parameter in the space time quasi taper (see \bold{Details}).}
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
}

\details{
 Best linear unbiased predictor and associated mean square error is computed 
 for Gaussian and some non Gaussian cases.
  Specifically, for a spatial or spatio-temporal or spatial bivariate dataset, given a set of  spatial locations and
  temporal istants and a correlation model
\code{corrmodel} with some fixed parameters and given  the type of RF (\code{model}) the function computes
  simple  kriging, for the specified spatial locations
  \code{loc} and temporal instants \code{time}, providing also the respective mean square error.
  For the choice of the spatial or spatio temporal correlation model see details in \code{\link{GeoCovmatrix}} function.
  The list \code{param} specifies mean and covariance parameters, see \code{\link{CorrParam}} and \code{\link{GeoCovmatrix}} for details. The \code{type_krig} parameter indicates the type of kriging. In the
case of simple kriging, the known mean can be specified by the parameter
  \code{mean} in the  list \code{param} (See examples).
}

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


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)

################################################################
########### Examples of spatial kriging ############
################################################################

################################################################
###
### Example 1. Spatial  kriging of a
### Gaussian random fields with Gen wendland correlation.
###
################################################################

model="Gaussian"
set.seed(79)
x = runif(400, 0, 1)
y = runif(400, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "GenWend"
mean=0; sill=5; nugget=0
scale=0.2;smooth=0;power2=4

param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth,power2=power2)

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,
              param=param)$data

## estimation with pairwise likelihood
fixed=list(nugget=nugget,smooth=0,power2=power2)
start=list(mean=0,scale=scale,sill=1)
I=Inf
lower=list(mean=-I,scale=0,sill=0)
upper=list(mean= I,scale=I,sill=I)
# Maximum pairwise likelihood fitting :
fit = GeoFit(data, coordx=coords, corrmodel=corrmodel,model=model,
                    likelihood='Marginal', type='Pairwise',neighb=3,
                    optimizer="nlminb", lower=lower,upper=upper,
                    start=start,fixed=fixed)

# locations to predict
xx=seq(0,1,0.03)
loc_to_pred=as.matrix(expand.grid(xx,xx))

pr=GeoKrig(loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,
       model=model,param= param, data=data,mse=TRUE)

colour = rainbow(100)

opar=par(no.readonly = TRUE)
par(mfrow=c(1,3))
quilt.plot(coords,data,col=colour)
# simple kriging map prediction
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=colour,
           xlab="",ylab="",
           main=" Kriging ")

# simple kriging MSE map prediction variance
image.plot(xx, xx, matrix(pr$mse,ncol=length(xx)),col=colour,
           xlab="",ylab="",main="Std error")
par(opar)

################################################################
###
### Example 2. Spatial  kriging of a Skew
### Gaussian random fields with Matern correlation.
###
################################################################
model="SkewGaussian"
set.seed(79)
x = runif(400, 0, 1)
y = runif(400, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "Matern"
mean=0
sill=2
nugget=0
scale=0.1
smooth=0.5
skew=3
param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth,skew=skew)

# Simulation of the spatial skew Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,
              param=param)$data

fixed=list(nugget=nugget,smooth=smooth)
start=list(mean=0,scale=scale,sill=1,skew=skew)
I=Inf
lower=list(mean=-I,scale=0,sill=0,skew=-I)
upper=list(mean= I,scale=I,sill=I,skew=I)
# Maximum pairwise likelihood fitting :
fit = GeoFit2(data, coordx=coords, corrmodel=corrmodel,model=model,
                    likelihood='Marginal', type='Pairwise',neighb=3,
                    optimizer="nlminb", lower=lower,upper=upper,
                    start=start,fixed=fixed)

# locations to predict
xx=seq(0,1,0.03)
loc_to_pred=as.matrix(expand.grid(xx,xx))
## optimal linear kriging
pr=GeoKrig(loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,
       model=model,param= param, data=data,mse=TRUE)

colour = rainbow(100)

opar=par(no.readonly = TRUE)
par(mfrow=c(1,3))
quilt.plot(coords,data,col=colour)
# simple kriging map prediction
image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=colour,
           xlab="",ylab="",
           main=" Kriging ")

# simple kriging MSE map prediction variance
image.plot(xx, xx, matrix(pr$mse,ncol=length(xx)),col=colour,
           xlab="",ylab="",main="Std error")
par(opar)

################################################################
###
### Example 3. Spatial  kriging  of a 
###   Gamma random field with mean spatial regression
###
###############################################################
set.seed(312)
model="Gamma"
corrmodel = "GenWend"  
# Define the spatial-coordinates of the points:
NN=400
coords=cbind(runif(NN),runif(NN))
## matrix covariates
a0=rep(1,NN)
a1=runif(NN,0,1)
X=cbind(a0,a1)
 ##Set model parameters
shape=2
## regression parameters
mean = 1;mean1= -0.2
# correlation parameters
nugget = 0;power2=4
scale = 0.3;smooth=0    

## simulation
param=list(shape=shape,nugget=nugget,mean=mean,mean1=mean1, 
  scale=scale,power2=power2,smooth=smooth)
data = GeoSim(coordx=coords,corrmodel=corrmodel, param=param,
                model=model,X=X)$data

#####starting and fixed parameters
fixed=list(nugget=nugget,power2=power2,smooth=smooth)
start=list(mean=mean,mean1=mean1, scale=scale,shape=shape)

## estimation with pairwise likelihood
fit2 = GeoFit(data=data,coordx=coords,corrmodel=corrmodel,X=X,
              neighb=3,likelihood="Conditional",type="Pairwise",
              start=start,fixed=fixed, model = model)
             
# locations to predict with associated covariates
xx=seq(0,1,0.03)
loc_to_pred=as.matrix(expand.grid(xx,xx))
NP=nrow(loc_to_pred)
a0=rep(1,NP)
a1=runif(NP,0,1)
Xloc=cbind(a0,a1)

#optimal linear  kriging 
pr=GeoKrig(data=data, coordx=coords,loc=loc_to_pred,corrmodel=corrmodel,model=model,mse=TRUE,
        X=X,Xloc=Xloc,
       sparse=TRUE,param=append(fit2$param,fit2$fixed))

## map 
opar=par(no.readonly = TRUE)
par(mfrow=c(1,3))
quilt.plot(coords,data,main="Data")
map=matrix(pr$pred,ncol=length(xx))
mapmse=matrix(pr$mse,ncol=length(xx))
image.plot(xx, xx, map,
           xlab="",ylab="",main="Kriging ")

image.plot(xx, xx, mapmse,
           xlab="",ylab="",main="MSE")
par(opar)


################################################################
########### Examples of spatio temporal kriging ############
################################################################

################################################################
###
### Example 3. Spatio temporal simple kriging of n locations
### sites and m temporal instants for a Gaussian random fields
### with estimated double Wendland correlation.
###
###############################################################
model="Gaussian"
# Define the spatial-coordinates of the points:
x = runif(30, 0, 1)
y = runif(30, 0, 1)
#coords=cbind(x,y)
times=1:5

# Define model correlation modek and associated parameters
corrmodel="Wend0_Wend0"
param=list(nugget=0,mean=0,power2_s=4,power2_t=4,
        scale_s=0.1,scale_t=2,sill=1)

# Simulation of the space time Gaussian random field:
set.seed(31)
data=GeoSim(coordx=coords,coordt=times,corrmodel=corrmodel,sparse=TRUE,
           param=param)$data

# Maximum pairwise likelihood fitting of the space time random field:
start = list(scale_s=0.15,scale_t=2,sill=1,mean=0)
fixed = list(nugget=0,power2_s=4,power2_t=4)
I=Inf
lower=list(scale_s=0,scale_t=0,sill=0,mean=-I)
upper=list(scale_s=I,scale_t=I,sill=I,mean=I)
fit = GeoFit(data, coordx=coords, coordt=times, model=model, corrmodel=corrmodel, 
             likelihood='Conditional', type='Pairwise',start=start,fixed=fixed,
              neighb=3,maxtime=1)


param=append(fit$param,fit$fixed)

# locations to predict
xx=seq(0,1,0.05)
loc_to_pred=as.matrix(expand.grid(xx,xx))
#  Define the times to predict
times_to_pred=1:2

pr=GeoKrig(loc=loc_to_pred,time=times_to_pred,coordx=coords,coordt=times,
      sparse=TRUE,corrmodel=corrmodel, param=param,data=data,mse=TRUE)

opar=par(no.readonly = TRUE)
par(mfrow=c(2,3))
zlim=c(-2.5,2.5)
colour = rainbow(100)

for(i in 1:2) {
quilt.plot(coords,data[i,] ,col=colour,main = paste(" data  at Time=" , i))   
image.plot(xx, xx, matrix(pr$pred[i,],ncol=length(xx)),col=colour,
        main = paste(" Kriging Time=" , i),ylab="")
image.plot(xx, xx, matrix(pr$mse[i,],ncol=length(xx)),col=colour,
           main = paste("Std err Time=" , i),ylab="")
}
par(opar)
}
\keyword{Composite}
