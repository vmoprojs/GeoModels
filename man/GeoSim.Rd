\name{GeoSim}
\alias{GeoSim}
\alias{print.GeoSim}
\encoding{UTF-8}
\title{Simulation of Gaussian and non Gaussian Random Fields.}
\description{
  Simulation of Gaussian and some non Gaussian spatial, spatio-temporal and spatial bivariate 
 random fields.
  The function return a realization of a Random Field for a given covariance model and covariance parameters.
  Simulation is based on Cholesky decomposition.
}
\usage{
GeoSim(coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl",
      GPU=NULL, grid=FALSE, local=c(1,1),method="cholesky", model='Gaussian', n=1, param,
      radius=6371, sparse=FALSE, X=NULL)
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
     \item{method}{String; the type of matrix decomposition used in the simulation. Default is cholesky. The other possible choices is \code{svd}.} 
  \item{model}{String; the type of RF and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{n}{Numeric; the number of trials for binomial RFs. The number of successes in the negative Binomial RFs. Default is \eqn{1}.} 
  \item{param}{A list of parameter values required in the
    simulation procedure of RFs, see \bold{Examples}.}
    \item{radius}{Numeric; a value indicating  the radius of the sphere when using the great 
    circle distance. Default value is the radius of the earth in Km (i.e. 6371)}  
  \item{sparse}{Logical; if \code{TRUE} then  cholesky decomposition is performed
  using sparse matrices algorithms (spam packake). It should be used
  with compactly supported covariance models.\code{FALSE} is the default.}
  \item{X}{Numeric; Matrix of space-time covariates.}
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
  \item{model}{The type of RF, see \code{\link{GeoFit}}.}
  \item{n}{The number of trial for Binomial RFs;the number of successes in a negative Binomial RFs;} 
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number the temporal realisations of the RF;}
  \item{param}{The vector of  parameters' estimates;}
  \item{radius}{The radius of the sphere if coordinates are passed in lon/lat format;}
  \item{randseed}{The seed used for the random simulation;}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if
  spatial RF;}
}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}


\examples{
library(GeoModels)
library(mapproj)
library(fields)

################################################################
###
### Example 1. Simulation of a spatial Gaussian RF on a regular grid
###
###############################################################

# Define the spatial-coordinates of the points:
x <- seq(0,1,0.045)
y <- seq(0,1,0.045)
set.seed(261)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(x,y,grid=TRUE, corrmodel="Matern", param=list(smooth=0.5,
             mean=0,sill=1,scale=0.4/3,nugget=0))$data

# Simulation of a spatial Gaussian RF with Generalized Wendland correlation function
# using sparse alghorithm matrices
set.seed(261)
data2 <- GeoSim(x,y,grid=TRUE, corrmodel="GenWend", param=list(smooth=0,
              power2=4,mean=0,sill=1,scale=0.4,nugget=0))$data
par(mfrow=c(1,2))
image.plot(x,y,data1,col=terrain.colors(100),main="Matern",xlab="",ylab="")
image.plot(x,y,data2,col=terrain.colors(100),main="Wendland",xlab="",ylab="")   

             

################################################################
###
### Example 2. Simulation of a spatial binomial RF 
### with  underlying exponential correlation
### on a regular grid
###
################################################################

# Define the spatial-coordinates of the points:
x <- seq(0, 1, 0.022)
y <- seq(0, 1, 0.022)
coords <- cbind(x,y)
set.seed(251)
n=5
# Simulation of a spatial Binomial RF:
sim <- GeoSim(x,y,grid=TRUE, corrmodel="Wend0",
             model="Binomial",n=n,sparse=TRUE,
             param=list(nugget=0,mean=0,scale=.2,sill=1,power2=4))

image.plot(x,y,sim$data,nlevel=n+1,col=terrain.colors(n+1),zlim=c(0,n)) 


################################################################
###
### Example 3. Simulation of a spatial Weibull  RF
### with  underlying exponential correlation 
###
###############################################################
# Define the spatial-coordinates of the points:
x <- seq(0,1,0.032)
y <- seq(0,1,0.032)
set.seed(261)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(x,y,grid=TRUE, corrmodel="Exponential",model="Weibull", 
         param=list(shape=1.2,mean=0,sill=1,scale=0.3/3,nugget=0))$data
image.plot(x,y,data1,col=terrain.colors(200),main="Weibull RF",xlab="",ylab="")



################################################################
###
### Example 4. Simulation of a spatial t  RF
### with  with  underlying exponential correlation 
###
###############################################################
# Define the spatial-coordinates of the points:
x <- seq(0,1,0.03)
y <- seq(0,1,0.03)
set.seed(268)
# Simulation of a spatial Gaussian RF with Matern correlation function
data1 <- GeoSim(x,y,grid=TRUE, corrmodel="GenWend",model="StudentT", sparse=TRUE,
         param=list(df=1/4,mean=0,sill=1,scale=0.3,nugget=0,smooth=1,power2=5))$data
image.plot(x,y,data1,col=terrain.colors(100),main="Student-t RF",xlab="",ylab="")




################################################################
###
### Example 5. Simulation of a sinhasinh RF
###   with  underlying Wend0 correlation.
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(800, 0, 2)
y <- runif(800, 0, 2)
coords <- cbind(x,y)
set.seed(261)
corrmodel="Wend0"
# Simulation of a spatial Gaussian RF:
param=list(power2=4,skew=0,tail=1,
             mean=0,sill=1,scale=0.2,nugget=0)  ## gaussian case
data0 <- GeoSim(coordx=coords, corrmodel=corrmodel,
               model="SinhAsinh", param=param,sparse=TRUE)$data
plot(density(data0),xlim=c(-7,7))

param=list(power2=4,skew=0,tail=0.7,
             mean=0,sill=1,scale=0.2,nugget=0) ## heavy tails
data1 <- GeoSim(coordx=coords, corrmodel=corrmodel,
               model="SinhAsinh", param=param,sparse=TRUE)$data
lines(density(data1),lty=2)

param=list(power2=4,skew=0.5,tail=1,
             mean=0,sill=1,scale=0.2,nugget=0)  ## asymmetry
data2 <- GeoSim(coordx=coords, corrmodel=corrmodel,
               model="SinhAsinh", param=param,sparse=TRUE)$data
lines(density(data2),lty=3)



################################################################
###
### Example 6. Simulation of a bivariate Gaussian RF
### with separable bivariate exponential correlation model
### on a regular grid.
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)
coords <- cbind(x,y)

# Simulation of a bivariate spatial Gaussian RF:
# with a separable Bivariate Matern
set.seed(12)
param=list(mean_1=4,mean_2=2,scale=0.12,smooth=0.5,
           sill_1=1,sill_2=1,nugget_1=0,nugget_2=0,pcol=0.5)
data <- GeoSim(coordx=coords,corrmodel="Bi_matern_sep",
              param=param)$data
par(mfrow=c(1,2))

quilt.plot(coords,data[1,],col=terrain.colors(100),main="1",xlab="",ylab="")
quilt.plot(coords,data[2,],col=terrain.colors(100),main="2",xlab="",ylab="")





################################################################
###
### Example 7. Simulation of a  spatio temporal Gaussian RF.
### observed on  dynamical location sites with double exponential correlation 
###
###############################################################

# Define the dynamical spatial-coordinates of the points:

coordt=1:5
coordx_dyn=list()
maxN=40
set.seed(8)
for(k in 1:length(coordt))
{
NN=sample(1:maxN,size=1)
x <- runif(NN, 0, 1)
y <- runif(NN, 0, 1)
coordx_dyn[[k]]=cbind(x,y)
}
coordx_dyn

param<-list(nugget=0,mean=0,scale_s=0.2/3,scale_t=2/3,sill=1)
data <- GeoSim(coordx_dyn=coordx_dyn, coordt=coordt, corrmodel="Exp_Exp",
                     param=param)$data
## spatial realization at first temporal instants
data[[1]]
## spatial realization at third temporal instants
data[[3]]




################################################################
###
### Example 8. Simulation of a Gaussian RF 
###  with a Wend0 correlation in the north emisphere of the planet earth
###
###############################################################
distance="Geod";radius=6371

NN=1000 ## total point on the sphere on lon/lat format
set.seed(80)
coords=cbind(runif(NN,-180,180),runif(NN,0,90))
## Set the wendland parameters
corrmodel <- "Wend0"
param<-list(mean=0,sill=1,nugget=0,scale=radius*0.1,power2=3)
# Simulation of a spatial Gaussian RF on the sphere
#set.seed(2)
data <- GeoSim(coordx=coords,corrmodel=corrmodel,sparse=TRUE,
        distance=distance,radius=radius,param=param)$data



################################################################
###
### Example 9. Simulation of a Gaussian RF 
###  with Wend0 model  on USA
###
###############################################################

distance="Geod";radius=6378.88
NN=40
x=seq(-125,-64,length.out=NN)
y=seq(27,50, length.out =NN)
nrow(expand.grid(x,y))
## Set the wendland parameters
corrmodel <- "Wend0"
param<-list(mean=0,sill=1,nugget=0,scale=radius*0.3,power2=3)
# Simulation of a spatial Gaussian RF on the sphere
#set.seed(2)
data <- GeoSim(x,y,grid=TRUE,corrmodel=corrmodel,sparse=TRUE,
        distance=distance,radius=radius,param=param)$data
image.plot(x,y,data,col=terrain.colors(100),xlab="",ylab="")
map("usa", add = TRUE)




################################################################
###
### Example 10. Simulation of a Wrapped  RF 
### with   underlying exponential correlation
### on a regular grid
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(200,0, 1)
y <- runif(200,0, 1)
coords <- cbind(x,y)
set.seed(251)

# Simulation of a spatial wrapped RF:
sim <- GeoSim(coordx=coords, corrmodel="Exp",
             model="Wrapped", 
             param=list(nugget=0,mean=0,scale=.1,sill=1))$data

long <- 0.08;          
x1 <- coords[,1] + long*cos(sim)
y1 <- coords[,2] + long*sin(sim)
eps <- 0.1
plot(0,xlim=c(0-eps,1+eps),ylim=c(0-eps,1+eps));   
require(shape)
Arrows(coords[,1], coords[,2], x1, y1, arr.length = 0.2, code = 2,  arr.type = "triangle",col=1)


}
\keyword{Simulation}
