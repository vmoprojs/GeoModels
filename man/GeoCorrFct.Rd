\name{GeoCorrFct}
\alias{GeoCorrFct}
\encoding{UTF-8}
\title{Spatial and Spatio-temporal correlation or  covariance of (non) Gaussian random fields}
\description{
  The function computes the  correlations 
  a spatial or spatio-temporal or a bivariate spatial Gaussian or non Gaussian randomm field
  witha  given  covariance model and a  set of spatial (temporal) distances.
}
\usage{
GeoCorrFct(x,t=NULL,corrmodel, model="Gaussian",distance="Eucl", param, radius=6371,n=1,covariance=FALSE)
}

\arguments{
  \item{x}{A set of spatial distances.}
  \item{t}{A set of (optional) temporal distances.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{model}{String; the type of RF. See \code{\link{GeoFit}}.}
  \item{distance}{String; the name of the spatial distance. The default
  is \code{Eucl}, the euclidean distance. See \code{\link{GeoFit}}.}
  \item{param}{A list of parameter values required  for the covariance
  model.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere 
     when using covariance models valid using the great circle distance. 
     Default value is the radius of the earth in Km (i.e. 6371)} 
  \item{n}{Numeric; the number of trials in a (negative) binomial random fields.
         Default is \eqn{1}.}
  \item{covariance}{Logic; if TRUE then the covariance is returned. Default is FALSE}
   }


\value{Returns  correlations or covariances values associated to a given parametric spatial and temporal correlation models.}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\examples{
library(GeoModels)

################################################################
###
### Example 1. Matern correlation model
###
###############################################################

# Define the spatial distances
x = seq(0,1,0.005)

# Correlation Parameters for Matern model 
CorrParam("Matern")
NuisParam("Gaussian")
# Matern Parameters 
param=list(smooth=0.5,sill=1,scale=0.2/3,nugget=0,mean=0)
corr= GeoCorrFct(x=x, corrmodel="Matern", param=param,model="Gaussian")
plot(x,corr,ylim=c(0,1),type="l")

################################################################
###
### Example 2. Generalized Wendland-Matern correlation model
###
###############################################################
CorrParam("GenWend_Matern")
NuisParam("Gaussian")
# GenWend Matern Parameters 
param=list(smooth=1,sill=1,scale=0.1,nugget=0,power2=1/4,mean=0)
corr= GeoCorrFct(x=x, corrmodel="GenWend_Matern", param=param,model="Gaussian")
plot(x,corr,ylim=c(0,1),type="l")


################################################################
###
### Example 3. Student t correlation with underlying 
### Generalized Wendland-Matern correlation model
###
###############################################################
CorrParam("GenWend_Matern")
NuisParam("StudentT")

x = seq(0,1,0.005)
param=list(smooth=1,sill=1,scale=0.1,nugget=0,power2=1/4,df=1/6,mean=0)
corr= GeoCorrFct(x=x, corrmodel="GenWend_Matern", param=param,model="StudentT")
plot(x,corr,ylim=c(0,1),type="l")

################################################################
###
### Example 4. Weibull correlation with underlying 
### Generalized Cauchy correlation model
###
###############################################################
CorrParam("Matern")
NuisParam("Weibull")

x = seq(0,1,0.005)
param=list(sill=1,scale=0.1,nugget=0,smooth=1.5,shape=4,mean=0)
corr= GeoCorrFct(x=x, corrmodel="Matern", param=param,model="Weibull")
plot(x,corr,ylim=c(0,1),type="l")
################################################################
###
### Example 4. Weibull correlation with underlying 
### Generalized Cauchy correlation model
###
###############################################################
CorrParam("Matern")
NuisParam("Poisson")

x = seq(0,1,0.005)
param=list(sill=1,scale=0.2/3,nugget=0,smooth=0.5,mean=2)
corr= GeoCorrFct(x=x, corrmodel="Matern", param=param,model="Poisson",covariance=TRUE)
plot(x,corr,ylim=c(0,exp(2)),type="l")
}

