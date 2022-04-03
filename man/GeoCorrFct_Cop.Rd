\name{GeoCorrFct_Cop}
\alias{GeoCorrFct_Cop}
\encoding{UTF-8}
\title{Spatial and Spatio-temporal correlation or  covariance of (non) Gaussian random fields (copula models)}
\description{
  The function computes the  correlations 
  a spatial or spatio-temporal or a bivariate spatial Gaussian or non Gaussian randomm field
  witha  given  covariance model and a  set of spatial (temporal) distances.
}
\usage{
GeoCorrFct_Cop(x,t=NULL,corrmodel, model="Gaussian",copula="Gaussian",distance="Eucl", param, radius=6371,n=1,covariance=FALSE,variogram=FALSE)
}

\arguments{
  \item{x}{A set of spatial distances.}
  \item{t}{A set of (optional) temporal distances.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{model}{String; the type of RF. See \code{\link{GeoFit}}.}
  \item{copula}{String; the type of copula. The two options are Gaussian and Clayton.}
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
     \item{variogram}{Logic; if FALSE then the covariance/coorelation is returned. Otherwise the associated
   semivariogram is returned}
   }


\value{Returns  correlations or covariances values associated to a given parametric spatial and temporal correlation models.}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\examples{
library(GeoModels)

################################################################
###
### Example 1. Correlation of a (mean reparametrized) beta random field with underlying 
### Matern correlation model using Gaussian and Clayton copulas
###
###############################################################

# Define the spatial distances
x = seq(0,0.4,0.005)

# Correlation Parameters for Matern model 
CorrParam("Matern")
NuisParam("Beta2")
# corr Gaussian copula
param=list(smooth=0.5,sill=1,scale=0.2/3,nugget=0,mean=0,min=0,max=1,shape=0.5)
corr1= GeoCorrFct_Cop(x=x, corrmodel="Matern", param=param,copula="Gaussian",model="Beta2")

# corr Clayton copula
param=list(smooth=0.5,sill=1,scale=0.2/3,nugget=0,mean=0,min=0,max=1,shape=0.5,nu=2)
corr2= GeoCorrFct_Cop(x=x, corrmodel="Matern", param=param,copula="Clayton",model="Beta2")

plot(x,corr1,ylim=c(0,1),type="l",xlab="distance",ylab="corr")
lines(x,corr2,ylim=c(0,1),lty=2)




}

