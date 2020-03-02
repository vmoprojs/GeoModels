\name{GeoQQ}
\alias{GeoQQ}
\encoding{UTF-8}
\title{Quantile-quantile plot of residuals}
\description{
   The procedure  plots  a quantile-quantile plot for the residuals associated to a fitted model
}
\usage{GeoQQ(fit)}
\arguments{
  \item{fit}{A GeoFit object obtained from 
    \code{\link{GeoResiduals}}.}
}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}


\examples{
library(GeoModels)
set.seed(211)

model="Weibull";shape=4
N=700 # number of location sites
# Set the coordinates of the points:
x = runif(N, 0, 1)
y = runif(N, 0, 1)
coords=cbind(x,y)


# regression parameters
mean = 5
mean1=0.8

X=cbind(rep(1,N),runif(N))
# correlation parameters:
corrmodel = "Wend0"
sill = 1
nugget = 0
scale = 0.3
power2=4


param=list(mean=mean,mean1=mean1, sill=sill, nugget=nugget, 
	           scale=scale,shape=shape,power2=power2)
# Simulation of the Gaussian RF:
data = GeoSim(coordx=coords, corrmodel=corrmodel, X=X,model=model,param=param)$data

start=list(mean=mean,mean1=mean1, scale=scale,shape=shape)
fixed=list(nugget=nugget,sill=sill,power2=power2)
# Maximum composite-likelihood fitting 
fit = GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,X=X,
                    likelihood="Marginal",type='Pairwise',start=start,
                    fixed=fixed,maxdist=0.1)

res=GeoResiduals(fit)

GeoQQ(res)
}

\keyword{Composite}
