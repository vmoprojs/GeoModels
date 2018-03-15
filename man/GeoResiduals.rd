\name{GeoResiduals}
\alias{GeoResiduals}
\encoding{UTF-8}
\title{Computes fitted covariance and/or  variogram}
\description{
  The procedure return a GeoFit object associated to  the estimated residuals
}
\usage{GeoResiduals(fit)}
\arguments{
  \item{fit}{A fitted object obtained from the
    \code{\link{GeoFit}}.}
}

\value{
  A GeoFit object with the estimated residuals  
}


\seealso{\code{\link{GeoFit}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales-Oñate, \email{victor.morales@uv.cl}
}



\examples{
library(GeoModels)
set.seed(211)

model="Weibull";shape=4
# Set the coordinates of the points:
x <- runif(700, 0, 1)
y <- runif(700, 0, 1)
coords=cbind(x,y)

# Set the model's parameters:
corrmodel <- "Wend0"
mean <- 5
sill <- 1
nugget <- 0
scale <- 0.3
power2=4


param=list(mean=mean,sill=sill, nugget=nugget, scale=scale,shape=shape,power2=power2)
# Simulation of the Gaussian RF:
data <- GeoSim(coordx=coords, corrmodel=corrmodel, model=model,param=param)$data

start=list(mean=mean,scale=scale,shape=shape)
fixed=list(nugget=nugget,sill=sill,power2=power2)
# Maximum composite-likelihood fitting of the BinomGaussian random field:
fit <- GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,
                    likelihood="Marginal",type='Pairwise',start=start,
                    fixed=fixed,maxdist=0.1)

res=GeoResiduals(fit)
# Empirical estimation of the variogram for the residuals:
vario <- GeoVariogram(res$data,coordx=coords,maxdist=0.5)

# Plot of covariance and variogram functions:
GeoCovariogram(res, show.vario=TRUE, vario=vario,pch=20)

}

\keyword{Composite}
