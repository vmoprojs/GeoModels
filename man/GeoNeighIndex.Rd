\name{GeoNeighIndex}
\alias{GeoNeighIndex}
\encoding{UTF-8}
\title{Spatial or spatiotemporal  near neighbour indices.}
\description{
  The function return  the spatial or spatiotemporal  near neighbour indices for a given neighbour  and temporal distance. 
}
\usage{
GeoNeighIndex(coordx,coordy=NULL,coordx_dyn=NULL,coordt=NULL,distance="Eucl",neighb=5,maxtime=1,radius=6371,bivariate=FALSE)
}
\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) assigning 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector assigning
    1-dimension of spatial coordinates.  Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL} then a spatial random field is expected.}
   \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }   
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details} of \code{\link{GeoFit}}.}
\item{neighb}{Numeric; an optional (vector of) positive integer indicating the 
    order of neighborhood. See the Section \bold{Details} for more information.}
   \item{maxtime}{A numeric value denoting the temporal distance \bold{Details}.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere when 
  using the great circle distance. Default value is the radius of the earth in Km (i.e. 6371)} 
    \item{bivariate}{Logical; if \code{FALSE} (the default) the data
    are interpreted as univariate spatial or spatial-temporal realisations. 
    Otherwise they are intrepreted as a a realization from a bivariate field.}
}

\details{
 The function return the  spatial or spatiotemporal indices of the pairs tha are neighboord of a certain order
}


\value{
  Returns a list containing the following components:

  \item{colidx}{First vector of indices}
  \item{rowidx}{Second vector of indices}
  \item{lags}{Vector of spatial distances}
  \item{lagt}{Vector of temporal distances}
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}


\examples{
library(GeoModels)

################################################################
###############################################################
set.seed(514)

require(GeoModels)
NN = 800
coords = cbind(runif(NN),runif(NN))
scale=0.5/3
param = list(mean=0,sill=1,nugget=0,scale=0.5/3,smooth=0.5)
corrmodel = "Matern"; 

param = list(mean=0,sill=1,nugget=0,scale=scale,smooth=0.5)
data = GeoSim(coordx = coords,corrmodel = corrmodel,
                  model = "Gaussian",param = param)$data

sel=GeoNeighIndex(coordx=coords,neighb=10)        

data1=data[sel$colidx]; data2=data[sel$rowidx]
## plotting pairs  that are neighboord of order 5
plot(data1,data2,xlab="",ylab="",main="h-scatterplot, neighb=10")

}
\keyword{Scatterplot}
