\name{GeoNA}
\alias{GeoNA}
\encoding{UTF-8}
\title{Deleting  NA values (missing values) from a spatial or spatio-temporal dataset.}
\description{The function deletes  NA  values from a spatial or spatio-temporal dataset  }
\usage{GeoNA(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, grid=FALSE, X=NULL, setting="spatial")}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)  or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) giving the data.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates. qn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates.   Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates ; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates; the default is \code{NULL}
    then a spatial random field is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL} }
  \item{grid}{Logical; if \code{FALSE} (the default) the data 
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
  \item{X}{Numeric; Matrix of spatio(temporal) covariates in the linear mean specification.}
  \item{setting}{String; are data spatial, spatio-temporal or spatial bivariate
                     (respectively  \code{spatial}, \code{spacetime}, \code{bivariate}) }

}

\value{
  Returns a list containing the following components:
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{data}{The data without NAvalues }
  \item{grid}{\code{TRUE} if the spatial data  are observed in a regular grid, otherwise  \code{FALSE};}
  \item{perc}{The percentage of NA values .}
  \item{setting}{Are data of spatial or spatio-temporal or spatial bivariate type}
  \item{X}{Covariates matrix } 
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\examples{

library(GeoModels)

# Define the spatial-coordinates of the points:
set.seed(79)
x = runif(200, 0, 1)
y = runif(200, 0, 1)
coords=cbind(x,y)
# Set the exponential cov parameters:
corrmodel = "Matern"
mean=0
sill=1
nugget=0
scale=0.3/3
smooth=0.5
param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,
              param=param)$data


data[1:100]=NA

a=GeoNA(data,coordx=coords)
a$perc # percentage of NA values 
#a$coordx# spatial coordinates without missing values
#a$data # data without missinng values


}
\keyword{Composite}
