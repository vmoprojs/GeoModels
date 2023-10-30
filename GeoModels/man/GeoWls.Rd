\name{GeoWLS}
\alias{GeoWLS}
\alias{print.GeoWLS}
\encoding{UTF-8}
\title{WLS of Random Fields}
\description{
  the function returns the parameters' estimates of a random field obtained by 
  the weigthed least squares estimator.
}
\usage{
GeoWLS(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, 
             distance="Eucl", fixed=NULL, grid=FALSE, maxdist=NULL,neighb=NULL,
             maxtime=NULL,  model='Gaussian', optimizer='Nelder-Mead',
             numbins=NULL, radius=6371, start=NULL, weighted=FALSE)
}
\arguments{
 \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation)  or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or an (\eqn{d \times d \times n}{d x d x n})-array (\eqn{n} iid spatial realisations on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation) or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid). See \code{\link{GeoFit}} for details.}
 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL} then a spatial random field is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}
    }  
  \item{corrmodel}{String; the name of a correlation model, for the
    description (see \code{\link{GeoFit}}).}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl},
                  the euclidean distance. See the Section \bold{Details}
		  of \code{\link{GeoFit}}.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{maxdist}{A numeric value denoting the maximum distance, see
    \bold{Details} in \code{\link{GeoFit}}.}
      \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhood. See \bold{Details} and \code{\link{GeoFit}}}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal lag considered.See \bold{Details} and \code{\link{GeoFit}}.}
   \item{model}{String; the type of random field. \code{Gaussian} is the default, see
    \code{\link{GeoFit}} for the different types.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the
    Section \bold{Details}}
       \item{radius}{Numeric; a value indicating  the radius of the sphere when using the great circle distance.
          Default value is the radius of the earth in Km (i.e. 6371)} 
  \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \code{\link{GeoFit}}).}
  \item{weighted}{Logical; if \code{TRUE} then the weighted least square
    estimator is considered. If \code{FALSE} (the
    default) then the classic least square is used.}
}

\details{
  The  \code{numbins} parameter indicates the number of adjacent
  intervals to consider in order to grouped distances with which to
  compute the (weighted) lest squares.

  The \code{maxdist} parameter indicates the maximum distance below which
  the shorter distances will be considered in the calculation of
  the (weigthed) least squares.
}

\value{
   Returns an object of class \code{WLS}.
  An object of class \code{WLS} is a list containing
  at most the following components:

  \item{bins}{Adjacent intervals of grouped distances;}
  \item{bint}{Adjacent intervals of grouped temporal separations}
  \item{centers}{The centers of the bins;}
  \item{coordx}{The vector or matrix of spatial coordinates;}
  \item{coordy}{The vector of spatial coordinates;}
  \item{coordt}{The vector of temporal coordinates;}
  \item{convergence}{A string that denotes if convergence is reached;}
  \item{corrmodel}{The correlation model;}
  \item{data}{The vector or matrix of data;}
  \item{distance}{The type of spatial distance;}
  \item{fixed}{The vector of fixed parameters;}
  \item{iterations}{The number of iteration used by the numerical routine;}
  \item{maxdist}{The maximum spatial distance used for the calculation of the variogram used in least square estimation.
                 If no spatial distance is specified then it is NULL;}
  \item{maxtime}{The maximum temporal distance used for the calculation of the variogram used in least square estimation.
                 If no temporal distance is specified then it is NULL;}
  \item{message}{Extra message passed from the numerical routines;}
  \item{model}{The type of random fields;}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number the temporal realisations of the random field;}
  \item{param}{The vector of  parameters' estimates;}
  \item{variograms}{The empirical spatial variogram;}
  \item{variogramt}{The empirical temporal variogram;}
  \item{variogramst}{The empirical spatial-temporal variogram;}
  \item{weighted}{A logical value indicating if its the weighted
    method;}
  \item{wls}{The value of the least squares at the minimum.}

}

\references{

  Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York:
  Wiley.

  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.

}

\seealso{\code{\link{GeoFit}}, \code{\link{optim}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

  
\examples{
library(GeoModels)


# Set the coordinates of the sites:

set.seed(211)
x <- runif(200, 0, 1)
set.seed(98)
y <- runif(200, 0, 1)
coords <- cbind(x,y)

################################################################
###
### Example 1. Least square fitting of a Gaussian random field
### with exponential correlation.
###
###############################################################

# Set the model's parameters:
corrmodel <- "Exponential"
mean <- 0
sill <- 1
nugget <- 0
scale <- 0.15/3
param <- list(mean=0,sill=sill, nugget=nugget, scale=scale)
# Simulation of the Gaussian random field:
set.seed(2)
data <- GeoSim(coordx=coords, corrmodel=corrmodel, param=param)$data

fixed=list(nugget=0,mean=mean)
start=list(scale=scale,sill=sill)
# Least square fitting of the random field:
fit <- GeoWLS(data=data,coordx=coords, corrmodel=corrmodel,
         fixed=fixed,start=start,maxdist=0.5)

# Results:
print(fit)

}

\keyword{LeastSquare}
