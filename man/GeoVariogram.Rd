\name{GeoVariogram}
\alias{GeoVariogram}
\encoding{UTF-8}
\title{Empirical Variogram estimation}
\description{
  The function returns an empirical estimate of the  variogram  for spatio (temporal)  and bivariate 
  random fields. 
}
\usage{
GeoVariogram(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, cloud=FALSE,
           distance='Eucl', grid=FALSE, maxdist=NULL,neighb=NULL,
           maxtime=NULL, numbins=NULL, radius=6371,  
           type='variogram',bivariate=FALSE)
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{n \times d}{n x d})-matrix
   (\eqn{n} iid spatial realisations) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or an (\eqn{d \times d \times n}{d x d x n})-array (\eqn{n} iid spatial realisations on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation) or an (\eqn{t \times d \times n }{t x d x n})-array
   (\eqn{n} iid spatial-temporal realisations) or or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) or an (\eqn{d \times d \times t \times n }{d x d x t x n})-array
   (\eqn{n} iid spatial-temporal realisations on regular grid). See \code{\link{GeoFit}} for details.}
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
  \item{cloud}{Logical; if \code{TRUE} the variogram cloud is computed,
  otherwise if \code{FALSE} (the default) the empirical (binned)
  variogram is returned.}
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details} of \code{\link{GeoFit}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set of 
    non-equispaced spatial sites.}
  \item{maxdist}{A numeric value denoting the spatial maximum distance,
   see the Section \bold{Details}.}
   \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhood  in the composite likelihood
    computation. See the Section \bold{Details} for more information.}
   \item{maxtime}{A numeric value denoting the temporal maximum distance,
   see the Section \bold{Details}.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the
    Section \bold{Details}.}
  \item{radius}{Numeric; a value indicating  the radius of the sphere when 
  using the great circle distance. Default value is the radius of the earth in Km (i.e. 6371)} 
  \item{type}{A String denoting the type of variogram. The option
     available is : \code{variogram}.}
    \item{bivariate}{Logical; if \code{FALSE} (the default) the data
    are interpreted as univariate spatial or spatial-temporal realisations. 
    Otherwise they are intrepreted as a a realization from a bivariate field.}
}

\details{

  We briefly report the definitions of semi-variogram used in this function.
  In the case of a spatial Gaussian random field
  the sample \code{variogram} estimator is defined by
  \deqn{\hat{\gamma}(h) = 0.5 \sum_{x_i, x_j \in N(h)} (Z(x_i) - Z(x_j))^2 / |N(h)|}
  where \eqn{N(h)} is the set of all the sample pairs whose distances fall into a tolerance region with size \eqn{h}
  (equispaced intervalls are considered).
  
  In the case of a spatio-temporal Gaussian random field the sample
  \code{variogram} estimator is defined by
  \deqn{\hat{\gamma}(h, u) = 0.5 \sum_{(x_i, l), (x_j, k) \in N(h, u)} (Z(x_i, l) - Z(x_j, k))^2 / |N(h, u)|}
  where \eqn{N(h, u)}{N(h, u)} is the set of all the sample pairs whose
  spatial distances fall into a tolerance region with size \eqn{h}{h}
  and \eqn{|k - l| = u}{\|k-l\|=u}. Note, that \eqn{Z(x_i, l)}{Z(x_i,l)}
  is the observation at site \eqn{x_i}{x_i} and time
  \eqn{l}{l}. 

  The  \code{numbins} parameter indicates the number of adjacent
  intervals to consider in order to grouped distances with which to
  compute the (weighted) lest squares.

  The \code{maxdist} parameter indicates the maximum spatial distance below which
  the shorter distances will be considered in the calculation of
  the (weigthed) least squares.

  The \code{maxtime} parameter indicates the maximum temporal distance below which
  the shorter distances will be considered in the calculation of
  the (weigthed) least squares.
}

\value{
  Returns an object of class \code{Variogram}.
  An object of class \code{Variogram} is a list containing
  at most the following components:

  \item{bins}{Adjacent intervals of grouped spatial distances if
  \code{cloud=FALSE}. Otherwise if \code{cloud=TRUE} all the spatial pairwise distances;}
  \item{bint}{Adjacent intervals of grouped temporal distances if
  \code{cloud=FALSE}. Otherwise if \code{cloud=TRUE} all the temporal pairwise distances;}
  \item{cloud}{If the variogram cloud is returned (\code{TRUE}) or the
  empirical variogram (\code{FALSE});}
  \item{centers}{The centers of the spatial bins;}
  \item{distance}{The type of spatial distance;}
  \item{lenbins}{The number of pairs in each spatial bin;}
  \item{lenbinst}{The number of pairs in each spatial-temporal bin;}
  \item{lenbint}{The number of pairs in each temporal bin;}
  \item{maxdist}{The maximum spatial distance used for the calculation of the variogram.
                 If no spatial distance is specified then it is NULL;}
  \item{maxtime}{The maximum temporal distance used for the calculation of the variogram.
                If no temporal distance is specified then it is NULL;}
  \item{spacetime_dyn}{If the space-time variogram is obtained using dynamical coordinates
  then it is(\code{TRUE}).   }
  \item{variograms}{The empirical spatial variogram;}
  \item{variogramst}{The empirical spatial-temporal variogram;}
  \item{variogramt}{The empirical temporal variogram;}
  \item{type}{The type of estimated variogram}
}

\references{


Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York: Wiley.

  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.
}

\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

 \examples{
library(GeoModels)

################################################################
###
### Example 1. Empirical estimation of the semi-variogram from a
### spatial Gaussian random field with exponential correlation.
###
###############################################################
set.seed(514)
# Set the coordinates of the sites:
x = runif(200, 0, 1)
y = runif(200, 0, 1)
coords = cbind(x,y)
# Set the model's parameters:
corrmodel = "exponential"
mean = 0
sill = 1
nugget = 0
scale = 0.3/3

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel, param=list(mean=mean,
              sill=sill, nugget=nugget, scale=scale))$data

# Empirical spatial semi-variogram estimation:
vario = GeoVariogram(coordx=coords,data=data,maxdist=0.6)

plot(vario,pch=20,ylim=c(0,1))


################################################################
###
### Example 2. Empirical estimation of the variogram from a
### spatio-temporal Gaussian random fields with Gneiting
### correlation function.
###
###############################################################

set.seed(331)
# Define the temporal sequence:
# Set the coordinates of the sites:
x = runif(200, 0, 1)
y = runif(200, 0, 1)
coords = cbind(x,y)
times = seq(1,10,1)

# Simulation of a spatio-temporal Gaussian random field:
data = GeoSim(coordx=coords, coordt=times, corrmodel="gneiting",
              param=list(mean=0,scale_s=0.08,scale_t=0.4,sill=1,
              nugget=0,power_s=1,power_t=1,sep=0.5))$data

# Empirical spatio-temporal semi-variogram estimation:
vario_st = GeoVariogram(data=data, coordx=coords, coordt=times, maxtime=7,maxdist=0.5)

plot(vario_st)
      
################################################################
###
### Example 3. Empirical estimation of the (cross) semivariograms 
### from a bivariate Gaussian random fields with Matern
### correlation function.
###
###############################################################
# Simulation of a bivariate spatial Gaussian random field:
set.seed(293)
# Define the spatial-coordinates of the points:
x = runif(400, 0, 1)
y = runif(400, 0, 1)
coords=cbind(x,y)

# Simulation of a bivariate Gaussian Random field 
# with matern (cross)  covariance function
param=list(mean_1=0,mean_2=0,scale_1=0.1/3,scale_2=0.15/3,scale_12=0.15/3,
           sill_1=1,sill_2=1,nugget_1=0,nugget_2=0,
           smooth_1=0.5,smooth_12=0.5,smooth_2=0.5,pcol=0.3)  
data = GeoSim(coordx=coords, corrmodel="Bi_matern", param=param)$data

# Empirical  semi-(cross)variogram estimation:
biv_vario=GeoVariogram(data,coordx=coords, bivariate=TRUE,maxdist=c(0.5,0.5,0.5))  

plot(biv_vario,pch=20)
      
}

\keyword{Variogram}
