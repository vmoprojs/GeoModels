\name{GeoVarestbootstrap}  
\alias{GeoVarestbootstrap}
\encoding{UTF-8}
\title{Update a \code{GeoFit} object   using parametric bootstrap for std error estimation}
\description{
  The procedure update a \code{GeoFit} object  computing stderr estimation, confidence intervals
  and p-values  using parametric bootstrap.}
\usage{GeoVarestbootstrap(fit,K=100,sparse=FALSE,
GPU=NULL,local=c(1,1),optimizer="Nelder-Mead",
                     lower=NULL, upper=NULL, 
                     method="cholesky",alpha=0.95, M=30,L=3000)}
\arguments{
  \item{fit}{A fitted object obtained from the
    \code{\link{GeoFit}}.}
     \item{K}{The number of simulations in the parametric bootstrap.}
       \item{sparse}{Logical; if \code{TRUE} then  cholesky decomposition is performed
  using sparse matrices algorithms (spam packake).}
       \item{GPU}{Numeric; if \code{NULL} (the default) 
      no OpenCL computation is performed. The user can choose the device to be used. Use \code{DeviceInfo()} function to see available devices, only double precision devices are allowed} 
        \item{local}{Numeric; number of local work-items of the OpenCL setup}
        \item{optimizer}{The type of optimization algorithm. See  \code{\link{GeoFit}} for details.}
        \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize}.}
       \item{upper}{An optional named list giving the values  for the upper bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize}.}
       \item{method}{String; The method of simulation. Default is \code{cholesky}. For large data set three options are  \code{TB} or \code{CE}  or \code{Vecchia} (see the  \code{GeoSimapprox}) function. } 
         \item{alpha}{ Numeric; The level of the confidence interval. }
         \item{M}{Numeric; the number of neighboords in the Vecchia method.} 
      \item{L}{Numeric; the number of lines in the  turning band method.} 
     }
\details{ The function update a \code{GeoFit} object  estimating stderr estimation 
and confidence interval using parametric bootstrap.}
\value{  
  Returns an (updated) object of class \code{GeoFit}.
}


\seealso{\code{\link{GeoFit}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{


library(GeoModels)

################################################################
###
### Example 1. Test on the parameter
### of a regression model using conditional composite likelihood
###
###############################################################
set.seed(342)
model="Gaussian" 
# Define the spatial-coordinates of the points:
NN=3500
x = runif(NN, 0, 1)
y = runif(NN, 0, 1)
coords = cbind(x,y)
# Parameters
mean=1; mean1=-1.25;  # regression parameters
 sill=1 # variance

# matrix covariates
X=cbind(rep(1,nrow(coords)),runif(nrow(coords)))

# model correlation 
corrmodel="Matern"
smooth=0.5;scale=0.1; nugget=0;

# simulation
param=list(smooth=smooth,mean=mean,mean1=mean1,
              sill=sill,scale=scale,nugget=nugget)
data = GeoSim(coordx=coords, corrmodel=corrmodel,
                model=model, param=param,X=X)$data

I=Inf

fixed=list(nugget=nugget,smooth=smooth)
start=list(mean=mean,mean1=mean1,scale=scale,sill=sill)

lower=list(mean=-I,mean1=-I,scale=0,sill=0)
upper=list(mean=I,mean1=I,scale=I,sill=I)
# Maximum pairwise composite-likelihood fitting of the RF:
fit = GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
              likelihood="Conditional",type="Pairwise",sensitivity=TRUE,
                   lower=lower,upper=upper,neighb=3,
                   optimizer="nlminb",X=X,
                    start=start,fixed=fixed)

unlist(fit$param)


#fit_update=GeoVarestbootstrap(fit,K=100,method='TB')
#fit_update$stderr
#fit_update$conf.int
#fit_update$pvalues

}


\keyword{Composite}
