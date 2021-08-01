\name{GeoCovariogram}
\alias{GeoCovariogram}
\encoding{UTF-8}
\title{Computes the fitted  variogram model.}
\description{
  The procedure computes and plots covariance  or variogram estimated fitting a Gaussian,
  and non Gaussian spatio (temporal) bivariate random  fields. Allows to add 
  the empirical estimates in order to compare them with the fitted model.
}
\usage{
GeoCovariogram(fitted, distance="Eucl",answer.cov=FALSE,
            answer.vario=FALSE, answer.range=FALSE, fix.lags=NULL, 
             fix.lagt=NULL, show.cov=FALSE, show.vario=FALSE, 
            show.range=FALSE, add.cov=FALSE, add.vario=FALSE,
            pract.range=95, vario, \dots)
}
\arguments{
  \item{fitted}{A fitted object obtained from the
    \code{\link{GeoFit}} or \code{\link{GeoWLS}} procedures.}
  \item{distance}{String; the name of the spatial distance. The default
  is \code{Eucl}, the euclidean distance. See \code{\link{GeoFit}}.}
  \item{answer.cov}{Logical; if \code{TRUE} a vector with the
    estimated covariance function is returned; if \code{FALSE} (the
    default) the covariance is not returned.}
  \item{answer.vario}{Logical; if \code{TRUE} a vector with the
    estimated variogram is returned; if \code{FALSE} (the
    default) the variogram is not returned.}
  \item{answer.range}{Logical; if \code{TRUE} the estimated
    pratical range is returned; if \code{FALSE} (the
    default) the pratical range is not returned.}
  \item{fix.lags}{Integer; a positive value denoting the spatial lag to
    consider for the plot of the temporal profile.}
  \item{fix.lagt}{Integer; a positive value denoting the temporal lag to
    consider for the plot of the spatial profile.}
  \item{show.cov}{Logical; if \code{TRUE} the estimated
    covariance function is plotted; if \code{FALSE} (the
    default) the covariance function is not plotted.}
  \item{show.vario}{Logical; if \code{TRUE} the estimated
    variogram is plotted; if \code{FALSE} (the
    default) the variogram is not plotted.}
  \item{show.range}{Logical; if \code{TRUE} the estimated
    pratical range is added on the plot; if \code{FALSE} (the
    default) the pratical range is not added.}
  \item{add.cov}{Logical; if \code{TRUE} the vector of the
    estimated covariance function is added on the current plot; if \code{FALSE} (the
    default) the covariance is not added.}
  \item{add.vario}{Logical; if \code{TRUE} the vector with the
    estimated variogram is added on the current plot; if \code{FALSE} (the
    default) the correlation is not added.}
  \item{pract.range}{Numeric; the percent of the sill to be reached.}
  \item{vario}{A \code{Variogram} object obtained from the
    \code{\link{GeoVariogram}} procedure.}
  \item{\dots}{other optional parameters which are passed to plot functions.}
}

\details{The function computes the fitted variogram model}
\value{
  The returned object is eventually a list with:

    \item{covariance}{The vector of the estimated covariance function;}
    \item{variogram}{The vector of the estimated variogram function;}
   
    
}

\references{
Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York:
Wiley.

Gaetan, C. and Guyon, X. (2010) Spatial Statistics and Modelling.
  \emph{Spring Verlang, New York}.
}

\seealso{\code{\link{GeoFit}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}



\examples{
library(GeoModels)
library(scatterplot3d)

################################################################
###
### Example 1. Plot of fitted covariance and fitted 
### and empirical variograms from a  Gaussian RF 
### with Matern correlation. 
###
###############################################################
set.seed(21)
# Set the coordinates of the points:
x = runif(300, 0, 1)
y = runif(300, 0, 1)
coords=cbind(x,y)

# Set the model's parameters:
corrmodel = "Matern"
model = "Gaussian"
mean = 0
sill = 1
nugget = 0
scale = 0.2/3
smooth=0.5

param=list(mean=mean,sill=sill, nugget=nugget, scale=scale, smooth=smooth)
# Simulation of the Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel, model=model,param=param)$data

start=list(mean=0,scale=scale,sill=sill)
fixed=list(nugget=nugget,smooth=smooth)
# Maximum composite-likelihood fitting of the Gaussian random field:
fit = GeoFit(data=data,coordx=coords, corrmodel=corrmodel,model=model,
            likelihood="Marginal",type='Pairwise',start=start,
            optimizer="BFGS", fixed=fixed,neighb=4)

# Empirical estimation of the variogram:
vario = GeoVariogram(data=data,coordx=coords,maxdist=0.5)

# Plot of covariance and variogram functions:
GeoCovariogram(fit, show.cov=TRUE,show.vario=TRUE, vario=vario,pch=20)

################################################################
###
### Example 2. Plot of fitted covariance and fitted 
### and empirical variograms from a  Binomial  
###  RF with exponential correlation. 
###
###############################################################
set.seed(2111)

model="Binomial";n=20
# Set the coordinates of the points:
x = runif(500, 0, 1)
y = runif(500, 0, 1)
coords=cbind(x,y)

# Set the model's parameters:
corrmodel = "exponential"
mean = 0
sill = 1
nugget = 0
scale = 0.2/3

param=list(mean=mean,sill=sill, nugget=nugget, scale=scale)
# Simulation of the Gaussian RF:
data = GeoSim(coordx=coords, corrmodel=corrmodel, model=model,param=param,n=n)$data

start=list(mean=0,scale=scale,sill=sill)
fixed=list(nugget=nugget)
# Maximum composite-likelihood fitting of the BinomGaussian random field:
fit = GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,
             likelihood="Marginal",type='Pairwise',start=start,n=n,
             optimizer="BFGS", fixed=fixed,neighb=4)

# Empirical estimation of the variogram:
vario = GeoVariogram(data,coordx=coords,maxdist=0.5)

# Plot of covariance and variogram functions:
GeoCovariogram(fit, show.cov=TRUE,show.vario=TRUE, vario=vario,pch=20)


################################################################
###
### Example 3. Plot of fitted covariance and fitted 
### and empirical variograms from a  RF
###  RF with Wend0 correlation. 
###
###############################################################
set.seed(211)

model="Gamma";shape=4
# Set the coordinates of the points:
x = runif(700, 0, 1)
y = runif(700, 0, 1)
coords=cbind(x,y)

# Set the model's parameters:
corrmodel = "Wend0"
mean = 0
sill = 1
nugget = 0
scale = 0.3
power2=4


param=list(mean=mean,sill=sill, nugget=nugget, scale=scale,shape=shape,power2=power2)
# Simulation of the Gaussian RF:
data = GeoSim(coordx=coords, corrmodel=corrmodel, model=model,param=param)$data

start=list(mean=0,scale=scale,shape=shape)
fixed=list(nugget=nugget,sill=sill,power2=power2)
# Maximum composite-likelihood fitting of the BinomGaussian random field:
fit = GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,
             likelihood="Marginal",type='Pairwise',start=start,
             optimizer="BFGS", fixed=fixed,neighb=4)

# Empirical estimation of the variogram:
vario = GeoVariogram(data,coordx=coords,maxdist=0.5)

# Plot of covariance and variogram functions:
GeoCovariogram(fit, show.cov=TRUE,show.vario=TRUE, vario=vario,pch=20)

################################################################
###
### Example 4. Plot of  fitted  and empirical variograms
### from a space time Gaussian random fields  
### with double exponential correlation. 
###
###############################################################
set.seed(92)
# Define the spatial-coordinates of the points:
x = runif(50, 0, 1)
y = runif(50, 0, 1)
coords=cbind(x,y)
# Define the temporal sequence:
time = seq(0, 15, 1)


# Simulation of the spatio-temporal Gaussian random field:
data = GeoSim(coordx=coords, coordt=time, corrmodel="Exp_Exp",param=list(mean=mean,
              nugget=nugget,scale_s=0.5/3,scale_t=2/2,sill=sill))$data

fixed=list(nugget=0, mean=0)
start=list(scale_s=0.2, scale_t=0.5, sill=1)
# Maximum composite-likelihood fitting of the space-time Gaussian random field:
fit = GeoFit(data, coordx=coords, coordt=time, corrmodel="Exp_Exp", maxtime=2,
             neighb=4, likelihood="Marginal", type="Pairwise", 
             optimizer="BFGS",fixed=fixed, start=start)

# Empirical estimation of spatio-temporal covariance:
vario = GeoVariogram(data,coordx=coords, coordt=time, maxtime=5,maxdist=0.5)

# Plot of the fitted space-time variogram
GeoCovariogram(fit,vario=vario,show.vario=TRUE)

# Plot of covariance, variogram and spatio and temporal profiles:
GeoCovariogram(fit,vario=vario,fix.lagt=1,fix.lags=1,show.vario=TRUE,pch=20)




################################################################
###
### Example 5. Plot of parametric and empirical variograms
### estimated from a Bivariate Gaussian random fields with 
### Matern correlation. 
###
###############################################################

# Simulation of a bivariate spatial Gaussian random field:
set.seed(892)
# Define the spatial-coordinates of the points:
x = runif(200, -1, 1)
y = runif(200, -1, 1)
coords=cbind(x,y)

# Simulation of a bivariate Gaussian Random field 
# with matern (cross)  covariance function
scale_1 = 0.25/3
scale_2 = 0.2/3
scale_12 = 0.15/3
sill_1=1
sill_2=1
smooth=0.5
pcol=0.3
param=list(mean_1=0,mean_2=0,scale_1=scale_1,scale_2=scale_2,scale_12=scale_12,
           sill_1=sill_1,sill_2=sill_2,nugget_1=0,nugget_2=0,
           smooth_1=smooth,smooth_12=smooth,smooth_2=smooth,pcol=pcol)  
data = GeoSim(coordx=coords, corrmodel="Bi_Matern", param=param)$data

# Empirical bivariate variogram estimation:
biv_vario=GeoVariogram(data,coordx=coords, bivariate=TRUE,maxdist=c(0.5,0.5,0.5)) 

# selecting fixed and estimating parameters
fixed=list(mean_1=0,mean_2=0,nugget_1=0,nugget_2=0,
          smooth_1=smooth,smooth_12=smooth,smooth_2=smooth)
start=list(sill_1=var(data[1,]),sill_2=var(data[2,]),
          scale_1=scale_1,scale_2=scale_2,scale_12=scale_12,
          pcol=cor(data[1,],data[2,]))

# Maximum likelihood fitting of the bivariate random field:
fit= GeoFit(data, coordx=coords, corrmodel="Bi_Matern",likelihood="Marginal",
            optimizer="BFGS", type="Pairwise",
            start=start,fixed=fixed,neighb=c(3,3,3))
            
GeoCovariogram(fit, vario=biv_vario,show.vario=TRUE,pch=20)


}

\keyword{Composite}
