\name{GeoScores}
\alias{GeoScores}
\encoding{UTF-8}
\title{Computation of predictive scores}
\description{
  The function computes some predictive scores for a spatial, spatiotemporal and bivariate Gaussian RFs}
\usage{GeoScores(data_to_pred, probject=NULL,pred=NULL,mse=NULL,
   score=c("brie","crps","lscore","pit","pe"))
}

\arguments{
  \item{data_to_pred}{A numeric vector of data to predict about a response} 
  \item{probject}{A Geokrig object obtained using the function Geokrig}
  \item{pred}{A numeric vector with predictions for the response.}
  \item{mse}{a numeric vector with prediction variances.}
  \item{score}{A  character  defining what statistic of the prediction errors should be
computed. Possible values are lscore, crps, brie  and pe. In the latter case 
scores based on prediction errors such as rmse, mae, mad are computed. Finally, the character pit allows to compute 
the probability integral transform for each value} 
}
\details{
GeoScores computes the items required to evaluate  the diagnostic criteria
proposed by Gneiting et al. (2007) for assessing the calibration and the sharpness of probabilistic predictions of (cross-) validation data. To this aim, GeoScores uses the assumption that the prediction errors  are Gaussian with zero mean and standard
deviations equal to the Kriging standard errors. This assumption is an approximation if the errors are not Gaussian.}

\value{
  Returns a list  containing the following informations:
  \item{LSCORE}{Logarithmic predictive score}
  \item{CRPS}{Continuous ranked probability predictive score}
  \item{RMSE}{Root mean squared error}
  \item{MAE}{Mean absolute error}
  \item{MAD}{Median absolute error}
  \item{PIT}{A vector of  probability integral transformation }
}


\references{
Gneiting T. and Raftery A. \emph{Strictly Proper Scoring Rules, 
   Prediction, and Estimation}.
  Journal of the American Statistical Association, \bold{102}
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{

library(GeoModels)

################################################################
######### Examples of predictive score computation  ############
################################################################

library(GeoModels)
model="Gaussian"
set.seed(79)
N=1000
x = runif(N, 0, 1)
y = runif(N, 0, 1)
coords=cbind(x,y)

# Set the exponential cov parameters:
corrmodel = "GenWend"
mean=0; sill=5; nugget=0
scale=0.2;smooth=0;power2=4

param=list(mean=mean,sill=sill,nugget=nugget,scale=scale,smooth=smooth,power2=power2)

# Simulation of the spatial Gaussian random field:
data = GeoSim(coordx=coords, corrmodel=corrmodel,
              param=param)$data


sel=sample(1:N,N*0.8)
coords_est=coords[sel,]; coords_to_pred=coords[-sel,]
data_est=data[sel]; data_to_pred=data[-sel]

## estimation with pairwise likelihood
fixed=list(nugget=nugget,smooth=0,power2=power2)
start=list(mean=0,scale=scale,sill=1)
I=Inf
lower=list(mean=-I,scale=0,sill=0)
upper=list(mean= I,scale=I,sill=I)
# Maximum pairwise likelihood fitting :
fit = GeoFit(data_est, coordx=coords_est, corrmodel=corrmodel,model=model,
                    likelihood='Marginal', type='Pairwise',neighb=3,
                    optimizer="nlminb", lower=lower,upper=upper,
                    start=start,fixed=fixed)
# locations to predict
xx=seq(0,1,0.03)
loc_to_pred=as.matrix(expand.grid(xx,xx))

pr=GeoKrig(loc=coords_to_pred,coordx=coords_est,corrmodel=corrmodel,
       model=model,param= param, data=data_est,mse=TRUE)

Pr_scores =GeoScores(data_to_pred,pred=pr$pred,mse=pr$mse)
Pr_scores$rmse;Pr_scores$brie
hist(Pr_scores$pit,freq=FALSE)
}
\keyword{Predictive scores}
