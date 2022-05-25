\name{GeoQQ}
\alias{GeoQQ}
\encoding{UTF-8}
\title{Quantile-quantile plot }
\description{
   Based on a GeoFit object, the procedure  plots  a quantile-quantile plot or compares the fitted  density with the histogram of the data. It is useful as diagnostic tool.
}
\usage{GeoQQ(fit,type="Q",add=FALSE,ylim=c(0,1),breaks=10,\dots)}

\arguments{
  \item{fit}{A GeoFit object possibly obtained from \code{\link{GeoResiduals}}.}
\item{type}{The type of plot. If Q then a qq-plot (default) is performed. If D then  a comparison between histrogram and the estimated marginal density is performed}
\item{add}{Logical; if TRUE the the estimated density ia added over an existing one} 
\item{ylim}{Numeric; a vector of length 2 used for the ylab parameter of the histogram plot.} 
\item{breaks}{Numeric; an integer number specifyng  the number of cells ofthe histogram plot if the option type=D is chosen.} 
\item{\dots}{Optional parameters passed to the plot function.}
}

\value{ 
Produces a plot. No values are returned.
}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{
library(GeoModels)


##################
### Example 1
##################
set.seed(21)
model="Tukeyh";tail=0.1
N=400 # number of location sites
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
	           scale=scale,tail=tail,power2=power2)
# Simulation of the Gaussian RF:
data = GeoSim(coordx=coords, corrmodel=corrmodel, X=X,model=model,param=param)$data

start=list(mean=mean,mean1=mean1, scale=scale,tail=tail)
fixed=list(nugget=nugget,sill=sill,power2=power2)
# Maximum composite-likelihood fitting 
fit = GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,X=X,
                    likelihood="Conditional",type='Pairwise',start=start,
                    fixed=fixed,neighb=4)

res=GeoResiduals(fit)
GeoQQ(res,type="Q")
GeoQQ(res,type="D",col="red",lwd=2,ylim=c(0,0.5),breaks=20)


##################
### Example 2
##################
set.seed(21)
model="Weibull";shape=1.5
N=600 # number of location sites
# Set the coordinates of the points:
x = runif(N, 0, 1)
y = runif(N, 0, 1)
coords=cbind(x,y)


# regression parameters
mean = 0

# correlation parameters:
corrmodel = "Matern"
smooth=0.5
nugget = 0
scale = 0.2/3


param=list(mean=mean, sill=1, nugget=nugget, 
             scale=scale,smooth=smooth, shape=shape)
# Simulation of the Gaussian RF:
data = GeoSim(coordx=coords, corrmodel=corrmodel,model=model,param=param)$data

start=list(mean=mean, scale=scale,shape=shape)
I=Inf
lower=list(mean=-I, scale=0,shape=0)
upper=list(mean= I, scale=I,shape=I)
I=Inf
fixed=list(nugget=nugget,sill=1,smooth=smooth)
# Maximum composite-likelihood fitting 
fit = GeoFit(data,coordx=coords, corrmodel=corrmodel,model=model,
                    likelihood="Conditional",type='Pairwise',start=start,
                    optimizer="nlminb",lower=lower,upper=upper,
                    fixed=fixed,neighb=3)
GeoQQ(fit,type="Q")
GeoQQ(fit,type="D",col="red",lwd=2,ylim=c(0,1),breaks=20)
}

\keyword{Composite}
