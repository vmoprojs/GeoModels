\name{GeoPit}
\alias{GeoPit}
\encoding{UTF-8}
\title{Probability integral  or normal score tranformation}
\description{
   The procedure  for a given GeoFit object applies  the probability integral tranformation or the normal score transformation to the data
}
\usage{GeoPit(fit,type="Uniform")}
\arguments{
  \item{fit}{A GeoFit object}.
\item{type}{The type of transformation. If "Uniform" then  the probability integral tranformation is performed.
If "Gaussian" then  the normal score transformation is performed. }
}

\value{Returns an (updated) object of class \code{GeoFit}}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\examples{

library(GeoModels)

model="Beta2"
copula="Clayton"

set.seed(221)
NN=800
x <- runif(NN);y <- runif(NN)
coords=cbind(x,y)


shape=1.5
scale=0.2;power2=4
smooth=0
nugget=0
nu=8

corrmodel="GenWend"

min=-2;max=1
mean=0


param=list(smooth=smooth,power2=power2, min=min,max=max,
             mean=mean, nu=nu,
             scale=scale,nugget=nugget,shape=shape)

optimizer="nlminb"

data <- GeoSimCopula(coordx=coords, corrmodel=corrmodel, 
model=model,param=param,copula=copula)$data

I=50
fixed<-list(nugget=nugget,sill=1,scale=scale,smooth=smooth,power2=power2,min=min,max=max,nu=nu)
start<-list(shape=shape,mean=mean)
lower<-list(shape=0,mean=-I)
upper<-list(shape=10,mean=I)

#### maximum independence likelihood
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
model=model,likelihood="Marginal",type="Independence",
                      optimizer=optimizer,lower=lower,
                      upper=upper,copula=copula,
                    start=start,fixed=fixed)

## PIT transformation
aa=GeoPit(fit1,type="Uniform")
hist(aa$data,freq=FALSE)
GeoScatterplot(aa$data,coords,neighb=c(1,2))
## Normal score transformation
bb=GeoPit(fit1,type="Gaussian")
hist(bb$data,freq=FALSE)
GeoScatterplot(bb$data,coords,neighb=c(1,2))


}

\keyword{Composite}
