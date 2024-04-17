\name{sp2Geo}
\alias{sp2Geo}
\encoding{UTF-8}
\title{Extracting information from an sp or spacetime object}
\description{Extracting information from an sp or spacetime object}
\usage{sp2Geo(spobj,spdata = NULL)}
\arguments{
    \item{spobj}{An object of class sp or spacetime}
   \item{spdata}{Character: The name of data in the sp or spacetime object}

}
\value{A list with spatio-temporal informations}


\details{
The function accepts as input a sp or spacetime object and the name of the  data of interest in the object  and it returns a list with some  useful informations  for Geomodels functions.
}
\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}
\examples{
# Define the spatial-coordinates of the points:
set.seed(3)
N <- 30  # number of location sites
x <- runif(N, 0, 1)
set.seed(6)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

# Define spatial matrix covariates and regression parameters
X <- cbind(rep(1,N),runif(N))
# Define spatial matrix dependent variable
Y <- rnorm(nrow(X))

obj1 <- sp::SpatialPoints(coords)
obj2 <- sp::SpatialPointsDataFrame(coords,data = data.frame(X,Y))

# sp2Geo info extraction
b <- sp2Geo(obj2,spdata  = "Y")
class(b)
b
}

\keyword{Composite}
