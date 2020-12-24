\name{NuisParam}
\alias{NuisParam}
\encoding{UTF-8}
\title{Lists the Nuisance Parameters of a Random Field}
\description{
 The procedure returns a list with the nuisance parameters of a given
  random field model.
}
\usage{
NuisParam(model, bivariate=FALSE,num_betas=c(1,1))
}
\arguments{
  \item{model}{String; the name of a random field.}
  \item{bivariate}{Logical; if \code{FALSE} (the default) the correlation  model is univariate spatial or spatial-temporal.  
       Otherwise is bivariate.}
  \item{num_betas}{Numerical; the nunber of mean parameters in the linear specification (default is 1)     }
}

\details{The function returns a list with the nuisance parameters of a given
  random field model.}
\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}


\examples{
library(GeoModels)

NuisParam("Gaussian")

NuisParam("Binomial")

NuisParam("Weibull",num_betas=2)

NuisParam("SkewGaussian", num_betas=3)

NuisParam("SinhAsinh")

NuisParam("StudentT")
## note that in the bivariate case sill and nugget are considered as correlation parameteres
NuisParam("Gaussian", bivariate=TRUE)

}
\keyword{Composite}
