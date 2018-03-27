\name{NuisParam}
\alias{NuisParam}
\encoding{UTF-8}
\title{Lists the Nuisance Parameters of a Random Field}
\description{
  Subroutine called by InitParam and other procedures.
  The procedure returns a list with the nuisance parameters of a given
  random field model.
}
\usage{
NuisParam(model, bivariate,num_betas)
}
\arguments{
  \item{model}{String; the name of a random field.}
  \item{bivariate}{Logical; if \code{FALSE} (the default) the correlation  model is univariate spatial or spatial-temporal.  
       Otherwise is bivariate.}
  \item{num_betas}{Numerical; the nunber of mean parameters in the linear specification.     }
}

\seealso{\code{\link{GeoFit}}}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}
}


\examples{
library(GeoModels)

NuisParam("Gaussian", FALSE,1)
## note that in the bivariate case sill and nugget are considered as correlation parameteres
NuisParam("Gaussian", TRUE,1)

NuisParam("BinomialGaussian", FALSE,1)

NuisParam("Chisq", FALSE,2)

NuisParam("SkewGaussian", FALSE,3)

NuisParam("SinhAsinhGaussian",FALSE,1)
}
\keyword{Composite}
