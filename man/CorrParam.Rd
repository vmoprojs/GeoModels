\name{CorrParam}
\alias{CorrParam}
\encoding{UTF-8}
\title{Lists the Parameters of a Correlation Model}
\description{
  The procedure returns a list with the parameters of a given
  correlation model.
}
\usage{
CorrParam(corrmodel)
}

\arguments{
  \item{corrmodel}{String: the name associated  to a given correlation model.}
}

\value{Return a vector string of correlation parameters.}

\details{The function return a list with the Parameters of a Correlation Model}

\seealso{\code{\link{GeoCovmatrix}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}

\examples{
require(GeoModels)
################################################################
###
### Example 1. Parameters of the  Matern  model
###
###############################################################

CorrParam("Matern")


################################################################
###
### Example 2. Parameters of the  Generalized Wendland  model
###
###############################################################

CorrParam("GenWend")


################################################################
###
### Example 3. Parameters of the  Generalized Cauchy  model
###
###############################################################

CorrParam("GenCauchy")


################################################################
###
### Example 4. Parameters of the  space time Gneiting  model
###
###############################################################

CorrParam("Gneiting")


################################################################
###
### Example 5. Parameters of the bi-Matern separable  model.
###            Note that in the bivariate case variance paramters are
###            included
###############################################################

CorrParam("Bi_Matern_sep")
}
\keyword{Composite}
