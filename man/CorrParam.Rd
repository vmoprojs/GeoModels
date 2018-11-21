\name{CorrParam}
\alias{CorrParam}
\encoding{UTF-8}
\title{Lists the Parameters of a Correlation Model}
\description{
  Subroutine called by InitParam and other procedures.
  The procedure returns a list with the parameters of a given
  correlation model.
}
\usage{
CorrParam(corrmodel)
}
\arguments{
  \item{corrmodel}{String; the name of a correlation model.}
}

\seealso{\code{\link{GeoCovmatrix}}}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
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
### Example 3. Parameters of the  Generalized Wendland  model
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
### Example 5. Parameters of the bi-Matern separable  model
###
###############################################################

CorrParam("Bi_Matern_sep")
}
\keyword{Composite}
