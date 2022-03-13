\name{GeoAniso}
\alias{GeoAniso}
\encoding{UTF-8}
\title{Spatial Anisotropy correction}
\description{Transforms or back-transforms a set of coordinates according to the geometric anisotropy parameters. }
\usage{GeoAniso(coords, anisopars=c(0,1), inverse = FALSE)}
\arguments{
  \item{coords}{An n x 2 matrix with the coordinates to be transformed.}
    \item{anisopars}{ A bivariate vector with the the anisotropy angle and the anisotropy ratio, respectively. The angle must be given in radians in [0,pi] and the anisotropy ratio must be greater than 1.}
  \item{inverse}{Logical: Default to FALSE. If TRUE the reverse transformation is performed.}

}

\value{Returns a  matrix of transformed coordinates}


\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}


\keyword{Composite}
