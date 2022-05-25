\name{GeoAniso}
\alias{GeoAniso}
\encoding{UTF-8}
\title{Spatial Anisotropy correction}
\description{Transforms or back-transforms a set of coordinates according to the geometric anisotropy parameters. }
\usage{GeoAniso(coords, anisopars=c(0,1), inverse = FALSE)}
\arguments{
  \item{coords}{An n x 2 matrix with the coordinates to be transformed.}
    \item{anisopars}{ A bivariate vector with the the anisotropy angle and the anisotropy ratio, respectively. The angle must be given in radians in [0,pi] and the anisotropy ratio must be greater or equal than 1.}
  \item{inverse}{Logical: Default to FALSE. If TRUE the reverse transformation is performed.}

}

\value{Returns a  matrix of transformed coordinates}



\details{
Geometric anisotropy is defined by a linear  tranformation from the anisotropic space to the isotropic space  that is  \deqn{Y = X R S}
 where  X is a matrix with original coordinates (anisotropic space), and Y is a matrix with  transformed coordinates (isotropic space). 
 Here R  is  a  rotation matrix  with associated anisotropy angle parameter (in [0,pi]) and a S is a shrinking matrix with associated anisotropy ratio 
 parameter (greeater or equal than one).
The two parameters are specified in the anisopars argument as a bivariate numeric vector. The case (0,1) corresponds to the isotropic case.
}
\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}


\keyword{Composite}
