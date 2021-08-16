\name{DeviceInfo}
\alias{DeviceInfo}
\encoding{UTF-8}
\title{Prints Device Information}
\description{
  Prints the device details available in your computer. Device name, Max compute units, whether it supports double precision, among others. 
}
\usage{
DeviceInfo()
}


\details{
  The user can take this information into account so that the \code{local} parameter is set up in \code{GeoFit} when \code{GPU} computation is chosen.
}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\examples{

library(GeoModels)
DeviceInfo()

}

\keyword{Devices}