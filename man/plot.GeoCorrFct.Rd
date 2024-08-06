\name{plot.GeoCorrFct}
\alias{plot.GeoCorrFct}
\encoding{UTF-8}
\title{Plot Spatial and Spatio-temporal correlation or  covariance of (non) Gaussian random fields}

\description{
  Plot Spatial and Spatio-temporal correlation or  covariance of (non) Gaussian random fields
  for a given set of spatial or spatiotemporal distances
    \code{\link{GeoCorrFct}}.
}

\usage{
\method{plot}{GeoCorrFct}(x,type="p", \dots)
}

\arguments{
  \item{x}{an object of the class \code{"GeoCorrFct"} }
    \item{type}{The type of graphic.  The possible options are  "p" and "l". If "p" then a point type graphic is displayed. 
    Otherwise a lines type graphic displayed. }
    \item{\dots}{Other graphical options arguments.
    \code{\link{plot}}   }
}

\details{
Plot Spatial and Spatio-temporal correlation or  covariance of (non) Gaussian random fields
}

\value{
  Produces a plot.
  No values are returned.
}

\seealso{
  \code{\link{GeoCorrFct}} for examples.
}

\keyword{Correlation}
