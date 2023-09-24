\name{plot.GeoVariogram}
\alias{plot.GeoVariogram}
\encoding{UTF-8}
\title{Plot empirical spatial, spatio-temporal and spatial bivariate semi-Variogram}

\description{
  Plot empirical spatial, spatio-temporal and spatial bivariate semi-Variogram using 
   on object  \code{\link{GeoVariogram}}.
}

\usage{
\method{plot}{GeoVariogram}(x, \dots)
}

\arguments{
  \item{x}{an object of the class \code{"GeoVariogram"} }
    \item{\dots}{other arguments to be passed to the function
    \code{\link{plot}}   }
}

\details{
  This function plots empirical semi   variogram in the spatial, spatio-temporal and spatial
  bivariate case
}

\value{
  Produces a plot.
  No values are returned.
}

\seealso{
  \code{\link{GeoVariogram}} for variogram computation and examples.
}

\keyword{Variogram}
