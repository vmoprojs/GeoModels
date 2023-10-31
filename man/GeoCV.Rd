\name{GeoCV}  
\alias{GeoCV}
\encoding{UTF-8}
\title{n-fold  kriging Cross-validation}
\description{The procedure use the \code{\link{GeoKrig}} or \code{\link{GeoKrigloc}} function to compute n-fold  kriging cross-validation using informations from a \code{\link{GeoFit}} object.  The function returns some prediction scores.}
\usage{GeoCV(fit, K=100, estimation=FALSE, optimizer="Nelder-Mead",
     lower=NULL, upper=NULL, n.fold=0.05,local=FALSE,
    neighb=NULL, maxdist=NULL,maxtime=NULL,sparse=FALSE,
    type_krig="Simple", which=1,seed=1)}

\arguments{
  \item{fit}{An object of class
    \code{\link{GeoFit}}.}
     \item{K}{The number of iterations in cross-validation.}
       \item{estimation}{Logical; if \code{TRUE} then  an estimation is performed at each iteration and the estimates are used in the prediction. 
       Otherwise the estimates in the object fit are used.}
       \item{optimizer}{The type of optimization algorithm if estimation is  \code{TRUE}. See  \code{\link{GeoFit}} for details.}
        \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize} if estimation is  \code{TRUE}.}
       \item{upper}{An optional named list giving the values  for the upper bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize} if estimation is  \code{TRUE}.}
       \item{n.fold}{Numeric; the percentage of data to be deleted (and predicted) in the cross-validation procedure.}
       \item{local}{Logical; If local is TRUE, then local kriging is performed. The default is FALSE.}
    \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhood if local kriging is performed.}
         \item{maxdist}{Numeric; an optional positive value indicating the distance in the spatial neighborhood
         if local kriging is performed.}
       \item{maxtime}{Numeric; an optional positive value indicating the distance in the temporal neighborhood
       if local kriging is performed.}
       \item{sparse}{Logical; if \code{TRUE} kriging and simulation are computed with sparse matrices algorithms 
          using spam package. Default is FALSE. It should be used with compactly supported covariances.} 
         \item{type_krig}{String; the type of kriging. If \code{Simple} (the
    default) then simple kriging is performed. If \code{Optim} then optimal kriging is performed for some non-Gaussian RFs}
      \item{which}{Numeric; In the case of bivariate  cokriging it indicates which variable to predict.
           It can be 1 or 2}
     \item{seed}{Numeric; The seed used in the  n-fold  kriging cross-validation. Default is 1. Comparison between
     different models in terms  of n-fold  kriging cross-validation must be performed using the same seed}
}

\value{
  Returns an object  containing the following informations:
  \item{predicted}{A list  of  the predicted values   in the CV procedure;}
  \item{data_to_pred}{A list  of  the data to predict  in the CV procedure;}
   \item{mae}{The vector of mean  absolute error in the CV procedure;}
  \item{rmse}{The vector of root mean  squared error in the CV procedure;}
    \item{lscore}{The vector of log-score in the CV procedure;}
      \item{crps}{The vector of continuous ranked probability score  in the CV procedure;}
}


\seealso{\code{\link{GeoKrig}}.}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua89@gmail.com},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/},
Christian", Caamaño-Carrillo, \email{chcaaman@ubiobio.cl},\url{https://www.researchgate.net/profile/Christian-Caamano}
}




\keyword{Composite}
