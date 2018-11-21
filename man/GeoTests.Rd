\name{GeoTests}
\alias{GeoTests}
\encoding{UTF-8}
\title{Statistical Hypothesis Tests for Nested Models}
\description{
  The function performs statistical hypothesis tests for nested models
  based on composite or standard likelihood versions of: Wald-type and
  Wilks-type (likelihood ratio) statistics.
}

\usage{
GeoTests(object1, object2, \dots, statistic)
}

\arguments{
  \item{object1}{An object of class \code{GeoFit}.}
  \item{object2}{An object of class \code{GeoFit} that is a nested
    model within \code{object1}.}
  \item{\dots}{Further successively nested objects.}
  \item{statistic}{String; the name of the statistic used within the
  hypothesis test (see \bold{Details}).}
}


\details{
  The implemented hypothesis tests for nested models are based on the
  following statistics:
  \enumerate{
     \item Wald-type (\code{Wald});
     \item Likelihood ratio or Wilks-type (\code{Wilks} under standard likelihood); 
     For composite likelihood available variants of the basic version are:
     \itemize{
       \item Rotnitzky and Jewell adjustment (\code{WilksRJ});
       \item Satterhwaite adjustment (\code{WilksS});
       \item Chandler and Bate adjustment (\code{WilksCB});
       \item Pace, Salvan and Sartori adjustment (\code{WilksPSS});
     }
   }

 More specifically,
 consider an \eqn{p}-dimensional random vector \eqn{\mathbf{Y}}{Y} with
 probability density function \eqn{f(\mathbf{y};\mathbf{\theta})}{f(y;theta)},
where \eqn{\mathbf{\theta} \in \Theta}{theta in Theta} is a
\eqn{q}-dimensional vector of parameters. Suppose that
\eqn{\mathbf{\theta}=(\mathbf{\psi},\mathbf{\tau})}{theta=(psi, tau)}
 can be partitioned in a \eqn{q'}-dimensional subvector \eqn{\psi}{psi}
 and \eqn{q''}-dimensional subvector \eqn{\tau}{tau}. Assume also to be
 interested in testing the specific values of the vector
 \eqn{\psi}{psi}. Then, one can use some statistical
 hypothesis tests for testing the null hypothesis \eqn{H_0:
 \psi=\psi_0}{H_0: psi=psi_0} against the alternative \eqn{H_1:
 \psi \neq \psi_0}{H_1: psi <> psi_0}. Composite likelihood versions
of 'Wald' and 'score' statistics have the usual asymptotic
chi-square distribution with \eqn{q'} degree of freedom. The Wald-type statistic is
\deqn{W=(\hat{\psi}-\psi_0)^T (G^{\psi \psi})^{-1}(\hat{\theta})(\hat{\psi}-\psi_0),}{%
W=(hat{psi}-psi_0)^T (G^{psi psi})^{-1} (hat{theta})(hat{psi}-psi_0),}
where \eqn{G_{\psi \psi}}{G_{psi psi}} is the \eqn{q' \times q'}{q' x q'}
 submatrix of the Godambe or Fisher information pertaining to \eqn{\psi}{psi} and
 \eqn{\hat{\theta}}{hat{theta}} is the maximum likelihood estimator from
 the full model. 
 This statistic can be called from the
 routine \code{GeoTests} assigning at the argument \code{statistic}
 the value: \code{Wald}.

 Alternatively to the Wald-type  statistic one can use the
 composite version of the Wilks-type or likelihood ratio statistic, given by
 \deqn{W=2[C \ell(\hat{\mathbf{\theta}};\mathbf{y}) - C \ell\{\mathbf{\psi}_0,
   \hat{\mathbf{\tau}}(\mathbf{\psi}_0);\mathbf{y}\}].}{%
 W=2[Cl(hat{theta};y) - Cl{psi_0, hat{tau}(psi_0);y}].}
 In the composite likelihood case, the asymptotic distribution of the composite 
 likelihood ratio statistic is given by
 \deqn{W \dot{\sim} \sum_{i} \lambda_i \chi^2,}{%
   W ~ sum_i lambda_i Chi^2_i,}
 for \eqn{i=1,\ldots,q'}{i=1,...,q'}, where \eqn{\chi^2_i}{Chi^2_i} are
 \eqn{q'} iid copies of a chi-square one random variable and
 \eqn{\lambda_1,\ldots,\lambda_{q'}}{lambda_1,...,lambda_{q'}}
 are the eigenvalues of the matrix \eqn{(H^{\psi \psi})^{-1} G^{\psi
 \psi}}{(H^{psi psi})^-1 G^{psi psi}}. There exist several adjustments
 to the composite likelihood ratio statistic in order to get an
 approximated \eqn{\chi^2_{q'}}{Chi^2_{q'}}. For example, Rotnitzky and Jewell
 (1990) proposed the adjustment \eqn{W'= W / \bar{\lambda}}{W'= W /
   bar{lambda}} where \eqn{\bar{\lambda}}{bar{lambda}} is the average
 of the eigenvalues \eqn{\lambda_i}{lambda_i}. This statistic can be
 called within the routine by the value: \code{WilksRJ}. A better
 solution is proposed by Satterhwaite (1946) defining \eqn{W''=\nu W /
 (q' \bar{\lambda})}{W''= nu W / (q' bar{lambda})}, where \eqn{\nu=(\sum_i
 \lambda)^2 / \sum_i \lambda^2_i}{nu = sum_i lambda / sum_i lambda^2_i} for
\eqn{i=1\ldots,q'}{i=1...,q'}, is the effective number of the degree of
freedom. Note that in this case the distribution of the likelihood ratio
 statistic is a chi-square random variable with \eqn{\nu}{nu} degree of
 freedom. This statistic can be called from the routine assigning the
 value: \code{WilksS}. For the adjustments suggested by Chandler and
 Bate (2007) and Pace, Salvan and Sartori (2011) we refere to the articles (see
 \bold{References}), these versions can be called from the routine assigning
 respectively the values: \code{WilksCB} and \code{WilksPSS}.

}

\value{
 An object of class \code{c("data.frame")}. The object contain a table
 with the results of the tested models. The rows represent the
 responses for each model and the columns the following results:
 \item{Num.Par}{The number of the model's parameters.}
 \item{Diff.Par}{The difference between the number of parameters of the model in the previous row
   and those in the actual row.}
 \item{Df}{The effective number of degree of freedom of the chi-square distribution.}
 \item{Chisq}{The observed value of the statistic.}
 \item{Pr(>chisq)}{The p-value of the quantile
 \code{Chisq} computed using a chi-squared distribution with \code{Df} degrees
    of freedom.}
}

\references{

Chandler, R. E., and Bate, S. (2007). Inference for Clustered Data Using the
Independence log-likelihood. \emph{Biometrika}, \bold{94}, 167--183.

Pace, L., Salvan, A. and Sartori, N. (2011). Adjusting Composite
Likelihood Ratio Statistics. \emph{Statistica Sinica}, \bold{21}, 129--148.

Rotnitzky, A. and Jewell, N. P. (1990). Hypothesis Testing of Regression
Parameters in Semiparametric Generalized Linear Models for Cluster
Correlated Data. \emph{Biometrika}, \bold{77}, 485--497.

Satterthwaite, F. E. (1946). An Approximate Distribution of Estimates of
Variance Components. \emph{Biometrics Bulletin}, \bold{2}, 110--114.

Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

}

\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}

\seealso{
  \command{\link{GeoFit}}.
}
\examples{


library(GeoModels)
################################################################
###
### Example 1. Parametric test of Gaussianity
### using SinhAsinh random fields using standard likelihood
###
###############################################################
set.seed(314)
model="SinhAsinh" 
# Define the spatial-coordinates of the points:
NN=500
x <- runif(NN, 0, 1)
y <- runif(NN, 0, 1)
coords <- cbind(x,y)
# Parameters
mean=0; nugget=0; sill=1
### skew and tail parameters
skew=0;tail=1   ## H0 i.e Gaussianity
# underlying model correlation 
corrmodel="Wend0"
power2=4;c_supp=0.2

# simulation
param=list(power2=power2,skew=skew,tail=tail,
             mean=mean,sill=sill,scale=c_supp,nugget=nugget)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model, param=param)$data


##### H1
fixed=list(power2=power2,nugget=nugget,mean=mean)
start=list(scale=c_supp,skew=skew,tail=tail,sill=sill)
# Maximum pairwise composite-likelihood fitting of the RF:
fitH1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                   likelihood="Full",type="Standard",varest=TRUE,
                    start=start,fixed=fixed)

fitH1$param

##### H0: Gaussianity (i.e tail1=1, skew=0 fixed)
fixed=list(power2=power2,nugget=nugget,mean=mean,tail=1,skew=0)
start=list(scale=c_supp,sill=sill)
# Maximum pairwise composite-likelihood fitting of the RF:
fitH0 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    likelihood="Full",type="Standard",varest=TRUE,
                    start=start,fixed=fixed)


#  Wald statistic test
GeoTests(fitH1, fitH0 ,statistic='Wald')
# likelihood ratio statistic test
GeoTests(fitH1, fitH0 , statistic='Wilks')

################################################################
###
### Example 2. Composite likelihood-based hypothesis testing.
### Testing significance of a covariate parameter 
###
###############################################################
set.seed(31)
# Define the spatial-coordinates of the points:
N=1500
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
X=cbind(rep(1,N),runif(N))
coords=cbind(x,y)

# Set the model's parameters:
corrmodel <- "Exp"
mean <- 1;  mean1=0.3 ## H0: mean1=0
sill <- 1
nugget <- 0
scale <- 0.15/3

# Simulation of the spatial Gaussian random field:
data <- GeoSim(coordx=coords, corrmodel=corrmodel, param=list(mean=mean,mean1=mean1,
              sill=sill,nugget=nugget,scale=scale),X=X)$data

# Pairwise-likelihood fitting of the random field, full model:
start=list(mean=mean,mean1=mean1,sill=sill,scale=scale)
fixed=list(nugget=nugget)
fitH1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, maxdist=0.01,
                    varest=TRUE,likelihood="Marginal",type="Pairwise",
                    fixed=fixed,start=start,X=X)

# Pairwise-likelihood fitting of the random field, with a nested model:
start=list(mean=mean,sill=sill,scale=scale)
fixed=list(nugget=nugget,mean1=0) # setting H0 restriction
fitH0 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, maxdist=0.01,
                    varest=TRUE,likelihood="Marginal",type="Pairwise",
                    fixed=fixed,start=start,X=X)


# Hypothesis testing results:
# composite  Wald statistic test:
GeoTests(fitH1, fitH0 ,statistic='Wald')
# composite likelihood ratio statistic test with PSS adjustment:
GeoTests(fitH1, fitH0, statistic='WilksPSS')


}
\keyword{spatial}
