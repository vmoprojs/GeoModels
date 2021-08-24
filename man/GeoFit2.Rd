\name{GeoFit2}
\alias{GeoFit2}
\encoding{UTF-8}
\title{Max-Likelihood-Based Fitting of Gaussian and non Gaussian RFs.}
\description{
  Maximum weighted composite-likelihood fitting for Gaussian  and some Non-Gaussian  
  univariate spatial, spatio-temporal and bivariate spatial RFs.
  A first preliminary estimation is performed using independence composite-likelihood for the 
  marginal parameters of the model. The estimates are then used as starting values in the second 
  final estimation step.
  The function allows to fix any of the parameters and setting upper/lower bound in the optimization.
}
\usage{
GeoFit2(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,copula=NULL,corrmodel,distance="Eucl",fixed=NULL,GPU=NULL, grid=FALSE, likelihood='Marginal',
     local=c(1,1),  lower=NULL,maxdist=Inf,neighb=NULL,
      maxtime=Inf, memdist=TRUE,method="cholesky", model='Gaussian',n=1, onlyvar=FALSE ,
      optimizer='Nelder-Mead', parallel=FALSE, 
      radius=6371,  sensitivity=FALSE,sparse=FALSE, start=NULL, taper=NULL, tapsep=NULL, 
      type='Pairwise', upper=NULL, varest=FALSE, vartype='SubSamp', weighted=FALSE, winconst=NULL, winstp=NULL, 
      winconst_t=NULL, winstp_t=NULL,X=NULL)
}

\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation)  or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)   or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid).
   For the description see the Section \bold{Details}.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) assigning 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector assigning
    1-dimension of spatial coordinates.  Coordinates on a sphere for a  fixed radius \code{radius} 
    are passed in lon/lat format expressed in decimal degrees.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates.  Optional argument, the default is \code{NULL}
    then a spatial RF is expected.}
  \item{coordx_dyn}{A list of \eqn{m} numeric (\eqn{d_t \times 2}{d x 2})-matrices
       containing dynamical (in time) spatial coordinates. Optional argument, the default is \code{NULL}}
  \item{copula}{String; the type of copula. It can be "Clayton" or "Gaussian"}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl}, the euclidean distance.
    See the Section \bold{Details}.}
  \item{fixed}{An optional named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated.}
    \item{GPU}{Numeric; if \code{NULL} (the default) 
      no OpenCL computation is performed. The user can choose the device to be used. Use \code{DeviceInfo()} function to see available devices, only double precision devices are allowed} 
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
     of non-equispaced spatial sites (irregular grid).}
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default, see the Section \bold{Details}.}
    \item{local}{Numeric; number of local work-items of the OpenCL setup}
   \item{lower}{An optional named list giving the values  for the lower bound of the space parameter
    when the optimizer is  \code{L-BFGS-B} or \code{nlminb}  or \code{optimize}. The 
    names of the list must be the same of the names in the   \code{start} list.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite or tapered likelihood
    computation. See the Section \bold{Details} for more information.}
      \item{neighb}{Numeric; an optional positive integer indicating the 
    order of neighborhood  in the composite likelihood
    computation. See the Section \bold{Details} for more information.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal separation considered in the composite or tapered likelihood
    computation (see \bold{Details}).}
        \item{memdist}{ Logical; if \code{TRUE} then  all  the distances  useful in the  composite likelihood 
      estimation are computed before  the optimization. \code{FALSE} is deprecated. }
   \item{method}{String; the type of matrix decomposition used in the simulation. Default is cholesky.
    The other possible choices is \code{svd}.}   
  \item{model}{String; the type of RF and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
    \item{n}{Numeric; number of trials in a binomial RF; number of successes in a negative binomial RF }  
  \item{onlyvar}{Logical; if \code{TRUE}  (and varest is  \code{TRUE}) only the variance covariance matrix is computed
  without optimizing. \code{FALSE} is the default.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). \code{Nelder-Mead} is the
    default. Other possible choices are \code{ucminf},\code{nlm}, \code{BFGS}  \code{L-BFGS-B} and \code{nlminb}. 
    In these last two cases upper and lower bounds can be passed by the user. In the case
    of one-dimensional optimization, the function \code{optimize} is used.  Two option  for global searching optimization using mcGlobaloptim package are possible with \code{multinlminb} and  \code{multiNelder-Mead} options.}  
  \item{parallel}{Logical; if \code{TRUE} optmization is performed using optimParallel using the maximum number of cores,  when optimizer is  \code{L-BFGS-B}.\code{FALSE} is the default.}  
   \item{radius}{Numeric; the radius of the sphere in the case of lon-lat coordinates.
   The default is 6371, the radius of the earth.}
   \item{sensitivity}{Logical; if \code{TRUE} then  the sensitivy  matrix is computed}
  \item{sparse}{Logical; if \code{TRUE} then  maximum likelihood is computed
  using sparse matrices algorithms (spam packake).It should be used
  with compactly supported covariance models.\code{FALSE} is the default.}
  \item{start}{An optional named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \bold{Details}).}
  \item{taper}{String; the name of the type of covariance matrix.
     It can be \code{Standard} (the default value) or
    \code{Tapering} for taperd covariance matrix.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the 
     space time adaptive taper (see \bold{Details}).}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods (see \bold{Details}).}
  \item{upper}{An optional named list giving the values  for the upper bound 
      of the space parameter when the optimizer is or \code{L-BFGS-B} or  \code{nlminb} or \code{optimize}. 
       The names of the list must be the same of the names in the   \code{start} list.}  
  \item{varest}{Logical; if \code{TRUE} the estimates' variances and
    standard errors are returned. For composite likelihood estimation it is deprecated.
    Use sensitivity TRUE and update the object using the function \code{GeoVarestbootstrap}
    \code{FALSE} is the default.}
  \item{vartype}{String; (\code{SubSamp} the default) the type of
    method used for computing
  the estimates' variances, see the Section \bold{Details}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see the Section \bold{Details}. If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
  \item{winconst}{Numeric; a bivariate positive vector for computing the spatial  sub-window
      in the sub-sampling procedure. See \bold{Details} for more information.}
 \item{winstp}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the spatial  sub-sampling procedure. The case \eqn{1} 
    correspond to no overlapping. See \bold{Details} for more information.}
   \item{winconst_t}{Numeric; a positive value for computing the temporal sub-window
         in the sub-sampling procedure. See \bold{Details} for more information.}
 \item{winstp_t}{Numeric; a value in \eqn{(0,1]} for defining the the proportion of 
      overlapping in the temporal  sub-sampling procedure. The case \eqn{1} 
    correspond to no overlapping. See \bold{Details} for more information.}
  \item{X}{Numeric; Matrix of spatio(temporal)covariates in the linear mean specification.}
}

\details{
  \code{GeoFit} provides standard maximum likelihood fitting for Gaussian models. 
  For Gaussian and non Gaussian models, only composite likelihood estimation based on pairs. Specifically marginal and conditional  pairwise 
  likelihood is considered for each type of random field (Gaussian and not Gaussian).
  The optimization method is specified using \code{optimizer}. The default method is \code{Nelder-mead} and other available methods are \code{ucminf}, \code{nlm}, \code{BFGS}. \code{L-BFGS-B} and  \code{nlminb}. In the last two cases upper and lower bounds constraints in the optimization can be specified using \code{lower} and \code{upper} parameters.

  Depending on the dimension of \code{data} and on  the name of the correlation model, 
  the observations   are assumed as a realization of
  a spatial, spatio-temporal or bivariate RF.
  Specifically, with \code{data}, \code{coordx}, \code{coordy}, \code{coordt} parameters:
  \itemize{
    \item If \code{data} is a numeric \eqn{d}{d}-dimensional vector, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}), then the data are interpreted as a single spatial
    realisation observed on \eqn{d}{d} spatial sites;
    \item If \code{data} is a numeric (\eqn{t \times d}{t x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}),
    \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data
    are interpreted as a single spatial-temporal realisation of a RF observed on \code{d} spatial sites and for \code{t} times.
      \item If \code{data} is a numeric (\eqn{2 \times d}{2 x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}), then the data
    are interpreted as a single spatial realisation of a bivariate RF observed on \code{d} spatial sites. 
    \item If \code{data} is a list,  \code{coordxdyn}  is a list and   \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, then the data
    are interpreted as a single spatial-temporal realisation of a RF observed on dynamical spatial sites (different locations sites for each temporal instants) and for \code{t} times.
 }

Is is also possible to specify a matrix of covariates using \code{X}.
Specifically:
  \itemize{
    \item  In the spatial case  \code{X} must be a (\eqn{d \times k}{d x k}) covariates matrix associated to  \code{data}  a numeric \eqn{d}{d}-dimensional vector;
     \item  In the spatiotemporal case  \code{X} must be a (\eqn{N \times k}{N x k}) covariates matrix associated to  \code{data}  a numeric (\eqn{t \times d}{t x d})-matrix, where \eqn{N=t\times d};
     \item  In the spatiotemporal case  \code{X} must be a (\eqn{N \times k}{N x k}) covariates matrix associated to  \code{data}  a numeric (\eqn{t \times d}{t x d})-matrix, where \eqn{N=2\times d};
 }

  The \code{corrmodel} parameter allows to select a specific correlation
  function for the RF. (See   \code{\link{GeoCovmatrix}} ).

  The \code{distance} parameter allows to consider differents kinds of spatial distances.
  The settings alternatives are:
  \enumerate{
    \item \code{Eucl}, the euclidean distance (default value);
    \item \code{Chor}, the chordal distance;
    \item \code{Geod}, the geodesic distance;
  }

  The \code{likelihood} parameter represents the composite-likelihood
  configurations. The settings alternatives are:
  \enumerate{
    \item \code{Conditional}, the composite-likelihood is formed by
    conditionals likelihoods;
    \item \code{Marginal}, the composite-likelihood is formed by
  marginals likelihoods (default value);
    \item \code{Full}, the composite-likelihood turns out to be the standard likelihood;
  }

  The \code{model} parameter indicates the type of RF
  considered. The available options are:

  RF with marginal symmetric distribution:
  \itemize{
    \item \code{Gaussian}, for a Gaussian RF.   
    \item \code{StudentT}, for a StudentT RF (see  Bevilacqua M., Caamaño C., Arellano Valle R.B., Morales-Oñate V., 2020).
    \item \code{Tukeyh}, for a Tukeyh RF.
    \item \code{Logistic}, for a Logistic RF.
 }   
   RF with positive values and  right skewed marginal distribution:
  \itemize{
    \item \code{Gamma} for a Gamma  RF (see  Bevilacqua M., Caamano C.,  Gaetan, 2020)
    \item \code{Weibull} for a Weibull  RF (see  Bevilacqua M., Caamano C.,  Gaetan, 2020)
    \item \code{LogGaussian} for a  LogGaussian  RF (see  Bevilacqua M., Caamano C.,  Gaetan, 2020)
    \item \code{LogLogistics} for a  LogLogistic RF. 
}
RF with with possibly asymmetric marginal distribution:
  \itemize{
    \item \code{SkewGaussian} for a skew Gaussian  RF (see   Alegrıa et al. (2017))
    \item \code{SinhAsinh} for a Sinh-arcsinh  RF.
 }
 RF with for directional data
 \itemize{   
    \item \code{Wrapped} for a wrapped Gaussian  RF (see  Alegria A., Bevilacqua, M., Porcu, E. (2016))
    }
RF with marginal counts data
 \itemize{  
     \item \code{Poisson} for a Poisson  RF.
     \item \code{PoissonZIP} for a zero inflated Poisson  RF.
     \item \code{Binomial} for a Binomial  RF.
     \item \code{BinomialNeg} for a negative Binomial RF.
     \item \code{BinomialNegZINB} for a zero inflated negative Binomial RF.
  }
RF using Gaussian and Clayton copula with the following marginal distribution:
 \itemize{   
    \item \code{Gaussian} 
    \item \code{Beta} for a Beta   RF, \code{Beta2} using mean parametrization.
    \item \code{Kumaraswamy} for a  Kumaraswamy  RF,\code{Kumaraswamy2} using median parametrization.
    }

  For  a given \code{model} the associated parameters are given by nuisance and covariance parameters. In order to obtain the nuisance parameter associated to a specific model use \code{\link{NuisParam}}. 
  In order to obtain the covariance parameter associated to a given covariance model use 
  \code{\link{CorrParam}}.

   All the nuisance and covariance parameters must be specified
   by the user using  the  \code{start}  and  the \code{fixed} parameter.
   Specifically:

   The \code{start} parameter allows to specify (as starting values for the optimization)  the parameters  to be estimated.
   The \code{fixed} parameter allows to fix  some of the parameters.

   Regression parameters in the linear specfication must be specified as \code{mean,mean1,..meank} (see \code{\link{NuisParam}}).
   In this case a matrix of covariates with suitable dimension can be specified using the parameter  \code{X}. 
   In the case of a single mean then \code{X} should not be specified and it is interpreted as a vector of ones.

   The \code{taper} parameter, optional in case that
   \code{type=Tapering}, indicates the type of taper
   correlation model. (See \code{\link{GeoCovmatrix}})


  If a Gaussian or (any) non Gaussian RF is considered then the possible combination is
  marginal pairwise likelihoods (\code{likelihood=Marginal}) and \code{type="Pairwise"})
  or (\code{likelihood=Conditional}) and \code{type="Pairwise"})


  If a Gaussian RF is considered  (\code{model=Gaussian}) then:
    \itemize{
      \item If the composite is formed by marginal likelihoods (\code{likelihood=Marginal}):
      \itemize{
	\item \code{Pairwise}, the composite-likelihood is defined by
	the pairwise likelihoods;
	\item \code{Difference}, the composite-likelihood is defined by
	likelihoods
	which are obtained as difference of the pairwise likelihoods.}
      \item If the composite is formed by conditional likelihoods
      (\code{likelihood=Conditional})
      \itemize{
	\item \code{Pairwise}, the composite-likelihood is defined by
	the pairwise conditional likelihoods.}
      \item If the composite is formed by a full likelihood
      (\code{likelihood=Full}):
      \itemize{
	\item \code{Standard}, the objective function is the classical
        multivariate likelihood;
	\item \code{Restricted}, the objective function is the
        restricted version of the full likelihood (e.g. Harville 1977, see \bold{References});
	\item \code{Tapering}, the objective function is the tapered 2 (unbiased estimating equation)
        version of the full likelihood (e.g. Kaufman et al. 2008, see \bold{References});
    \item \code{Tapering1}, the objective function is the tapered 1 (biased estimating equation)
        version of the full likelihood (e.g. Kaufman et al. 2008, see \bold{References}); 
   \item \code{CV}, the objective function  is the cross validation estimation method ;           
      }
      }



  Two type  of binary weights can be used in  the composite likelihood estimation,
  one based on neighboords and one based on distances.
   
  In the first  case  binary weights are set to 1 and 0 otherwise
  depending if the pairs are neighboords of a certain order (1, 2, 3, ..) specified by the parameter (\code{neighb}). This weighting scheme is effecient for large-data sets since the computation of the 'useful' pairs  in based on  the package RANN that  provides
  fast nearest neighbour search.

  In the second case, the \code{maxdist} parameter set the maximum
  spatial distance below which pairs of sites with inferior distances
  Pairs with distance less than \code{maxdist} have weight 1
  and are included in the likelihood computation, instead those with
  greater distance have weight 0 and then excluded.
  This weighting scheme is less ineffecient for large data.
  The same arguments of \code{maxdist} are valid for \code{maxtime} but
  here the weigthed composite-likelihood regards the case of
  spatial-temporal field.

  The \code{varest=TRUE} parameter specifies if the standard error estimation of the estimated parameters must be computed.
  For Gaussian RFs and standard  likelihood estimation, standard  errors  are computed as square root of the diagonal elements of the Fisher
  Information matrix (asymptotic covariance matrix of the estimates under increasing domain).
  For Gaussian and non Gaussian  RFs and composite  likelihood estimation, standard errors estimate are computed as square root of the diagonal elements of the
  Godambe Information matrix. (asymptotic covariance matrix of the estimates under increasing domain  (see  Bevilacqua et. al. (2012) , Bevilacqua and Gaetan (2013)).

  For standard error estimation  of composite likelihood estimatin
  use the option \code{sensitivity=TRUE} and then update the object with the function \code{GeoVarestbootstrap}.
  The option  \code{varest=TRUE}  is deprecated and the comments below should not be considered.
  The \code{varest=TRUE} option allows std error estimation trough a sub-sampling procedure.
  In the  the sub-sampling procedure,\code{winconst}
   and \code{winstp} parameters represent respectively a positive constant used to
  determine the sub-window size and the the step with which the
  sub-window moves. 
  In the spatial case (subset of \eqn{R^2}{R^2}), the domain is seen as
  a rectangle \eqn{B \times H}{BxH}, therefore the size of the
  sub-window side \eqn{b}{b} is given by \eqn{b=winconst \times
  \sqrt(B)}{b=winconst x sqrt(B)} (similar is of \eqn{H}{H}).
  For a complete description see Lee and Lahiri (2002).
  By default \code{winconst} is set \eqn{B / (4 \times \sqrt(B))}{B / (4
  x sqrt(B))}.  The \code{winstp} parameter is used to determine the sub-window step. The latter is given by the
  proportion of the sub-window size, so that when \code{winstp=1} there
  is not overlapping between contiguous sub-windows.
  In the spatial case by default \code{winstp=0.5}.
  The sub-window is moved
  by successive steps in order to cover the entire spatial
  domain.
  Observations, that fall in disjoint or overlapping windows
  are considered indipendent samples.

  In the spatio-temporal case  \code{winconst_t} represents
  the lenght of the temporal sub-window. By default the size of the
  sub-window is computed following the rule established in Li et al. (2007).
  By default \code{winstp} is the time step.
}

\value{
  Returns an object of class \code{GeoFit}.
  An object of class \code{GeoFit} is a list containing
  at most the following components:
  \item{bivariate}{Logical:\code{TRUE} if the Gaussian RF is bivariate, otherwise  \code{FALSE};}
  \item{clic}{The composite information criterion, if the full
  likelihood is considered then it coincides with the Akaike information
  criterion;}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{coordx_dyn}{A  list of dynamical (in time) spatial coordinates;}
  \item{convergence}{A string that denotes if convergence is reached;}
  \item{copula}{The type of copula;}
  \item{corrmodel}{The correlation model;}
  \item{data}{The vector or matrix or array of data;}
  \item{distance}{The type of spatial distance;}
  \item{fixed}{The vector of fixed parameters;}
  \item{iterations}{The number of iteration used by the numerical routine;}
  \item{likelihood}{The configuration of the composite likelihood;}
  \item{logCompLik}{The value of the log composite-likelihood at the maximum;}
  \item{maxdist}{The maximum spatial distance used in the weigthed composite likelihood. 
                       If no spatial distance is specified then it is NULL;}
  \item{maxtime}{The maximum temporal distance  used in the weigthed composite likelihood.
                          If no spatial distance is specified then it is NULL;}
  \item{message}{Extra message passed from the numerical routines;}
  \item{model}{The density associated to the likelihood objects;}
   \item{missp}{True if a misspecified  Gaussian model is ued in  the composite likelihhod;}
  \item{n}{The number of trials in a binominal RF;the number of successes in a  negative Binomial RFs;} 
  \item{ns}{The number of (different) location sites in the bivariate case;} 
   \item{nozero}{In the case of tapered likelihood  the percentage of non zero values in the covariance matrix.
  Otherwise is NULL.}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number of the temporal realisations of the RF;}
  \item{param}{The vector of  parameters' estimates;}
  \item{radius}{The radius of the sphere in the case of great  circle distance;}
   \item{stderr}{The vector of standard errors;}
  \item{sensmat}{The sensitivity matrix;}
  \item{varcov}{The matrix of the variance-covariance of the estimates;}
  \item{varimat}{The variability matrix;}
  \item{vartype}{The method used to compute the variance of the estimates;}
  \item{type}{The type of the likelihood objects.}
  \item{winconst}{The constant used to compute the window size in the
  spatial sub-sampling;}
  \item{winstp}{The step used for moving the window in the spatial sub-sampling;}
    \item{winconst_t}{The constant used to compute the window size in the
  spatio-temporal sub-sampling;}
  \item{winstp_}{The step used for moving the window in the spatio-temporal sub-sampling;}
  \item{X}{The matrix of covariates;}
}


\references{

  Maximum Restricted Likelihood Estimator:

  Harville, D. A. (1977)
  Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems.
  \emph{Journal of the American Statistical Association}, \bold{72},
  320--338.

  Tapered likelihood:

    Kaufman, C. G., Schervish, M. J. and Nychka, D. W. (2008)
  Covariance Tapering for Likelihood-Based Estimation in Large Spatial Dataset.
  \emph{Journal of the American Statistical Association}, \bold{103},
  1545--1555.

  Composite-likelihood:

  Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
  Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

  Varin, C. and Vidoni, P. (2005)
  A Note on Composite Likelihood Inference and Model Selection.
  \emph{Biometrika}, \bold{92}, 519--528.


  Weighted Composite-likelihood for non Gaussian RF: 

  Alegrıa A., Caro S., Bevilacqua M., Porcu E., Clarke J. (2017) \emph{Estimating covariance functions of multivariate skew-Gaussian random fields on the sphere}. Spatial Statistics \bold{22} 388-402

  Alegria A., Bevilacqua, M., Porcu, E. (2016) 
  Likelihood-based inference for multivariate space-time wrapped-Gaussian fields. 
    \emph{Journal of Statistical Computation and Simulation}. \bold{86(13)}, 2583--2597.

  Bevilacqua M., Caamano C.,  Gaetan C. (2020) On modeling positive continuous data with spatio-temporal dependence.  \emph{Environmetrics}   \bold{31(7)}

  Bevilacqua M., Caamaño C., Arellano Valle R.B., Morales-Oñate V.  (2020) Non-Gaussian Geostatistical Modeling using (skew) t Processes. \emph{Scandinavian Journal of Statistics}.

  Weighted Composite-likelihood for Gaussian RFs:

    Bevilacqua, M. Gaetan, C., Mateu, J. and Porcu, E. (2012)
  Estimating space and space-time covariance functions for large data
  sets: a weighted composite likelihood approach.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{107}, 268--280.

  Bevilacqua, M., Gaetan, C. (2015) Comparing composite likelihood methods based on pairs for spatial Gaussian random fields.  \emph{Statistics and Computing}, \bold{25(5)},  877-892.


  Sub-sampling estimation:

  Carlstein, E. (1986) The Use of Subseries Values for Estimating the
  Variance.
  \emph{The Annals of Statistics},
  \bold{14}, 1171--1179.

  Heagerty, P. J. and Lumley T. (2000) Window Subsampling of Estimating
  Functions with Application to Regression Models.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{95}, 197--211.

  Lee, Y. D. and Lahiri S. N. (2002) Variogram Fitting by Spatial
  Subsampling.
  \emph{Journal of the Royal Statistical Society. Series B},
  \bold{64}, 837--854.

}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/view/moreno-bevilacqua/home},
Víctor Morales Oñate, \email{victor.morales@uv.cl}, \url{https://sites.google.com/site/moralesonatevictor/}
}
  
\examples{
library(GeoModels)
library(fields)

###############################################################
############ Examples of spatial Gaussian RFs ################
###############################################################

################################################################
###
### Example 1, 2 :  Estimation of a spatial Gaussian RF with 
### Matern correlation using  pairwise likelihood and 
### maximum likelihood with  BGGS and nlminb optimization
###
###############################################################

# Define the spatial-coordinates of the points:
set.seed(3)
N=400  # number of location sites
x <- runif(N, 0, 1)
set.seed(6)
y <- runif(N, 0, 1)
coords <- cbind(x,y)

# Define spatial matrix covariates
X=cbind(rep(1,N),runif(N))

# Set the covariance model's parameters:
corrmodel <- "Matern"
mean <- 0.2
mean1 <- -0.5
sill <- 1
nugget <- 0
scale <- 0.2/3
smooth=0.5
param<-list(mean=mean,mean1=mean1,sill=sill,nugget=nugget,scale=scale,smooth=smooth)

# Simulation of the spatial Gaussian RF:
data <- GeoSim(coordx=coords,corrmodel=corrmodel, param=param,X=X)$data

fixed<-list(nugget=nugget,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)



################################################################
###
### Example 1. Maximum pairwise likelihood fitting of
### Gaussian RFs with exponential correlation.
### 
###############################################################
fit1 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel, 
                    optimizer="BFGS",neighb=3,likelihood="Marginal",
                    type="Pairwise", start=start,fixed=fixed,X=X)
print(fit1)

################################################################
###
### Example 2. Standard Maximum likelihood fitting of
### Gaussian RFs with exponential correlation.
###
###############################################################
I=Inf
lower<-list(mean=-I,mean1=-I,scale=0,sill=0)
upper<-list(mean=I,mean1=I,scale=I,sill=I)
fit2 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel,
                    optimizer="nlminb",upper=upper,lower=lower,
                    likelihood="Full",type="Standard", 
                    start=start,fixed=fixed,X=X)
print(fit2)


###############################################################
############ Examples of spatial non-Gaussian RFs #############
###############################################################


################################################################
###
### Example 3. Maximum pairwise likelihood fitting of spatial
### Gamma and Weibull  RFs with Generalized Wendland correlation
### using Nelder-Mead and BFGS 
###############################################################
set.seed(524)
# Define the spatial-coordinates of the points:
N=500
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
X=cbind(rep(1,N),runif(N))
mean=1; mean1=2 # regression parameters
nugget=0
shape=2
scale=0.2
smooth=0

model="Weibull"
corrmodel="GenWend"
param=list(mean=mean,mean1=mean1,sill=1-nugget,scale=scale,
                     shape=shape,nugget=nugget,power2=4,smooth=smooth)
# Simulation of a  non stationary weibull RF:
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model,X=X,
           param=param)$data

fixed<-list(nugget=nugget,power2=4,sill=1-nugget,smooth=smooth)
start<-list(mean=mean,mean1=mean1,scale=scale,shape=shape)

# Maximum pairwise composite-likelihood fitting of the RF:
fit <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    neighb=3,likelihood="Marginal",type="Pairwise",X=X,
                    start=start,fixed=fixed)
print(fit$param)

model="Gamma"
start<-list(mean=mean,mean1=mean1,scale=scale)
fixed<-list(nugget=nugget,power2=4,sill=1-nugget,shape=6)
# Maximum pairwise composite-likelihood fitting of the RF:
fit <- GeoFit2(data=data,coordx=coords,corrmodel="Wend0", model=model,
                    neighb=3,likelihood="Marginal",type="Pairwise",X=X,
                    start=start,fixed=fixed)

print(fit$param)



################################################################
###
### Example 5. Maximum pairwise likelihood fitting of
### SinhAsinh-Gaussian spatial  RFs with Wendland correlation
###
###############################################################
set.seed(261)
model="SinhAsinh"
# Define the spatial-coordinates of the points:
x <- runif(500, 0, 1)
y <- runif(500, 0, 1)
coords <- cbind(x,y)

corrmodel="Wend0"
mean=0;nugget=0
sill=1
skew=-0.5
tail=1.5
power2=4
c_supp=0.2

# model parameters
param=list(power2=power2,skew=skew,tail=tail,
             mean=mean,sill=sill,scale=c_supp,nugget=nugget)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model, param=param)$data

plot(density(data))
fixed=list(power2=power2,nugget=nugget)
start=list(scale=c_supp,skew=skew,tail=tail,mean=mean,sill=sill)
# Maximum pairwise likelihood:
fit1 <- GeoFit2(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    neighb=3,likelihood="Marginal",type="Pairwise",
                    start=start,fixed=fixed)



}
\keyword{Composite}