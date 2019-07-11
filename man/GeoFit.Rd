\name{GeoFit}
\alias{GeoFit}
\alias{print.GeoFit}
\encoding{UTF-8}
\title{Max-Likelihood-Based Fitting of Gaussian and non Gaussian RFs.}
\description{
  Maximum weighted composite-likelihood fitting for Gaussian  and some Non-Gaussian  
  univariate spatial, spatio-temporal and bivariate spatial RFs
  The function returns the model parameters' estimates and the estimates'
  variances. Moreover the function allows to fix any of the parameters and setting upper/lower bound in the optimization.
}
\usage{
GeoFit(data, coordx, coordy=NULL, coordt=NULL,coordx_dyn=NULL, corrmodel, distance="Eucl",
            fixed=NULL,GPU=NULL, grid=FALSE, likelihood='Marginal',  local=c(1,1),lower=NULL, 
            maxdist=NULL, maxtime=NULL, method="cholesky", model='Gaussian',
            n=1, onlyvar=FALSE ,optimizer='Nelder-Mead', parallel=FALSE, radius=6371, sensitivity=FALSE,
            sparse=FALSE, start=NULL, taper=NULL, tapsep=NULL,type='Pairwise',
            upper=NULL, varest=FALSE, vartype='SubSamp', weighted=FALSE, winconst=NULL,
            winstp=NULL,winconst_t=NULL, winstp_t=NULL,X=NULL)
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
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal separation considered in the composite or tapered likelihood
    computation (see \bold{Details}).}
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
    of one-dimensional optimization, the function \code{optimize} is used.}  
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
    standard errors are returned.
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
  Note, that the standard likelihood may be seen as particular case of the
  composite likelihood. In this respect \code{GeoFit} provides maximum likelihood
  fitting. Only composite likelihood estimation based on pairs are considered. Specifically marginal pairwise
  likelihood is considered for each type of random field (Gaussian and not Gaussian).
  In the Gaussian case  other types of estimation are available: conditional  and difference pairwise likelihood, covariance tapering and  cross-validation method.


  The optimization method is specified using \code{optimizer}. The default method is \code{Nelder-mead}
  and other available methods are \code{ucminf}, \code{nlm}, \code{BFGS}. \code{L-BFGS-B} and  \code{nlminb}. In the last two cases upper and lower bounds constraints in the optimization can be specified using \code{lower} and \code{upper}.




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

 It is possible consider data on a regular grid  using option \code{grid=TRUE}. Specifically:
  \itemize{
   \item If \code{data} is a numeric (\eqn{d \times d}{d x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors,  \code{grid=TRUE}, then the data are interpreted as
    a single spatial RF realisation observed on \code{d} equispaced spatial sites (named regular grid).

    \item If \code{data} is a numeric (\eqn{d \times d \times t}{d x d x t})-array, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors, \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, \code{grid=TRUE}, then the data are interpreted as a single spatial-temporal realisation of a RF observed on
    \code{d} equispaced spatial sites and for \code{t} times.

    \item If \code{data} is a numeric (\eqn{d \times d \times 2}{d x d x 2})-array, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors, \code{grid=TRUE}, then the data are interpreted as a single 
    realisation of a bivariate  RF observed on
    \code{d} equispaced spatial sites.
}



When  \code{grid=FALSE}  is is possible to specify a matrix of covariates using \code{X}.
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
  marginals likelihoods;
    \item \code{Full}, the composite-likelihood turns out to be the standard likelihood;
  }





  The \code{model} parameter indicates the type of RF
  considered. The available options are:

  RF with marginal symmetric distribution:
  \itemize{
    \item \code{Gaussian}, for a Gaussian RF 
    \item \code{Logistic}, for a Logistic RF  (see  Bevilacqua, M. Caamano C., Gaetan , C. 2018)
    \item \code{StudentT}, for a StudentT RF  (see  Bevilacqua, M. Caamano C., 2018)
 }   
   RF with positive right skewed marginal distribution:
  \itemize{
    \item \code{Gamma} for a Gamma  RF (see  Bevilacqua, M. Caamano C., 2018)
     \item \code{Gamma} for a Weibull  RF (see  Bevilacqua, M. Caamano C., 2018)
    \item \code{LogGaussian} for a  LogGaussian  RF
    \item \code{LogLogistics} for a  LogLogistic RF   (see  Bevilacqua, M. Caamano C., Gaetan , C. 2018)
}
RF with with possibly asymmetric marginal distribution:
  \itemize{
    \item \code{SkewGaussian} for a skew Gaussian  RF (see Hao Zhang, H. and El-Shaarawi A. (2010))
    \item \code{SinhAsinh} for Sinh-arcsinh  RF (see  Bevilacqua, M. Caamano C., 2018)
 }
 RF with for directional data
 \itemize{   
    \item \code{Wrapped} for a wrapped Gaussian  RF (see  Alegria A., Bevilacqua, M., Porcu, E. (2016))
    }
Rf with marginal counts data
 \itemize{  
     \item \code{Binomial} for a Binomial  RF (see Bevilacqua, M. Caamano C., Gaetan C. (2018))
    \item \code{BinomialNeg} for a negative Binomial RF(see Bevilacqua, M. Caamano C., Gaetan C. (2018))
  }

  For  a given \code{model} the associated parameters are given by nuisance and covariance parameters.
  In order to obtain the nuisance parameter use \code{\link{NuisParam}}. 
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

   The \code{type} parameter represents the type of likelihood used in the
  composite-likelihood definition. The possible alternatives are listed
  in the following scheme.

  If a Gaussian or (any) non Gaussian RF is considered then the possible combination is
  marginal pairwise likelihoods (\code{likelihood=Marginal}) and \code{type="Pairwise"}).

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


    The \code{maxdist} parameter set the maximum
  spatial distance below which pairs of sites with inferior distances
  are considered in the composite-likelihood. This can be
  lower of the  maximum spatial distance. \bold{Note} that
  this corresponds to use a weighted composite-likelihood with binary
  weights. Pairs with distance less than \code{maxdist} have weight 1
  and are included in the likelihood computation, instead those with
  greater distance have weight 0 and then excluded.
  The default is \code{NULL}, in this case the effective maximum spatial distance
  between sites is considered.


  The same arguments of \code{maxdist} are valid for \code{maxtime} but
  here the weigthed composite-likelihood regards the case of
  spatial-temporal field. At the moment is
  implemented only for Gaussian RFs. The default
  is \code{NULL}, in this case the effective maximum temporal lag
  between pairs of observations is considered.

  In the case of tapering likelihood \code{maxdist} and \code{maxtime} describes
   the spatial and temporal compact support of the taper (see \code{\link{GeoCovmatrix}}).
  If they are not specified then the maximum spatial and temporal distances are considered.
  In the case of space time adaptive taper the \code{tapsep} parameter allows to specify 
  the spatio temporal compact support (see \code{\link{GeoCovmatrix}}).    
      
  The \code{varest} parameter specifies if the standard error estimation of the estimated parameters must be computed.
  For Gaussian RF and standard (restricted) likelihood estimation, standard  errors  are computed as square root of the diagonal elements of the Fisher
  Information matrix (asymptotic covariance matrix of the estimates under increasing domain).
  For Gaussian RF and tapered and composite  likelihood estimation, standard errors estimate are computed as square root of the diagonal elements of the
  Godambe Information matrix. (asymptotic covariance matrix of the estimates under increasing domain  (see Shaby, B. and D. Ruppert (2012) for tapering and
  Bevilacqua et. al. (2012) , Bevilacqua and Gaetan (2013) for weighted composite likelihood)).
  The \code{vartype} parameter specifies the method used to compute the estimates' variances in the composite likelihood case. In particular
  for estimating the variability matrix \code{J} in the Godambe expression matrix.
  This parameter is considered if \code{varest=TRUE}. The options are:
  \itemize{
    \item \code{SubSamp} (the default), indicates the Sub-Sampling method;
  }

  The \code{weighted} parameter specifies if the likelihoods forming the
  composite-likelihood must be weighted. If \code{TRUE} the weights are
  selected by opportune procedures that improve the efficient of the
  maximum composite-likelihood estimator (not implemented yet). If
  \code{FALSE} the efficient improvement procedure is not used.

  For computing the standard errors by the sub-sampling procedure,
  \code{winconst} and \code{winstp} parameters represent respectively a positive constant used to
  determine the sub-window size and the the step with which the
  sub-window moves.

  In the spatial case (subset of \eqn{R^2}{R^2}), the domain is seen as
  a rectangle \eqn{B \times H}{BxH}, therefore the size of the
  sub-window side \eqn{b}{b} is given by \eqn{b=winconst \times
  \sqrt(B)}{b=winconst x sqrt(B)} (similar is of \eqn{h}{h}).
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
  \item{n}{The number of trials in a binominal RF;the number of successes in a  negative Binomial RFs;} 
   \item{nozero}{In the case of tapered likelihood  the percentage of non zero values in the covariance matrix.
  Otherwise is NULL.}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number of the temporal realisations of the RF;}
  \item{param}{The vector of  parameters' estimates;}
  \item{param}{The radius of the sphere in the case of great  circle distance;}
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

    Shaby, B. and D. Ruppert (2012).
  Tapered covariance: Bayesian estimation and asymptotics. \emph{J. Comp. Graph.
Stat.}, \bold{21-2}, 433--452.


  Composite-likelihood:

  Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
  Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

  Varin, C. and Vidoni, P. (2005)
  A Note on Composite Likelihood Inference and Model Selection.
  \emph{Biometrika}, \bold{92}, 519--528.


  Weighted Composite-likelihood for binary RFs:

   Patrick, J. H. and Subhash, R. L. (1998)
  A Composite Likelihood Approach to Binary Spatial Data.
  \emph{Journal of the American Statistical Association, Theory & Methods}, \bold{93}, 1099--1111.


  Weighted Composite-likelihood for wrapped directional RFs:

   Alegria A., Bevilacqua, M., Porcu, E. (2016) 
   Likelihood-based inference for multivariate space-time wrapped-Gaussian fields. 
    \emph{Journal of Statistical Computation and Simulation}. \bold{86(13)}, 2583--2597.


  Weighted Composite-likelihood for Gaussian RFs:

    Bevilacqua, M. Gaetan, C., Mateu, J. and Porcu, E. (2012)
  Estimating space and space-time covariance functions for large data
  sets: a weighted composite likelihood approach.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{107}, 268--280.

  Bevilacqua, M. Gaetan, C. (2013)
  On composite likelihood inference based on pairs for
spatial Gaussian RFs
  \emph{Techical Report, Department of Statistics, de Valparaiso University}.


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

  Li, B., Genton, M. G. and Sherman, M. (2007). A nonparametric assessment of properties of space-time
  covariance functions. \emph{Journal of the American Statistical
    Association}, \bold{102}, 736--744
}



\author{Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},\url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home},
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
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, 
                    optimizer="BFGS",maxdist=0.05,likelihood="Marginal",
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
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel,
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
fit <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    maxdist=.05,likelihood="Marginal",type="Pairwise",X=X,
                    optimizer="Nelder-Mead",
                    start=start,fixed=fixed)
print(fit$param)

model="Gamma"
start<-list(mean=mean,mean1=mean1,scale=scale)
fixed<-list(nugget=nugget,power2=4,sill=1-nugget,shape=6)
# Maximum pairwise composite-likelihood fitting of the RF:
fit <- GeoFit(data=data,coordx=coords,corrmodel="Wend0", model=model,
                    maxdist=.05,likelihood="Marginal",type="Pairwise",X=X,
                      optimizer="BFGS",
                    start=start,fixed=fixed)

print(fit$param)


################################################################
###
### Example 4. Maximum pairwise likelihood fitting of
### StudendT spatial  RFs with Wendland correlation
###
###############################################################

set.seed(15274)
# Define the spatial-coordinates of the points:
N=300
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
X=cbind(rep(1,N),runif(N))
mean=1; mean1=2 # regression parameters
nugget=0
sill=0.5
scale=0.2
df=4   ## degrees of freedom
model="StudentT"
corrmodel="Wend0"

# Simulation of a studentT RF:
param=list(mean=mean,mean1=mean1,sill=sill,scale=scale,df=1/df,nugget=nugget,power2=4)
data <- GeoSim(coordx=coords, corrmodel=corrmodel,model=model,X=X,
           param=param)$data
##  estimation assuming df unknown
fixed<-list(nugget=nugget,power2=4)
start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill,df=1/df)



# Maximum pairwise composite-likelihood fitting of the RF:

fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
               maxdist=.05,likelihood="Marginal",type="Pairwise",X=X,
               optimizer="BFGS",start=start,fixed=fixed)
print(fit1$param)





##  df must be rounded and fixed
df=round(1/(as.numeric(fit1$param['df'])))
 fixed<-list(nugget=nugget,power2=4,df=1/df)
 start<-list(mean=mean,mean1=mean1,scale=scale,sill=sill)
# Maximum pairwise composite-likelihood fitting of the RF:
fit <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    maxdist=.05,likelihood="Marginal",type="Pairwise",X=X,
                    optimizer="BFGS" , start=start,fixed=fixed)
print(fit$param)


################################################################
###
### Example 5. Maximum pairwise likelihood fitting of
### SinhAsinh-Gaussian spatial  RFs with Wendland correlation
### using  BFGS
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
fit1 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    maxdist=0.05,likelihood="Marginal",type="Pairwise",
                   optimizer="BFGS" , start=start,fixed=fixed)
print(fit1$param)
# Maximum likelihood:
fit2 <- GeoFit(data=data,coordx=coords,corrmodel=corrmodel, model=model,
                    likelihood="Full",type="Standard",
                   optimizer="BFGS" , start=start,fixed=fixed)
print(fit2$param)


################################################################
###
### Example 6. Maximum pairwise likelihood fitting of 
### Binomial and negative Binomial   RFs 
### with exponential correlation using  BFGS.
###
###############################################################

set.seed(422)
N=350
x <- runif(N, 0, 1)
y <- runif(N, 0, 1)
coords <- cbind(x,y)
mean=0.1; mean1=0.8; mean2=-0.5 # regression parameters
X=cbind(rep(1,N),runif(N),runif(N)) # marix covariates
corrmodel <- "Wend0"
param=list(mean=mean,mean1=mean1,mean2=mean2,sill=1,nugget=0,scale=0.2,power2=4)
# Simulation of the spatial Binomial-Gaussian RF:
data <- GeoSim(coordx=coords, corrmodel=corrmodel, model="Binomial", n=10,X=X,
             param=param)$data

fixed <- list(nugget=nugget,power2=4,sill=1)
start <- list(scale=0.2,mean=mean,mean1=mean1,mean2=mean2)
# Maximum pairwise likelihood:
fit1 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel,n=10, X=X, start=start,
                    maxdist=0.05,model="Binomial",  fixed=fixed, optimizer="BFGS")
                  
print(fit1)

set.seed(220)
# Simulation of the spatial Negative Binomial-Gaussian RF:
data <- GeoSim(coordx=coords, corrmodel=corrmodel, model="BinomialNeg", n=5,X=X,
              param=param)$data
# Maximum pairwise likelihood:
fit2 <- GeoFit(data=data, coordx=coords, corrmodel=corrmodel, n=5,X=X, start=start,
              maxdist=0.05, model="BinomialNeg", fixed=fixed, optimizer="BFGS")
                

print(fit2)


###############################################################
######### Examples of spatio-temporal RFs ###########
###############################################################
set.seed(52)
# Define the temporal sequence:
time <- seq(1, 10, 1)

# Define the spatial-coordinates of the points:
x <- runif(20, 0, 1)
set.seed(42)
y <- runif(20, 0, 1)
coords=cbind(x,y)

# Set the covariance model's parameters:
corrmodel="Exp_Exp"
scale_s=0.2/3
scale_t=1
sill=1
nugget=0
mean=0

param<-list(mean=0,scale_s=scale_s,scale_t=scale_t,
            sill=sill,nugget=nugget)

# Simulation of the spatial-temporal Gaussian RF:
data <- GeoSim(coordx=coords,coordt=time,corrmodel=corrmodel,
              param=param)$data

################################################################
###
### Example 7. Maximum pairwise likelihood fitting of a
### space time Gaussian RF with double-exponential correlation
### using  BFGS
###
###############################################################
# Fixed parameters
fixed<-list(nugget=nugget)
# Starting value for the estimated parameters
start<-list(mean=mean,scale_s=scale_s,scale_t=scale_t,sill=sill)

# Maximum composite-likelihood fitting of the RF:
fit <- GeoFit(data=data,coordx=coords,coordt=time,
                    corrmodel="Exp_Exp",maxtime=1,maxdist=0.3,
                    likelihood="Marginal",type="Pairwise",
                    optimizer="BFGS" , start=start,fixed=fixed)
print(fit)


################################################################
###
### Example 8. Maximum standard likelihood fitting of a
### space time Gaussian RF observed on  dynamical spatial coordinates
### with double-exponential correlation using  BFGS
### 
###############################################################

maxN=50
coordx_dyn=list()
set.seed(31)
for(k in 1:length(time))
{
NN=sample(1:maxN,size=1)
x <- runif(NN, 0, 1)
y <- runif(NN, 0, 1)
coordx_dyn[[k]]=cbind(x,y)
}

data <- GeoSim(coordx_dyn=coordx_dyn, coordt=time, corrmodel="Exp_Exp",
                param=param)$data
I=Inf
lower<-list(mean=-I,scale_s=0,scale_t=0,sill=0)
upper<-list(mean=-I,scale_s=I,scale_t=I,sill=I)
fit <- GeoFit(data=data,coordx_dyn=coordx_dyn,coordt=time,
              corrmodel="Exp_Exp",likelihood="Full",type="Standard",
              optimizer="BFGS", lower=lower,upper=upper,
              start=start,fixed=fixed)
print(fit)




###############################################################
######### Examples of  spatial bivariate RFs ###########
###############################################################


################################################################
###
### Example 9. Maximum, and pairwise  likelihood fitting of a
### bivariate Gaussian RF with separable Bivariate  matern 
### (cross) correlation model using BFGS and a parallel 
### implementation of L-BFGS-B
###############################################################

# Define the spatial-coordinates of the points:
set.seed(5)
x <- runif(250, 0, 1)
y <- runif(250, 0, 1)
coords=cbind(x,y)
# parameters
param=list(mean_1=0,mean_2=0,scale=0.1,smooth=0.5,sill_1=1,sill_2=1,
           nugget_1=0,nugget_2=0,pcol=0.2)
# Simulation of a spatial Gaussian RF:
data <- GeoSim(coordx=coords, corrmodel="Bi_Matern_sep", 
              param=param)$data

# selecting fixed and estimated parameters
fixed=list(nugget_1=0,nugget_2=0,smooth=0.5)
start=list(mean_1=0,mean_2=0,sill_1=var(data[1,]),sill_2=var(data[2,]),
           scale=0.1,pcol=cor(data[1,],data[2,]))


# Maximum pairwise likelihood

fitcl<- GeoFit(data=data, coordx=coords, corrmodel="Bi_Matern_sep",
                     likelihood="Marginal",type="Pairwise",
                     optimizer="BFGS" , start=start,fixed=fixed,
                     maxdist=c(0.05,0.05,0.05))
print(fitcl)

I=Inf
lower=list(mean_1=-I,mean_2=-I,sill_1=0,sill_2=0,
           scale=0,pcol=-1)
upper=list(mean_1=I,mean_2=I,sill_1=I,sill_2=I,
           scale=I,pcol=1)
# Maximum likelihood :
require(optimParallel)
fitml<- GeoFit(data=data, coordx=coords, likelihood="Full",
               corrmodel="Bi_Matern_sep", type="Standard",
               upper=upper,lower=lower,optimizer="L-BFGS-B", 
               parallel=TRUE,
               start=start,fixed=fixed)
print(fitml)

}
\keyword{Composite}
