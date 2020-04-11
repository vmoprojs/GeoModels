####################################################
### Emails: moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: Fitting.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 27/01/2020.
####################################################


### Procedures are in alphabetical order.

### Fitting procedure:

GeoFit <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl",
                         fixed=NULL,GPU=NULL, grid=FALSE, likelihood='Marginal', local=c(1,1),
                         lower=NULL,maxdist=NULL,
                          maxtime=NULL, memdist=TRUE,method="cholesky", model='Gaussian',n=1, onlyvar=FALSE ,
                          optimizer='Nelder-Mead', parallel=FALSE,
                         radius=6371,  sensitivity=FALSE,sparse=FALSE, start=NULL, taper=NULL, tapsep=NULL, 
                         type='Pairwise', upper=NULL, varest=FALSE, vartype='SubSamp', weighted=FALSE, winconst=NULL, winstp=NULL, 
                         winconst_t=NULL, winstp_t=NULL,X=NULL)
{
    call <- match.call()

    #CMdl<-CkCorrModel(corrmodel)
    #Stime <- CheckST(CMdl)  
    
    ### Check the parameters given in input:
    corrmodel=gsub("[[:blank:]]", "",corrmodel)
    model=gsub("[[:blank:]]", "",model)
    distance=gsub("[[:blank:]]", "",distance)
    optimizer=gsub("[[:blank:]]", "",optimizer)
    likelihood=gsub("[[:blank:]]", "",likelihood)
    type=gsub("[[:blank:]]", "",type)
    if(!is.logical(memdist)) memdist=FALSE
 
    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting",
                             fixed, grid, likelihood, maxdist, maxtime, model, n,
                              optimizer, NULL, radius, start, taper, tapsep, 
                             type, varest, vartype, weighted, X)
   
    if(!is.null(checkinput$error))
      stop(checkinput$error)
    ### Initialization global variables:
    GeoFit <- NULL
    score <- sensmat <- varcov <- varimat <- parscale <- NULL
    ### Initialization parameters:
    unname(coordt);
    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy)}

    initparam <- WlsStart(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting", fixed, grid,#10
                         likelihood, maxdist, maxtime,  model, n, NULL,#16
                         parscale, optimizer=='L-BFGS-B', radius, start, taper, tapsep,#22
                         type, varest, vartype, weighted, winconst, winstp,winconst_t, winstp_t, X,memdist)#32


    if(!is.null(initparam$error))   stop(initparam$error)
    ## checking for upper and lower bound for method 'L-BFGS-B' and optimize method

    if(optimizer %in% c('L-BFGS-B','nlminb','nmkb') || length(initparam$param)==1){
   # if(optimizer=='L-BFGS-B'|| optimizer=='nlminb' || optimizer=='nmkb' ||optimizer=='lbfgsb3c' || length(initparam$param)==1){

    if(!is.null(lower)||!is.null(upper)){
       if(!is.list(lower)||!is.list(upper))  stop("lower and upper bound must be a list\n")
    #setting alphabetic order
      lower=lower[order(names(lower))]
      upper=upper[order(names(upper))] 
      npar<-length(initparam$param) 
      ll<-as.numeric(lower);uu<-as.numeric(upper)
      if(length(ll)!=npar||length(uu)!=npar)
           stop("lower and upper bound must be of the same length of starting values\n") 

      #if(sum(uu<=initparam$upper)<npar||sum(ll>=initparam$lower)>npar)
      #     stop("one or more values of the lower and upper bounds are out of the valid  range\n")  
      if(sum(names(initparam$param)==names(upper))<npar || sum(names(initparam$param)==names(lower))<npar){
           stop("the names of  parameters in the lower and/or  upper bounds do not match with starting parameters names .\n") }
      ll[ll==0]=.Machine$double.eps ## when 0 we don't want exactly zero
      uu[uu==Inf]=1e+12
      initparam$upper <- uu;initparam$lower <- ll
     }}
      #if(length(param)==1)  {initparam$upper=1000}

   # Full likelihood:
    if(likelihood=='Full')
          # Fitting by log-likelihood maximization:
         fitted <- Lik(initparam$bivariate,initparam$coordx,initparam$coordy,initparam$coordt,coordx_dyn, initparam$corrmodel,
                               unname(initparam$data),initparam$fixed,initparam$flagcorr,
                               initparam$flagnuis,grid,initparam$lower,method,initparam$model,initparam$namescorr,
                               initparam$namesnuis,initparam$namesparam,initparam$numcoord,initparam$numpairs,
                               initparam$numparamcorr,initparam$numtime,optimizer,onlyvar,parallel,
                               initparam$param,initparam$radius,initparam$setup,initparam$spacetime,sparse,varest,taper,initparam$type,
                               initparam$upper,initparam$ns,unname(initparam$X))

    # Composite likelihood:
    if(likelihood=='Marginal' || likelihood=='Conditional' || likelihood=='Marginal_2'){
  
    if(!memdist)
          fitted <- CompLik(initparam$bivariate,initparam$coordx,initparam$coordy,initparam$coordt,coordx_dyn,initparam$corrmodel,unname(initparam$data), #6
                                   initparam$distance,initparam$flagcorr,initparam$flagnuis,initparam$fixed,GPU,grid, #12
                                   initparam$likelihood,local, initparam$lower,initparam$model,initparam$n,#17
                                   initparam$namescorr,initparam$namesnuis,#19
                                   initparam$namesparam,initparam$numparam,initparam$numparamcorr,optimizer,onlyvar,parallel,
                                   initparam$param,initparam$spacetime,initparam$type,#27
                                   initparam$upper,varest,initparam$vartype,initparam$weighted,initparam$winconst,initparam$winstp,#33
                                   initparam$winconst_t,initparam$winstp_t,initparam$ns,
                                   unname(initparam$X),sensitivity)
    if(memdist)
        fitted <- CompLik2(initparam$bivariate,initparam$coordx,initparam$coordy,initparam$coordt,
                                   coordx_dyn,initparam$corrmodel,unname(initparam$data), #6
                                   initparam$distance,initparam$flagcorr,initparam$flagnuis,initparam$fixed,GPU,grid, #12
                                   initparam$likelihood,local, initparam$lower,initparam$model,initparam$n,#17
                                   initparam$namescorr,initparam$namesnuis,#19
                                   initparam$namesparam,initparam$numparam,initparam$numparamcorr,optimizer,onlyvar,parallel,
                                   initparam$param,initparam$spacetime,initparam$type,#27
                                   initparam$upper,varest,initparam$vartype,initparam$weighted,initparam$winconst,initparam$winstp,#33
                                   initparam$winconst_t,initparam$winstp_t,initparam$ns,
                                   unname(initparam$X),sensitivity,initparam$colidx,initparam$rowidx)
      }

    numtime=1
    if(initparam$spacetime) numtime=length(coordt)
    if(initparam$bivariate) numtime=2
    dimat <- initparam$numcoord*numtime#
    if(is.null(dim(initparam$X)))  initparam$X=as.matrix(rep(1,dimat))
    # Delete the global variables:
    .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE) #if(mem)
    ### Set the output object:
    GeoFit <- list(bivariate=initparam$bivariate,
                         claic = fitted$claic,
                         clbic = fitted$clbic,
                         coordx = initparam$coordx,
                         coordy = initparam$coordy,
                         coordt = initparam$coordt,
                         coordx_dyn=coordx_dyn,
                         convergence = fitted$convergence,
                         corrmodel = corrmodel,
                         data = initparam$data,
                         distance = distance,
                         fixed = initparam$fixed,
                         grid = grid,
                         iterations = fitted$counts,
                         likelihood = likelihood,
                         logCompLik = fitted$value,
                         message = fitted$message,
                         model = model,
                         n=initparam$n,
                         ns=initparam$ns,
                         numbetas=initparam$num_betas,
                         numcoord=initparam$numcoord,
                         numtime=initparam$numtime,
                         optimizer=optimizer,
                         param = fitted$par,
                         nozero = initparam$setup$nozero,
                         score = fitted$score,
                         maxdist =maxdist,
                         maxtime = maxtime,
                         radius = radius,
                         spacetime = initparam$spacetime,
                         stderr = fitted$stderr,
                         sensmat = fitted$sensmat,
                         varcov = fitted$varcov,
                         varimat = fitted$varimat,
                         vartype = vartype,
                         type = type,
                         weighted=initparam$weighted,
                         winconst = initparam$winconst,
                         winstp = initparam$winstp,
                         winconst_t = initparam$winconst_t,
                         winstp_t = initparam$winstp_t,
                         X = initparam$X)
    structure(c(GeoFit, call = call), class = c("GeoFit"))
  }

print.GeoFit <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(x$likelihood=='Full'){
        method <- 'Likelihood'
        if(x$type=="Tapering") {claic <- "CLAIC";clbic <- "CLBIC";}
        else { claic <- "AIC";clbic <- "BIC";    }
      }
    else{
        method <- 'Composite-Likelihood'; claic <- 'CLAIC';clbic <- 'CLBIC';}
  if(x$model=='Gaussian'||x$model=='Gauss'){ process <- 'Gaussian';model <- 'Gaussian'}
  if(x$model=='Gamma') { process <- 'Gamma'; model <- 'Gamma'}
    if(x$model=='TwoPieceBimodal') { process <- 'TwoPieceBimodal'; model <- 'TwoPieceBimodal'}
  if(x$model=='LogLogistic') { process <- 'LogLogistic'; model <- 'LogLogistic'}
    if(x$model=='Gaussian_misp_Poisson') { process <- 'Poisson'; model <- 'Misspecified Gaussian '}
     if(x$model=='Poisson') { process <- 'Poisson'; model <- 'Poisson'}
  if(x$model=='Gaussian_misp_StudentT') { process <- 'StudentT'; model <- 'Misspecified Gaussian '}
    if(x$model=='Gaussian_misp_SkewStudentT') { process <- 'SkewStudentT'; model <- 'Misspecified Gaussian '}
  if(x$model=='Logistic') { process <- 'Logistic'; model <- 'Logistic'}
    if(x$model=='Tukeyh') { process <- 'Tukeyh'; model <- 'Tukeyh'}
  if(x$model=='Gamma2'){ process <- 'Gamma2'; model <- 'Gamma2'}
  if(x$model=='LogGauss'||x$model=='LogGaussian'){ process <- 'Log Gaussian'; model <- 'LogGaussian'}
  if(x$model=='SkewGauss'||x$model=='SkewGaussian'){ process <- 'Skew Gaussian';model <- 'SkewGaussian'}
  if(x$model=='StudentT'){ process <- 'StudentT';model <- 'StudentT'}
  if(x$model=='TwoPieceStudentT'){ process <- 'TwoPiece StudentT';model <- 'TwoPieceStudentT'}
  if(x$model=='TwoPieceTukeyh'){ process <- 'TwoPiece Tukeyh';model <- 'TwoPieceTukeyh'}
  if(x$model=='TwoPieceGaussian'||x$model=='TwoPieceGauss'){ process <- 'TwoPiece Gaussian';model <- 'TwoPieceGaussian'}
  if(x$model=='SinhAsinh'){ process <- 'SinhAsinh'; model <- 'SinhAsinh'}    
  if(x$model=='Wrapped'){ process <- 'Wrapped'; model <- 'Wrapped'}
    if(x$model=='Weibull'){ process <- 'Weibull'; model <- 'Weibull'}
  if(x$model=='Binomial'){ process <- 'Binomial';model <- 'Binomial'}
  if(x$model=='Kumaraswamy'){ process <- 'Kumaraswamy';model <- 'Kumaraswamy'}
  if(x$model=='Binomial_TwoPieceGaussian'||x$model=='Binomial_TwoPieceGauss'){ process <- 'Binomial TwoPiece Gaussian';model <- 'Binomial_TwoPieceGauss'}
  if(x$model=='BinomialNeg_TwoPieceGaussian'||x$model=='BinomialNeg_TwoPieceGauss'){ process <- 'Negative Binomial TwoPiece Gaussian';model <- 'BinomialNeg_TwoPieceGauss'}
  if(x$model=='Binomial2'){ process <- 'Binomial';model <- 'Binomial2'}     
  if(x$model=='BinomialNeg'){ process <- 'BinomialNeg'; model <- 'BinomialNeg'}
  if(x$model=='Geom'||x$model=='Geometric'){ process <- 'Geometric';model <- 'Geometric'}
  if(x$model=='PoisBin'){ process <- 'Poisson Binomial';model <- 'PoisBin'}
  if(x$model=='PoisBinNeg'){ process <- 'Poisson NegBinomial';model <- 'PoisBinNeg'}    
  if(x$bivariate){ biv <- 'bivariate';x$numtime=1}
  else { biv <- 'univariate'}                       

    cat('\n##################################################################')
    cat('\nMaximum', method, 'Fitting of', process, 'Random Fields\n')
    cat('\nSetting:', x$likelihood, method, '\n')
    cat('\nModel:', model, '\n')
    cat('\nType of the likelihood objects:', x$type, x$method,'\n')
    cat('\nCovariance model:', x$corrmodel, '\n')
    cat('\nOptimizer:', x$optimizer, '\n')
    cat('\nNumber of spatial coordinates:', x$numcoord, '\n')
    cat('Number of dependent temporal realisations:', x$numtime, '\n')
    cat('Type of the random field:', biv, '\n')
    cat('Number of estimated parameters:', length(x$param), '\n')
    cat('\nType of convergence:', x$convergence, '')
    cat('\nMaximum log-', method, ' value: ',
        format(x$logCompLik, digits = digits, nsmall = 2), '\n', sep='')

    if(!is.null(x$claic))
      cat(claic,':', format(x$claic, digits = digits),'\n')

    if(!is.null(x$clbic))
      cat(clbic,':', format(x$clbic, digits = digits),'\n')  

    cat('\nEstimated parameters:\n')
    print.default(x$param, digits = digits, print.gap = 2,
                  quote = FALSE)

    if(!is.null(x$stderr))
      {
        cat('\nStandard errors:\n')
        print.default(x$stderr, digits = digits, print.gap = 2,
                      quote = FALSE)
      }

    #if(!is.null(x$varcov))
      #{
      #  cat('\nVariance-covariance matrix of the estimates:\n')
      #  print.default(x$varcov, digits = digits, print.gap = 3,
      #                quote = FALSE)
     # }

    cat('\n##################################################################\n')
    invisible(x)
  }

