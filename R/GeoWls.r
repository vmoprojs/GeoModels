####################################################
### File name: WeightedLeastSquare.r
####################################################


### Procedures are in alphabetical order.

print.GeoWLS <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    if(x$model=='Gaussian'||x$model=='Gauss'){process <- x$model
                                              model <- x$model}
    if(x$weighted)
      method <- 'Weighted Least Squares'
    else
      method <- method <- 'Least Squares'

    cat('\n##############################################################')
    cat('\nResults:', method,'Fitting of', process, 'Random Fields.\n')
    cat('\nModel used from the', method, ':', model, '\n')
    cat('\nCovariance model:', x$corrmodel, '\n')
    cat('Number of spatial coordinates:', x$numcoord, '\n')
    cat('Number of dependent temporal realisations:', x$numtime, '\n')
    cat('Number of estimated parameters:', length(x$param), '\n')
    cat('The value of the', method, 'at the minimum:',-x$wls,'\n')
    cat('Number of spatial bins', length(x$bins),'\n')
    cat('Number of temporal bins', length(x$bint),'\n')
    #cat('Min and max spatial distances:', x$srange,'\n')
    #if(length(x$coordt)>1) cat('Min and max temporal interval:', x$trange,'\n')

    cat('\nEstimated parameters:\n')
    print.default(x$param, digits = digits, print.gap = 2,
                  quote = FALSE)

    cat('\n##############################################################\n')
    invisible(x)
  }


WlsStart <- function(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, fcall, fixed, grid,
                    likelihood, maxdist, neighb,maxtime, model, n, param, parscale,
                    paramrange, radius, start, taper, tapsep, type, varest, vartype,
                    weighted, winconst,winconst_t, winstp_t, winstp,copula,X,memdist,nosym)
  {
    # Determines the range of the parameters for a given correlation
    SetRangeParam <- function(namesparam, numparam)
    {
        low <- 0#.Machine$double.eps
        big<- Inf
        lower <- NULL
        upper <- NULL
        # Check for the param set:
        for(i in 1:numparam){ 
            if(namesparam[i]=='mean'||namesparam[i]=='mean_1'||namesparam[i]=='mean_2'){
                lower <- c(lower, -big)
                upper <- c(upper, big)}
              if(namesparam[i]==paste('mean',1,sep="")||
                 namesparam[i]==paste('mean',2,sep="")||
                 namesparam[i]==paste('mean',3,sep="")||
                 namesparam[i]==paste('mean',4,sep="")||
                 namesparam[i]==paste('mean',5,sep="")||
                 namesparam[i]==paste('mean',6,sep="")||
                 namesparam[i]==paste('mean',7,sep="")||
                 namesparam[i]==paste('mean',8,sep="")||
                 namesparam[i]==paste('mean',9,sep="")||
                 namesparam[i]==paste('mean',10,sep="")||
                 namesparam[i]==paste('mean',11,sep="")||
                 namesparam[i]==paste('mean',12,sep="")||
                 namesparam[i]==paste('mean',13,sep="")||
                 namesparam[i]==paste('mean',14,sep="")||
                 namesparam[i]==paste('mean',15,sep="")||
                 namesparam[i]==paste('mean',16,sep="")||
                 namesparam[i]==paste('mean',17,sep="")||
                 namesparam[i]==paste('mean',18,sep="")||
                 namesparam[i]==paste('mean',19,sep="")||
                 namesparam[i]==paste('mean',20,sep="")||
                 namesparam[i]==paste('mean',21,sep="")||
                 namesparam[i]==paste('mean',22,sep="")||
                 namesparam[i]==paste('mean',23,sep="")||
                 namesparam[i]==paste('mean',24,sep="")||
                 namesparam[i]==paste('mean',25,sep=""))
              {
                lower <- c(lower, -big)
                upper <- c(upper, big)}   
            if(namesparam[i]=='skew'){
                lower <- c(lower, -big)
                upper <- c(upper, big)}
            if(namesparam[i]=='tail'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
           if(namesparam[i]=='df'){
                lower <- c(lower, 0)
                upper <- c(upper, big)}
            if(namesparam[i]=='nugget'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
             if(namesparam[i]=='nugget1'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
          if(namesparam[i]=='nugget2'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i]=='power'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power_s'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power_t'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='power1'){
                lower <- c(lower, low)
                upper <- c(upper, 2)}
            if(namesparam[i]=='scale'||namesparam[i]=='scale_11'||namesparam[i]=='scale_12'||namesparam[i]=='scale_22'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
                  if(namesparam[i]=='smooth_11'||namesparam[i]=='smooth_12'||namesparam[i]=='smooth_22'){
                lower <- c(lower, low)
                upper <- c(upper, big)}  
            if(namesparam[i]=='power_2'||namesparam[i]=='power2_11'||namesparam[i]=='power2_12'||namesparam[i]=='power2_22'){
                lower <- c(lower, low)
                upper <- c(upper, big)}       
            if(namesparam[i]=='scale_s'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i]=='scale_t'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i]=='sep'){
                lower <- c(lower, low)
                upper <- c(upper, 1)}
            if(namesparam[i]=='pcol'){
                lower <- c(lower, -1)
                upper <- c(upper,  1)}
            if(namesparam[i]=='sill'||namesparam[i]=='sill_1'||namesparam[i]=='sill_2'||namesparam[i]=='nugget_1'||namesparam[i]=='nugget2_2'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
            if(namesparam[i]=='smooth'){
                lower <- c(lower, low)
                upper <- c(upper, big)}
              }
              
        return(list(lower=lower, upper=upper))
    }

    ### Initialization parameters:
    initparam <- StartParam(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, fcall, fixed,
                           grid, likelihood, maxdist,neighb, maxtime, model, n, 
                           param, parscale, paramrange, radius,  start, taper, tapsep,
                           "GeoWLS", type, varest, vartype,
                           weighted, winconst,winconst_t, winstp_t, winstp,copula, X, memdist, nosym)

       
  
    if(!is.null(initparam$error))     stop(initparam$error)
    if(length(coordt)>0&&is.list(X)) X=X[[1]]

    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    
    if(!bivariate) {if(is.null(X))  {X=1;num_betas=1} 
                    else num_betas=ncol(X)  }

    if( bivariate) {if(is.null(X))  {X=1;num_betas=c(1,1)} 
                    else 
                   { if(is.list(X)) num_betas=c(ncol(X[[1]]),ncol(X[[2]])) 
                     else num_betas=c(ncol(X),ncol(X))
                    } } 

    
    ### Set the initial type of likelihood objects:
    initparam$type <- CkType(type)
   # if(substr(model,1,6)=='Binary'||substr(model,1,8)=='Binomial'||substr(model,1,11)=='BinomialNeg'||substr(model,1,4)=='Geom'||substr(model,1,4)=='Pois') return(initparam)
    if(is.null(start)) start <- NA else start <- unlist(start)
    if(is.null(fixed)) fixed <- NA else fixed <- unlist(fixed)
    ### Checks if all the starting values have been passed by the user:



    if(initparam$numstart==initparam$numparam) {
   
        if((model %in% c('Gaussian','Gauss','Chisq','LogLogistic','Logistic','Gamma','Gamma2','Beta','Beta2','LogGaussian','LogGauss','Binomial_TwoPieceGaussian','Binomial_TwoPieceGauss',
          'Tukeygh','Tukeyh','Tukeyh2','Kumaraswamy','Kumaraswamy2','Weibull','SkewGaussian','SkewGauss','SinhAsinh','StudentT','SkewStudentT',
          "Gaussian_misp_StudentT","Gaussian_misp_Poisson","Gaussian_misp_Tukeygh",
          "Gaussian_misp_SkewStudentT","PoissonGamma","PoissonWeibull","Gaussian_misp_PoissonGamma",
          "TwoPieceStudentT",'Wrapped',"TwoPieceGaussian","TwoPieceGauss","TwoPieceTukeyh","TwoPieceBimodal")) & 
          (type %in% c('Standard','Pairwise','Tapering','Tapering1','Independence')))
        {
 
##########################################################        
        if(!initparam$bivariate){  ###spatial or temporal univariate case
          if(is.na(fixed["mean"])&is.na(fixed["mean2"])){
  
              if(is.na(start["mean"])) {initparam$param <- c(initparam$fixed["mean"], initparam$param)}
              else {initparam$param <- c(start["mean"], initparam$param)}
            initparam$namesparam <- sort(names(initparam$param))
            initparam$param <- initparam$param[initparam$namesparam]
            initparam$numparam <- initparam$numparam+1
            initparam$flagnuis['mean'] <- 1
            initparam$numfixed <- initparam$numfixed-1
          if(initparam$numfixed > 0) {initparam$fixed <- fixed}
          else {initparam$fixed <- NULL}
          }
          else {initparam$fixed['mean'] <- fixed["mean"]} 
          if(num_betas>1){
          for(i in 1:(num_betas-1)) {
            if(is.na(fixed[paste("mean",i,sep="")]))
          {
              if(is.na(start[paste("mean",i,sep="")])) {initparam$param <- c(initparam$fixed[paste("mean",i,sep="")], initparam$param)}
              else {initparam$param <- c(start[paste("mean",i,sep="")], initparam$param)}
            initparam$namesparam <- sort(names(initparam$param))
            initparam$param <- initparam$param[initparam$namesparam]
            initparam$numparam <- initparam$numparam+1
            initparam$flagnuis[paste("mean",i,sep="")] <- 1
            initparam$numfixed <- initparam$numfixed-1}
            else {initparam$fixed[paste("mean",i,sep="")] <- fixed[paste("mean",i,sep="")]} 
         }
         if(initparam$numfixed > 0) {initparam$fixed <- fixed}
          else {initparam$fixed <- NULL}
          }

      }   

################################################################################################################################################
     if(initparam$bivariate) {           ###bivariate case
   
              if(is.na(fixed["mean_1"])){
              if(is.na(start["mean_1"])) {initparam$param <- c(initparam$fixed["mean_1"], initparam$param)}
              else {initparam$param <- c(start["mean_1"], initparam$param)}
              initparam$namesparam <- sort(names(initparam$param))
              initparam$param <- initparam$param[initparam$namesparam]
              initparam$numparam <- initparam$numparam+1
              initparam$flagnuis['mean_1'] <- 1
              initparam$numfixed <- initparam$numfixed-1
              if(initparam$numfixed > 0) {initparam$fixed <- fixed}
              else {initparam$fixed <- NULL}}
              else { initparam$fixed['mean_1'] <- fixed["mean_1"] }

              if(is.na(fixed["mean_2"])){
              initparam$namesparam<-names(initparam$namesparam)
              if(is.na(start["mean_2"])) {initparam$param <- c(initparam$fixed["mean_2"], initparam$param)}
              else {initparam$param <- c(start["mean_2"], initparam$param)}
              initparam$namesparam <- sort(names(initparam$param))
              initparam$param <- initparam$param[initparam$namesparam]
              initparam$numparam <- initparam$numparam+1
              initparam$flagnuis['mean_2'] <- 1
              initparam$numfixed <- initparam$numfixed-1
              if(initparam$numfixed > 0) {initparam$fixed <- fixed}
              else {initparam$fixed <- NULL}}
              else {initparam$fixed['mean_2'] <- fixed["mean_2"]}  
        if(num_betas[1]>1){
          for(i in 1:(num_betas[1]-1)) {
            if(is.na(fixed[paste("mean_1",i,sep="")]))
          {
              if(is.na(start[paste("mean_1",i,sep="")])) {initparam$param <- c(initparam$fixed[paste("mean_1",i,sep="")], initparam$param)}
              else {initparam$param <- c(start[paste("mean_1",i,sep="")], initparam$param)}
            initparam$namesparam <- sort(names(initparam$param))
            initparam$param <- initparam$param[initparam$namesparam]
            initparam$numparam <- initparam$numparam+1
            initparam$flagnuis[paste("mean_1",i,sep="")] <- 1
            initparam$numfixed <- initparam$numfixed-1}
            else {initparam$fixed[paste("mean_1",i,sep="")] <- fixed[paste("mean_1",i,sep="")]} 
         }}
           if(num_betas[2]>1){
          for(i in 1:(num_betas[2]-1)) {
            if(is.na(fixed[paste("mean_2",i,sep="")]))
          {
              if(is.na(start[paste("mean_2",i,sep="")])) {initparam$param <- c(initparam$fixed[paste("mean_2",i,sep="")], initparam$param)}
              else {initparam$param <- c(start[paste("mean_2",i,sep="")], initparam$param)}
            initparam$namesparam <- sort(names(initparam$param))
            initparam$param <- initparam$param[initparam$namesparam]
            initparam$numparam <- initparam$numparam+1
            initparam$flagnuis[paste("mean_2",i,sep="")] <- 1
            initparam$numfixed <- initparam$numfixed-1}
            else {initparam$fixed[paste("mean_2",i,sep="")] <- fixed[paste("mean_2",i,sep="")]} 
         }}}
###########################
              
        }
        paramrange=TRUE
        if(paramrange) paramrange <- SetRangeParam(names(initparam$param), length(initparam$param))
        else  paramrange <- list(lower=NULL, upper=NULL)
        initparam$lower<-paramrange$lower
        initparam$upper<-paramrange$upper
        return(initparam)
      }
        initparam$error="Some starting and/or fixed parameters are missing. (All the covariance and nuisance  parameters must be included) "
        return(initparam)
  }

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
GeoWLS <- function(data, coordx, coordy=NULL, coordt=NULL,  coordx_dyn=NULL, corrmodel, distance="Eucl",
                         fixed=NULL,grid=FALSE, maxdist=NULL, neighb=NULL,maxtime=NULL, model='Gaussian',
                         optimizer='Nelder-Mead', numbins=NULL, radius=6371,  start=NULL,
                         weighted=FALSE)
  {
    ### Check first if the model is not binary:
    if(substr(model,1,6)=='Binary'||substr(model,1,6)=='Binomial') stop("The weighted least squares method can not be used with binary data")

    call <- match.call()
    ### Check the parameters given in input:
    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance,"Fitting", fixed, grid, 'None',
                              maxdist, maxtime,  model,NULL,  optimizer, NULL, radius, start, NULL,
                             NULL, 'GeoWLS', FALSE, 'SubSamp', weighted, NULL,NULL)
    

    if(!is.null(checkinput$error))
      stop(checkinput$error)
    # check the number of bins:
    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')
    # set the default number of spatial bins:
    if(is.null(numbins))
      numbins <- 13
    ### Define the object function for the weighted least squares method:
    WLsquare <- function(bins, bint, corrmodel, fixed, fun, lenbins, moments,
                         namescorr, namesnuis, numbins, numbint, param)
      {
        param <- c(param, fixed)#set the parameters set:
        paramcorr <- param[namescorr]#set the correlation parameters:
        nuisance <- param[namesnuis]#set the nuisance parameters:
        #computes the weighted least squares:
        result <- .C(fun, as.double(bins), as.double(bint), as.integer(corrmodel),
                     as.double(lenbins), as.double(moments), as.integer(numbins),
                     as.integer(numbint), as.double(nuisance), as.double(paramcorr),
                     res=double(1), PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)$res
       return(result)
      }
    ### Initializes global variables:
    GeoWLS <- NULL
    fname <- NULL
    variogramt <- NULL
    variogramst <- NULL
    ### Initializes the parameter values:
    parscale <- NULL
    initparam <- StartParam(coordx, coordy, coordt, coordx_dyn,corrmodel, data, distance, "Fitting", fixed, grid,
    'None', maxdist, neighb,maxtime,  model, NULL,  NULL,
                           parscale, optimizer=='L-BFGS-B', radius, start,NULL,  NULL,
                           'GeoWLS', 'GeoWLS', FALSE, 'SubSamp', FALSE, 1, 1,1,1,NULL, NULL,0)
    if(!is.null(initparam$error))
      stop(initparam$error)
     coordx=initparam$coordx
     coordy=initparam$coordy
    ###### ----------- START Estimation of the empirical variogram ---------- #####
    numvario <- numbins-1
    bins <- double(numbins) # vector of spatial bins
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    ### Checks the type of variogram:
    fname <- 'Binned_Variogram'


    if(initparam$spacetime){### Computes the spatial-temporal variogram:
      spacetime_dyn=FALSE
      if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
      ns=initparam$ns
      NS=cumsum(ns)
      numbint <- initparam$numtime-1 # number of temporal bins
      bint <- double(numbint)        # vector temporal bins
      momentt <- double(numbint)     # vector of temporal moments
      lenbint <- integer(numbint)    # vector of temporal bin sizes
      numbinst <- numvario*numbint    # number of spatial-temporal bins
      binst <- double(numbinst)      # spatial-temporal bins
      momentst <- double(numbinst)   # vector of spatial-temporal moments
      lenbinst <- integer(numbinst)  # vector of spatial-temporal bin sizes
      if(!spacetime_dyn){
                                  data=c(t(data))
                                  coordx=rep(coordx,length(coordt))
                                  coordy=rep(coordy,length(coordt))
                         }
      if(spacetime_dyn) data=unlist(data)
         NS=c(0,NS)[-(length(ns)+1)]
      fname <- 'Binned_Variogram_st';fname <- paste(fname,"2",sep="") 
      # Compute the spatial-temporal moments:
      EV=.C(fname, bins=bins, bint=bint, as.double(coordx),as.double(coordy),as.double(coordt),as.double(initparam$data),
           lenbins=lenbins,lenbinst=lenbinst,lenbint=lenbint,moments= moments, momentst=momentst, momentt=momentt, as.integer(numbins), as.integer(numbint),
           as.integer(ns),as.integer(NS), PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
      bins=EV$bins
      bint=EV$bint
      lenbins=EV$lenbins
      lenbint=EV$lenbint
      lenbinst=EV$lenbinst
      moments=EV$moments
      momentst=EV$momentst
      momentt=EV$momentt
      indbin <- lenbins>0
      centers <- bins[1:numvario]+diff(bins)/2
      bins <- bins[indbin]
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      # Computes the spatial marginal variogram:
      variograms <- moments/lenbins
      numbins <- sum(indbin)
      indbint <- lenbint>0
      bint <- bint[indbint]
      momentt <- momentt[indbint]
      lenbint <- lenbint[indbint]
      numbint <- sum(indbint)
      # Computes the temporal marginal variogram:
      variogramt <- momentt/lenbint
      indbinst <- lenbinst>0
      momentst <- momentst[indbinst]
      lenbinst <- lenbinst[indbinst]
      numbinst <- sum(indbinst)
      # Computes the spatial-temporal variogram:
      variogramst <- momentst/lenbinst
      # Set the moment vectors and their sizes:
      moment <- matrix(momentst,nrow=numbins,ncol=numbint,byrow=TRUE)
      lenbin <- matrix(lenbinst,nrow=numbins,ncol=numbint,byrow=TRUE)
      moment <- rbind(momentt, moment)
      moment <- cbind(c(0,moments),moment)
      lenbin <- rbind(lenbint, lenbin)
      lenbin <- cbind(c(1,lenbins),lenbin)
      moments <- moment
      lenbins <- lenbin
      bins <- c(-bins[1],bins)
      bint <- c(0,bint)
      numbins <- numbins+1
      numbint <- numbint+1}
      # Set an initial value for the scale parameter:
    #if(!is.null(initparam$param['scale_s']))
          #initparam$param['scale_s'] <- bins[max(variograms)==variograms]}
    else{### Computes the spatial variogram:
      numbint <- 1 # number of temporal bins
      bint <- double(numbint) # vector temporal bins
      momentt <- double(1) # vector of temporal moments
      momentst <- double(1)   # vector of spatial-temporal moments
      lenbint <- integer(1) # vector of temporal bin sizes
      lenbinst <- integer(1)  # vector of spatial-temporal bin sizes
      fname <- paste(fname,"2",sep="")
      EV=.C(fname, bins=bins, as.double(coordx),as.double(coordy),as.double(coordt),as.double(initparam$data), lenbins=lenbins,
         moments=moments, as.integer(numbins),PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
      bins=EV$bins
      lenbins=EV$lenbins
      moments=EV$moments
     
      indbin <- lenbins>0
      centers <- bins[1:numvario]+diff(bins)/2
      bins <- bins[indbin]
      centers <- centers[indbin]
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      numbins <- sum(indbin)
      variograms <- moments/lenbins
 
      #plot(centers,variograms)
      # Set an initial value for the scale parameter:
      if(!is.null(initparam$param['scale']))
        initparam$param['scale'] <- bins[max(variograms)==variograms]}
    ###### ----------- END Estimation of the empirical variogram ---------- #####

    ###### ---------- START model fitting by weighted least squares method ----------######
    if(model=='Gaussian'||model=='Gauss') # Gaussian random field:
      if(weighted) fname <- 'GeoWLS_G'
      else fname <- 'LeastSquare_G'

    ### Computes estimates by the weighted least squares method:
    if(optimizer=='L-BFGS-B')
      fitted <- optim(initparam$param, WLsquare, bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                      numbins=numbins, numbint=numbint, control=list(fnscale=-1, factr=1, pgtol=1e-14,
                      maxit = 1e8), lower=initparam$lower, upper=initparam$upper, hessian=FALSE)
    else
      fitted <- optim(initparam$param, WLsquare, bins=bins, bint=bint, corrmodel=initparam$corrmodel,
                      fixed=initparam$fixed, fun=fname, lenbins=lenbins, method=optimizer,
                      moments=moments, namescorr=initparam$namescorr, namesnuis=initparam$namesnuis,
                      numbins=numbins, numbint=numbint, control=list(fnscale=-1, reltol=1e-14, maxit=1e8),
                      hessian=FALSE)
    ###### ---------- END model fitting by weighted least squares method ----------######
    ### Removes the global variobales:
     .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)  
    ### Set the output:
    GeoWLS <- list(bins=bins,
                         bint=bint,
                         centers=centers,
                         coordx = initparam$coordx,
                         coordy = initparam$coordy,
                         coordt = initparam$coordt,
                         coordx_dyn=coordx_dyn,
                         convergence = fitted$convergence,
                         corrmodel = corrmodel,
                         data = initparam$data,
                         fixed = initparam$fixed,
                         grid = grid,
                         iterations = fitted$counts,
                         maxdist =maxdist,
                         maxtime = maxtime,
                         message = fitted$message,
                         model=model,
                         numcoord=initparam$numcoord,
                         numtime=initparam$numtime,
                         param = fitted$par,
                         variograms = variograms,
                         variogramt = variogramt,
                         variogramst = variogramst,
                         weighted = weighted,
                         wls = fitted$value)
    structure(c(GeoWLS, call = call), class = c("GeoWLS"))
  }
