####################################################
### File name: GeoFit.r
####################################################


GeoFit2 <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,copula=NULL,corrmodel, distance="Eucl",
                         fixed=NULL,anisopars=NULL,est.aniso=c(FALSE,FALSE),GPU=NULL, grid=FALSE, likelihood='Marginal', local=c(1,1),
                         lower=NULL,maxdist=Inf,neighb=NULL,
                          maxtime=Inf, memdist=TRUE,method="cholesky", model='Gaussian',n=1, onlyvar=FALSE ,
                          optimizer='Nelder-Mead', parallel=FALSE,
                         radius=6371,  sensitivity=FALSE,sparse=FALSE, start=NULL, taper=NULL, tapsep=NULL, 
                         type='Pairwise', upper=NULL, varest=FALSE, vartype='SubSamp', weighted=FALSE, winconst=NULL, winstp=NULL, 
                         winconst_t=NULL, winstp_t=NULL,X=NULL,nosym=FALSE)
{
    call <- match.call()

  
    if(!is.null(copula))
     { if((copula!="Clayton")&&(copula!="Gaussian")) stop("the type of copula is wrong")}

    if(type=='Independence'&&likelihood!='Marginal') stop("Independence likelihood must be coupled with 
        Marginal likelihood")
    ### Check the parameters given in input:
      if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")
    corrmodel=gsub("[[:blank:]]", "",corrmodel)
    model=gsub("[[:blank:]]", "",model)
    distance=gsub("[[:blank:]]", "",distance)
    optimizer=gsub("[[:blank:]]", "",optimizer)
    likelihood=gsub("[[:blank:]]", "",likelihood)
    type=gsub("[[:blank:]]", "",type)
    if(!is.logical(memdist)) memdist=FALSE
    if(!is.null(X)) X=as.matrix(X)
    if(is.numeric(neighb)) {
            neighb=round(neighb)
            if(all(neighb<1))  stop("neighb must be an integer >=1")
          }
    if(!is.null(anisopars)) {if(!is.list(anisopars)) stop("anisopars must be a list with two elements")}
    
    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting",
                             fixed, grid, likelihood, maxdist, maxtime, model, n,
                              optimizer, NULL, radius, start, taper, tapsep, 
                             type, varest, vartype, weighted,copula, X)
   
    if(!is.null(checkinput$error))
      stop(checkinput$error)
    ### Initialization global variables:
    GeoFit <- NULL
    score <- sensmat <- varcov <- varimat <- parscale <- NULL
    ### Initialization parameters:
    cooordt=unname(coordt);
    if(is.null(coordx_dyn)){
    coordx=unname(coordx);coordy=unname(coordy)}

    initparam <- WlsStart(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting", fixed, grid,#10
                         likelihood, maxdist,neighb,maxtime,  model, n, NULL,#16
                         parscale, optimizer=='L-BFGS-B', radius, start, taper, tapsep,#22
                         type, varest, vartype, weighted, winconst, winstp,winconst_t, winstp_t, copula,X,memdist,nosym)#32

   


   ## in the case on external fixed mean

  MM=NULL
  if(is.na(initparam$fixed['mean'])&length(c(initparam$X))==1) {MM=fixed$mean;initparam$mean=0}
  


    ## moving sill from starting to fixed parameters if necessary
        if(sum(initparam$namesparam=='sill')==1){
    if(initparam$model %in%  c(2,14,16,21,42,50,26,24,30,46,43,11)) 
    {initparam$param=initparam$param[initparam$namesparam!='sill'];initparam$namesparam=names(initparam$param)
    a=1; names(a)="sill";initparam$fixed=c(initparam$fixed,a)}}


    if(!is.null(initparam$error))   stop(initparam$error)
    ## checking for upper and lower bound for method 'L-BFGS-B' and optimize method

      if(!(optimizer %in% c('L-BFGS-B','nlminb','nlm','nmkb','nmk','multiNelder-Mead','multinlminb',"BFGS","Nelder-Mead","ucminf","optimize","SANN")))
             stop("optimizer is not correct\n")
     ####        
    if(optimizer %in% c('L-BFGS-B','nlminb','nmkb','multinlminb','multiNelder-Mead') || length(initparam$param)==1){
   
    if(!is.null(lower)||!is.null(upper)){
       if(!is.list(lower)||!is.list(upper))  stop("lower and upper bound must be a list\n")

       if(sum(unlist(lower)>unlist((upper)))>0) stop("some values of the lower bound is greater of the upper bound \n")
    #setting alphabetic order

      if(sum(names(lower)=='sill')==1){
          if(initparam$model %in%  c(2,14,16,21,42,50,26,24,30,46,43,11)) 
            {lower=lower[names(lower)!='sill'];upper=upper[names(upper)!='sill']; }}

      lower=lower[order(names(lower))]
      upper=upper[order(names(upper))] 
      npar<-length(initparam$param) 
      ll<-as.numeric(lower);uu<-as.numeric(upper)
      if(length(ll)!=npar||length(uu)!=npar)
           stop("lower and upper bound must be of the same length of starting values\n")   
      if(sum(names(initparam$param)==names(upper))<npar || sum(names(initparam$param)==names(lower))<npar){
           stop("the names of  parameters in the lower and/or  upper bounds do not match with starting parameters names .\n") }
      ll[ll==0]=.Machine$double.eps ## when 0 we don't want exactly zero
      uu[uu==Inf]=1e+12

      initparam$upper <- uu;initparam$lower <- ll
     }}


###############################################################################################
    fitted_ini<-CompIndLik2(initparam$bivariate,initparam$coordx,initparam$coordy,initparam$coordt,
                                   coordx_dyn,unname(initparam$data), 
                                   initparam$flagcorr,initparam$flagnuis,initparam$fixed,grid,
                                    initparam$lower,initparam$model,initparam$n ,
                                     initparam$namescorr,initparam$namesnuis,
                                   initparam$namesparam,initparam$numparam,optimizer,onlyvar,parallel, initparam$param,initparam$spacetime,initparam$type,#27
                                   initparam$upper,names(upper),varest, initparam$ns, unname(initparam$X),sensitivity,copula,MM)


######################################################
######updating starting and names  parameters 
######################################################
namespp=names(fitted_ini$par) # names of the parameters estimaded with Ind cl
aa=append(initparam$param,initparam$fixed) ## all the parameters
sel=match(namespp,names(aa));sel=sel[!is.na(sel)]  # indices to replace
aa[sel]=fitted_ini$par      #replacing
#nn=names(initparam$param)  ## selecting new starting parameters
sel=match(names(aa),initparam$namesparam);sel=sel[!is.na(sel)]  
initparam$param=aa[sel]   
######################################################



#updating with aniso parameters
update.aniso=function(param,namesparam,fixed,namesfixed,lower,upper,anisopars,estimate_aniso)
{
 un_anisopars=unlist(anisopars); namesaniso=names(un_anisopars)  
 kk=unlist(anisopars)*estimate_aniso
 ll=c(0,1)
 uu=c(pi,9999);
 lwr=c(lower,ll[estimate_aniso])
 upr=c(upper,uu[estimate_aniso])
 anisostart=kk[kk>0]
 anisofixed=kk[kk==0]
 param=c(param,anisostart)
 fixed=c(fixed,anisofixed)
 namesparam=names(param);namesfixed=names(fixed)
 if(sum(!is.na(fixed[namesaniso]))){ # updating fixed values
  if(!estimate_aniso[2]& estimate_aniso[1]) fixed["ratio"]=un_anisopars['ratio']
  if(!estimate_aniso[1]& estimate_aniso[2]) fixed["angle"]=un_anisopars['angle']
  if(!estimate_aniso[1]&!estimate_aniso[2]) {fixed["angle"]=un_anisopars['angle'];fixed["ratio"]=un_anisopars['ratio']}
    }
a=list(param=param,fixed=fixed,namesparam=namesparam,namesfixed=namesfixed,lower=lwr,upper=upr)
return(a)
}



aniso=FALSE
if(!is.null(anisopars)) {
                 aniso=TRUE;namesaniso=c("angle","ratio")
                 qq=update.aniso(initparam$param,initparam$namesparam,initparam$fixed,initparam$namesfixed,initparam$lower,initparam$upper,
                 anisopars,est.aniso)
                  initparam$param=qq$param ; initparam$fixed=qq$fixed
                  initparam$namesparam=qq$namesparam; initparam$namesfixed=qq$namesfixed
                  initparam$lower=qq$lower; initparam$upper=qq$upper
                       }
   # Full likelihood:
    if(likelihood=='Full')
          # Fitting by log-likelihood maximization:
         fitted <- Lik(copula,initparam$bivariate,initparam$coordx,initparam$coordy,initparam$coordt,coordx_dyn, initparam$corrmodel,
                               unname(initparam$data),initparam$fixed,initparam$flagcorr,
                               initparam$flagnuis,grid,initparam$lower,method,initparam$model,initparam$namescorr,
                               initparam$namesnuis,initparam$namesparam,initparam$numcoord,initparam$numpairs,
                               initparam$numparamcorr,initparam$numtime,optimizer,onlyvar,parallel,
                               initparam$param,initparam$radius,initparam$setup,initparam$spacetime,sparse,varest,taper,initparam$type,
                               initparam$upper,initparam$ns,unname(initparam$X),initparam$neighb,MM,aniso)

    # Composite likelihood:
    if((likelihood=='Marginal' || likelihood=='Conditional' || likelihood=='Marginal_2')&&type=="Pairwise"){


    if(!memdist)
          fitted <- CompLik(copula,initparam$bivariate,initparam$coordx,initparam$coordy,initparam$coordt,coordx_dyn,initparam$corrmodel,unname(initparam$data), #6
                                   initparam$distance,initparam$flagcorr,initparam$flagnuis,initparam$fixed,GPU,grid, #12
                                   initparam$likelihood,local, initparam$lower,initparam$model,initparam$n,#17
                                   initparam$namescorr,initparam$namesnuis,#19
                                   initparam$namesparam,initparam$numparam,initparam$numparamcorr,optimizer,onlyvar,parallel,
                                   initparam$param,initparam$spacetime,initparam$type,#27
                                   initparam$upper,varest,initparam$vartype,initparam$weighted,initparam$winconst,initparam$winstp,#33
                                   initparam$winconst_t,initparam$winstp_t,initparam$ns,
                                   unname(initparam$X),sensitivity,MM,aniso)
    if(memdist)
        fitted <- CompLik2(copula,initparam$bivariate,initparam$coordx,initparam$coordy,initparam$coordt,
                                   coordx_dyn,initparam$corrmodel,unname(initparam$data), #6
                                   initparam$distance,initparam$flagcorr,initparam$flagnuis,initparam$fixed,GPU,grid, #12
                                   initparam$likelihood,local, initparam$lower,initparam$model,initparam$n,#17
                                   initparam$namescorr,initparam$namesnuis,#19
                                   initparam$namesparam,initparam$numparam,initparam$numparamcorr,optimizer,onlyvar,parallel,
                                   initparam$param,initparam$spacetime,initparam$type,#27
                                   initparam$upper,varest,initparam$vartype,initparam$weighted,initparam$winconst,initparam$winstp,#33
                                   initparam$winconst_t,initparam$winstp_t,initparam$ns,
                                   unname(initparam$X),sensitivity,initparam$colidx,initparam$rowidx,initparam$neighb,MM,aniso)
      }



 


     ##misspecified models
    missp=FALSE 
    if(model=="Gaussian_misp_Tukeygh"){model="Tukeygh";missp=TRUE}
    if(model=="Gaussian_misp_Poisson"){model="Poisson";missp=TRUE}
    if(model=="Gaussian_misp_Binomial"){model="Binomial";missp=TRUE}
    if(model=="Gaussian_misp_PoissonGamma"){model="PoissonGamma";missp=TRUE}
    if(model=="Gaussian_misp_PoissonZIP"){model="PoissonZIP";missp=TRUE}
    if(model=="Gaussian_misp_StudentT"){model="StudentT";missp=TRUE}
    if(model=="Gaussian_misp_SkewStudentT"){model="SkewStudentT";missp=TRUE}
    ##################
    numtime=1
    if(initparam$spacetime) numtime=length(coordt)
    if(initparam$bivariate) numtime=2
    dimat <- initparam$numcoord*numtime#
    if(is.null(dim(initparam$X)))  initparam$X=as.matrix(rep(1,dimat))
    # Delete the global variables:

    #if(is.null(neighb)) .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
    #if(is.numeric(neighb))    
    if( !(likelihood=='Marginal'&&type=="Independence"))
    {             
     if(memdist) .C('DeleteGlobalVar2', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
     else        .C('DeleteGlobalVar' , PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
    }
    #if(is.null(neighb)&is.numeric(maxdist)) .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)

if(is.null(neighb)&is.numeric(maxdist))  fitted$value=2*fitted$value ##!!

ff=as.list(initparam$fixed)
if(!is.null(MM)) ff$mean=MM

if(length(initparam$param)==1) optimizer="optimize"

if(aniso) anisopars=as.list(c(fitted$par,ff)[namesaniso])
    ### Set the output object:
    GeoFit <- list(      anisopars=anisopars,
                         bivariate=initparam$bivariate,
                         claic = fitted$claic,
                         clbic = fitted$clbic,
                         coordx = initparam$coordx,
                         coordy = initparam$coordy,
                         coordt = initparam$coordt,
                         coordx_dyn=coordx_dyn,
                         convergence = fitted$convergence,
                         copula=copula,
                         corrmodel = corrmodel,
                         data = initparam$data,
                         distance = distance,
                         est.aniso=est.aniso,
                         fixed = ff,
                         GPU=GPU,
                         grid = grid,
                         iterations = fitted$counts,
                         local=local,
                         likelihood = likelihood,
                         logCompLik = fitted$value,
                         lower=lower,
                         message = fitted$message,
                         model = model,
                         n=initparam$n,
                         ns=initparam$ns,
                         numbetas=initparam$num_betas,
                         numcoord=initparam$numcoord,
                         numtime=initparam$numtime,
                         optimizer=optimizer,
                         param = as.list(fitted$par),
                         nozero = initparam$setup$nozero,
                         score = fitted$score,
                         maxdist =maxdist,
                         maxtime = maxtime,
                         neighb=initparam$neighb,
                         numpairs=initparam$numpairs,
                         missp=missp,
                         radius = radius,
                         spacetime = initparam$spacetime,
                         stderr = fitted$stderr,
                         sensmat = fitted$sensmat,
                         upper=upper,
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
  missp=""
  if(x$missp) missp="misspecified"
  if(x$model=='Gaussian'||x$model=='Gauss'){ process <- 'Gaussian';model <- 'Gaussian'}
  if(x$model=='Gamma') { process <- 'Gamma'; model <- 'Gamma'}
  if(x$model=='TwoPieceBimodal') { process <- 'TwoPieceBimodal'; model <- 'TwoPieceBimodal'}
  if(x$model=='LogLogistic') { process <- 'LogLogistic'; model <- 'LogLogistic'}
  if(x$model=='Gaussian_misp_Poisson') { process <- 'Poisson'; model <- 'Misspecified Gaussian Poisson '}
  if(x$model=='Gaussian_misp_Binomial') { process <- 'Binomial'; model <- 'Misspecified Gaussian Binomial '}
  if(x$model=='Gaussian_misp_PoissonZIP') { process <- 'PoissonZIP'; model <- 'Misspecified Gaussian Poisson Inflated'}
  if(x$model=='Poisson') { process <- 'Poisson'; model <- 'Poisson'}
  if(x$model=='PoissonGamma') { process <- 'PoissonGamma'; model <- 'PoissonGamma'}
  if(x$model=='Gaussian_misp_PoissonGamma') { process <- 'PoissonGamma'; model <- 'Misspecified Gaussian PoissonGamma'}
  if(x$model=='PoissonZIP') { process <- 'PoissonZIP'; model <- 'PoissonZIP'}
  if(x$model=='Beta2') { process <- 'Beta2'; model <- 'Beta2'}
  if(x$model=='Gaussian_misp_StudentT') { process <- 'StudentT'; model <- 'Misspecified Gaussian  StudentT '}
  if(x$model=='StudentT'){ process <- 'StudentT';model <- 'StudentT'}
  if(x$model=='Gaussian_misp_Tukeygh') { process <- 'Tukeygh'; model <- 'Misspecified Gaussian Tukeygh '}
  if(x$model=='Tukeygh') { process <- 'Tukeygh'; model <- 'Tukeygh '}
  if(x$model=='Gaussian_misp_SkewStudentT') { process <- 'SkewStudentT'; model <- 'Misspecified Gaussian   SkewStudentT '}
  if(x$model=='SkewStudentT') { process <- 'SkewStudentT'; model <- 'SkewStudentT'}
  if(x$model=='Logistic') { process <- 'Logistic'; model <- 'Logistic'}
  if(x$model=='Tukeyh') { process <- 'Tukeyh'; model <- 'Tukeyh'}
  if(x$model=='Tukeyh2') { process <- 'Tukeyh2'; model <- 'Tukeyh2'}
  if(x$model=='Gamma2'){ process <- 'Gamma2'; model <- 'Gamma2'}
  if(x$model=='LogGauss'||x$model=='LogGaussian'){ process <- 'Log Gaussian'; model <- 'LogGaussian'}
  if(x$model=='SkewGauss'||x$model=='SkewGaussian'){ process <- 'Skew Gaussian';model <- 'SkewGaussian'}
  if(x$model=='TwoPieceStudentT'){ process <- 'TwoPiece StudentT';model <- 'TwoPieceStudentT'}
  if(x$model=='TwoPieceTukeyh'){ process <- 'TwoPiece Tukeyh';model <- 'TwoPieceTukeyh'}
  if(x$model=='TwoPieceGaussian'||x$model=='TwoPieceGauss'){ process <- 'TwoPiece Gaussian';model <- 'TwoPieceGaussian'}
  if(x$model=='SinhAsinh'){ process <- 'SinhAsinh'; model <- 'SinhAsinh'}    
  if(x$model=='Wrapped'){ process <- 'Wrapped'; model <- 'Wrapped'}
  if(x$model=='Weibull'){ process <- 'Weibull'; model <- 'Weibull'}
  if(x$model=='Binomial'){ process <- 'Binomial';model <- 'Binomial'}
  if(x$model=='BinomialLogistic'){ process <- 'BinomialLogistic';model <- 'BinomialLogistic'}
  if(x$model=='Kumaraswamy'){ process <- 'Kumaraswamy';model <- 'Kumaraswamy'}
  if(x$model=='Kumaraswamy2'){ process <- 'Kumaraswamy2';model <- 'Kumaraswamy2'}
  if(x$model=='Beta'){ process <- 'Beta';model <- 'Beta'}
  if(x$model=='Binomial_TwoPieceGaussian'||x$model=='Binomial_TwoPieceGauss'){ process <- 'Binomial TwoPiece Gaussian';model <- 'Binomial_TwoPieceGauss'}
  if(x$model=='BinomialNeg_TwoPieceGaussian'||x$model=='BinomialNeg_TwoPieceGauss'){ process <- 'Negative Binomial TwoPiece Gaussian';model <- 'BinomialNeg_TwoPieceGauss'}
  if(x$model=='Binomial2'){ process <- 'Binomial';model <- 'Binomial2'}     
  if(x$model=='BinomialNeg'){ process <- 'BinomialNeg'; model <- 'BinomialNeg'}
  if(x$model=='BinomialNegLogistic'){ process <- 'BinomialNegLogistic'; model <- 'BinomialNegLogistic'}
  if(x$model=='BinomialNegZINB'){ process <- 'BinomialNegZINB'; model <- 'BinomialNegZINB'}
  if(x$model=='Geom'||x$model=='Geometric'){ process <- 'Geometric';model <- 'Geometric'}
  if(x$model=='PoisBin'){ process <- 'Poisson Binomial';model <- 'PoisBin'}
  if(x$model=='PoisBinNeg'){ process <- 'Poisson NegBinomial';model <- 'PoisBinNeg'}    
  if(x$bivariate){ biv <- 'bivariate';x$numtime=1}
  else { biv <- 'univariate'}                       
    cat('\n##################################################################')
    cat('\nMaximum', missp, method, 'Fitting of', process, 'Random Fields\n')
    if(!is.null(x$copula)) {cat('\nCopula:', x$copula,'\n')}
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
    print.default(unlist(x$param), digits = digits, print.gap = 2,
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

