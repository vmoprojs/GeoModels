####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate.
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Institutions: 
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: GeoCovariogram.r
### Description:
### This file contains a set of procedures
### to compute and plot the estimated covariance
### function and the variogram after fitting a
### random field by composite-likelihood.
### Last change: 28/03/2018.
####################################################
   
 
### Procedures are in alphabetical order.

### Compute and plot the (estimated) covariance function and the variogram
### from a fitted model obtain from the GeoFit or the WLeastSquare procedure
GeoCovariogram <- function(fitted, distance="Eucl", answer.cov=FALSE, answer.vario=FALSE,
                         answer.range=FALSE, fix.lags=NULL, fix.lagt=NULL,
                        show.cov=FALSE, show.vario=FALSE,  show.range=FALSE,
                        add.cov=FALSE, add.vario=FALSE, pract.range=95,
                        vario=NULL, ...)
  {
    result <- NULL
    CheckDistance<- function(distance)
    {
        CheckDistance <- NULL
        CheckDistance <- switch(distance,
                                eucl=0,
                                Eucl=0,
                                chor=1,  
                                Chor=1,
                                geod=2,
                                Geod=2,
                                proj=3,
                                Proj=3)
        return(CheckDistance)
    }
    # define the bivariate Gaussian distribution:
    vpbnorm <- function(corrmodel,lags,lagt,nuisance,numlags,numlagt,param,threshold)
    {
        p=.C("vpbnorm",  as.integer(corrmodel), as.double(lags), as.double(lagt),
           as.integer(numlags), as.integer(numlagt), as.double(nuisance),
           as.double(param), rho=double(numlags*numlagt), as.double(threshold), PACKAGE='GeoModels',
           NAOK=TRUE)
        return(p$rho)
    }
               
    CorrelationFct <- function(bivariate,corrmodel, lags, lagt, numlags, numlagt, mu,model, nuisance,param)
    {
       if(!bivariate) {
                             p=.C('VectCorrelation', corr=double(numlags*numlagt), as.integer(corrmodel), as.double(lags),
                             as.integer(numlags), as.integer(numlagt), as.double(mu),as.integer(model),as.double(nuisance),as.double(param),
                             as.double(lagt), PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
                             cc=p$corr
                    }
        else    {
                             p=.C('VectCorrelation_biv', corr=double(numlags*4),vario=double(numlags*4), as.integer(corrmodel), as.double(lags),
                             as.integer(numlags), as.integer(numlagt),  as.double(mu),as.integer(model),as.double(nuisance), as.double(param),
                             as.double(lagt), PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
                             cc=c(p$corr,p$vario)   

                    }
        return(cc)
    }
    if(show.range)  {
    # Pratical range in the Gaussian case:
    PracRangeNorm <- function(corrmodel, lags, lagt, nuisance, numlags, numlagt, param, pract.range)
    { 
        return(nuisance["nugget"]+nuisance["sill"]*CorrelationFct(bivariate,corrmodel, lags, lagt, numlags, numlagt, mu,CkModel(fitted$model),nuisance, param)
            -(nuisance["nugget"]+nuisance["sill"]*(1-pract.range)))
    
    }
}
################    
### starting ###
################
    dyn=FALSE
    isvario <- !is.null(vario) # is empirical variogram is passed?
    bivariate <- fitted$bivariate
    if(bivariate) fitted$numtime=1
    ispatim <- fitted$spacetime
    dyn<- is.list(fitted$coordx_dyn)
    
par(mfrow=c(1,1))
if(!ispatim && !bivariate){ if( (show.cov && show.vario) || (show.cov)) par(mfrow=c(1,2))}
if(show.vario && ispatim) par(mfrow=c(1,2))
if(show.vario && ispatim && !is.null(fix.lags) && !is.null(fix.lagt)) par(mfrow=c(2,2))
if(show.cov && ispatim && !is.null(fix.lags) && !is.null(fix.lagt)) par(mfrow=c(2,2))
if(show.cov && ispatim && is.null(fix.lags) && is.null(fix.lagt)) par(mfrow=c(1,1))
if(show.vario && bivariate) {par(mfrow=c(2,2))}
if(bivariate&&dyn) par(mfrow=c(1,2))


    
    # START ---- check input --- #
    if(!class(fitted)=='GeoFit' & !class(fitted)=='GeoWLS')
        stop("Enter an object obtained GeoFit or WLeastSquare")

    if(isvario & !class(vario)=='GeoVariogram')
        stop("Enter an object obtained from the function GeoVariogram")

    if(!is.numeric(pract.range) & answer.range)
        stop("Enter a number for the parameter % of sill")
    else{
        if(pract.range < 0 || pract.range > 100)
            stop("Entered an incorrect value for the % of sill")
        else
          pract.range <- pract.range / 100}
    # Set the type of process:
    model <- CkModel(fitted$model)
    gaussian <- model==1
    skewgausssian<- model==10
    gamma<- model==21
    studentT<- model==12
    weibull<- model==26
    twopieceT<- model==27
    twopieceGauss<- model==29
    loggauss<- model==22
    binary <- model==2
    binomial <- model==11
    binomial2 <- model==19
    geom <- model==14
    studentT <- model ==12
    loglogistic <- model==24
    zero <- 0;slow=1e-2
    # lags associated to empirical variogram estimation
    if(isvario){
    lags <- c(0,vario$centers);numlags <- length(lags)
    if(ispatim) lagt <-c(0,vario$bint) else lagt=0
    numlagt <- length(lagt)
    }

    # check the fixed lags:
    if(!is.null(fix.lags)){
        if(!is.numeric(fix.lags) || fix.lags<0 || fix.lags>numlags){
            cat("invalid spatial fixed lag\n")
            return(result)}}
    if(!is.null(fix.lagt))
        if(!is.numeric(fix.lagt) || fix.lagt<0 || fix.lagt>numlagt){
            cat("invalid temporal fixed lag\n")
            return(result)}
    # END ---- check input --- #
    # set range (for practical range)
    if(fitted$corrmodel=="matern"){
        lower <- 1e-10
        upper <- 1e20}
    else{
        lower <- 0  
        upper <- 1e100}
    num_betas=fitted$numbetas
    mu=0;nuisance=0
    mm=0
    ## selecting nuisance mean annd corr parameters
      if(!bivariate){
       param <- c(fitted$fixed,fitted$param)[CorrelationPar(CkCorrModel(fitted$corrmodel))]
        nuisance <- c(fitted$fixed,fitted$param)[NuisParam(fitted$model,FALSE,num_betas)]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=nuisance[sel]
        nuisance=nuisance[!sel]
        }
      if(bivariate){
        param <- c(fitted$fixed,fitted$param)[CorrelationPar(CkCorrModel(fitted$corrmodel))]
        nuisance <- c(fitted$fixed,fitted$param)[NuisParam(fitted$model,FALSE,num_betas)]
    }

     #############################################
     # computing the spatio-temporal distances where to compute the fitted model
     #############################################
    if(isvario) {
    lags_m <- seq(slow,max(vario$centers),length.out =150)
    if (ispatim) lagt_m <-seq(slow,max(vario$bint),length.out =150)
    else         lagt_m<-0
        }
    else{
        mmm <- double(2)
        type_dist <- CheckDistance(distance)
        p=.C("Maxima_Minima_dist",mmm=as.double(mmm),as.double(fitted$coordx),as.double(fitted$coordy)
        ,as.integer(fitted$numcoord),as.integer(type_dist),as.double(fitted$radius),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
       # if(type_dist==0) mmx=max(c(dist(cbind(fitted$coordx,fitted$coordy))))
        lags_m <- seq(slow,p$mmm[2],length=150)
        if (ispatim) {
            tt <- double(2)
            #p=.C("Maxima_Minima_time",tt=as.double(tt),as.double(fitted$coordt),as.integer(fitted$numtime),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
            lagt_m <- seq(slow,max(c(dist(fitted$coordt))) ,length=150)
        }
        else lagt_m<-0
    }
    numlags_m <- length(lags_m)
    numlagt_m <- length(lagt_m)
if(!bivariate) {
    if(ncol(fitted$X)==1)  X=as.matrix(rep(1,numlags_m*numlagt_m))
    mu=X%*%as.numeric(mm)}
  else{}  
   ###############
    corrmodel <- CkCorrModel(fitted$corrmodel)

     nui=nuisance
     nui['sill']=1-nui['nugget']
    #nui['nugget']=1-nui['sill']
  #   nui=nuisance
  #  if(gamma||weibull||studentT||loglogistic) {nui['sill']=1;nui['nugget']=1-nui['sill']}
    correlation <- CorrelationFct(bivariate,corrmodel, lags_m, lagt_m, numlags_m, numlagt_m,mu,
                                     CkModel(fitted$model), nui,param)

    # Gaussian random field:
    if(gaussian){
    
        if(bivariate){  covariance11 <- correlation[(0*length(lags_m)+1):(1*length(lags_m))]
                        covariance12 <- correlation[(1*length(lags_m)+1):(2*length(lags_m))]
                        covariance22 <- correlation[(3*length(lags_m)+1):(4*length(lags_m))]
                        variogram11  <- correlation[(4*length(lags_m)+1):(5*length(lags_m))]
                        variogram12  <- correlation[(5*length(lags_m)+1):(6*length(lags_m))]
                        variogram22  <- correlation[(7*length(lags_m)+1):(8*length(lags_m))]
                           }
        else { 
        covariance <- nuisance["nugget"]+nuisance["sill"]*correlation
        variogram <- nuisance["nugget"]+nuisance["sill"]*(1-correlation)
        }
    }
##########################################
 if(skewgausssian) {    
            if(bivariate) {}
              else {
              vv=as.numeric(nuisance['sill']);sk=nuisance['skew'];sk2=sk^2;
              vs=(vv+sk2*(1-2/pi))
              correlation=(1-nuisance['nugget'] )*correlation
              corr2=correlation^2;  
              cc=((2*sk2/pi)*(sqrt(1-corr2) + correlation*asin(correlation)-1) + correlation*vv)/(vv+sk2*(1-2/pi))
              covariance=vs*cc;variogram=vs*(1-cc) }
                   }
##########################################
   if(twopieceT)        { if(bivariate) {}
                        else {

                              nu=1/as.numeric(nuisance['df']); sk=nuisance['skew'];sill=nuisance['sill']
                              vs= sill* nu*(3*sk^2+1)/(nu-2)-(4*sk^2*nu*gamma((nu-1)/2)^2)/(pi*gamma(nu/2)^2)
                              ##
                              corr2=correlation^2;sk2=sk^2
                              a1=Re(hypergeo::hypergeo(0.5, 0.5, nu/2,corr2))
                              a2=correlation*asin(correlation) + (1-corr2)^(0.5)
                              ll=qnorm((1-sk)/2)
                              p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                              a3=3*sk2 + 2*sk + 4*p11 - 1
                              KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
                              cc= KK*(a1*a2*a3-4*sk2);
                              ##
                              covariance=vs*cc;variogram=vs*(1-cc)  }
                  }  
##########################################
   if(twopieceGauss)        { if(bivariate) {}
                        else {
                              sk=nuisance['skew'];sill=nuisance['sill'];sk2=sk^2
                              vs= sill*(1+3*sk2)-8*sk2/pi
                              corr2=sqrt(1-correlation^2)
                              ll=qnorm((1-sk)/2)
                              p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                              KK=3*sk2+2*sk+ 4*p11 - 1
                              cc=(2*((corr2 + correlation*atan(correlation/corr2))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   ) 
                              covariance=vs*cc;variogram=vs*(1-cc)  }
                  } 


##########################################
   if(studentT)        { if(bivariate) {}
                        else {
                              nu=1/as.numeric(nuisance['df']);sill=as.numeric(nuisance['sill'])
                              vs=sill*(nu)/(nu-2)
                              
                              cc=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,correlation^2))*correlation)/(2*gamma(nu/2)^2)
                              covariance=vs*cc;variogram=vs*(1-cc)  }
                  }     
##########################################
 if(gamma)        { if(bivariate) {}
                        else {
                              vs=2*exp(mm)^2/nuisance['shape']
                              cc=((1-nuisance['nugget'] )*correlation)^2
                              covariance=vs*cc;variogram=vs*(1-cc)  }
                  }
##########################################
 if(weibull)        { if(bivariate) {} 
                        else {
                        vs=exp(mm)^2*(gamma(1+2/nuisance["shape"])/gamma(1+1/nuisance["shape"])^2-1)
                        auxcorr= (gamma(1+1/nuisance['shape']))^2/((gamma(1+2/nuisance['shape']))-(gamma(1+1/nuisance['shape']))^2)
                        cc=auxcorr*(Re(hypergeo::hypergeo(-1/nuisance['shape'], -1/nuisance['shape'], 1,((1-nuisance['nugget'] )*correlation)^2)) -1)
                        covariance=vs*cc;variogram=vs*(1-cc)  }
                    }
##########################################
  if(loglogistic)    { if(bivariate) {}  
                      else { 
                     sh=nuisance["shape"]
                     vs=exp(mm)^2*(2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)
                     cc=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                                    (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,(1-nuisance['nugget'] )*correlation^2))*
                                     Re(hypergeo::hypergeo( 1/sh,  1/sh, 1,(1-nuisance['nugget'] )*correlation^2)) -1)
                      covariance=vs*cc;variogram=vs*(1-cc)   }
                    }
##########################################
  if(loggauss)    { if(bivariate) {}  
                      else {    
                    vvar=as.numeric(nuisance["sill"]+nuisance["nugget"])
                    yy=nuisance["sill"]*correlation
                    cc<-(exp(nuisance["sill"]*correlation)-1)/(exp(vvar)-1)
                    vs<-(exp(vvar)-1)*exp(2*mm+vvar)
                    covariance=vs*cc;variogram=vs*(1-cc)   
                     }
                  }

   #                   
   if(binary||binomial||binomial2||geom) {
                    if(bivariate) {}
                    if(!bivariate) {          
                           pp=pnorm(mu)
                           if(binary||binomial||binomial2) vv=min(fitted$n)*pp*(1-pp)
                           if(geom)             vv=(1-pp)/pp^2;
                           covariance=vv*correlation
                           variogram=vv*(1-correlation)}
                   }
##########################################
##########################################
##########################################
      vario.main <- "Spatial semi-variogram"
      vario.ylab <- "Semi-Variogram"
        if(ispatim){
            dim(covariance) <- c(numlags_m, numlagt_m)
            dim(variogram) <- c(numlags_m, numlagt_m)
            vario.main <- "Space-time semi-variogram"
            vario.zlab <- "Semi-Variogram" }
        # compute the practical range:
        if(show.range || answer.range){
            Range <- uniroot(PracRangeNorm, c(lower, upper), corrmodel=corrmodel,
                             nuisance=nuisance, numlags=1, numlagt=1, lagt=lagt,
                             param=param, pract.range=pract.range)$root
        }
    # binary random field:
    #if(binary){
    #    covariance <- nuisance["nugget"]+nuisance["sill"]*correlation
    #    p <- pnorm(nuisance["mean"])
    #    q <- vpbnorm(corrmodel, lags_m, lagt_m, nuisance,
    #                 numlags_m, numlagt_m, param, 0)
    #    variogram <- log(q*(1-2*p+q)/(p-q)^2)
    #    vario.main <- "Spatial lorelogram"
    #    vario.ylab <- "Lorelogram"
    #    if(ispatim){
    #        dim(covariance) <- c(numlags_m, numlagt_m)
    #        dim(variogram) <- c(numlags_m, numlagt_m)
    #        vario.main <- "Space-time lorelogram"
    #        vario.zlab <- "Lorelogram"}
    #      }


    # display the covariance function
    if(show.cov){
        if(bivariate&&!dyn){
            #par(mfrow=c(2,2))
       plot(lags_m, covariance11, type='l', ylim=c(min(covariance11),
                     max(covariance11)), main="First covariance",
                     xlab="Distance", ylab="Covariance",...)
       plot(lags_m, covariance12, type='l', ylim=c(min(covariance12),
                     max(covariance12)), main="Cross covariance",
                     xlab="Distance", ylab="Covariance",...)      
       plot(lags_m, covariance12, type='l', ylim=c(min(covariance12),
                     max(covariance12)), main="Cross covariance",
                     xlab="Distance", ylab="Covariance",...)   
       plot(lags_m, covariance22, type='l', ylim=c(min(covariance22),
                     max(covariance22)), main="Second covariance",
                     xlab="Distance", ylab="Covariance",...)          
         }
          if(bivariate&&dyn){
            #par(mfrow=c(2,2))
       plot(lags_m, covariance11, type='l', ylim=c(min(covariance11),
                     max(covariance11)), main="First covariance",
                     xlab="Distance", ylab="Covariance",...)
       plot(lags_m, covariance22, type='l', ylim=c(min(covariance22),
                     max(covariance22)), main="Second covariance",
                     xlab="Distance", ylab="Covariance",...)          
         }
        if(ispatim){# spatio-temporal case:
            # build the covariance matrix:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            numplot <- 1+plags+plagt
            # par(mfrow=c(1,numplot))
            # temporal section
            par(mai=c(.2,.2,.2,.2))
            persp(lags_m, lagt_m, covariance, xlab="Distance", ylab="Time",
                  zlab="Covariance", ltheta=90,
                  shade=0.75, ticktype="detailed", phi=30,
                  theta=30,main="Fitted pace-time covariance",
                  , cex.axis=.8, cex.lab=.8,zlim=c(0,max(covariance))) #
            if(plagt){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lagt_m, covariance[fix.lags,], xlab="Time",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main="Space-time cov: temporal profile",...)}
            if(plags){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lags_m, covariance[,fix.lagt], xlab="Distance",
                     ylab="Covariance", type="l",cex.axis=.8,cex.lab=.8,
                     main="Space-time cov: spatial profile",...)}
        }
       if(!ispatim && !bivariate){# spatial case:
            if(add.cov & dev.cur()!=1){
                lines(lags_m, covariance,...)
                if(show.range) abline(v=Range)}
            else{
                plot(lags_m, covariance, type='l', ylim=c(min(covariance),
                     max(covariance)), main="Spatial covariance",
                     xlab="Distance", ylab="Covariance",...)
                if(show.range) abline(v=Range)}}
         }

    # display the variogram function
    if(show.vario){
      if(bivariate&&!dyn){
          #par(mfrow=c(2,2))
       plot(vario$centers,vario$variograms[1,], main="First semi-variogram",ylim=c(0,max(vario$variograms[1,])),
           xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram11, type='l',...)
       if(min(vario$variogramst)>0) {ll1=0;ll=max(vario$variogramst)}
       if(min(vario$variogramst)<0) {ll1=min(vario$variogramst);ll=-min(vario$variogramst)}
       plot(vario$centers,vario$variogramst, main="Cross semi-variogram",ylim=c(ll1,ll),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram12, type='l',...)
       plot(vario$centers,vario$variogramst, main="Cross semivariogram",ylim=c(ll1,ll),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram12, type='l',...)
       plot(vario$centers,vario$variograms[2,], main="Second semi-variogram",ylim=c(0,max(vario$variograms[2,])),
         xlim=c(0,max(vario$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram22, type='l',...)  }
 
   if(bivariate&&dyn){
          #par(mfrow=c(2,2))
       plot(vario$centers,vario$variograms[1,], main="First semi-variogram",ylim=c(0,max(vario$variograms[1,])),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram11, type='l',...)
       plot(vario$centers,vario$variograms[2,], main="Second semi-variogram",ylim=c(0,max(vario$variograms[2,])),
                     xlab="Distance", ylab="Semi-Variogram",...)
       lines(lags_m, variogram22, type='l',...)  }

    
        if(ispatim){# spatio-temporal case:
            plagt <- !is.null(fix.lags)
            plags <- !is.null(fix.lagt)
            nrowp <- 1
            ncolp <- 1
            if(isvario){
                ncolp <- ncolp+1
                if(plags || plagt) nrowp <- nrowp+1}
            else ncolp <- ncolp+plags+plagt
            #par(mfrow=c(nrowp,ncolp))
            sup <- 0
            tup <- 0
            if(isvario){
                nbins <- length(vario$centers)
                nbint <- length(vario$bint)
                evario <- matrix(vario$variogramst,nrow=nbins,
                                 ncol=nbint,byrow=TRUE)
                evario <- rbind(c(zero,vario$variogramt),cbind(vario$variograms,evario))
                evario.grid <- expand.grid(c(0,vario$centers),c(0,vario$bint))
                scatterplot3d::scatterplot3d(evario.grid[,1],evario.grid[,2], c(evario),
                              type="h",highlight.3d=TRUE,cex.axis=.7,cex.lab=.7,
                              main=paste("Empirical",vario.main),xlab="Distance",
                              ylab="Time",zlab=vario.zlab,mar=c(2,2,2,2),mgp=c(0,0,0))
                if(plagt) tup <- max(evario[fix.lags,],na.rm=TRUE)
                if(plags) sup <- max(evario[,fix.lagt],na.rm=TRUE)
            } 
            par(mai=c(.2,.2,.2,.2))
            persp(lags_m, lagt_m, variogram, xlab="Distance",
                  ylab="Time", zlab=vario.zlab, ltheta=90,
                  shade=0.75, ticktype="detailed", phi=30,
                  theta=30,main=vario.main, cex.axis=.8,
                   cex.lab=.8)  #zlim=c(0,max(variogram))
            vvv=nuisance["nugget"]+nuisance["sill"]
            ########
            if(gamma)    vvv=2*exp(mm["mean"])^2/nuisance["shape"]
            if(weibull)  vvv=exp(mm["mean"])^2*(gamma(1+2/nuisance["shape"])/gamma(1+1/nuisance["shape"])^2-1)
            if(loglogistic)  vvv=exp(mm["mean"])^2*
                               (2*nuisance['shape']*sin(pi/nuisance['shape'])^2/(pi*sin(2*pi/nuisance['shape']))-1)
            if(loggauss) vvv=(exp(nuisance["sill"])-1)*(exp(mm['mean']+0.5*nuisance["sill"]))^2
            if(binomial) vvv=fitted$N*pnorm(mm['mean'])*(1-pnorm(mm['mean']))
            if(geom)     vvv= (1-pnorm(mm['mean']))/pnorm(mm['mean'])^2
            if(skewgausssian) vvv=(nuisance["sill"]+nuisance["skew"])^2*(1-2/pi)
            if(studentT)      vvv=nuisance["df"]/(nuisance["df"]-2)
            ########
            if(plagt){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lagt_m, variogram[fix.lags,], xlab="Time",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(min(0,min(variogram)),max(vvv,tup)), main=paste(vario.ylab,": temporal profile",
                     sep=""),...)
                if(isvario) points(lagt, evario[fix.lags,],...)}
            if(plags){
                par(mai=c(.5,.5,.3,.3),mgp=c(1.6,.6,0))
                plot(lags_m, variogram[,fix.lagt], xlab="Distance",cex.axis=.8,cex.lab=.8,
                     ylab=vario.ylab, type="l", ylim=c(min(0,min(variogram)),max(vvv,sup)), main=paste(vario.ylab,": spatial profile",
                     sep=""),...)
                if(isvario) points(lags, evario[,fix.lagt],...)}}
        if(!ispatim && !bivariate){# spatial case:
            if(add.vario & dev.cur()!=1){
                points(vario$centers, vario$variograms,...)
                lines(lags_m, variogram,...)
                if(show.range) abline(v=Range)}
            else{
                bnds <- range(variogram)
                bnds[1] <- min(bnds[1], min(vario$variograms))
                bnds[2] <- max(bnds[2], max(vario$variograms))
                plot(lags_m, variogram, type='l',  ylim=c(bnds[1],bnds[2]),
                     main=vario.main,xlab="Distance",
                     ylab=vario.ylab,...)
                points(vario$centers, vario$variograms,...)
                if(show.range) abline(v=Range)}
        }}
        


    if(ispatim) par(mai=c(1.02 ,0.82 ,0.82 ,0.42),mgp=c(3,1,0))
    # return the estimated covariance function
    if(answer.cov) {result <- list(lags=lags_m,lagt=lagt_m, covariance=covariance)}
    # return the estimated variogram/lorelogram function
    if(answer.vario) {
        if(!is.list(result)) {if(!bivariate) if(gaussian) result <- list(lags=lags_m,lagt=lagt_m, variogram=variogram)
                              if(bivariate)  if(gaussian) result <- list(lags=lags_m,lagt=lagt_m, variogram11=variogram11,
                                                                        variogram12=variogram12,variogram22=variogram22)          
                              }
        else {
            if(!bivariate){if(gaussian) result$variogram <- variogram}
            if(bivariate){
                if(gaussian) {result$variogram11 <- variogram11;result$variogram12 <- variogram12;result$variogram22 <- variogram22}
                }}}
    if(!is.null(result))
    #par(mfrow=c(1,1))
    return(result)
  }
