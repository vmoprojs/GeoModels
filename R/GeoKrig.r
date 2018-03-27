####################################################
### Authors: Moreno Bevilacqua Víctor Morales Oñate.
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: GeoKrig.r  
### Description:  
### This file contains a set of procedures
### for computing simple (tapered) and ordinary kriging
### predictor  at an unknown space (time) locations.
### Last change: 28/02/2017.
#################################################### 

GeoKrig<- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, maxdist=NULL,
               maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE, lin_opt=TRUE, param, radius=6378.388, sparse=FALSE, 
               taper=NULL, tapsep=NULL, time=NULL, type="Standard",type_mse=NULL, type_krig="Simple",weigthed=TRUE, 
               which=1, X=NULL,Xloc=NULL)

{ 
  ##################################
  getInv<-function(covmatrix){
     if(!covmatrix$sparse){
              decompvarcov <- MatDecomp(covmatrix$covmatrix,method)
              if(is.logical(decompvarcov)){print(" Covariance matrix is not positive definite");stop()}      
               invcov <- MatInv(decompvarcov,method)}
        else{  decompvarcov  <- try(spam::as.spam(covmatrix$covmatrix),silent=TRUE)
               if(class(decompvarcov)=="try-error") {print(" Covariance matrix is not positive definite");stop()}
               invcov<-spam::solve.spam(decompvarcov)
        }
    }
######################################
###################################### 
######################################  
    if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
    if(!is.matrix(loc))   stop("loc parameter must be a matrix")
    if(!(ncol(loc)==2))   stop("loc parameter must be a matrix  N X 2")
    if(!is.null(Xloc)) Xloc=as.matrix(Xloc)
    if(is.matrix(X) &&is.null(Xloc))  stop("Covariates for locations to predict are missing ")
    if(is.null(X) &&is.matrix(Xloc))  stop("Covariates  are missing ")
###################################### 
###################################### 
###################################### 
    #### number of points to predict
    numloc <- nrow(loc); tloc <- length(time);
    if(!tloc)  tloc <- 1
    locx <- loc[,1];locy <- loc[,2]
    #######################################################
    ############ standard (tapered) kriging ###############
    #######################################################
    if(type=="Standard"||type=="standard"||type=="Tapering"||type=="tapering")  {
    #################################################
    ##### computing covariance  matrix ##############
    #################################################
    logGausstemp=FALSE
    if(model=="LogGaussian") {model="Gaussian"    # we need a "Gaussian" covariance matrix
                             logGausstemp=TRUE}

    if(model=="Weibull"||model=="Gamma")          # we need a x covariane matrix with with mean=0   x=gamma,weibull
     {
     paramtemp=param 
     sel=substr(names(param),1,4)=="mean";  ## selecting mean values 
     meantemp=names(param[sel])             ## saving  mean values
     param=param[!sel];param$mean=0;        ## mean must be =0 when calling covariance matrix
     Xtemp=X;X=NULL                         ## saving X and setting X=NULL
     }
    covmatrix <- GeoCovmatrix(coordx, coordy, coordt, coordx_dyn, corrmodel, distance, grid, maxdist, maxtime, model, n, param, 
      radius, sparse, taper, tapsep, type, X) 


    ###########
    spacetime_dyn=FALSE
    if(!is.null(covmatrix$coordx_dyn)) spacetime_dyn=TRUE
    ##############
    if(!spacetime_dyn) dimat=covmatrix$numcoord*covmatrix$numtime
    if(spacetime_dyn)  dimat =sum(covmatrix$ns)
    dimat2=numloc*tloc
    ###############
    ###############
    if(model=="Weibull"||model=="Gamma"){
          if(is.null(Xtemp)) X=matrix(rep(1,dimat))
          else               X=Xtemp    
          param=paramtemp
          covmatrix$namesnuis=unique(c(meantemp,covmatrix$namesnuis))
    }
    else 
        X=covmatrix$X 
    ###############
    ###############
    num_betas=ncol(X)
    
    if(is.null(Xloc))   Xloc=as.matrix(rep(1,dimat2))  

    nuisance <- param[covmatrix$namesnuis]
    sel=substr(names(nuisance),1,4)=="mean"
    betas=as.numeric(nuisance[sel])   ## mean paramteres
    other_nuis=as.numeric(nuisance[!sel]) 
    ################################################
    ################################################
    if(type=="Tapering"||type=="tapering") {
      covmatrix_true <-  GeoCovmatrix(coordx, coordy, coordt, coordx_dyn, corrmodel, distance, grid, maxdist, maxtime, model, n, param, 
      radius, sparse, NULL, NULL, "Standard",X) }
    ############
    tapmod<-NULL                                                
    cmodel<-corrmodel
    cdistance<-distance
    corrmodel <- CkCorrModel(covmatrix$corrmodel)
    distance <- CheckDistance(covmatrix$distance)
    corrparam <- unlist(covmatrix$param[covmatrix$namescorr])# selecting the correlation parametrs
    bivariate <- covmatrix$bivariate;   
    if(bivariate) tloc=1
    spacetime <- covmatrix$spacetime; 
    if(bivariate) if(!(which==1 || which==2) ) stop("which  parameter must be 1 or 2")
    pred <- NULL
    varpred<-varpred2<-vv<-vv2<-NULL
    k <- 0 
    if(grid) ccc=expand.grid(covmatrix$coordx,covmatrix$coordy)
    else     ccc=cbind(covmatrix$coordx,covmatrix$coordy)
    ###############################################################
    if((spacetime||bivariate)&&spacetime_dyn) dataT=t(unlist(data)) 
    else dataT=t(data)
########################################################################################
########################################################################################
########################################################################################

if(covmatrix$model %in% c(1,10,21,12,26,24,27))
{  
## gaussian=1 
## skew gaussian=10   
## student =12
## gamma=21 
## weibull=26
## loglogistic=24
## loggaussian=22 ojo
## twopieceStudentT=27
    ################################
    ## standard kriging  ##############
    ################################      
       mu=X%*%betas
       muloc=Xloc%*%betas

    if((type=="Standard"||type=="standard")) {
         ## Computing CORRELATIONS between the locations to predict and the locations observed
        cc=.C('Corr_c',corri=double(dimat*dimat2), as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
        as.integer(corrmodel),as.integer(FALSE),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
        as.integer(numloc),as.integer(tloc),as.integer(covmatrix$ns),as.integer(covmatrix$numtime),as.double(corrparam),
        as.integer(covmatrix$spacetime),
        as.integer(covmatrix$bivariate),as.double(time),as.integer(distance),as.integer(which-1),
        as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
        if(covmatrix$model==1)    #gaussian
                        {
                         vv=covmatrix$param['sill'];
                         corri=cc$corri  }
        if(covmatrix$model==10) {    #skew gaussian
                        corr2=cc$corri^2   
                        sk=covmatrix$param['skew'];sk2=sk^2
                        vv=covmatrix$param['sill']
                        corri=((2*sk2/pi)*(sqrt(1-corr2) + cc$corri*asin(cc$corri)-1) + cc$corri*vv)/(vv+sk2*(1-2/pi));
                        }        
        if(covmatrix$model==21)  # gamma
                        corri=((1-covmatrix$param["nugget"])*cc$corri)^2
        if(covmatrix$model==12) # student T
                         {
                        vv=covmatrix$param['sill']  
                        nu=1/covmatrix$param['df']
                        corri=((nu-2)*gamma((nu-1)/2)^2*gsl::hyperg_2F1(0.5,0.5 ,nu/2 ,cc$corri^2)*cc$corri)/(2*gamma(nu/2)^2)
                      }
        if(covmatrix$model==26) {  # weibull 
                        sh=covmatrix$param['shape']
                        bcorr=    (gamma(1+1/sh))^2/((gamma(1+2/sh))-(gamma(1+1/sh))^2)
                        cc1=(1-covmatrix$param["nugget"])*cc$corri
                        corri=bcorr*((1-cc1^2)^(1+2/sh)*gsl::hyperg_2F1(1+1/sh, 1+1/sh, 1,cc1^2) -1)              
         }
         if(covmatrix$model==27) {  # two piece StudenT
                        nu=1/covmatrix$param['df'];sk=covmatrix$param['skew']
                        vv=covmatrix$param['sill']
                        corr2=cc$corri^2;sk2=sk^2
                        a1=gsl::hyperg_2F1(0.5, 0.5, nu/2,corr2)
                        a2=cc$corri*asin(cc$corri) + (1-corr2)^(0.5)
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = cc$corri, recycle = TRUE)
                        a3=3*sk2 + 2*sk + 4*p11 - 1
                        KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma((nu)/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
                        corri= KK*(a1*a2*a3-4*sk2);                     
                      }
        if(mse){
             if(bivariate)  { 
                          if(which==1)    vvar=covmatrix$param["sill_1"] 
                          if(which==2)   vvar=covmatrix$param["sill_2"]}
             else    {if(covmatrix$model==1)   vvar=vv     #gaussian
                      if(covmatrix$model==10)  vvar= (vv+sk^2*(1-2/pi))   ## skewgaus
                      if(covmatrix$model==12)  vvar=nu/(nu-2)              ## studentT
                      if(covmatrix$model==27)  vvar=nu*(3*sk^2+1)/(nu-2)-
                                                     (4*sk^2*nu*gamma((nu-1)/2)^2)/(pi*gamma(nu/2)^2) # two pieceT
                      #if(covmatrix$model==21)  vvar=2*exp(muloc)/covmatrix$param['shape']
                      #if(covmatrix$model==22) {kk=exp(2*(muloc)+covmatrix$param['sill']);vvar=kk*(exp(covmatrix$param['sill'])-1)}
                      }
           }
########################################################################################
########################################################################################
########################################################################################
 #### inverse of var covar  ##################################   

        invcov <- getInv(covmatrix)  ### invserse of cov matrix
        #############################################################
        ##### multiplying the correlations for the variance
        cc <- t(matrix(corri,nrow=dimat,ncol=dimat2)) 
        if(!bivariate){
          if(covmatrix$model==1)  
                          cc=cc* vv
          if(covmatrix$model==10)  #skewgaussian
                          cc=cc* (vv+(sk)^2*(1-2/pi))
          if(covmatrix$model==12)  #studentT
                          cc=cc* vv *(nu/(nu-2))
          if(covmatrix$model==27)   ##two piece studentT
                          cc=cc* vv *(nu*(3*sk^2+1)/(nu-2)-(4*sk^2*nu*gamma((nu-1)/2)^2)/(pi*gamma(nu/2)^2) )       
          if(covmatrix$model==21)  { # gamma
                                    emuloc=exp(muloc) 
                                    emu=exp(mu)
                                   # V0=NULL;for(i in 1:dimat) V0=cbind(V0,emuloc^2) 
                                   # V0=2*V0/covmatrix$param['shape']
                                   # V1=2*diag(c(emuloc^2))/covmatrix$param['shape']
                                   # cc=V0*cc
                                    cc=2*cc/covmatrix$param['shape']
                                    }
           if(covmatrix$model==26)  { # weibull
                                    emuloc=exp(muloc)
                                    emu=exp(mu)
                                    #V0=NULL;for(i in 1:dimat) V0=cbind(V0,emuloc^2) 
                                    #V0=V0*(gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1) 
                                    #V1=diag(c( emuloc^2))*(gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1) 
                                    cc=cc*(gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1)
                                    #cc=V0%*%cc
                                  }
          if(covmatrix$model==24)  {}

         }
         else{}
##################################################################
#################kriging weights##################################
##################################################################
        krig_weights <- cc%*%invcov
##################################################################
################# simple kriging #################################
################################################################## 
if(type_krig=='Simple'||type_krig=='simple')  {  

      if(!bivariate) {  ## space and spacetime simple kringing
      
               ####gaussian  and StudenT  two piece simple kriging
               if(covmatrix$model %in% c(1,12,27,10))
               {
                             pp <- c(muloc)      +  krig_weights %*% (c(dataT)-c(mu))   
              }
                      
               ####gamma weibbull loglogistic simple kriging
               if(covmatrix$model %in% c(21,24,26))
                      {       ones=rep(1,length(c(dataT)))
                              one=rep(1,length(c(muloc)))
                              pp <- c(emuloc) * ( one+krig_weights %*% (c(dataT)/emu-ones) )    
                      }
                      
               ####log gaussian   simple kriging
               if(covmatrix$model==1&&logGausstemp)   {
                 pp <- c(muloc)      +  krig_weights %*% (c(log(dataT))-c(mu)) 
                QQ=diag(as.matrix(diag(covmatrix$param['nugget']+covmatrix$param['sill'],dimat2) - krig_weights%*%t(cc)))
                pp=exp(pp+QQ/2)
              }
               #pp <- (c(emuloc)+covmatrix$param['sill']/2) + 
                #                          krig_weights %*% (c(dataT)-exp(c(mu)+covmatrix$param['sill']/2)) 
        }     ####simple kriging
      else  {   ## bivariate  case   cokriging
                      dat <- c(dataT) - 
                                 as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$numcoord), rep(covmatrix$param['mean_2'],covmatrix$numcoord)))
                      if(which==1) pp <- param$mean_1 + krig_weights %*% dat
                      if(which==2) pp <- param$mean_2 + krig_weights %*% dat
            } 
   #######################################         
   #### MSE COMPUTATION ##################
   #######################################
      if(mse) {
if(covmatrix$model %in% c(1,12,27,10))  vv=diag(as.matrix(diag(vvar,dimat2) - krig_weights%*%t(cc)))  ## simple variance  kriging predictor variance
if(covmatrix$model %in% c(21)) vv=emuloc^2*diag(as.matrix(diag(2/covmatrix$param['shape'],dimat2)   
                                                - krig_weights%*%t(cc)))
if(covmatrix$model %in% c(26)) vv=emuloc^2*diag(as.matrix(diag( gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1,dimat2)   
                                                - krig_weights%*%t(cc)))
if(covmatrix$model==1&&logGausstemp)
       vv <-    exp(muloc+covmatrix$param['nugget']+covmatrix$param['sill']/2)^2 *  
                               diag(as.matrix(diag(exp(vvar),dimat2) - exp(krig_weights%*%t(cc)))) 
               }        
            }               
##################################################################
################# ordinary kriging ###############################
################################################################## 
if(type_krig=='Ordinary'||type_krig=='ordinary')  {
                 if(!bivariate) { 
                         betas=  solve(t(X)%*%invcov%*%X)%*%t(X)%*%invcov%*%data    # GLS estomator of Beta
                         pp <- c(Xloc%*%betas) + krig_weights %*% (c(dataT)-c(X%*%betas)) } 
                  else{}      ##todo
                  if(mse) {ss <- (Xloc-krig_weights%*%X)%*%(solve(t(X)%*%invcov%*%X)%*%t(X))%*%invcov   + krig_weights
                           vv <-  diag(as.matrix(diag(vvar,dimat2)+ krig_weights %*% t(cc) -2*t(ss)%*%cc)) }  ## ordinary kriging predictor variance
                      }
########################## 

######################################################
####formatting data ##################################
######################################################
if(mse)
{
     if(spacetime||bivariate)  varpred=matrix(c(vv),nrow=tloc,ncol=numloc)
     else                      varpred=c(vv)
}

if(spacetime||bivariate)  pred=matrix(t(pp),nrow=tloc,ncol=numloc)
else pred=c(pp)
}#end gaussian standard kriging
###################################
##### taper kriging  ##############
###################################
if(type=="Tapering"||type=="tapering")  {
     if(!is.null(covmatrix$tapmod)) tapmod <-CkCorrModel(covmatrix$tapmod)
    else    tapmod <-NULL
    tp=.C('Corr_c_tap',corri=double(covmatrix$numcoord * covmatrix$numtime* numloc * tloc), 
                       corri_tap=double(covmatrix$numcoord * covmatrix$numtime* numloc * tloc),
        as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
        as.integer(corrmodel),as.integer(tapmod),as.integer(FALSE),as.double(locx),as.double(locy),
        as.double(c(covmatrix$maxdist,covmatrix$tapsep)),as.double(covmatrix$maxtime),
        as.integer(covmatrix$numcoord),
        as.integer(numloc),as.integer(covmatrix$ns),as.integer(tloc),as.integer(covmatrix$numtime),as.double(corrparam),as.integer(covmatrix$spacetime),
        as.integer(covmatrix$bivariate),as.double(time),as.integer(distance),as.integer(which-1),
        as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
    corri_tap=tp$corri_tap;corri=tp$corri
     if(mse){
             if(bivariate)  { if(which==1)   vvar=covmatrix$param["sill_1"]; if(which==2)   vvar=covmatrix$param["sill_2"]}
             else   vvar=covmatrix$param["sill"]
           }
      #### inverse of var covar  ##################################   
        invcov <- getInv(covmatrix)  ### invserse of cov matrix
        #############################################################
        if(!bivariate) cc <- t(matrix(corri*vvar,nrow=dimat,ncol=dimat2))
        else           cc <- t(matrix(corri,nrow=dimat,ncol=dimat2))
        krig_weights <- cc%*%invcov

    cc_tap <- t(matrix(corri_tap,nrow=dimat,ncol=dimat2))
    krig_weights_tap1 <- cc_tap%*%as.matrix(invcov)
   
    if(type_krig=='Simple'||type_krig=='simple')  {
    if(!bivariate) pp <- c(Xloc%*%betas) + krig_weights_tap1 %*% (c(dataT)-c(X%*%betas)) 
     else           { 
             #dat <- as.numeric((c(dataT)-rep(c(covmatrix$param['mean_1'],covmatrix$param['mean_2']),covmatrix$numcoord))  ) 
              #dat=as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$numcoord), rep(covmatrix$param['mean_2'],covmatrix$numcoord)))
               dat <- c(dataT) - 
                                 as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$numcoord), rep(covmatrix$param['mean_2'],covmatrix$numcoord)))
             ww1<- krig_weights_tap1  %*% dat
             if(which==1) pp <- param$mean_1 + ww1  
             if(which==2) pp <- param$mean_2 + ww1    
                     }    
     if(mse) {vv <- diag(as.matrix(diag(vvar,dimat2) - 2*krig_weights_tap1%*%t(cc))
                            +krig_weights_tap1%*%covmatrix_true$covmatrix%*%t(krig_weights_tap1) )       ## simple variance kriging tapering predictor variance
              vv2 <- diag(as.matrix(diag(vvar,dimat2) - krig_weights_tap1%*%t(cc_tap)))}
     
       if(spacetime||bivariate) {varpred=matrix(c(vv),nrow=tloc,ncol=numloc); varpred2=matrix(c(vv2),nrow=tloc,ncol=numloc);} 
        else                    {varpred=c(vv);varpred2=c(vv2)} 
     }
           if(spacetime||bivariate) pred=matrix(t(pp),nrow=tloc,ncol=numloc)
           else pred=c(pp)
    }     ##### end tapering
} #### 
####################################################################################################################################
###################### binomial and geometric case (only simple kriging) #####################################
####################################################################################################################################

if(covmatrix$model %in% c(2,11,14,19))
{  
         if(type=="Standard"||type=="standard") {

     mu0 = Xloc%*%betas; mu  = X%*%betas 
     p0=pnorm(mu0); pmu=pnorm(mu)   
     kk=0
    if(covmatrix$model==2||covmatrix$model==11) kk=min(n)
    if(covmatrix$model==19) kk=min(nloc)
    ## Computing correlation between the locations to predict and the locations observed
    ccorr=.C('Corr_c_bin',corri=double(dimat*dimat2), as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
    as.integer(corrmodel),as.integer(FALSE),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
    as.integer(numloc),as.integer(covmatrix$model),as.integer(tloc),
    as.double(kk),as.integer(covmatrix$ns),as.integer(covmatrix$numtime),as.double(c(mu0)),as.double(other_nuis),as.double(corrparam),as.integer(covmatrix$spacetime),
    as.integer(bivariate),as.double(time),as.integer(distance),as.integer(which-1),
    as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
    corri=ccorr$corri
    ###inverse of cov matrix
    invcov <- getInv(covmatrix)  
        #############################################################
        if(!bivariate) cc <- t(matrix(corri,nrow=dimat,ncol=dimat2))
        else           {}
        krig_weights <- cc%*%invcov
       if(type_krig=='Simple'||type_krig=='simple')  {
       ##########################################################
       if(covmatrix$model==2||covmatrix$model==11){  ### binomial
            if(!bivariate) pp <- n*c(p0) + krig_weights %*% (c(dataT)-n*c(pmu))  ## simple kriging
            else{} #todo
           if(mse){
                   vvar=c(n*p0*(1-p0))  ### variance (possibly no stationary)
                   vv <- diag(as.matrix(diag(vvar,dimat2)   - krig_weights%*%t(cc)))}
          } 
      ###########################################################    
           if(covmatrix$model==19){  ### binomial2
            if(!bivariate)   pp <- nloc*c(p0) + krig_weights %*% (c(dataT)-c(n*pmu))  ## simple kriging
            else{} #todo
          if(mse) {vvar=c(nloc*p0*(1-p0)) ### variance (possibly no stationary)
                   vv <- diag(as.matrix(diag(vvar,dimat2)   - krig_weights%*%t(cc))) }
          }
         ##########################################################
       if(covmatrix$model==14){    ###geometric
            if(!bivariate) ## space and spacetime
            { k1=c(pnorm(mu0));k2=c(pnorm(mu)); 
              pp <- (1-k1)/k1 + krig_weights %*% (c(dataT)-(1-k2)/k2) }
            else{}   #tood
            if(mse) { vvar=(1-k1)/k1^2   ### variance (possibly no stationary)
                      vv <- diag(as.matrix(diag(vvar,dimat2)   - krig_weights%*%t(cc))) }
          }   
        }

     if(type_krig=='Ordinary'||type_krig=='ordinary')  {
                 if(!bivariate) { 
                          betas=  solve(t(X)%*%invcov%*%X)%*%t(X)%*%invcov%*%dataT    # GLS estomator of Beta
                         if(covmatrix$model==2||covmatrix$model==11||covmatrix$model==19)
                          pp <- nloc*c(pnorm(Xloc%*%betas)) + krig_weights %*% (c(dataT)-c(n*pnorm(X%*%betas)))   
                          if(covmatrix$model==14){
                          k1=c(pnorm(Xloc%*%betas));k2=c(pnorm(X%*%betas))
                          pp <- (1-k1)/k1 + krig_weights %*% (c(dataT)-(1-k2)/k2)  }
                        }
                  else{}     ### todo   
                  if(mse) {ss <- (Xloc-krig_weights%*%X)%*%(solve(t(X)%*%invcov%*%X)%*%t(X))%*%invcov   + krig_weights
                           vv <-  diag(as.matrix(diag(vvar,dimat2)+ krig_weights %*% t(cc) -2*t(ss)%*%cc)) }  ## ordinary kriging predictor variance
                      }    
      if(spacetime||bivariate) {
            pred=matrix(t(pp),nrow=tloc,ncol=numloc);
            varpred=matrix(c(vv),nrow=tloc,ncol=numloc);
            } 
          else{pred=c(pp);varpred=c(vv)}

    }}  ##end  binary or binomial or geometric kriging 
########################################################################################
######################  to do!!!!!##################################################################
########################################################################################
if(covmatrix$model==9) {   ### ### optimal linear  tukey gaussian case 
         if(type=="Standard"||type=="standard") {
          me=as.numeric(param$mean)
          nug=as.numeric(covmatrix$param["nugget"]); v=as.numeric(covmatrix$param["sill"])
          sk=as.numeric(covmatrix$param["skew"]); tl= as.numeric(covmatrix$param["tail"])
          vvar=2 #####
        nuis<-c(me,nug,v,sk,tl)
    ## Computing covariance between the locations to predict and the locations observed
    cc=.C('Corr_c_tukey',corri=double(dimat*dimat2), as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
    as.integer(corrmodel),as.integer(covmatrix$grid),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
    as.integer(numloc),as.integer(covmatrix$model),as.integer(tloc),
    as.integer(n),as.integer(covmatrix$ns),as.integer(covmatrix$numtime),as.double(nuis),as.double(corrparam),as.integer(covmatrix$spacetime),
    as.integer(bivariate),as.double(time),as.integer(distance),as.integer(which-1),
    as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)

     t1=1-tl;t2=t1^2-tl^2;sk2=sk^2
     if(sk==0) tm=0                         ###means
     else      tm=(exp(sk2/(2*t1))-1)/(sk*sqrt(t1))
     mm=me + sqrt(v) * tm
       ###inverse of cov matrix
    invcov <- getInv(covmatrix)
       if(!bivariate) cc <- t(matrix(cc$corri*v,nrow=dimat,ncol=dimat2))
       krig_weights <- cc%*%invcov
       if(type_krig=='Simple'||type_krig=='simple')  {
            if(!bivariate)  pp <- mm+ krig_weights %*% (c(dataT)-mm)   ## simple kriging predictor for skew data
            else {} ## todo
         }
       if(mse) vv <- diag(as.matrix(diag(vvar,dimat2) - krig_weights%*%t(cc)))  ## simple mse  kriging predictor variance for skew data           
        
       if(type_krig=='Ordinary'||type_krig=='ordinary')  {}    

          if(spacetime||bivariate) {
            pred=matrix(t(pp),nrow=tloc,ncol=numloc);
            varpred=matrix(c(vv),nrow=tloc,ncol=numloc);
            } 
          else{pred=c(pp);varpred=c(vv)}
    }}  ##end  tukey case
########################################################################################
########################################################################################
########################################################################################

if(tloc==1)  {c(pred);c(varpred);c(varpred2)}
    # Return the objects list:
    Kg <- list(    bivariate=bivariate,
                   coordx = covmatrix$coordx,
                   coordy = covmatrix$coordy,
                   coordt = covmatrix$coordt,
                   covmatrix=covmatrix$covmatrix,
                   corrmodel = corrmodel,
                   data=data,
                   distance = distance,
                   grid=covmatrix$grid,
                   loc=loc,
                   nozero=covmatrix$nozero,
                   numcoord = covmatrix$numcoord,
                   numloc= numloc,
                   numtime = covmatrix$numtime,
                   numt = tloc,
                   maxdist=maxdist,
                   maxtime=maxtime,
                   model=covmatrix$model,
                   param = covmatrix$param,
                   pred=pred,
                   radius=radius,
                   spacetime = covmatrix$spacetime,
                   tapmod=tapmod,
                   time=time,
                   type=type,
                   type_krig=type_krig,
                   mse=varpred,
                   mse2=varpred2)
    structure(c(Kg, call = call), class = c("Kg"))
}



















































































#####################################################
#####################################################
#####################################################
#####################################################
############# pairwise kriging  #####################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
 if(type=="Pairwise"||type=="pairwise"||type=="Univariate"||type=="univariate") {  ## pairwise kriging

     ### control on input data


  

    checkinput <- CkInput(coordx, coordy, coordt, corrmodel, data, distance, "Fitting",
                             NULL, grid, 'Marginal', maxdist, maxtime, model, n, 
                             "Nelder-Mead", NULL, radius,  NULL, NULL, NULL,
                             type, FALSE, 'SubSamp', FALSE)
     ### control on input kriging
    checkinput <- CkInput(loc, coordy, coordt, corrmodel, data, distance, "Kriging",
                             NULL, grid, NULL,  maxdist, maxtime, model=model, n,  
                             NULL, param, radius, NULL, taper, tapsep, type_krig, 
                             NULL, NULL, NULL)
    if(!is.null(checkinput$error)) stop(checkinput$error)
    # Initialising the parameters:


                             initparam <- StartParam(coordx, coordy, coordt, corrmodel, data, distance, "Simulation",
                             NULL, grid, NULL, maxdist, maxtime, model, n, 
                             param, NULL, NULL, radius, NULL, taper, tapsep, type, type,
                             NULL, NULL, FALSE, NULL, NULL,NULL,NULL,NULL)                                          
    cmodel<-corrmodel
    cdistance<-distance
    corrmodel <- initparam$corrmodel
    distance <- initparam$distance
    corrparam <- unlist(initparam$param[initparam$namescorr])# selecting the correlation parametrs
    nuisparam <- unlist(initparam$param[initparam$namesnuis])# selecting the correlation parametrs

    bivariate <- initparam$bivariate
    if(bivariate) if(!(which==1 || which==2) ) stop("which  parameter must be 1 or 2")
    spt=initparam$spacetime; if(!spt) times=1
    pred <- NULL
    varpred <- NULL;
    k <- 0
    ### if data on a grid 
    if(grid) ccc=expand.grid(initparam$coordx,initparam$coordy)
    else     ccc=cbind(initparam$coordx,initparam$coordy)
    if(type_krig=='Simple'||type_krig=='simple')  type_krig=0
    if(type_krig=='Ordinary'||type_krig=='ordinary')  type_krig=1  ## to be done
    
    ###########################################################
    ## Standard simple pairwise (co)-kriging  ##############
    ###########################################################
    # gaussian case lin_opt and pair kriging are equal
    #if(initparam$model==1)  {if(lin_opt==TRUE)  lin_opt=FALSE}  
    if(type=="Pairwise"||type=="Pairwise") pair=1 
    if(type=="Univariate"||type=="univariate") pair=0

 

    res=double(numloc*tloc)
     
        cc=.C("pair_k",as.double(ccc[,1]),as.double(ccc[,2]),as.double(coordt),as.integer(corrmodel),
          as.double(data),as.double(locx),as.double(locy),as.double(time),as.integer(n),as.integer(initparam$numcoord),as.integer(numloc),
          as.integer(tloc),as.integer(initparam$numtime),as.double(maxdist),as.double(maxtime),as.double(nuisparam),
          as.integer(initparam$model),as.double(corrparam),as.integer(pair),as.double(initparam$radius),as.integer(type_krig),as.integer(initparam$distance) ,
          as.integer(spt),as.integer(lin_opt),as.integer(weigthed),res=as.double(res),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
    pred=matrix(cc$res,ncol=numloc,nrow=tloc,byrow=T)    

   if(mse){
       if(type_mse=="Theoretical"){
         res1=double(numloc*tloc)
          cc1=.C("mse_pair_k",as.double(ccc[,1]),as.double(ccc[,2]),as.double(coordt),as.integer(corrmodel),
          as.double(locx),as.double(locy),as.double(time),as.integer(n),as.integer(initparam$numcoord),as.integer(numloc),
          as.integer(tloc),as.integer(initparam$numtime),as.double(maxdist),as.double(maxtime),as.double(nuisparam),
          as.integer(initparam$model),as.double(corrparam),as.double(initparam$radius),as.integer(type_krig),as.integer(initparam$distance) ,
          as.integer(spt),as.integer(lin_opt),as.integer(weigthed),res1=as.double(res1),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
         }   
       if(type_mse=="SubSamp"){a=0}
      
  varpred=matrix(cc1$res1,ncol=numloc,nrow=tloc,byrow=T)
   }
    
  
    #if(tloc==1)  {c(pred);c(varpred)}
    # Return the objects list:
    Kg <- list(    bivariate=bivariate,
                   coordx = initparam$coordx,
                   coordy = initparam$coordy,
                   coordt = initparam$coordt,
                   covmatrix=NULL,
                   corrmodel = corrmodel,
                   data=data,
                   distance = distance,
                   grid=initparam$grid,
                   loc=loc,
                   nozero=NULL,
                   numcoord = initparam$numcoord,
                   numloc= numloc,
                   numtime = initparam$numtime,
                   numt = tloc,
                   maxdist=maxdist,
                   maxtime=maxtime,
                   model=initparam$model,
                   param = initparam$param,
                   pred=pred,
                   radius=radius,
                   spacetime = initparam$spacetime,
                   tapmod=NULL,
                   time=time,
                   type=type,
                   type_krig=type_krig,
                   mse=varpred)
    structure(c(Kg, call = call), class = c("Kg"))
}
return(Kg)
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
Prscores<-function(data,method="cholesky",matrix)   {

varcov=matrix$covmatrix
data=c(data)
if(!matrix$sparse){
decompvarcov <- MatDecomp(varcov,method)
inv <- MatInv(decompvarcov,method)
}
if(matrix$sparse){ 
decompvarcov <- chol(varcov)
inv <- spam::solve.spam(decompvarcov)
}
vv=diag(inv)                                                                              
###########################
nsites=length(matrix$coordx)
ntime=1
if(matrix$spacetime) ntime=length(matrix$coordt)
if(matrix$bivariate) ntime=2
dime <- nsites*ntime
###########################
D=diag(1/vv,dime,dime)
DD=sqrt(D)
temp=inv%*%data
z=D%*%temp
zz=DD%*%temp
RMSE=sqrt((1/dime)*(t(z)%*%z))
LSCORE=(1/(2*dime))*(sum(log(2*pi/vv))+sum(zz^2))
CRPS=(1/dime)*(sum((1/vv)^0.5*zz*(2*pnorm(zz)-1))+2*sum((1/vv)^0.5*pnorm(zz))+sum((1/vv)^0.5)/sqrt(pi))
###########################
scores <- list(RMSE = RMSE,
               LSCORE = LSCORE,
               CRPS = CRPS)
return(scores)
}

