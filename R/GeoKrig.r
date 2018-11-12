####################################################
### Authors: Moreno Bevilacqua Víctor Morales Oñate.
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: GeoKrig.r  
### Description:  
### This file contains a set of procedures
### for computing simple (tapered)  kriging
### predictor  at an unknown space (time) locations.
### Last change: 28/08/2018.
#################################################### 

GeoKrig<- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, maxdist=NULL,
               maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE, lin_opt=TRUE, param, radius=6371, sparse=FALSE, 
               taper=NULL, tapsep=NULL, time=NULL, type="Standard",type_mse=NULL, type_krig="Simple",weigthed=TRUE, 
               which=1, X=NULL,Xloc=NULL)

{ 
######################################
    getInv<-function(covmatrix,b){
     if(!covmatrix$sparse){
               U <- MatDecomp(covmatrix$covmatrix,method)
               if(is.logical(U)){print(" Covariance matrix is not positive definite");stop()}      
               return(backsolve(U, backsolve(U, b, transpose = TRUE)))
             }
        else{  
               if(spam::is.spam(covmatrix))  U = try(spam::chol.spam(covmatrix$covmatrix),silent=TRUE)
               else                    U = try(spam::chol.spam(spam::as.spam(covmatrix$covmatrix)),silent=TRUE)
               
               if(class(U)=="try-error") {print(" Covariance matrix is not positive definite");stop()}
              return(spam::backsolve(U, spam::forwardsolve(U, b)))
        }
    }
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
    #### number of points to predict
    numloc <- nrow(loc); tloc <- length(time);
    if(!tloc)  tloc <- 1
    locx <- loc[,1];locy <- loc[,2]
    #######################################################
    ############ standard (tapered) kriging ###############
    #######################################################
    if(type %in% c("Standard","standard","Tapering","tapering")) {
    #################################################
    ##### computing covariance  matrix ##############
    #################################################
    logGausstemp=FALSE
    if(model=="LogGaussian") {model="Gaussian"    # we need a "Gaussian" covariance matrix
                             logGausstemp=TRUE}
     if(model=="SinhAsinh") {model="Gaussian"    # we need a "Gaussian" covariance matrix
                             vv=param['sill']; sk=param['skew']; tail=param['tail']
                             param['skew']=NULL; param['tail']=NULL; param['sill']=1
                             }
    if(model %in% c("Weibull","Gamma","LogLogistic"))          # we need a x covariane matrix with with mean=0   x=gamma,weibull,loglogistic
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
     if(model %in% c("Weibull","Gamma","LogLogistic")) {
          if(is.null(Xtemp)) X=matrix(rep(1,dimat))
          else               X=Xtemp    
          param=paramtemp
          covmatrix$namesnuis=unique(c(meantemp,covmatrix$namesnuis))
    }
    else { X=covmatrix$X }
    ###############
    ###############
    num_betas=ncol(X)
    if(is.null(Xloc))   Xloc=as.matrix(rep(1,dimat2)) 
    if(spacetime_dyn)
    { 
      if(!ncol(X)==1)
      {
      if(!is.list(Xloc)) {stop("covariates must be given as a list")}
      else               {env <- new.env();Xloc=do.call(rbind,args=c(Xloc),envir = env)
                          Xloc=as.matrix(Xloc)}
                        }
    } 
    nuisance <- param[covmatrix$namesnuis]
    sel=substr(names(nuisance),1,4)=="mean"
    betas=as.numeric(nuisance[sel])   ## mean paramteres
    other_nuis=as.numeric(nuisance[!sel]) 
    ################################################
    ################################################
    if(type %in% c("Tapering","tapering")) {
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

if(covmatrix$model %in% c(1,10,21,12,26,24,27,29,20))
{  
## gaussian=1 
## skew gaussian=10   
## student =12
## gamma=21 
## weibull=26
## loglogistic=24
## loggaussian=22 ojo
## twopieceStudentT=27
## twopieceGaussian=29
## sihasin=20
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
        if(covmatrix$model==1) { #gaussian 
                         vv=as.numeric(covmatrix$param['sill'])
                         corri=cc$corri 
                          }
        if(covmatrix$model==20) corri=cc$corri # sinh
                         
        if(covmatrix$model==10) {    #skew gaussian
                        corr2=(as.numeric((1-covmatrix$param["nugget"])*cc$corri))^2
                        sk=as.numeric(covmatrix$param['skew']);sk2=sk^2
                        vv=as.numeric(covmatrix$param['sill'])
                        corri=((2*sk2/pi)*(sqrt(1-corr2) + cc$corri*asin(cc$corri)-1) + cc$corri*vv)/(vv+sk2*(1-2/pi));
                        }        
        if(covmatrix$model==21)  # gamma
                        corri=((1-as.numeric(covmatrix$param["nugget"]))*cc$corri)^2
        if(covmatrix$model==12) # student T
                         {
                        cc=as.numeric((1-as.numeric(covmatrix$param["nugget"]))*cc$corri ) 
                        vv=as.numeric(covmatrix$param['sill']) 
                        nu=1/as.numeric(covmatrix$param['df'])
                        corri=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,cc^2))*cc)/(2*gamma(nu/2)^2)
                      }
        if(covmatrix$model==26) {  # weibull 
                        sh=as.numeric(covmatrix$param['shape'])
                        bcorr=    (gamma(1+1/sh))^2/((gamma(1+2/sh))-(gamma(1+1/sh))^2)
                        cc1=as.numeric((1-covmatrix$param["nugget"])*cc$corri)
                        corri=bcorr*((1-cc1^2)^(1+2/sh)*Re(hypergeo::hypergeo(1+1/sh,1+1/sh ,1 ,cc1^2)) -1) #ojo es una tranformada de la 1F2             
         }
          if(covmatrix$model==24) {  # loglogistic
                        sh=as.numeric(covmatrix$param['shape'])
                        cc1=(1-as.numeric(covmatrix$param["nugget"]))*cc$corri
corri=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                        (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,cc1^2))*
                         Re(hypergeo::hypergeo(1/sh, 1/sh, 1,cc1^2)) -1)              
         }
         if(covmatrix$model==27) {  # two piece StudenT
                        nu=1/as.numeric(covmatrix$param['df']);sk=as.numeric(covmatrix$param['skew'])
                        vv=as.numeric(covmatrix$param['sill'])
                        corr2=cc$corri^2;sk2=sk^2
                        a1=Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,corr2))
                        a2=cc$corri*asin(cc$corri) + (1-corr2)^(0.5)
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = cc$corri, recycle = TRUE)
                        a3=3*sk2 + 2*sk + 4*p11 - 1
                        KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma((nu)/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
                        corri= KK*(a1*a2*a3-4*sk2);                     
                      }
        if(covmatrix$model==29) {  # two piece Gaussian 
                          corr2=sqrt(1-cc$corri^2)
                          vv=as.numeric(covmatrix$param['sill'])
                          sk=as.numeric(nuisance['skew']); sk2=sk^2
                          ll=qnorm((1-sk)/2)
                          p11=pbivnorm::pbivnorm(ll,ll, rho = cc$corri, recycle = TRUE)
                          KK=3*sk2+2*sk+ 4*p11 - 1
                          corri=(2*((corr2 + cc$corri*atan(cc$corri/corr2))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )                  
                      }
        if(mse){ 
             if(bivariate)  { 
                          if(which==1)    vvar=covmatrix$param["sill_1"]+covmatrix$param["nugget_1"]
                          if(which==2)    vvar=covmatrix$param["sill_2"]+covmatrix$param["nugget_2"]}
             else    {if(covmatrix$model==1)   vvar=vv     #gaussian
                      if(covmatrix$model==10)  vvar= (vv+sk^2*(1-2/pi))   ## skewgaus
                      if(covmatrix$model==12)  vvar=vv*nu/(nu-2)              ## studentT
                      if(covmatrix$model==27)  vvar=nu*(3*sk2+1)/(nu-2)-
                                                     (4*sk2*nu*gamma((nu-1)/2)^2)/(pi*gamma(nu/2)^2) # two pieceT
                      if(covmatrix$model==29)  vvar=(1+3*sk2)-8*sk2/pi                  # two piece Gaussian
                      #if(covmatrix$model==21)  vvar=2*exp(muloc)/covmatrix$param['shape']
                      #if(covmatrix$model==22) {kk=exp(2*(muloc)+covmatrix$param['sill']);vvar=kk*(exp(covmatrix$param['sill'])-1)}
                      }
           }
########################################################################################
########################################################################################
##### multiplying the correlations for the variance
        cc <- matrix(corri,nrow=dimat,ncol=dimat2)
        if(!bivariate){
          #gaussian
          if(covmatrix$model==1)  cc=cc* vv 
          #sinh
          if(covmatrix$model==20)  cc=cc*1            
          #skewgaussian
          if(covmatrix$model==10) cc=cc* (vv+sk^2*(1-2/pi)) 
          #studentT
          if(covmatrix$model==12) cc=cc* vv *(nu/(nu-2)) 
          ##two piece studentT
          if(covmatrix$model==27)   
                          cc=cc* vv *(nu*(3*sk2+1)/(nu-2)-(4*sk2*nu*gamma((nu-1)/2)^2)/(pi*gamma(nu/2)^2) )       
          ##two piece gaussian
          if(covmatrix$model==29)   
                          cc=cc* vv *((1+3*sk2)-8*sk2/pi) 
          # gamma                
          if(covmatrix$model==21)  { 
                                    emuloc=exp(muloc) 
                                    emu=exp(mu)
                                    cc=2*cc/covmatrix$param['shape']
                                    }
           # weibull                         
           if(covmatrix$model==26)  { 
                                    emuloc=exp(muloc)
                                    emu=exp(mu)
                                    cc=cc*(gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1)
                                  }
          #loglogistic                        
          if(covmatrix$model==24)  { 
                                    emuloc=exp(muloc)
                                    emu=exp(mu)
                                    cc=cc*(2*covmatrix$param['shape']*sin(pi/covmatrix$param['shape'])^2/(pi*sin(2*pi/covmatrix$param['shape']))-1)
                                  }
         }
         else{}
##################################################################
#################kriging weights##################################
##################################################################
krig_weights <- t(getInv(covmatrix,cc))
##################################################################
################# simple kriging #################################
################################################################## 
if(type_krig=='Simple'||type_krig=='simple')  {  
      if(!bivariate) {  ## space and spacetime simple kringing
               ####gaussian, StudenT  two piece  skew gaussian simple kriging
               if(covmatrix$model %in% c(1,12,27,29,10))
                     pp <- c(muloc)      +  krig_weights %*% (c(dataT)-c(mu))   
              ####sinh
               if(covmatrix$model %in% c(20))
               {
                    kk=krig_weights %*% (c(dataT))
                     pp <- c(muloc)      +  sqrt(vv)* sinh( (1/tail)*(asinh(kk)+sk))     
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
                QQ=diag(as.matrix(diag(covmatrix$param['nugget']+covmatrix$param['sill'],dimat2) - krig_weights%*%cc))
                pp=exp(pp+QQ/2)
              } 
               #pp <- (c(emuloc)+covmatrix$param['sill']/2) + 
                #                          krig_weights %*% (c(dataT)-exp(c(mu)+covmatrix$param['sill']/2)) 
        }     ####simple kriging
      else  {   ## bivariate  case   cokriging
          dat <- c(dataT) - as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$ns[1]), 
                            rep(covmatrix$param['mean_2'],covmatrix$ns[2])))
                      if(which==1) pp <- param$mean_1 + krig_weights %*% dat
                      if(which==2) pp <- param$mean_2 + krig_weights %*% dat
            } 
   #######################################         
   #### MSE COMPUTATION ##################
   #######################################

   ####### 
      if(mse) {
# Gaussian,StudentT,skew-Gaussian,two piece linear kriging     
if(covmatrix$model %in% c(1,12,27,29,10))  

vv=diag(as.matrix(diag(vvar,dimat2) - krig_weights%*%cc))  ## simple variance  kriging predictor variance
#gamma
if(covmatrix$model %in% c(21)) 

vv=emuloc^2*diag(as.matrix(diag(2/covmatrix$param['shape'],dimat2)- krig_weights%*%cc))
#weibull
if(covmatrix$model %in% c(26)) 
vv=emuloc^2*diag(as.matrix(diag( gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1,dimat2)   
                                                - krig_weights%*%cc))
#loglogistic          
if(covmatrix$model %in% c(24)) 
vv=emuloc^2*diag(as.matrix(diag((2*covmatrix$param['shape']*sin(pi/covmatrix$param['shape'])^2/
                       (pi*sin(2*pi/covmatrix$param['shape']))-1),dimat2) - krig_weights%*%cc))

if(covmatrix$model==1&&logGausstemp)
       vv <-    exp(muloc+covmatrix$param['nugget']+covmatrix$param['sill']/2)^2 *  
                               diag(as.matrix(diag(exp(vvar),dimat2) - exp(krig_weights%*%cc))) 
               }        
            }               
##################################################################
################# ordinary kriging ###############################
################################################################## 
#if(type_krig=='Ordinary'||type_krig=='ordinary')  {
#                 if(!bivariate) { 
#                         betas=  solve(t(X)%*%invcov%*%X)%*%t(X)%*%invcov%*%data    # GLS estomator of Beta
#                         pp <- c(Xloc%*%betas) + krig_weights %*% (c(dataT)-c(X%*%betas)) } 
#                  else{}      ##todo
#                  if(mse) {ss <- (Xloc-krig_weights%*%X)%*%(solve(t(X)%*%invcov%*%X)%*%t(X))%*%invcov   + krig_weights
#                           vv <-  diag(as.matrix(diag(vvar,dimat2)+ krig_weights %*% t(cc) -2*t(ss)%*%cc)) }  ## ordinary kriging predictor variance
#                      }
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
             if(bivariate)  { if(which==1)   vvar=covmatrix$param["sill_1"]+covmatrix$param["nugget_1"]; 
                              if(which==2)   vvar=covmatrix$param["sill_2"]+covmatrix$param["nugget_2"]; }
             else   vvar=covmatrix$param["sill"]+covmatrix$param["nugget"];
           }
        #############################################################
        if(!bivariate) cc <- matrix(corri*(covmatrix$param["sill"]+covmatrix$param["nugget"]),nrow=dimat,ncol=dimat2)
        else           cc <- matrix(corri,nrow=dimat,ncol=dimat2)

    cc_tap <- matrix(corri_tap,nrow=dimat,ncol=dimat2)
    krig_weights_tap1 <- t(getInv(covmatrix,cc_tap))    # cc_tap%*%as.matrix(invcov)
    if(type_krig=='Simple'||type_krig=='simple')  {
    if(!bivariate) pp <- c(Xloc%*%betas) + krig_weights_tap1 %*% (c(dataT)-c(X%*%betas)) 
     else           { 
          dat <- c(dataT) - as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$numcoord), rep(covmatrix$param['mean_2'],covmatrix$numcoord)))
             ww1<- krig_weights_tap1  %*% dat
             if(which==1) pp <- param$mean_1 + ww1  
             if(which==2) pp <- param$mean_2 + ww1    
                     }    
     if(mse) {
         vv <- diag(as.matrix(diag(vvar,dimat2) - 2*krig_weights_tap1%*%cc)
                            +krig_weights_tap1%*%covmatrix_true$covmatrix%*%t(krig_weights_tap1 ) )      ## simple variance kriging tapering predictor variance

              vv2 <- diag(as.matrix(diag(vvar,dimat2) - krig_weights_tap1%*%cc_tap))}
     
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
    
        #############################################################
        if(!bivariate) cc <- matrix(corri,nrow=dimat,ncol=dimat2)
        else           {}

        krig_weights <- t(getInv(covmatrix,cc))

       if(type_krig=='Simple'||type_krig=='simple')  {
       ##########################################################
       if(covmatrix$model==2||covmatrix$model==11){  ### binomial
            if(!bivariate) pp <- n*c(p0) + krig_weights %*% (c(dataT)-n*c(pmu))  ## simple kriging
            else{} #todo
           if(mse){
                   vvar=c(n*p0*(1-p0))  ### variance (possibly no stationary)
                   vv <- diag(as.matrix(diag(vvar,dimat2)   - krig_weights%*%cc))}
          } 
      ###########################################################    
           if(covmatrix$model==19){  ### binomial2
            if(!bivariate)   pp <- nloc*c(p0) + krig_weights %*% (c(dataT)-c(n*pmu))  ## simple kriging
            else{} #todo
          if(mse) {vvar=c(nloc*p0*(1-p0)) ### variance (possibly no stationary)
                   vv <- diag(as.matrix(diag(vvar,dimat2)   - krig_weights%*%cc)) }
          }
         ##########################################################
       if(covmatrix$model==14){    ###geometric
            if(!bivariate) ## space and spacetime
            { k1=c(pnorm(mu0));k2=c(pnorm(mu)); 
              pp <- (1-k1)/k1 + krig_weights %*% (c(dataT)-(1-k2)/k2) }
            else{}   #tood
            if(mse) { vvar=(1-k1)/k1^2   ### variance (possibly no stationary)
                      vv <- diag(as.matrix(diag(vvar,dimat2)   - krig_weights%*%cc)) }
          }   
        }

    # if(type_krig=='Ordinary'||type_krig=='ordinary')  {
            #     if(!bivariate) { 
            #              betas=  solve(t(X)%*%invcov%*%X)%*%t(X)%*%invcov%*%dataT    # GLS estomator of Beta
            #             if(covmatrix$model==2||covmatrix$model==11||covmatrix$model==19)
            #              pp <- nloc*c(pnorm(Xloc%*%betas)) + krig_weights %*% (c(dataT)-c(n*pnorm(X%*%betas)))   
            #              if(covmatrix$model==14){
            #              k1=c(pnorm(Xloc%*%betas));k2=c(pnorm(X%*%betas))
            #              pp <- (1-k1)/k1 + krig_weights %*% (c(dataT)-(1-k2)/k2)  }
            #            }
            #      else{}     ### todo   
            #      if(mse) {ss <- (Xloc-krig_weights%*%X)%*%(solve(t(X)%*%invcov%*%X)%*%t(X))%*%invcov   + krig_weights
            #               vv <-  diag(as.matrix(diag(vvar,dimat2)+ krig_weights %*% t(cc) -2*t(ss)%*%cc)) }  ## ordinary kriging predictor variance
           #           }    
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
    
       if(!bivariate) cc <- matrix(cc$corri*v,nrow=dimat,ncol=dimat2)

       krig_weights <- t(getInv(covmatrix,cc))
       
       if(type_krig=='Simple'||type_krig=='simple')  {
            if(!bivariate)  pp <- mm+ krig_weights %*% (c(dataT)-mm)   ## simple kriging predictor for skew data
            else {} ## todo
         }
       if(mse) vv <- diag(as.matrix(diag(vvar,dimat2) - krig_weights%*%cc))  ## simple mse  kriging predictor variance for skew data           
        
       #if(type_krig=='Ordinary'||type_krig=='ordinary')  {}    

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


return(Kg)
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
Prscores<-function(data,method="cholesky",matrix)   {
if(class(matrix)!="CovMat") stop("A CovMat object is needed as input\n")
varcov=matrix$covmatrix
rownames(varcov)=c();colnames(varcov)=c()
if(nrow(varcov)!=length(data)) stop("The dimension of the covariance  matrix and/or the vector data are not correct  \n")
data=c(unname(data))
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
print(dime)
###########################
D=diag(1/vv,dime,dime)
DD=sqrt(D)
temp=inv%*%data
z=D%*%temp
zz=DD%*%temp
RMSE=sqrt((1/dime)*(t(z)%*%z))
LSCORE=(1/(2*dime))*(sum(log(2*pi/vv))+sum(zz^2))
CRPS=(1/dime)*(sum((1/vv)^0.5*zz*(2*pnorm(zz)-1))+2*sum((1/vv)^0.5*pnorm(zz))+sum((1/vv)^0.5)/sqrt(pi))
MAE=(1/dime)*(sum(abs(z)))
###########################
scores <- list(RMSE = RMSE,
               LSCORE = LSCORE,
               CRPS = CRPS,
               MAE=MAE)
return(scores)
}

