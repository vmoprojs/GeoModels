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
### Last change: 28/04/2020.
#################################################### 
GeoKrig<- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, maxdist=NULL,
               maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE, lin_opt=TRUE, param, radius=6371, sparse=FALSE, 
               taper=NULL, tapsep=NULL, time=NULL, type="Standard",type_mse=NULL, type_krig="Simple",weigthed=TRUE, 
               which=1, X=NULL,Xloc=NULL)

{ 
######################################
   getInv<-function(covmatrix,b){
     if(!covmatrix$sparse){
               U =MatDecomp(covmatrix$covmatrix,method)
               Inv=MatInv(U,method)
               if(is.logical(U)){print(" Covariance matrix is not positive definite");stop()}      
               return(list(a=backsolve(U, backsolve(U, b, transpose = TRUE)),b=Inv))
             }
 if(covmatrix$sparse){ 
               if(spam::is.spam(covmatrix))  U = try(spam::chol(covmatrix$covmatrix),silent=TRUE)
               else                    U = try(spam::chol(spam::as.spam(covmatrix$covmatrix)),silent=TRUE)
               if(class(U)=="try-error") {print(" Covariance matrix is not positive definite");stop()}
               Inv=spam::solve.spam(U)
              return(list(a=spam::backsolve(U, spam::forwardsolve(U, b)),b=Inv))
        }
    }
###################################### 
######################################  
    corrmodel=gsub("[[:blank:]]", "",corrmodel)
    model=gsub("[[:blank:]]", "",model)
    distance=gsub("[[:blank:]]", "",distance)
    method=gsub("[[:blank:]]", "",method)
    type_krig=gsub("[[:blank:]]", "",type_krig)
    type=gsub("[[:blank:]]", "",type)
#####################################
    if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
    if(!is.matrix(loc))   stop("loc parameter must be a matrix")
    if(!(ncol(loc)==2))   stop("loc parameter must be a matrix  N X 2")
    if(!is.null(Xloc)) Xloc=as.matrix(Xloc)
    if(is.matrix(X) &&is.null(Xloc))  stop("Covariates for locations to predict are missing ")
    if(is.null(X) &&is.matrix(Xloc))  stop("Covariates  are missing ")
    if(CheckST(CkCorrModel(corrmodel))) if(is.null(time)) 
              stop("At least one temporal instants is needed for space-time kriging ")
###################################### 
###################################### 
    #### number of points to predict
    numloc <- nrow(loc); tloc <- length(time);
    if(!tloc)  tloc <- 1
    locx <- loc[,1];locy <- loc[,2]
    bb=0
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
     if(model=="SinhAsinh"||model=="Tukeygh") {
                             model="Gaussian"    # we need a "Gaussian" covariance matrix
                             vv=param['sill']; sk=param['skew']; tail=param['tail']
                             param['skew']=NULL; param['tail']=NULL; param['sill']=1
                             }
    if(model=="Tukeyh")    { model="Gaussian"    # we need a "Gaussian" covariance matrix
                             vv=param['sill']; tail=param['tail']
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

    bivariate <- covmatrix$bivariate;   
    if(bivariate) tloc=1
    spacetime <- covmatrix$spacetime; 

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
    NS=0
    if(spacetime||bivariate)
         { NS=cumsum(covmatrix$ns); NS=c(0,NS)[-(length(covmatrix$ns)+1)] }

    if(is.null(Xloc)) Xloc=as.matrix(rep(1,dimat2))
    else {
    if(spacetime_dyn) Xloc=as.matrix(Xloc)
    } 
    nuisance <- param[covmatrix$namesnuis]
    sel=substr(names(nuisance),1,4)=="mean"
    betas=as.numeric(nuisance[sel])   ## mean paramteres
  
    if(bivariate) {
                 sel1=substr(names(nuisance),1,6)=="mean_1"
                 betas1=as.numeric(nuisance[sel1])   ## mean1 paramteres
                 sel2=substr(names(nuisance),1,6)=="mean_2"
                 betas2=as.numeric(nuisance[sel2])   ## mean1 paramteres
               }
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
    if(bivariate) if(!(which==1 || which==2) ) stop("which  parameter must be 1 or 2")
    pred <- NULL
    varpred<-varpred2<-vv<-vv2<-NULL
    k <- 0 
    ccc=cbind(covmatrix$coordx,covmatrix$coordy)
    if(grid) ccc=expand.grid(covmatrix$coordx,covmatrix$coordy)
    else  { 
     if((spacetime||bivariate)&&(!spacetime_dyn)) ccc=cbind(rep(covmatrix$coordx,covmatrix$numtime),
                                                            rep(covmatrix$coordy,covmatrix$numtime))

     if((spacetime||bivariate)&&( spacetime_dyn)) ccc=do.call(rbind,args=c(coordx_dyn)) 
      }
    ###############################################################
    if((spacetime||bivariate)&&spacetime_dyn) dataT=t(unlist(data)) 
    else dataT=t(data)

    if(bivariate){ X11=X[1:covmatrix$ns[1],]
                   X22=X[(covmatrix$ns[1]+1):(covmatrix$ns[1]+covmatrix$ns[2]),]               
                   if(!is.null(Xloc))
                         {
                          X11_loc=Xloc[(1:(nrow(Xloc)/2)),]
                          X22_loc=Xloc[(nrow(Xloc)/2+1):nrow(Xloc),]}
                  }
########################################################################################
########################################################################################
########################################################################################

if(covmatrix$model %in% c(1,10,21,12,26,24,27,38,29,20,34,39))
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
## twopieceTukeyh=38
## sihasin=20
## tukey=34
    ################################
    ## standard kriging  ##############
    ################################   
       if(!bivariate) mu=X%*%betas
       if(bivariate)  mu=c(X11%*%matrix(betas1),X22%*%matrix(betas2))

   
       if(!bivariate) muloc=Xloc%*%betas
       if(bivariate)  {if(!is.null(Xloc)) muloc=c(X11_loc%*%matrix(betas1),X22_loc%*%matrix(betas2))}
    if((type=="Standard"||type=="standard")) {
         ## Computing CORRELATIONS between the locations to predict and the locations observed
        cc=.C('Corr_c',corri=double(dimat*dimat2), as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
        as.integer(corrmodel),as.integer(FALSE),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
        as.integer(numloc),as.integer(tloc),as.integer(covmatrix$ns),
        as.integer(NS),
        as.integer(covmatrix$numtime),as.double(corrparam),
        as.integer(covmatrix$spacetime),
        as.integer(covmatrix$bivariate),as.double(time),as.integer(distance),as.integer(which-1),
        as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
        if(covmatrix$model==1) { #gaussian 
                         vv=as.numeric(covmatrix$param['sill'])
                         corri=cc$corri 
                          }
        if(covmatrix$model==20||covmatrix$model==34)  # sinh tukeyh
                 { corri=((1-as.numeric(covmatrix$param["nugget"]))*cc$corri) }

        if(covmatrix$model==10) {    #skew gaussian
                        cc=((1-as.numeric(covmatrix$param["nugget"]))*cc$corri)
                        corr2=cc^2
                        sk=as.numeric(covmatrix$param['skew']);sk2=sk^2
                        vv=as.numeric(covmatrix$param['sill'])
                        corri=((2*sk2/pi)*(sqrt(1-corr2) + cc*asin(cc)-1) + cc*vv)/(vv+sk2*(1-2/pi));                  
                        }        
        if(covmatrix$model==21) { # gamma
                        cc=((1-as.numeric(covmatrix$param["nugget"]))*cc$corri)
                        corri=cc^2
                        sh=as.numeric(covmatrix$param["shape"])
                                }
        if(covmatrix$model==12) # student T
                         {
                        cc=as.numeric((1-as.numeric(covmatrix$param["nugget"]))*cc$corri ) 
                        vv=as.numeric(covmatrix$param['sill']) 
                        nu=1/as.numeric(covmatrix$param['df'])
                        if(nu<170) corri=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,cc^2))*cc)/(2*gamma(nu/2)^2)
                        else       corri=exp(log(nu-2)+2*lgamma(0.5*(nu-1))-log(2)-2*lgamma(nu/2)+log(Re(hypergeo::hypergeo(0.5,0.5, nu/2,cc^2)))+log(cc))
                      }
           if(covmatrix$model==18) # skew student T
                         {
                        cc=as.numeric((1-as.numeric(covmatrix$param["nugget"]))*cc$corri ) 
                        vv=as.numeric(covmatrix$param['sill']) 
                        nu=1/as.numeric(covmatrix$param['df'])
                        sk=as.numeric(covmatrix$param['skew']);sk2=sk^2
                        w=sqrt(1-sk2)
                        KK=2*sk2/pi
                        D1=(nu-1)/2;D2=nu/2;
                        corr2=(2*sk2/(pi*w^2+sk2*(pi-2)))*(sqrt(1-cc^2)+cc*asin(cc)-1)+w^2*cc/(w^2+sk2*(1-2/pi))  
                        corri=(pi*(nu-2)*gamma(D1)^2/(2*(pi*gamma(D2)^2-sk2*(nu-2)*gamma(D1)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,nu/2,cc^2))*((1-KK)*corr2+KK)-KK) 
                      }
        if(covmatrix$model==26) {  # weibull 
                        cc=as.numeric((1-covmatrix$param["nugget"])*cc$corri)
                        sh=as.numeric(covmatrix$param['shape'])
                        bcorr=    (gamma(1+1/sh))^2/((gamma(1+2/sh))-(gamma(1+1/sh))^2)
                        corri=bcorr*((1-cc^2)^(1+2/sh)*Re(hypergeo::hypergeo(1+1/sh,1+1/sh ,1 ,cc^2)) -1)                       
         }
          if(covmatrix$model==24) {  # loglogistic
                        cc1=(1-as.numeric(covmatrix$param["nugget"]))*cc$corri
                        sh=as.numeric(covmatrix$param['shape'])
                        corri=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                                  (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,cc1^2))*
                                        Re(hypergeo::hypergeo(1/sh, 1/sh, 1,cc1^2)) -1)              
         }
         if(covmatrix$model==27) {  # two piece StudenT
                        cc=(1-as.numeric(covmatrix$param["nugget"]))*cc$corri
                        nu=1/as.numeric(covmatrix$param['df']);sk=as.numeric(covmatrix$param['skew'])
                        vv=as.numeric(covmatrix$param['sill'])
                        corr2=cc^2;sk2=sk^2
                        a1=Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,corr2))
                        a2=cc*asin(cc) + (1-corr2)^(0.5)
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = cc, recycle = TRUE)
                        a3=3*sk2 + 2*sk + 4*p11 - 1
                        KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma((nu)/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
                        corri= KK*(a1*a2*a3-4*sk2);         

                      }
         if(covmatrix$model==38) {  # two piece Tukey h
                        cc=(1-as.numeric(covmatrix$param["nugget"]))*cc$corri
                        tail=as.numeric(covmatrix$param['tail']);sk=as.numeric(covmatrix$param['skew'])
                        vv=as.numeric(covmatrix$param['sill'])
                        corr2=cc^2;sk2=sk^2;
                        gg2=(1-(1-corr2)*tail)^2
                        xx=corr2/gg2
                        A=(asin(sqrt(xx))*sqrt(xx)+sqrt(1-xx))/(1-xx)^(1.5)
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = cc, recycle = TRUE)
                        a3=3*sk2 + 2*sk + 4*p11 - 1
                        mm=8*sk2/(pi*(1-tail)^2); 
                        ff=(1+3*sk2)/(1-2*tail)^(1.5)
                        M=(2*(1-corr2)^(3/2))/(pi*gg2)
                        corri=  (M*A*a3-mm)/( ff- mm)      
                      }
             if(covmatrix$model==39) {  # bimodal
                                 correlation=(1-as.numeric(covmatrix$param["nugget"]))*cc$corri
                                 nu=as.numeric(1/covmatrix$param['df']); sk=as.numeric(covmatrix$param['skew'])
                                 vv=as.numeric(covmatrix$param['sill'])
                                 ll=qnorm((1-sk)/2)
                                 p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                                 corr2=correlation^2;sk2=sk^2
                                 a1=Re(hypergeo::hypergeo(-0.5,-0.5,nu/2,corr2))
                                 a3=3*sk2 + 2*sk + 4*p11 - 1
                                 MM=nu*(1+3*sk2)*gamma(nu/2)^2-8*sk2*gamma(0.5*(nu+1))^2
                                 KK=2*gamma((nu+1)/2)^2 / MM
                                 corri= KK*(a1*a3-4*sk2)  
                      }
        if(covmatrix$model==29) {  # two piece Gaussian 
                          cc=(1-as.numeric(covmatrix$param["nugget"]))*cc$corri
                          corr2=sqrt(1-cc^2)
                          vv=as.numeric(covmatrix$param['sill'])
                          sk=as.numeric(nuisance['skew']); sk2=sk^2
                          ll=qnorm((1-sk)/2)
                          p11=pbivnorm::pbivnorm(ll,ll, rho = cc, recycle = TRUE)
                          KK=3*sk2+2*sk+ 4*p11 - 1
                          corri=(2*((corr2 + cc*asin(cc))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )                  
                      }
        if(mse){ 
             if(bivariate)  { 
                          if(which==1)    vvar=covmatrix$param["sill_1"]+covmatrix$param["nugget_1"]
                          if(which==2)    vvar=covmatrix$param["sill_2"]+covmatrix$param["nugget_2"]
                            }
             else    {
                  if(covmatrix$model==1)   vvar= vv     #gaussian
                  if(covmatrix$model==10)  vvar= vv+sk^2*(1-2/pi)   ## skewgaus
                  if(covmatrix$model==12)  vvar= vv*nu/(nu-2)              ## studentT
                  if(covmatrix$model==27)  vvar= vv*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*(gamma(0.5*(nu-1))/gamma(0.5*nu))^2) # two piece studentT
                  if(covmatrix$model==38)  vvar= vv*((1-2*tail)^(-1.5)* (1+3*(sk2)) - 4*(sk2)*2/(pi*(1-tail)^2)) # two piecetukeyh
                  if(covmatrix$model==39)  vvar= vv*MM/(gamma(0.5*nu)^2)      # bimodal                   
                  if(covmatrix$model==29)  vvar= vv*((1+3*sk2) - 8*sk2/pi )                 # two piece Gaussian
                  if(covmatrix$model==18)  vvar= vv*(nu/(nu-2) -  (nu*sk2/pi)*(gamma(D1)/gamma(D2))^2)

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
              #studentT
          if(covmatrix$model==12) cc=cc* vv*nu/(nu-2) 
          #sinh and Tukey ojo
         # if(covmatrix$model==20||covmatrix$model==34)  cc=cc*1  
          #skewgaussian
          if(covmatrix$model==10) {
                                 MM=sk*sqrt(2/pi)
                                 muloc=muloc + MM
                                 mu=mu +       MM
                                 cc=cc* (vv+sk^2*(1-2/pi) )
                                 }
          #skewstudendT
          if(covmatrix$model==18) 
          { 
              KK=gamma(0.5*(nu-1))/gamma(0.5*nu)
              MM=sqrt(vv*nu/pi)*sk*KK
              muloc=muloc + MM
              mu=mu +       MM
              cc=cc* vv*(nu/(nu-2) -  (nu*sk2/pi)*(gamma(D1)/gamma(D2))^2)
          }
          ##two piece gaussian
          if(covmatrix$model==29)   
                           {
                            MM=sqrt(vv)*2*sk*sqrt(2/pi)
                            muloc=muloc - MM
                            mu=mu  -  MM
                            cc=cc*  vv*((1+3*sk2) - 8*sk2/pi ) 
                           }
          ##two piece studentT
          if(covmatrix$model==27)  
          { 
             KK=gamma(0.5*(nu-1))/gamma(0.5*nu)
             MM=sqrt(vv)*2*sk*sqrt(nu/pi)*KK
             muloc=muloc - MM
             mu=mu  -      MM
             cc=cc*  vv*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*(gamma(0.5*(nu-1))/gamma(0.5*nu))^2)
          }
          ##two piece tukeyh
          if(covmatrix$model==38)   {
            MM=sqrt(vv)*2*sk*sqrt(2/pi)/(1-tail)
            muloc=muloc - MM
            mu=mu       - MM
            cc=cc* vv*((1-2*tail)^(-1.5)* (1+3*(sk2)) - 4*(sk2)*2/(pi*(1-tail)^2))
                        }
          ### bimodal
          if(covmatrix$model==39)   
             {
          MM=sk*sqrt(8)*gamma(0.5*(nu+1))/gamma(nu*0.5)
          muloc=muloc -MM
          mu=mu  -   MM
          KK=nu*(1+3*sk2)*gamma(nu/2)^2-8*sk2*gamma(0.5*(nu+1))^2
          cc=cc* vv*KK/(gamma(0.5*nu)^2)  
             }
          #####################################################################
          # gamma                
          if(covmatrix$model==21)   { 
                                    emuloc=exp(muloc) 
                                    emu=exp(mu)
                                    cc=2*cc/sh
                                    }
           # weibull                         
           if(covmatrix$model==26)  { 
                                    emuloc=exp(muloc)
                                    emu=exp(mu)
                                    cc=cc*(gamma(1+2/sh)/gamma(1+1/sh)^2-1)
                                    }
          #loglogistic                        
          if(covmatrix$model==24)  { 
                                    emuloc=exp(muloc)
                                    emu=exp(mu)
                                    cc=cc*(2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)
                                   }
         }
         else{}
##################################################################
#################kriging weights##################################
##################################################################
MM=getInv(covmatrix,cc)
krig_weights <- t(MM$a)
##################################################################
################# simple kriging #################################
################################################################## 
if(type_krig=='Simple'||type_krig=='simple')  {  
      if(!bivariate) {  ## space and spacetime simple kringing
               ###################################################
               if(covmatrix$model %in% c(1,12,27,38,29,10,18,39))   ####gaussian, StudenT, two piece  skew gaussian 
              {
                     pp <- c(muloc)      +  krig_weights %*% (c(dataT)-c(mu))   
              }
              ###################################################
               if(covmatrix$model %in% c(20)) # Sinh
               {
                    kk=krig_weights %*% (c(dataT))
                     pp <- c(muloc)      +  sqrt(vv)* sinh( (1/tail)*(asinh(kk)+sk))     
               }
              ###################################################
                 if(covmatrix$model %in% c(34)) # Tukeyh
               {
                    kk=krig_weights %*% (c(dataT))
                     pp <- c(muloc)      +  sqrt(vv)* kk*exp(tail*kk^2/2)    
               } 
               ###################################################
               ####gamma weibull loglogistic 
               if(covmatrix$model %in% c(21,24,26))
                      {       ones=rep(1,length(c(dataT)))
                              one=rep(1,length(c(muloc)))
                              pp <- c(emuloc) * ( one + krig_weights %*% (c(dataT)/emu-ones) )    
                      }
               ####log gaussian   simple kriging
               if(covmatrix$model==1&&logGausstemp)   {
                rp=as.numeric(covmatrix$param['sill'])

                 pp <- c(muloc-0.5*rp)      +  krig_weights %*% (c(log(dataT))-c(mu-0.5*rp)) 

                QQ=diag(as.matrix(diag(covmatrix$param['sill'],dimat2) - krig_weights%*%cc))
                pp=exp(pp+QQ/2) #/exp(covmatrix$param['sill']/2)
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
                aa=Xloc-krig_weights%*%X
                AA=chol2inv(chol(crossprod(X,(MM$b) %*% X)))
                bb=tcrossprod(aa%*%AA,aa)
# Gaussian,StudentT,skew-Gaussian,two piece linear kriging     
if(covmatrix$model %in% c(1,12,27,38,29,10,18,39))  
        {vv=diag(as.matrix(diag(vvar,dimat2) - krig_weights%*%cc  + bb)) } ## simple variance  kriging predictor variance
#gamma
if(covmatrix$model %in% c(21)) 
       { vv=emuloc^2*diag(as.matrix(diag(2/covmatrix$param['shape'],dimat2)- krig_weights%*%cc + bb))}
#weibull
if(covmatrix$model %in% c(26)) 
           {vv=emuloc^2*diag(as.matrix(diag( gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1,dimat2)   
                                                - krig_weights%*%cc+ bb))}
#loglogistic          
if(covmatrix$model %in% c(24)) 
           {vv=emuloc^2*diag(as.matrix(diag((2*covmatrix$param['shape']*sin(pi/covmatrix$param['shape'])^2/
                       (pi*sin(2*pi/covmatrix$param['shape']))-1),dimat2) - krig_weights%*%cc + bb))}

if(covmatrix$model==1&&logGausstemp)
       {vv <-    exp(muloc + covmatrix$param['sill']/2)^2 *  
                               diag(as.matrix(diag(exp(vvar),dimat2) - exp(krig_weights%*%cc+ bb))) }
               }     # end if(mse)   
}               


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
        as.integer(numloc),as.integer(covmatrix$ns),as.integer(NS),as.integer(tloc),as.integer(covmatrix$numtime),as.double(corrparam),as.integer(covmatrix$spacetime),
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
    MM=getInv(covmatrix,cc_tap)
    krig_weights_tap1 <- t(MM$a)    # cc_tap%*%as.matrix(invcov)
    if(type_krig=='Simple'||type_krig=='simple')  {
    if(!bivariate) 
                { pp <- c(Xloc%*%betas) + krig_weights_tap1 %*% (c(dataT)-c(X%*%betas)) }
     else           
                { 

                 dat <- c(dataT) - as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$numcoord), rep(covmatrix$param['mean_2'],covmatrix$numcoord)))
                 ww1<- krig_weights_tap1  %*% dat
                 if(which==1) pp <- param$mean_1 + ww1  
                 if(which==2) pp <- param$mean_2 + ww1    
                }    
     if(mse) {

         aa=Xloc-krig_weights_tap1%*%X
         AA=chol2inv(chol(crossprod(X,(MM$b) %*% X)))
         bb=tcrossprod(aa%*%AA,aa)

         vv <- diag(as.matrix(diag(vvar,dimat2) - 2*krig_weights_tap1%*%cc)
                            +krig_weights_tap1%*%covmatrix_true$covmatrix%*%t(krig_weights_tap1) + bb )      ## simple variance kriging tapering predictor variance

         vv2 <- diag(as.matrix(diag(vvar,dimat2) - krig_weights_tap1%*%cc_tap + bb))}
     
       if(spacetime||bivariate) {varpred=matrix(c(vv),nrow=tloc,ncol=numloc); varpred2=matrix(c(vv2),nrow=tloc,ncol=numloc);} 
        else                    {varpred=c(vv);varpred2=c(vv2)} 
     }
           if(spacetime||bivariate) pred=matrix(t(pp),nrow=tloc,ncol=numloc)
           else pred=c(pp)
    }     ##### end tapering

} #### 
####################################################################################################################################
###################### binomial  binomial negative and poisson #####################################
####################################################################################################################################

if(covmatrix$model %in% c(2,11,14,19,30))
{  
     if(type=="Standard"||type=="standard") {

     mu0 = Xloc%*%betas; 
     if(!bivariate) mu  = X%*%betas 
     if(bivariate)  mu  = c(X11%*%betas1,X22%*%betas2)
     kk=0
     if(covmatrix$model==2||covmatrix$model==11) kk=min(n)
     if(covmatrix$model==19) kk=min(nloc)
## ojo que es la covarianza
if(covmatrix$model %in% c(2,11,14,19))
    ## Computing correlation between the locations to predict and the locations observed
    ccorr=.C('Corr_c_bin',corri=double(dimat*dimat2), as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
    as.integer(corrmodel),as.integer(FALSE),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
    as.integer(numloc),as.integer(covmatrix$model),as.integer(tloc),
    as.double(kk),as.integer(covmatrix$ns),as.integer(NS),as.integer(covmatrix$numtime),as.double(c(mu0)),as.double(other_nuis),as.double(corrparam),as.integer(covmatrix$spacetime),
    as.integer(bivariate),as.double(time),as.integer(distance),as.integer(which-1),
    as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
    
if(covmatrix$model %in% c(30))
      ccorr=.C('Corr_c_poi',corri=double(dimat*dimat2), as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
    as.integer(corrmodel),as.integer(FALSE),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
    as.integer(numloc),as.integer(covmatrix$model),as.integer(tloc),
    as.double(kk),as.integer(covmatrix$ns),as.integer(NS),as.integer(covmatrix$numtime),as.double(c(mu0)),as.double(other_nuis),as.double(corrparam),as.integer(covmatrix$spacetime),
    as.integer(bivariate),as.double(time),as.integer(distance),as.integer(which-1),
    as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
 
    corri=ccorr$corri

        ##################inverse of cov matrix##############################################
        if(!bivariate) cc <- matrix(corri,nrow=dimat,ncol=dimat2)
        else           {}
        MM=getInv(covmatrix,cc)
        krig_weights <- t(MM$a)

       if(type_krig=='Simple'||type_krig=='simple')  {

          ##########################################################
       if(covmatrix$model==30){  ### poisson
        p0=exp(mu0); pmu=exp(mu) 
            if(!bivariate) 
                   {  pp <- c(p0) + krig_weights %*% (c(dataT)-c(pmu)) }  ## simple kriging
            else{} #todo
           if(mse)  vvar=p0  ### variance (possibly no stationary)     
          } 
       ##########################################################
       if(covmatrix$model==2||covmatrix$model==11){  ### binomial
        p0=pnorm(mu0); pmu=pnorm(mu) 
            if(!bivariate) 
                       { pp <- n*c(p0) + krig_weights %*% (c(dataT)-n*c(pmu)) }  ## simple kriging
            else{} #todo
           if(mse) vvar=n*p0*(1-p0)  ### variance (possibly no stationary
          } 
      ###########################################################    
           if(covmatrix$model==19){  ### binomial2
            p0=pnorm(mu0); pmu=pnorm(mu) 
            if(!bivariate)   pp <- nloc*c(p0) + krig_weights %*% (c(dataT)-c(n*pmu))  ## simple kriging
            else{} #todo
          if(mse)    vvar=nloc*p0*t(1-p0) ### variance (possibly no stationary)   
          }
         ##########################################################
       if(covmatrix$model==14){    ###geometric
         p0=pnorm(mu0); pmu=pnorm(mu) 
            if(!bivariate) ## space and spacetime
            { k1=c(p0);k2=c(pmu); 
              pp <- (1-k1)/k1 + krig_weights %*% (c(dataT)-(1-k2)/k2) }
            else{}   #tood
            if(mse) vvar=(1-k1)/k1^2   ### variance (possibly no stationary)
                
          }

        if(mse){
                  aa=Xloc-krig_weights%*%X
                  AA=chol2inv(chol(crossprod(X,(MM$b) %*% X)  ))
                  bb=tcrossprod(aa%*%AA,aa)
                  vv <- diag(sqrt(vvar%*%t(vvar))  - krig_weights%*%cc  + bb) 
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
########################################################################################
########################################################################################

if(tloc==1)  {c(pred);c(varpred);c(varpred2)}
    # Return the objects list:
    Kg <- list(    bivariate=bivariate,
                   coordx = covmatrix$coordx,
                   coordy = covmatrix$coordy,
                   coordt = covmatrix$coordt,
                   coordx_dyn=covmatrix$coordx_dyn,
                   covmatrix=covmatrix$covmatrix,
                   corrmodel = corrmodel,
                   data=data,
                   distance = distance,
                   grid=covmatrix$grid,
                   loc=loc,
                   nozero=covmatrix$nozero,
                   ns=covmatrix$ns,
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

###################################################################################################
###################################################################################################

GeoNeighborhood <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, distance="Eucl", grid=FALSE, 
                  loc, maxdist=NULL,maxtime=NULL, radius=6371, time=NULL, X=NULL)
{
  XX=NULL
  sel_ss=1
  sel_tt=1
  if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
  if(!is.matrix(loc))   stop("loc parameter must be a matrix")
    #if(!(ncol(loc)==2))   stop("loc parameter must be a matrix  N X 2")
  spacetime=FALSE
  corrmodel="Exponential"
  
  if(!is.null(coordt)) {spacetime=TRUE;corrmodel="Exp_Exp"}
  checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting",
                             NULL, grid, 'Marginal', maxdist, maxtime, 'Gaussian', 1,
                              'Nelder-Mead', NULL, radius, NULL, NULL, NULL, 
                          'Pairwise', FALSE, 'SubSamp', FALSE, X)
  if(!is.null(checkinput$error)) stop(checkinput$error)
  spacetime_dyn=FALSE
  if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
    
## handling spatial coordinates
    if(is.null(coordy)) coords=as.matrix(coordx)
    else{
    if(grid) coords=as.matrix(expand.grid(coordx,coordy))
    else     coords=cbind(coordx,coordy)  
    }

############### total dimension #####
NN=nrow(coords);TT=length(coordt);dimat2=NN*TT
#####################################
sel_tt=NULL
colnames(loc)=NULL;colnames(coords)=NULL;
a_s=rbind(loc,coords)
## type of dist
if(distance=="Eucl") dd_s=as.matrix(dist(a_s))
if(distance=="Geod") dd_s=fields::rdist.earth(a_s,miles=F,R=radius)
if(distance=="Chor") dd_s=2*sin(fields::rdist.earth(a_s,miles=F,R=radius)/2)

ss=c((dd_s<maxdist)[,1])[-1]
sel_ss=coords[ss,]
## spatial
if(!spacetime){ data_sel=data[ss]
                if(!is.null(X)) XX=X[ss,]
 #if(dim(as.matrix(sel_ss,nrow=dimat2))[1]==0) stop("spatial distance for local kringing is too small")
            }
## space-time
if(spacetime)
{
  a_t=c(time,coordt);
  dd_t=dist(a_t)
  tt=c((as.matrix(dd_t)<maxtime)[,1])[-1]
  sel_tt=coordt[tt]
  data_sel=data[tt,ss]
# if(dim(as.matrix(sel_ss,nrow=dimat2))[1]==0) stop("spatial distance for local kringing is too small")
if(!is.null(X)){
  XX=list()
  for(i in 0:(TT-1)) XX[[i+1]]=(X[(i*NN+1):((i*NN)+NN),])[ss,]
  RR=as.numeric(!tt)*seq(1:TT)
  RR=RR[RR>0]
  XX=XX[-(RR)]
  XX=do.call(rbind,args=c(XX))
  }}
return(list(data=data_sel,coordx=sel_ss,coordt=sel_tt,distance=distance,
      numcoord=nrow(sel_ss),numtime=length(sel_tt),radius=radius,spacetime=spacetime,X=XX))
} 
