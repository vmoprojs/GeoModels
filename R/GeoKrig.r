####################################################
### File name: GeoKrig.r
####################################################


GeoKrig= function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc, maxdist=NULL,
               maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE, lin_opt=TRUE, param, radius=6371, sparse=FALSE, 
               taper=NULL, tapsep=NULL, time=NULL, type="Standard",type_mse=NULL, type_krig="Simple",weigthed=TRUE, 
               which=1, copula=NULL,X=NULL,Xloc=NULL)

{ 
######################################
   getInv=function(covmatrix,b){
     if(!covmatrix$sparse){
               U =MatDecomp(covmatrix$covmatrix,method)
               #Inv=MatInv(U,method)
               Inv=0
               if(is.logical(U)){print(" Covariance matrix is not positive definite");stop()}  
               Invc=backsolve(U, backsolve(U, b, transpose = TRUE))    
               return(list(a=Invc,bb=Inv))
             }
 if(covmatrix$sparse){ 
          
               if(spam::is.spam(covmatrix$covmatrix))  U = try(spam::chol.spam(covmatrix$covmatrix),silent=TRUE)
               else                    U = try(spam::chol.spam(spam::as.spam(covmatrix$covmatrix)),silent=TRUE)
               if(class(U)=="try-error") {print(" Covariance matrix is not positive definite");stop()}
               #Inv=spam::chol2inv.spam(U)
               Inv=0
               Invc= spam::backsolve(U, spam::forwardsolve(U, b)) ## R^-1 %*% c
              return(list(a=Invc,bb=Inv))
        }
    }
###################################### 
########## START ##################### 
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
    if(!is.null(Xloc)) 
         { 
           if(is.vector(Xloc)) Xloc=matrix(Xloc,nrow=1)
           else                Xloc=as.matrix(Xloc)
         }
 
    if(is.matrix(X) &&is.null(Xloc))  stop("Covariates for locations to predict are missing ")
    if(is.null(X) &&is.matrix(Xloc))  stop("Covariates  are missing ")
    if(CheckST(CkCorrModel(corrmodel))) if(is.null(time)) 
              stop("At least one temporal instants is needed for space-time kriging ")
###################################### 
###################################### 

    #### number of points to predict
     if(is.null(time)) time=0
    numloc = nrow(loc); tloc = length(time);
    if(!tloc)  tloc = 1
    locx = loc[,1];locy = loc[,2]
    bb=0
    #######################################################
    ############ standard (tapered) kriging ###############
    #######################################################
    if(type %in% c("Standard","standard","Tapering","tapering")) {
    #################################################
    ##### computing covariance  matrix ##############
    #################################################
    logGausstemp=SinhAsinhtemp=FALSE  #Tukeyhtemp=Tukeyghtemp=

    #### cases for optimal predictor!
    if(model=="LogGaussian") {model="Gaussian"    # we need a "Gaussian" covariance matrix
                             logGausstemp=TRUE}
     #if(model=="SinhAsinh"||model=="Tukeygh") {
    if(model=="SinhAsinh") {
                              if(model=="SinhAsinh")  SinhAsinhtemp=TRUE
                              #if(model=="Tukeygh")    Tukeyghtemp=TRUE
                             model="Gaussian"    # we need a "Gaussian" covariance matrix
                             vv=as.numeric(param['sill']); sk=as.numeric(param['skew']); tail=as.numeric(param['tail'])
                             param['skew']=NULL; param['tail']=NULL; param['sill']=1
                             }
    #if(model=="Tukeyh")    { model="Gaussian"    # we need a "Gaussian" covariance matrix
     #                        vv=as.numeric(param['sill']); tail=as.numeric(param['tail'])
      #                       Tukeyhtemp=TRUE
       #                      param['skew']=NULL; param['tail']=NULL; param['sill']=1
        #                   }
   if(model %in% c("Weibull","Gamma","LogLogistic"))          # we need a x covariane matrix with with mean=0   x=gamma,weibull,loglogistic
     {
     paramtemp=param 
     sel=substr(names(param),1,4)=="mean";  ## selecting mean values 
     meantemp=names(param[sel])             ## saving  mean values
     param=param[!sel];param$mean=0;        ## mean must be =0 when calling covariance matrix
     Xtemp=X;X=NULL                         ## saving X and setting X=NULL
     }
     

    covmatrix = GeoCovmatrix(coordx=coordx, coordy=coordy, coordt=coordt, coordx_dyn=coordx_dyn, 
         corrmodel=corrmodel, distance= distance,grid=grid,maxdist= maxdist,maxtime=maxtime,model=model,n=n, 
          param=param,radius=radius,sparse=sparse,taper=taper,tapsep=tapsep,type=type,copula=copula,X=X) 
    ###########
 
    bivariate = covmatrix$bivariate;   
    if(bivariate) tloc=1
    spacetime = covmatrix$spacetime; 

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
    nuisance = param[covmatrix$namesnuis]

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
      covmatrix_true =  GeoCovmatrix(coordx, coordy, coordt, coordx_dyn, corrmodel, distance, grid, maxdist, maxtime, model, n, param, 
      radius, sparse, NULL, NULL, "Standard",X)

       }
    ############
    tapmod=NULL                                                
    cmodel=corrmodel
    cdistance=distance
    corrmodel = CkCorrModel(covmatrix$corrmodel)
    distance = CheckDistance(covmatrix$distance)
    corrparam = unlist(covmatrix$param[covmatrix$namescorr])# selecting the correlation parametrs
    if(bivariate) if(!(which==1 || which==2) ) stop("which  parameter must be 1 or 2")
    pred = NULL
    varpred=varpred2=vv=vv2=NULL
    k = 0 
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

if(covmatrix$model %in% c(1,10,18,21,12,26,24,27,38,29,20,34,39,28,40,9))    ## contnuos model 
{  
## gaussian=1 
## skew gaussian=10   
## student =12  
## skewstudent =18 
## gamma=21 
## weibull=26
## loglogistic=24
## loggaussian=22 ojo
## twopieceStudentT=27
## beta=28
## twopieceGaussian=29
## twopieceTukeyh=38
## twopiecebimodal=39
## sihasin=20
## tukeygh=9
## tukey=34
## tukeyyh2=40
    ################################
    ## standard kriging  ##############
    ################################   
       if(!bivariate) {
                      mu=X%*%betas;
                      muloc=Xloc%*%betas

                      }
       if(bivariate) { mu=c(X11%*%matrix(betas1),X22%*%matrix(betas2))
                       if(!is.null(Xloc)) 
                         muloc=c(X11_loc%*%matrix(betas1),X22_loc%*%matrix(betas2))
                     }

    if((type=="Standard"||type=="standard")) {
        
        corri=double(dimat*dimat2)
         #Computing gaussian CORRELATIONS between the locations to predict and the locations observed
     #  cc=.C('Corr_c',corri=corri, as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),as.integer(corrmodel),
     #   as.integer(FALSE),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
     #   as.integer(numloc),as.integer(tloc),as.integer(covmatrix$ns),as.integer(NS),
     #   as.integer(covmatrix$numtime),as.double(corrparam),as.integer(covmatrix$spacetime),
     #   as.integer(covmatrix$bivariate),as.double(time),as.integer(distance),as.integer(which-1),
    #  as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)

  cc=dotCall64::.C64('Corr_c',
    SIGNATURE = c("double","double","double","double", "integer", #5
                   "integer","double","double","integer","integer",
                   "integer","integer","integer","integer","double",
                   "integer","integer","double","integer","integer","double"), 
   corri=corri,ccc[,1] , ccc[,2] , covmatrix$coordt ,corrmodel , #5
         0, locx , locy , covmatrix$numcoord ,numloc ,
         tloc , covmatrix$ns , NS ,covmatrix$numtime , corrparam ,
         covmatrix$spacetime ,covmatrix$bivariate , time , distance , which-1 , covmatrix$radius ,
 INTENT = c("rw","r","r","r","r",
            "r","r","r","r","r",
            "r","r","r","r","r",
            "r","r","r","r","r","r"),
        PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

        ####  transforming gaussian correlations depending on the (non)Gaussian  model
   if(bivariate){ 
                 #  adding models.....
                corri=cc$corri
                }
   else {     
       ## for each model..
        rho=cc$corri
        cc=(1-as.numeric(covmatrix$param["nugget"]))*rho
       ##################################################
        if(covmatrix$model==1)   {  #Gaussian
                            vv=as.numeric(covmatrix$param['sill'])
                            corri=cc   } 
 ############################################################    
        if(covmatrix$model==10) {    #skew gaussian
                        corr2=rho^2
                        sk=as.numeric(covmatrix$param['skew']);sk2=sk^2
                        vv=as.numeric(covmatrix$param['sill'])
                        corri=(2*sk2)*(sqrt(1-corr2) + rho*asin(rho)-1)/(pi*vv+sk2*(pi-2)) + (cc*vv)/(vv+sk2*(1-2/pi))                
                        }    
############################################################
        if(covmatrix$model==12) # student T
                         { 
                        vv=as.numeric(covmatrix$param['sill']) 
                        nu=1/as.numeric(covmatrix$param['df'])
                        if(nu<170) corri=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,rho^2))*cc)/(2*gamma(nu/2)^2)
                        else       corri=cc
                      }
############################################################
        if(covmatrix$model==40) # tukeyh2
                         { 
                           vv=as.numeric(covmatrix$param['sill']);tail1=as.numeric(covmatrix$param['tail1']);tail2=as.numeric(covmatrix$param['tail2']) 
                           hr=tail1;hl=tail2
                           x1=1-(1-cc^2)*hr
                           x2=(1-hr)^2-(cc*hr)^2
                           y1=1-(1-cc^2)*hl
                           y2=(1-hl)^2-(cc*hl)^2
                           g=1-hl-hr+(1-cc^2)*hl*hr
                           h1=sqrt(1-cc^2/(x1^2))+(cc/x1)*asin(cc/x1)
                           h2=sqrt(1-cc^2/(y1^2))+(cc/y1)*asin(cc/y1)
                           h3=sqrt(1-cc^2/(x1*y1))+sqrt(cc^2/(x1*y1))*asin(sqrt(cc^2/(x1*y1)))
                           p1=x1*h1/(2*pi*(x2)^(3/2))+cc/(4*(x2)^(3/2))
                           p2=y1*h2/(2*pi*(y2)^(3/2))+cc/(4*(y2)^(3/2))
                           p3=-(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2))+cc/(4*(g)^(3/2))
                           mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
                           vv1=0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)-(mm)^2
                           corri=(p1+p2+2*p3-mm^2)/vv1
                         }
############################################################
    if(covmatrix$model==34) # tukeyh
                         { 
                          vv=as.numeric(covmatrix$param['sill']);
                          tail=as.numeric(covmatrix$param['tail'])
                          if(tail>0){
                              aa=(1-2*tail)^(-1.5) # variance
                              corri=(-cc/((1+tail*(cc-1))*(-1+tail+tail*cc)*(1+tail*(-2+tail-tail*cc^2))^0.5))/aa
                            }
                         if(!tail) corri=cc  
                         }
if(covmatrix$model==9) # tukeygh
                         { 
                         vv=as.numeric(covmatrix$param['sill']);
                         tail=as.numeric(covmatrix$param['tail']);
                         skew=as.numeric(covmatrix$param['skew']);
                         h=tail;g=skew
                         if(!g&&!h) { corri=cc } 
                         if(g&&!h){  
                             aa=( -exp(g^2)+exp(g^2*2))*g^(-2)
                             corri= (( -exp(g^2)+exp(g^2*(1+cc)))*g^(-2))/aa} 
                         if(!g&&h){ 
                             aa=(1-2*h)^(-1.5) 
                             corri = cc/(aa*( (1-h)^2-h^2*cc^2 )^(1.5)) } 
                        if(h&&g){ # ok
                             rho=cc; rho2=cc*cc;
                             h2=h*h; g2=g*g; u=1-h;a=1+rho;
                             A1=exp(a*g2/(1-h*a));
                             A2=2*exp(0.5*g2*  (1-h*(1-rho2))  / (u*u- h2*rho2)  );
                             A3=g2*sqrt(u*u- rho2*h2)
                             kk=(exp(g2/(2*u))-1)/(g*sqrt(u));
                             cova=vv*((A1-A2+1)/A3-kk*kk)
                             vari=vv*((exp(2*g2/(1-2*h))-2*exp(g2/(2*(1-2*h)))+1)/(g2*sqrt(1-2*h))-kk*kk)
                             corri=cova/vari

                             }
                          } 
############################################################
           if(covmatrix$model==18) # skew student T
                         {
                  vv=as.numeric(covmatrix$param['sill']) 
                  nu=1/as.numeric(covmatrix$param['df']); 
                  sk=as.numeric(covmatrix$param['skew'])

                  sk2=sk*sk;l=nu/2; f=(nu-1)/2; w=sqrt(1-sk2);y=rho;
                  CorSkew=(2*sk2/(pi*w*w+sk2*(pi-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*cc/(w*w+sk2*(1-2/pi)) ;
                  corri=(pi*(nu-2)*gamma(f)^2/(2*(pi*gamma(l)^2-sk2*(nu-2)*gamma(f)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,l,y*y))
                                *((1-2*sk2/pi)*CorSkew+2*sk2/pi)-2*sk2/pi);
                                       }
         if(covmatrix$model==27) {  # two piece StudenT
             nu=as.numeric(1/covmatrix$param['df']); 
             sk=as.numeric(covmatrix$param['skew']);
             vv=as.numeric(covmatrix$param['sill'])
             corr2=cc^2;sk2=sk^2
             a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
             a2=cc*asin(cc) + (1-corr2)^(0.5)
             ll=qnorm((1-sk)/2)
             p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
             a3=3*sk2 + 2*sk + 4*p11 - 1
             KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
             corri= KK*(a1*a2*a3-4*sk2);        
                      }
############################################################
         if(covmatrix$model==38) {  # two piece Tukey h
                        tail=as.numeric(covmatrix$param['tail']);sk=as.numeric(covmatrix$param['skew']);vv=as.numeric(covmatrix$param['sill'])
                        corr2=cc^2;sk2=sk^2;
                        gg2=(1-(1-corr2)*tail)^2
                        xx=corr2/gg2
                        A=(asin(sqrt(xx))*sqrt(xx)+sqrt(1-xx))/(1-xx)^(1.5)
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
                        a3=3*sk2 + 2*sk + 4*p11 - 1
                        mm=8*sk2/(pi*(1-tail)^2); 
                        ff=(1+3*sk2)/(1-2*tail)^(1.5)
                        M=(2*(1-corr2)^(3/2))/(pi*gg2)
                        corri=  (M*A*a3-mm)/( ff- mm)      
                      }
############################################################
         if(covmatrix$model==39) {  # bimodal
                                  nu=as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew'])
                                  vv=as.numeric(covmatrix$param['sill'])
                                  delta=as.numeric(nuisance['shape'])
                                  alpha=2*(delta+1)/nu
                                  nn=2^(1-alpha/2)
                                  ll=qnorm((1-sk)/2)
                                  p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
                                  corr2=cc^2;sk2=sk^2
                                  a1=Re(hypergeo::hypergeo(-1/alpha ,-1/alpha,nu/2,corr2))
                                  a3=3*sk2 + 2*sk + 4*p11 - 1
                                  MM=(2^(2/alpha)*(gamma(nu/2 + 1/alpha))^2) 
                                  vari=2^(2/alpha)*(gamma(nu/2 + 2/alpha))*gamma(nu/2)* (1+3*sk2) - sk2*2^(2/alpha+2)*gamma(nu/2+1/alpha)^2 
                                  corri= MM*(a1*a3-4*sk2)/vari
                      }
############################################################
        if(covmatrix$model==29) {  # two piece Gaussian 
                          corr2=sqrt(1-cc^2)
                          vv=as.numeric(covmatrix$param['sill'])
                          sk=as.numeric(nuisance['skew']); sk2=sk^2
                          ll=qnorm((1-sk)/2)
                          p11=pbivnorm::pbivnorm(ll,ll, rho = rho, recycle = TRUE)
                          KK=3*sk2+2*sk+ 4*p11 - 1
                          corri=(2*((corr2 + cc*asin(cc))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )                  
                      }
############################################################
       if(covmatrix$model==28)   ##  beta case
    {
         corr2=cc^2
         shape1=as.numeric(covmatrix$param['shape1']);     
         shape2=as.numeric(covmatrix$param['shape2']);  
         cc1=0.5*(shape1+shape2)
         vv=shape1*shape2/((cc1+1)*(shape1+shape2)^2)
         idx=which(abs(corr2)>1e-10);corr22=corr2[idx]
         nu2=shape1/2;alpha2=shape2/2
         res=0;ss=0;k=0
         while(k<=100){
             p1=2*(lgamma(cc1+k)-lgamma(cc1)+lgamma(nu2+1+k)-lgamma(nu2+1))
             p2=lgamma(k+1)+(lgamma(nu2+k)-lgamma(nu2))+2*(lgamma(cc1+1+k)-lgamma(cc1+1))
             b1=p1-p2
             b2=log(hypergeo::genhypergeo(U=c(cc1+k,cc1+k,alpha2), L=c(cc1+k+1,cc1+k+1), polynomial=TRUE,maxiter=1000, z=corr22))
             b3=k*log(corr22)
             sum=exp(b1+b2+b3)
             res=res+sum
             if (all(sum<1e-6)){ break} else{ A=res}
         k=k+1
        }
         cc[idx]=A
         corri=shape1*(cc1 + 1 ) * ((1-corr2)^(cc1) *cc -1)/shape2 ## correlation
         corri[-idx]=0
    }
############################################################
############################################################
          if(covmatrix$model==21) { # gamma
                        sh=as.numeric(covmatrix$param["shape"])
                        corri=cc^2
                                }
            if(covmatrix$model==26) {  # weibull 
                        sh=as.numeric(covmatrix$param['shape'])
                        bcorr=    (gamma(1+1/sh))^2/((gamma(1+2/sh))-(gamma(1+1/sh))^2)
                        corri=bcorr*((1-cc^2)^(1+2/sh)*Re(hypergeo::hypergeo(1+1/sh,1+1/sh ,1 ,cc^2)) -1)                       
         }

          if(covmatrix$model==24) {  # loglogistic
                        sh=as.numeric(covmatrix$param['shape'])
                        corri=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                                  (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,cc^2))*
                                        Re(hypergeo::hypergeo(1/sh, 1/sh, 1,cc^2)) -1)              
         }
    }
############################################################
############################################################
if(bivariate)  { 
                          if(which==1)    vvar=covmatrix$param["sill_1"]+covmatrix$param["nugget_1"]
                          if(which==2)    vvar=covmatrix$param["sill_2"]+covmatrix$param["nugget_2"]
               }
else    {
    #####  computing mean and variances for each model
     if(covmatrix$model==1)   {vvar= vv #gaussian
                               M=0 
                               }   
     if(covmatrix$model==10)  {vvar= vv+sk^2*(1-2/pi) ## skewgaus 
                               M=sk*sqrt(2/pi) 
                               }  
     if(covmatrix$model==12)  {vvar= vv*nu/(nu-2)    ## studentT
                               M=0      
                               } 
     if(covmatrix$model==34)  {vvar= vv*(1-2*tail)^(-1.5)            ## tukey h
                               M=0      
                               }  
      if(covmatrix$model==9)  {                                      #tukeygh
                               tail2=tail*tail;skew2=skew*skew; u=1-tail;
                               mm=(exp(skew2/(2*u))-1)/(skew*sqrt(u));
                               vvar=vv* ((exp(2*skew2/(1-2*tail))-2*exp(skew2/(2*(1-2*tail)))+1)/(skew2*sqrt(1-2*tail))-mm^2)
                               M=sqrt(vv)*mm
                               } 
     if(covmatrix$model==40) {  ## tukeyh2
                               mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
                               vvar= vv* (0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)- mm^2 )   
                               M=sqrt(vv)*mm      
                             }     
     if(covmatrix$model==18)  { #skew student T
                               D1=(nu-1)*0.5; D2=nu*0.5;
                               mm=sqrt(nu)*gamma(D1)*sk/(sqrt(pi)*gamma(D2));
                               vvar=vv*(nu/(nu-2)-mm^2)
                               M=sqrt(vv)*mm
                              }
     if(covmatrix$model==27)  { # two piece studentT
                              ttemp=gamma(0.5*(nu-1))/gamma(0.5*nu)
                              vvar= vv*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*ttemp^2)
                              M= - sqrt(vv)*2*sk*sqrt(nu/pi)*ttemp
                              }
     if(covmatrix$model==28)  { #beta
                              ssup=as.numeric((covmatrix$param["max"]-covmatrix$param["min"]))
                              vvar=(ssup^2)*(shape1*shape2/((cc1+1)*(shape1+shape2)^2))
                              M=ssup*shape1/(shape1+shape2)+as.numeric(covmatrix$param["min"])
                              muloc=0;mu=0
                              } 
     if(covmatrix$model==29)  {vvar= vv*((1+3*sk2) - 8*sk2/pi ) # two piece Gaussian
                               M=-sqrt(vv)*2*sk*sqrt(2/pi)  
                              }      
     if(covmatrix$model==38)  { # two piecetukeyh
                              vvar= vv*(ff - mm)
                              M= -sqrt(vv)*2*sk*sqrt(2/pi)/(1-tail)
                            
                              } 
     if(covmatrix$model==39)  { #twopiecebimodal
                               vvar= vv*vari/(nn^(2/alpha)*gamma(nu/2)^2)
                               M=-sqrt(vv)*sk*2^(1/alpha+1)*gamma(nu/2+1/alpha)/(nn^(1/alpha)*gamma(nu*0.5))  
                              }    

     if(covmatrix$model==21)  {vvar= 2/sh }          #gamma
     if(covmatrix$model==24)  {vvar= 2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1 }  ##loglogistic
     if(covmatrix$model==26)  {vvar= gamma(1+2/sh)/gamma(1+1/sh)^2-1   }    ## weibull            
     }
########################################################################################
##### multiplying the  correlations for the variance

         CC = matrix(corri*vvar,nrow=dimat,ncol=dimat2)
#### updating mean
        if(!bivariate){
          #### additive model on the real line 
         if(covmatrix$model %in% c(10,18,29,27,38,39,28,34,9))
                                {
                                 muloc=muloc + M
                                 mu=mu +       M
                                 }
          ### multiplicative model on the positive real line
          if(covmatrix$model %in% c(21,26,24))  {emuloc=exp(muloc);emu=exp(mu) }          
         }
         else{}    
##################################################################
##########computing kriging weights##################################
##################################################################
MM=getInv(covmatrix,CC)  #compute (\Sigma^-1) %*% cc
krig_weights = t(MM$a)
##################################################################
################# simple kriging #################################
################################################################## 
if(type_krig=='Simple'||type_krig=='simple')  {  
      if(!bivariate) {  ## space and spacetime simple kringing
               
      
   #### optimal predictorss ###    
    #if(SinhAsinhtemp||Tukeyhtemp||Tukeyghtemp){
      if(SinhAsinhtemp){
          ###################################################
               if(SinhAsinhtemp) # Sinh
               {
                   
                     kk=krig_weights %*% (c(dataT))
                     pp = c(muloc)      +  sqrt(vv)* sinh( (1/tail)*(asinh(kk)+sk))     
               }
              ###################################################
            #     if(Tukeyhtemp) # Tukeyh
            #   {
            #         kk=krig_weights %*% (c(dataT))
            #         pp = c(muloc)      +  sqrt(vv)* kk*exp(tail*kk^2/2)    
            #   } 
            #  ################################################### 
            #       if(Tukeyghtemp) # Tukeygh
            #   {
            #          kk=krig_weights %*% (c(dataT))
            #         pp = c(muloc)      +  sqrt(vv)* (exp(sk*kk)-1)*exp(0.5*tail*kk^2)/sk    
            #   } 
        }
  else {
               ############################ optimal linear predictors #######################
               if(covmatrix$model %in% c(1,12,27,38,29,10,18,39,37,28,40,34,9))   ####gaussian, StudenT, two piece  skew gaussian bimodal
              {
        
                     pp = c(muloc)      +  krig_weights %*% (c(dataT)-c(mu))  
              }
          } 
               ###################################################
               #### gamma weibull loglogistic 
               if(covmatrix$model %in% c(21,24,26))
                      {       ones=rep(1,length(c(dataT)))
                              one=rep(1,length(c(muloc)))
                              pp = c(emuloc) * ( one + krig_weights %*% (c(dataT)/emu-ones) )    
                      }
               ####log gaussian   simple kriging
               if(covmatrix$model==1&&logGausstemp)   {
                rp=as.numeric(covmatrix$param['sill'])
                pp = c(muloc-0.5*rp)      +  krig_weights %*% (c(log(dataT))-c(mu-0.5*rp)) 

                QQ=diag(as.matrix(diag(covmatrix$param['sill'],dimat2) - krig_weights%*%CC))
                pp=exp(pp+QQ/2) #/exp(covmatrix$param['sill']/2)
              } 
               #pp = (c(emuloc)+covmatrix$param['sill']/2) + 
                #                          krig_weights %*% (c(dataT)-exp(c(mu)+covmatrix$param['sill']/2)) 
        }     ####simple kriging
      else  {   ## bivariate  case   cokriging
          dat = c(dataT) - as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$ns[1]), 
                            rep(covmatrix$param['mean_2'],covmatrix$ns[2])))
                      if(which==1) pp = param$mean_1 + krig_weights %*% dat
                      if(which==2) pp = param$mean_2 + krig_weights %*% dat
            } 
   #######################################         
   #### MSE COMPUTATION ##################
   #######################################

   ####### 
      if(mse) {
                #aa=Xloc-krig_weights%*%X
                #AA=chol2inv(chol(crossprod(X,(MM$bb) %*% X)))
                #bb=tcrossprod(aa%*%AA,aa)
                bb=0

BB= krig_weights%*%CC

#BB=crossprod(t(krig_weights),CC)
# Gaussian,StudentT,skew-Gaussian,two piece linear kriging     
if(covmatrix$model %in% c(1,12,27,38,29,10,18,39,28,40,34,9))  
        {vv=diag(as.matrix(diag(vvar,dimat2) - BB  + bb)) } ## simple variance  kriging predictor variance

#gamma
if(covmatrix$model %in% c(21)) 
       { vv=emuloc^2*diag(as.matrix(diag(2/covmatrix$param['shape'],dimat2)- BB + bb))}
#weibull
if(covmatrix$model %in% c(26)) 
           {vv=emuloc^2*diag(as.matrix(diag( gamma(1+2/covmatrix$param["shape"])/gamma(1+1/covmatrix$param["shape"])^2-1,dimat2)   
                                                - BB+ bb))}
#loglogistic          
if(covmatrix$model %in% c(24)) 
           {vv=emuloc^2*diag(as.matrix(diag((2*covmatrix$param['shape']*sin(pi/covmatrix$param['shape'])^2/
                       (pi*sin(2*pi/covmatrix$param['shape']))-1),dimat2) - BB + bb))}

if(covmatrix$model==1&&logGausstemp)
       {vv =    exp(muloc + covmatrix$param['sill']/2)^2 *  
                               diag(as.matrix(diag(exp(vvar),dimat2) - exp(BB+ bb))) }
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
     if(!is.null(covmatrix$tapmod)) tapmod =CkCorrModel(covmatrix$tapmod)
    else    tapmod =NULL
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
             else   vvar=covmatrix$param["sill"]
           }
        #############################################################
        if(!bivariate) cc = matrix((1-(covmatrix$param["nugget"]))*corri*(covmatrix$param["sill"]),nrow=dimat,ncol=dimat2)
        else           cc = matrix(corri,nrow=dimat,ncol=dimat2)
    cc_tap = matrix(corri_tap,nrow=dimat,ncol=dimat2)
    MM=getInv(covmatrix,cc_tap)      #compute (\Sigma^-1) %*% cc
    krig_weights_tap1 = t(MM$a)    # cc_tap%*%as.matrix(invcov)
    if(type_krig=='Simple'||type_krig=='simple')  {
    if(!bivariate) 
                { pp = c(Xloc%*%betas) + krig_weights_tap1 %*% (c(dataT)-c(X%*%betas)) }
     else           
                { 
                 dat = c(dataT) - as.numeric(c(rep(covmatrix$param['mean_1'],covmatrix$numcoord), rep(covmatrix$param['mean_2'],covmatrix$numcoord)))
                 ww1= krig_weights_tap1  %*% dat
                 if(which==1) pp = param$mean_1 + ww1  
                 if(which==2) pp = param$mean_2 + ww1    
                }    
     if(mse) {
         #aa=Xloc-krig_weights_tap1%*%X
         #AA=chol2inv(chol(crossprod(X,(MM$bb) %*% X)))
         #bb=tcrossprod(aa%*%AA,aa)
         bb=0
         vv = diag(as.matrix(diag(vvar,dimat2) - 2*krig_weights_tap1%*%cc)
                            +krig_weights_tap1%*%covmatrix_true$covmatrix%*%t(krig_weights_tap1) + bb )      ## simple variance kriging tapering predictor variance
         vv2 = diag(as.matrix(diag(vvar,dimat2) - krig_weights_tap1%*%cc_tap + bb))}
     
       if(spacetime||bivariate) {varpred=matrix(c(vv),nrow=tloc,ncol=numloc); varpred2=matrix(c(vv2),nrow=tloc,ncol=numloc);} 
        else                    {varpred=c(vv);varpred2=c(vv2)} 
     }
           if(spacetime||bivariate) pred=matrix(t(pp),nrow=tloc,ncol=numloc)
           else pred=c(pp)
    }     ##### end tapering

} #### 
####################################################################################################################################
###################### binomial  binomial negative and poisson (inflated) #####################################
####################################################################################################################################


if(covmatrix$model %in% c(2,11,14,19,30,36,16,43,44,45))
{  
     if(type=="Standard"||type=="standard") {



     mu0 = Xloc%*%betas; 
     if(!bivariate) mu  = X%*%betas 
     if(bivariate)  mu  = c(X11%*%betas1,X22%*%betas2)
     kk=0
     if(covmatrix$model==2||covmatrix$model==11) kk=min(n)
     if(covmatrix$model==19) kk=min(nloc)
     if(covmatrix$model==16||covmatrix$model==45) kk=n
## ojo que es la covarianza
if(covmatrix$model %in% c(2,11,14,16,19,30,36,43,44,45))
{
  corri=double(dimat*dimat2)

    ## Computing correlation between the locations to predict and the locations observed
    ccorr=.C('Corr_c_bin',corri=corri, as.double(ccc[,1]),as.double(ccc[,2]),as.double(covmatrix$coordt),
    as.integer(corrmodel),as.integer(FALSE),as.double(locx),as.double(locy),as.integer(covmatrix$numcoord),
    as.integer(numloc),as.integer(covmatrix$model),as.integer(tloc),
    as.integer(kk),as.integer(covmatrix$ns),as.integer(NS),as.integer(covmatrix$numtime),
    as.double(rep(c(mu),dimat2)),as.double(other_nuis),as.double(corrparam),as.integer(covmatrix$spacetime),
    as.integer(bivariate),as.double(time),as.integer(distance),as.integer(which-1),
    as.double(covmatrix$radius),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)


      

#     ccorr=dotCall64::.C64('Corr_c_bin',
 #   SIGNATURE = c("double","double","double","double", "integer","integer",  #6
 #                  "double","double","integer","integer","integer",        #5
 #                  "integer","integer","integer","integer","integer","double", #6
 #                  "double", "double","integer","integer","double","integer","integer","double"),  #8
 #  corri=corri,  ccc[,1], ccc[,2], covmatrix$coordt,corrmodel,0,
 #   locx, locy,covmatrix$numcoord,numloc,covmatrix$model,
 #   tloc,kk,covmatrix$ns,NS,covmatrix$numtime, rep(c(mu),dimat2),
 #    other_nuis, corrparam,covmatrix$spacetime, bivariate, time,distance,which-1, covmatrix$radius,
#INTENT = c("rw","r","r","r","r","r",
     #       "r","r","r","r","r",
    #        "r","r","r","r","r","r",
   #        "r","r","r","r","r","r","r","r"),
   #   PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)
 }   
    corri=ccorr$corri

 ##################inverse of cov matrix##############################################
        if(!bivariate) cc = matrix(corri,nrow=dimat,ncol=dimat2)
        else           {}

        MM=getInv(covmatrix,cc)  #compute (\Sigma^-1) %*% cc
        krig_weights = t(MM$a)
######################################################################################



       if(type_krig=='Simple'||type_krig=='simple')  {
          ##########################################################
       if(covmatrix$model==30||covmatrix$model==36){  ### poisson
        p0=exp(mu0); pmu=exp(mu) 
        #print(p0[1]);print(pmu[1]);print(c(dataT[1]))
            if(!bivariate) 
                   {  pp = c(p0) + krig_weights %*% (c(dataT)-c(pmu)) }  ## simple kriging
            else{} #todo
           if(mse)  vvar=p0  ### variance (possibly no stationary)     
          }
            ##########################################################
       if(covmatrix$model==43||covmatrix$model==44){  ### poisson  inflated
        p=pnorm(covmatrix$param['pmu'])
        p0=exp(mu0); pmu=exp(mu) 
            if(!bivariate) 
                   {  pp = (1-p)*c(p0) + krig_weights %*% (c(dataT)-(1-p)*c(pmu)) }  ## simple kriging
            else{} #todo
           if(mse)  vvar=(1-p)*p0*(1+p*p0)  ### variance (possibly no stationary)  

          }  

       ##########################################################
       if(covmatrix$model==2||covmatrix$model==11){  ### binomial
        p0=pnorm(mu0); pmu=pnorm(mu) 
            if(!bivariate) 
                       { pp = n*c(p0) + krig_weights %*% (c(dataT)-n*c(pmu)) }  ## simple kriging
            else{} #todo
           if(mse) vvar=n*p0*(1-p0)  ### variance (possibly no stationary
          } 
      ###########################################################    
           if(covmatrix$model==19){  ### binomial2
            p0=pnorm(mu0); pmu=pnorm(mu) 
            if(!bivariate)   pp = nloc*c(p0) + krig_weights %*% (c(dataT)-c(n*pmu))  ## simple kriging
            else{} #todo
          if(mse)    vvar=nloc*p0*t(1-p0) ### variance (possibly no stationary)   
          }
         ##########################################################
       if(covmatrix$model==14||covmatrix$model==16){    ###geometric or negative binomial
         p0=pnorm(mu0); pmu=pnorm(mu) 
            if(!bivariate) ## space and spacetime
            { k1=c(p0);k2=c(pmu); 
              pp = n*(1-k1)/k1 + krig_weights %*% (c(dataT)-n*(1-k2)/k2) }
            else{}   #tood
            if(mse) vvar=n*(1-k1)/k1^2   ### variance (possibly no stationary)
                
          }
        if(covmatrix$model==45){    ###inflated negative binomial
            p0=pnorm(mu0); pmu=pnorm(mu)
            p=as.numeric(pnorm(covmatrix$param['pmu']))
            if(!bivariate) ## space and spacetime
            { k1=c(p0);k2=c(pmu);    
              pp = (1-p)*n*(1-k1)/k1 + krig_weights %*% (c(dataT)-(1-p)*n*(1-k2)/k2) 
              }
            else{}   #tood
            if(mse) vvar=n*(1-k1)*(1-p)*(1+n*p*(1-k1))/k1^2
          
          }
        if(mse){
                  ##aa=Xloc-krig_weights%*%X
                  ##AA=chol2inv(chol(crossprod(X,(MM$bb) %*% X)  ))
                  ##bb=tcrossprod(aa%*%AA,aa)
                 
                  bb=0
                  vv = diag(sqrt(tcrossprod(vvar))  - krig_weights%*%cc  + bb) 
                }  
        }

    # if(type_krig=='Ordinary'||type_krig=='ordinary')  {
            #     if(!bivariate) { 
            #              betas=  solve(t(X)%*%invcov%*%X)%*%t(X)%*%invcov%*%dataT    # GLS estomator of Beta
            #             if(covmatrix$model==2||covmatrix$model==11||covmatrix$model==19)
            #              pp = nloc*c(pnorm(Xloc%*%betas)) + krig_weights %*% (c(dataT)-c(n*pnorm(X%*%betas)))   
            #              if(covmatrix$model==14){
            #              k1=c(pnorm(Xloc%*%betas));k2=c(pnorm(X%*%betas))
            #              pp = (1-k1)/k1 + krig_weights %*% (c(dataT)-(1-k2)/k2)  }
            #            }
            #      else{}     ### todo   
            #      if(mse) {ss = (Xloc-krig_weights%*%X)%*%(solve(t(X)%*%invcov%*%X)%*%t(X))%*%invcov   + krig_weights
            #               vv =  diag(as.matrix(diag(vvar,dimat2)+ krig_weights %*% t(cc) -2*t(ss)%*%cc)) }  ## ordinary kriging predictor variance
           #           }   
           
      if(spacetime||bivariate) {
            pred=matrix(t(pp),nrow=tloc,ncol=numloc);
            varpred=matrix(c(vv),nrow=tloc,ncol=numloc);
            } 
          else{
           
            pred=c(pp);varpred=c(vv)}

    }}  ##end  binary or binomial or geometric kriging  poisson
########################################################################################
########################################################################################
########################################################################################

if(tloc==1)  {c(pred);c(varpred);c(varpred2)}
    # Return the objects list:
    Kg = list(    bivariate=bivariate,
                   coordx = covmatrix$coordx,
                   coordy = covmatrix$coordy,
                   coordt = covmatrix$coordt,
                   coordx_dyn=covmatrix$coordx_dyn,
                   covmatrix=covmatrix$covmatrix,
                   corrmodel = corrmodel,
                   copula=copula,
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

Prscores=function(data,method="cholesky",matrix)   {
if(class(matrix)!="CovMat") stop("A CovMat object is needed as input\n")
varcov=matrix$covmatrix
rownames(varcov)=c();colnames(varcov)=c()
if(nrow(varcov)!=length(data)) stop("The dimension of the covariance  matrix and/or the vector data are not correct  \n")
data=c(unname(data))
if(!matrix$sparse){
decompvarcov = MatDecomp(varcov,method)
inv = MatInv(decompvarcov,method)
}

if(matrix$sparse){ 
decompvarcov = spam::chol.spam(varcov)
inv = spam::solve.spam(decompvarcov)
}
vv=diag(inv)                                                                              
###########################
nsites=length(matrix$coordx)
ntime=1
if(matrix$spacetime) ntime=length(matrix$coordt)
if(matrix$bivariate) ntime=2
dime = nsites*ntime


MM=0
param=matrix$param
namesnuis=matrix$namesnuis
nuisance <- param[namesnuis]
sel=substr(names(nuisance),1,4)=="mean"
mm=as.numeric(nuisance[sel])
MM=(matrix$X)%*%mm
########
if(matrix$model %in% c(1,35,12,34,25))  data=data-MM#Gaussian #StudentT Tukey  Logistic
if(matrix$model %in% c(10))             #SkewGauussian
                 { kk=param['skew'];
                   data=data-(MM+kk*sqrt(2/pi))}
if(matrix$model %in% c(11))     data=data-(matrix$n)*pnorm(MM) #binomial
if(matrix$model %in% c(30,36))  data=data-exp(MM) #poisson
if(matrix$model %in% c(43,44))  {p=pnorm(param['pmu']);data=data-(1-p)*exp(MM)} #poisson inlated
if(matrix$model %in% c(27)) #two piece t models
     {kk=param['skew'];
      dd=param['df'];
      ss=param['sill'];
      data=data-(MM-(2*kk*sqrt(ss*dd)*gamma((dd-1)/2))/(gamma(dd/2)*sqrt(pi)))
     }
if(matrix$model %in% c(29)) #two piece gaussian
     {kk=param['skew'];
      ss=param['sill'];
      data=data-(MM-(2*kk*sqrt(2*ss/pi)))
     }
if(matrix$model %in% c(38)) #two piece tukeyh
     {kk=param['skew'];
      ss=param['sill'];
      tt=param['tail'];
      data=data-(MM-(2*kk*sqrt(2*ss/pi)/(1-tt)))
     }
if(matrix$model %in% c(40))      #tukeyh2
     {ss=param['sill'];
      t1=param['tail1'];
      t2=param['tail2'];
      data=data-(MM+sqrt(ss)*(t1-t2)/(sqrt(2*pi)*(1-t1)*(1-t2)))
     }
if(matrix$model %in% c(16))     data=data-(matrix$n)*(1-pnorm(MM))/pnorm(MM) #binomialnegative
if(matrix$model %in% c(45))     data=data-(1-pnorm(param['pmu']))*(matrix$n)*(1-pnorm(MM))/pnorm(MM) #binomialnegative inflated
if(matrix$model %in% c(20))      #sas  
     {ss=param['sill'];
      kk=param['skew'];tt=param['tail'];
      data=data-(MM+sqrt(ss)*sinh(kk/tt)*exp(0.25)*(besselK(.25,(tt+1)/(2*tt))+besselK(.25,(1-tt)/(2*tt)))/(sqrt(8*pi)))
     } 
if(matrix$model %in% c(39))      #twopiecebimodal
     {vv=param['sill'];
      sk=param['skew'];
      nu=param['df'];
      delta=param['shape'];
      alpha=2*(delta+1)/nu
      nn=2^(1-alpha/2)
      data=data-(MM-sqrt(vv)*sk*2^(1/alpha+1)*gamma(nu/2+1/alpha)/(nn^(1/alpha)*gamma(nu*0.5))) 
     } 


#######
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
scores = list(RMSE = RMSE,
               LSCORE = LSCORE,
               CRPS = CRPS,
               MAE=MAE)
return(scores)
}

###################################################################################################
###################################################################################################