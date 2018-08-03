####################################################
### Authors: Moreno Bevilacqua, Víctor Morales Oñate.
### Email: moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Universidad de Valparaiso, Departamento de Estad?stica
### File name: Simulation.r
### Description:
### This file contains a set of procedures
### for the simulation of Gaussian random fields and
### related functions.
### Last change: 28/03/2018.
####################################################
 

# Simulate spatial and spatio-temporal random felds:
GeoSim <- function(coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl",GPU=NULL, grid=FALSE, 
     local=c(1,1),method="cholesky",model='Gaussian', n=1, param, radius=6378.388, sparse=FALSE,X=NULL)
{
####################################################################
############ internal function #####################################
####################################################################
ddim<-function(coordx,coordy,coordt)  
{
dimt=1
if(is.null(coordy))  dims=dim(coordx)[1]
else                 dims=length(coordx)*length(coordy)
if(!is.null(coordt)) dimt=length(coordt)
return(dims*dimt)
} 
forGaussparam<-function(model,param,bivariate)
{
   if(model %in% c("SkewGaussian","SkewGauss","TwoPieceGaussian","TwoPieceGauss"))  {
     if(!bivariate) param[which(names(param) %in% c("skew"))] <- NULL
     if(bivariate)  param[which(names(param) %in% c("skew_1","skew_2"))] <- NULL
                  
   }

     if(model %in% c("SkewStudentT","TwoPieceStudentT")){
     if(!bivariate) param[which(names(param) %in% c("df","skew"))] <- NULL
     if(bivariate)  param[which(names(param) %in% c("df_1","df_2","skew_1","skew_2"))] <- NULL
   }

    if(model %in% c("Tukey","SinhAsinh")){
     if(!bivariate) param[which(names(param) %in% c("skew","tail"))] <- NULL
     if(bivariate)  param[which(names(param) %in% c("skew_1","skew_2","tail_1","tail_2"))] <- NULL
   }
    if(model %in% c("Gamma","LogLogistic","Weibull"))  {
     if(!bivariate) param[which(names(param) %in% c("shape"))] <- NULL
     if(bivariate)  param[which(names(param) %in% c("shape_1","shape_2"))] <- NULL
   }  
     if(model %in% c("StudentT"))  {
     if(!bivariate) param[which(names(param) %in% c("df"))] <- NULL
     if(bivariate)  param[which(names(param) %in% c("df_1","df_2"))] <- NULL
   }  
    if(model %in% c("Gamma2"))  {
     if(!bivariate) param[which(names(param) %in% c("shape1","shape2"))] <- NULL
    # if(bivariate)  param[which(names(param) %in% c("shape1_1","shape1_2","shape2_1","shape2_2"))] <- NULL
   }    
 return(param)   
}
##############################################################################
########### for Gaussian and non Gaussian RF obtained using Gaussian RF ######
##############################################################################
   RFfct1<- function(ccov,dime,nuisance,param,simd,X,ns)
    {
        numcoord=ccov$numcoord; numtime=ccov$numtime;grid=ccov$grid;

        if(is.null(dim(X))) X=as.matrix(rep(1,numcoord*numtime))  ## in the case of no covariates
        spacetime=ccov$spacetime;bivariate=ccov$bivariate
        if(grid){
            numcoordx=ccov$numcoordx; numcoordy=ccov$numcoordy
            sim <- array(double(dime), c(numcoordx, numcoordy, numtime, 1))
                if(!bivariate) da <- as.numeric(nuisance['mean'])+simd
                else da <- c(rep(nuisance$mean_1,numcoord),rep(nuisance$mean_2,numcoord))+simd
                l=0
                for(k in 1:(numtime)) {
                    sim[,,l+1,1]=da[seq(dime * l/numtime + 1,  dime * (l + 1)/numtime)]
                    l=l+1 }
              if(!spacetime&&!bivariate) {sim <- array(sim, c(numcoordx, numcoordy))}
                else  sim <- array(sim, c(numcoordx, numcoordy, numcoord))
        }
        else{ 
                if(!bivariate) {
                               sel=substr(names(nuisance),1,4)=="mean"; num_betas=sum(sel) ;mm=NULL
                               if(num_betas==1) mm=nuisance$mean
                               if(num_betas>1)  mm=c(mm,as.numeric((nuisance[sel])))
                               sim <- X%*%mm+simd  }
                else   { 
                   if(is.null(ns))  sim <- c(rep(as.numeric(c(nuisance['mean_1'])),numcoord),
                                     rep(as.numeric(c(nuisance['mean_2'])),numcoord)) +simd 
                    else            sim <- c(rep(as.numeric(c(nuisance['mean_1'])),ns[1]),
                                     rep(as.numeric(c(nuisance['mean_2'])),ns[2])) +simd 
                          }
            if(!spacetime&&!bivariate) sim <- c(sim)
            else sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
          } 
        return(sim)
    }
####################################################################
############# END internal functions ###############################
####################################################################
    checkinput <- CkInput(coordx, coordy, coordt,coordx_dyn, corrmodel, NULL, distance, "Simulation",
    NULL, grid, NULL, NULL, NULL, model, n,  NULL, param,radius,
     NULL, NULL, NULL, "Standard", NULL, NULL, NULL,X)

    if(!is.null(checkinput$error)) stop(checkinput$error)
    spacetime_dyn=FALSE
    ################################################################################################################
    ##############################################################################################################
    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    spacetime<-CheckST(CkCorrModel(corrmodel))
    if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
  ################################################################################ 
  ################ setting parameters for each model #############################
  ################################################################################
    sel=substr(names(param),1,4)=="mean"; 
    num_betas=sum(sel)   ## number of covariates
    k=1

  
#################################
    if(model %in% c("SkewGaussian","SkewGauss","Beta",
                    "StudentT","SkewStudentT",
                    "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss",
                    "Gamma","Gamma2","Weibull",
                    "LogLogistic","Logistic")) 
       {

        if(spacetime_dyn){
          env <- new.env()
          #coords=do.call(rbind,args=c(coordx_dyn),envir = env) 
          if(is.list(X))  X=do.call(rbind,args=c(X),envir = env)}

        if(!bivariate){
           if(num_betas==1)  mm<-param$mean
           if(num_betas>1)   mm<- X%*%as.numeric((param[sel]))
           param$mean=0;if(num_betas>1) {for(i in 1:(num_betas-1)) param[[paste("mean",i,sep="")]]=0}
        if((model %in% c("SkewGaussian","SkewGauss","TwoPieceGaussian","TwoPieceGauss",
                    "StudentT","SkewStudentT","TwoPieceStudentT"))) 
        {
          vv<-param$sill;
          param$sill=1-param$nugget
        }
        if(model%in% c("SkewGaussian","SkewGauss","SkewStudentT",
               "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss")) sk<-param$skew
        }
        else {
            mm1<-param$mean_1;param$mean_1=0; 
            mm2<-param$mean_2;param$mean_2=0;mm=c(mm1,mm2)
            vv1<-param$sill_1;param$sill_1=1-param$nugget_1;
            vv2<-param$sill_2;param$sill_2=1-param$nugget_2;;vv=c(vv1,vv2)
            sk1<-param$skew_1;sk2<-param$skew_2;sk=c(sk1,sk2)
        }}
#################################
    if(model %in% c("Tukey","SinhAsinh"))  {
         if(!bivariate){
          mm<-param$mean;param$mean=0
          vv<-param$sill;param$sill=1
          sk<-param$skew; tl<-param$tail}
         else {
            mm1<-param$mean_1;param$mean_1=0; mm2<-param$mean_2;param$mean_2=0;mm=c(mm1,mm2)
            vv1<-param$sill_1;param$sill_1=1;vv2<-param$sill_2;param$sill_2=1;vv=c(vv1,vv2)
            sk1<-param$skew_1;sk2<-param$skew_2;sk=c(sk1,sk2)
            tl1<-param$tail_1;tl2<-param$tail_2;sk=c(tl1,tl2)
        }}
#################################
    if(model %in% c("Wrapped"))  {
        k=2;
        if(!bivariate){
            if(num_betas==1) mm<-2*atan(param$mean)+pi;   
            if(num_betas>1)  mm<-2*atan(X%*%as.numeric((param[sel])))+pi;
            param$mean=0    
            if(num_betas>1) {for(i in 1:(num_betas-1)) param[[paste("mean",i,sep="")]]=0}}
        else {
            mm1<-2*atan(param$mean_1)+pi;param$mean_1=0;
            mm2<-2*atan(param$mean_2)+pi;param$mean_2=0;
            mm=c(mm1,mm2)
        }} 
################################# how many random fields ################
    if(model %in% c("SkewGaussian","SkewGauss","Weibull","TwoPieceGaussian","TwoPieceGauss")) k=2    
    if(model %in% c("LogLogistic","Logistic")) k=4
    if(model %in% c("Binomial"))   k=round(n)
    if(model %in% c("Geometric","BinomialNeg")){ k=99999; if(model %in% c("Geometric")) {model="BinomialNeg";n=1}} 
    if(model %in% c("Gamma"))  k=round(param$shape)
    if(model %in% c("StudentT"))  k=round(1/param$df)+1
    if(model %in% c("SkewStudentT","TwoPieceStudentT"))  k=round(1/param$df)+2
#################################
     if(model %in% c("Gamma2")) {    
             k=round(param$shape1)
             if(!bivariate) {  mm<-param$mean;param$mean=0
             vv<-param$sill;param$sill=1 }}
     #if(model %in% c("Beta")) {  k=round(param$shape1)+round(param$shape2)       
       #  if(!bivariate) {  mm<-param$mean;param$mean=0
        #     vv<-param$sill;param$sill=1 }} 
  ################################################################################ 
  ################################################################################ 
   ns=NULL
   if(spacetime_dyn) {
        coords=NULL
        if(bivariate) coordt=c(1,2)
       coords=do.call(rbind,args=c(coordx_dyn))      
       ns=lengths(coordx_dyn)/2 
       coordx <- coords[,1]; coordy <- coords[,2]
       dime=sum(ns)
   }
   else { dime=ddim(coordx,coordy,coordt) }
   dd=array(0,dim=c(dime,1,k))    
   cumu=NULL;#s=0 # for negative binomial  case
 #########################################
     
#### computing covariance matrix of the Gaussian random field
#print(forGaussparam(model,param,bivariate)) #pay attention to the parameter
    ccov = GeoCovmatrix(coordx, coordy, coordt,coordx_dyn, corrmodel, distance, grid,NULL,NULL, "Gaussian", n, 
                forGaussparam(model,param,bivariate), radius, FALSE,NULL,NULL,"Standard",X)
    if(spacetime_dyn) ccov$numtime=1
    numcoord=ccov$numcoord;numtime=ccov$numtime;
    dime<-numcoord*numtime
    varcov<-ccov$covmat;  ######covariance matrix!!

#########################################################    
  for(i in 1:k) {  
    ss=matrix(rnorm(dime) , nrow=dime, ncol = 1)
   

   #### simulating with cholesky decomposition using GPU
    if(!is.null(GPU)&&sparse) sparse=FALSE   ### if gpu no sparse 

    if(!is.null(GPU)) {  ## todo...
                         ## here we wave to set the context!!!
                         ## setContext(id=3L)  example
                         ##varcov=gpuR::vclMatrix(varcov, type="float")
                         ##ss=gpuR::vclMatrix(ss, type="float")
                       }
    #### simulating with matrix decomposition using sparse or dense matrices
    if(sparse) {  A=spam::as.spam(ccov$covmat);
                  simd=as.numeric(spam::rmvnorm.spam(1,mu=rep(0, dime), A) )
               }
    else
    {
        decompvarcov <- MatDecomp(varcov,method)
        if(is.logical(decompvarcov)){print(" Covariance matrix is not positive definite");stop()}
        sqrtvarcov <- MatSqrt(decompvarcov,method)
       if(!is.null(GPU)) simd=(gpuR::crossprod(sqrtvarcov,ss))[]
       else simd=crossprod(sqrtvarcov,ss)
    }
    #######################################################################
    nuisance<-param[ccov$namesnuis]
    if(i==1&&(model=="SkewGaussian"||model=="SkewGauss")&&bivariate) ccov$param["pcol"]=0
    ####################################
    #####formatting simulation #########
    sim<-RFfct1(ccov,dime,nuisance,param,simd,ccov$X,ns)
    ####################################
    ####### starting cases #############
    ####################################
      if(model %in% c("Binomial", "BinomialNeg")) {    
        simdim <- dim(sim)
        sim <- as.numeric(sim>0)
        dim(sim) <- simdim
        }
    ####################################    
    if(model %in% c("Weibull","SkewGaussian","SkewGauss","Binomial",
                "Gamma","Gamma2","LogLogistic","Logistic","StudentT",
                "SkewStudentT","TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss")) {
       if(!bivariate) dd[,,i]=t(sim)
       if(bivariate)  dd[,,i]=sim[i,]
      }
     ####################################     
    if(model %in% c("BinomialNeg")){ 
                # s=s+1 ;
                 cumu=rbind(cumu,c(sim));
                 if(sum(colSums(cumu)>=n)==dime) {break;}### checking if at least n success have ben achived
               }
    }
 ####### end for #########################  
############################################
### using  gausssian random fields  in order to generate non gaussausan random fiels
###########################################
 if(model %in% c("Tukey"))   { 
     t1=1-tl;   t2=t1^2-tl^2;   sk2=sk^2;   
     tm=(exp(sk^2/(2*t1))-1)/(sk*sqrt(t1))
     if(sk) {   vvar=((exp(sk2 * 2/(1-2*tl)) - 2* exp( (sk2 *0.5)/(t1-tl^2))    +1))/(sk2*t1 - tl^2) - tm^2
                trans=(exp(sk*sim)-1)*exp(0.5*tl*sim^2)/sk
                sim=mm+sqrt(vv/vvar)*trans }
     else {
             vvar=(1-2*tl)^(-1.5)
             if(tl)   { trans=sim*exp(0.5*tl*sim^2); sim=mm+sqrt(vv/vvar)*trans }
             else  { sim=mm+sqrt(vv)*sim }
            }
      if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numcoordx, numcoordy))
        else                        sim <- array(sim, c(numcoordx, numcoordy, numtime)) 
            }
        }
    #########################################
    if (model %in% c("SinhAsinh","SinhAsinh")) 
    { trans=sinh( (1/tl)*(asinh(sim)+sk))
      sim=mm+sqrt(vv)*trans
       if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numcoordx, numcoordy))
        else                        sim <- array(sim, c(numcoordx, numcoordy, numtime)) 
            }
    }
    #######################################
    if(model %in% c("SkewGaussian","SkewGauss"))   {
        sim=mm+sk*abs(dd[,,1])+sqrt(vv)*dd[,,2]

             if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numcoordx, numcoordy))
        else                        sim <- array(sim, c(numcoordx, numcoordy, numtime)) 
            }
        }
    #######################################
    #######################################    
    #######################################
    if(model %in% c("Binomial"))   { 
                  sim[sim==1]=0;
                  sim=c(sim)
                  for(i in 1:k) sim=sim+dd[,,i] 
          if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numcoordx, numcoordy))
        else                        sim <- array(sim, c(numcoordx, numcoordy, numtime)) 
            }
    }
     #######################################
    if(model %in% c("BinomialNeg"))   {
          sim_bn=NULL
          for(p in 1:dime) sim_bn=c(sim_bn,which(cumu[,p]>0,arr.ind=T)[n]-n)
          # RE-Formatting  output:
        if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim_bn)
                else                       sim <- matrix(sim_bn, nrow=numtime, ncol=numcoord)
        }
        else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim_bn, c(numcoordx, numcoordy)) 
        else                        sim <- array(sim_bn, c(numcoordx, numcoordy, numtime)) 
        }}   


################################################        
if(model %in% c("StudentT"))   { 
     sim=NULL
     for(i in 1:(k-1))  sim=cbind(sim,dd[,,i]^2)
        aa=mm+sqrt(vv)*(c(dd[,,k])/sqrt(rowSums(sim)/(k-1)))

            if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(aa)
                else                       sim <- matrix(aa, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(aa, c(numcoordx, numcoordy))
        else                        sim <- array(aa, c(numcoordx, numcoordy, numtime)) 
            }
        }
################################################
if(model %in% c("TwoPieceGaussian","TwoPieceGauss"))   { 
        sim=dd[,,1]
        discrete=dd[,,2] 
        pp=qnorm((1-sk)/2)
        sel=(discrete<=pp);discrete[sel]=1-sk;discrete[!sel]=-1-sk;
        aa=mm+sqrt(vv)*(abs(sim)*discrete)
            if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(aa)
                else                       sim <- matrix(aa, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(aa, c(numcoordx, numcoordy))
        else                        sim <- array(aa, c(numcoordx, numcoordy, numtime)) 
            }
        }
################################################

if(model %in% c("TwoPieceStudentT"))   { 
     sim=NULL
     for(i in 1:(k-2))  sim=cbind(sim,dd[,,i]^2)

        aa=(c(dd[,,k-1])/sqrt(rowSums(sim)/(k-2)))
        pp=qnorm((1-sk)/2)
        discrete=dd[,,k] 
        sel=(discrete<=pp);discrete[sel]=1-sk;discrete[!sel]=-1-sk;
        aa=mm+sqrt(vv)*(abs(aa)*discrete)

            if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(aa)
                else                       sim <- matrix(aa, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(aa, c(numcoordx, numcoordy))
        else                        sim <- array(aa, c(numcoordx, numcoordy, numtime)) 
            }
        }
################################################
if(model %in% c("SkewStudentT"))   { 
     sim=NULL
     print(k)
     print(mm)
     print(vv)
     print(sk)
     for(i in 1:(k-2))  sim=cbind(sim,dd[,,i]^2)

        aa=0+sk*abs(dd[,,k-1])+dd[,,k]*sqrt(1-sk^2)
        sim=mm+sqrt(vv)*(aa/sqrt(rowSums(sim)/(k-2)))
            if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numcoordx, numcoordy))
        else                        sim <- array(sim, c(numcoordx, numcoordy, numtime)) 
            }
        }    

#######################################
    if(model %in% c("LogLogistic","Logistic"))   { 
      sim1=sim2=NULL
    for(i in 1:2)  sim1=cbind(sim1,dd[,,i]^2)
    for(i in 3:4)  sim2=cbind(sim2,dd[,,i]^2)
     sim1=rowSums(sim1)/2; sim2=rowSums(sim2)/2;
     ######################################################
      if(model %in% c("LogLogistic"))   
       sim=exp(mm)*(sim1/sim2)^((1/param$shape))/(gamma(1+1/param$shape)*gamma(1-1/param$shape))
    if(model %in% c("Logistic"))   
       sim=mm+log(sim1/sim2)*(param$sill)^(0.5)   
  if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numcoordx, numcoordy))
        else                        sim <- array(sim, c(numcoordx, numcoordy, numtime)) 
            }
        }

#######################################
    if(model %in% c("Gamma","Gamma2","Weibull"))   { 
      sim=NULL
    for(i in 1:k)  sim=cbind(sim,dd[,,i]^2)
     ######################################################
      if(model %in% c("Weibull"))   
                sim=exp(mm)*(rowSums(sim)/2)^(1/param$shape)/(gamma(1+1/param$shape))
      if(model %in% c("Gamma","Gamma2"))  
      { 
      #print(rowSums(sim)/k)
      sim=exp(mm)*rowSums(sim)/k      ## gamma)
      if(model %in% c("Gamma2")) 
               sim=sim+rgamma(length(sim),shape=param$shape2/2)   ##gamma2
      }
     
         if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{numcoordx=length(coordx);numcoordy=length(coordy);
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numcoordx, numcoordy))
        else                        sim <- array(sim, c(numcoordx, numcoordy, numtime)) 
            }
        }
    #######################################
    if(model %in% c("Wrapped"))   {
        if(spacetime) mm=matrix(mm,nrow=nrow(sim),ncol=ncol(sim),byrow=TRUE)
        sim=(sim+mm)%%(2*pi)
      }
    #######################################   
     if(model %in% c("LogGaussian","LogGauss"))   {
        sim=exp(sim)}      ### 
    ###########. formatting data for space time dynamic case. #########
    if(spacetime_dyn) {
                    sim_temp=list()
                    for(k in 1:length(coordt))
                       { if(k==1) {indx=1:(sum(ns[1:k]))}
                         if(k>1)    {indx=(sum(ns[1:(k-1)])+1):(sum(ns[1:k]))}
                         sim_temp[[k]]=c(sim)[indx] }
    sim=sim_temp     
    }
    #######################################
    if(ccov$bivariate)   ccov$numtime=1

    # Delete the global variables:
    # Return the objects list:
    GeoSim <- list(bivariate = bivariate,
    coordx = ccov$coordx,
    coordy = ccov$coordy,
    coordt = ccov$coordt,
    coordx_dyn =coordx_dyn,
    corrmodel = corrmodel,
    data = sim,
    distance = distance,
    grid = grid,
    model = model,
    n=n,
    numcoord = ccov$numcoord,
    numtime = ccov$numtime,
    param = ccov$param,
    radius = radius,
    randseed=.Random.seed,
    spacetime = spacetime,
    sparse=ccov$sparse,
    X=X)
#}
##############################################
    structure(c(GeoSim, call = call), class = c("GeoSim"))
}

