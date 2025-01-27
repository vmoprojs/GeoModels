####################################################
### File name: GeoSimapprox.r
####################################################

# Simulate approximate spatial and spatio-temporal random felds:
GeoSimapprox <- function(coordx, coordy=NULL,coordz=NULL,coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl",GPU=NULL, grid=FALSE,
     local=c(1,1),max.ext=1,method="TB", L=1000,model='Gaussian', parallel=FALSE,ncores=NULL,
     n=1, param, anisopars=NULL,radius=6371,X=NULL,spobj=NULL,nrep=1,progress=TRUE)
{
####################################################################
############  starting internal function ###########################
####################################################################

###############################
ddim<-function(coordx,coordy,coordz,coordt)
{
dimt=1
if(is.null(coordz))
{
if(is.null(coordy))  dims=dim(coordx)[1]
else                 dims=length(coordx)*length(coordy)
}
else   dims=length(coordx)*length(coordy)*length(coordz)
if(!is.null(coordt)) dimt=length(coordt)
return(dims*dimt)
}
###############################
RFfct1<- function(numcoord,numtime,spacetime,bivariate,dime,nuisance,simd,X,ns)
    {
        if(!bivariate) {if(is.null(dim(X))) {X=as.matrix(rep(1,numcoord*numtime))}}  ## in the case of no covariates
        if( bivariate) {if(is.null(dim(X))) {X=as.matrix(rep(1,ns[1]+ns[2]))}}
        if(!bivariate) {
                              sel=substr(names(nuisance),1,4)=="mean";
                               num_betas=sum(sel);mm=NULL
                               if(num_betas==1) {mm=nuisance$mean;
                                                 sim = X*mm+simd
                                                }
                               if(num_betas>1)  { mm=c(mm,as.numeric((nuisance[sel])));
                                                  sim = X%*%mm+simd
                                                }            
                       }
        if(bivariate)  {
                  sel1=substr(names(nuisance),1,6)=="mean_1";
                  sel2=substr(names(nuisance),1,6)=="mean_2";
                  num_betas1=sum(sel1);mm1=NULL;
                  num_betas2=sum(sel2);mm2=NULL;
                  if(num_betas1==1) mm1=nuisance$mean_1; if(num_betas1>1)  mm1=c(mm1,as.numeric((nuisance[sel1])))
                  if(num_betas2==1) mm2=nuisance$mean_2; if(num_betas2>1)  mm2=c(mm2,as.numeric((nuisance[sel2])))
                  X11=as.matrix(X[1:ns[1],]);
                  X22=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]);
                  if(is.null(ns))  {sim <- c(X11%*%mm1,X22%*%mm2) + simd }
                  else            sim <- c(rep(as.numeric(nuisance['mean_1']),ns[1]),
                                         rep(as.numeric(nuisance['mean_2']),ns[2])) + simd
                  }


            if(!spacetime&&!bivariate) sim <- c(sim)

            else sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        return(sim)
    }
#########################################################

######################################################################
############## main  internal function  ##############################
######################################################################
simu_approx=function(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel,ncores)
{
##### spatial case 


if(!spacetime){    
   ## Turning Bands
 if(method=="TB") { 
        if(bivariate)
        { param$nugget1=0;par;param$nugget2=0;
        simu=tbm2d(coords,coordt, param, corrmodel,L,bivariate)
        }
        else
        {   tau2=param$nugget
            vvv=param$sill
            ##----------###
            param$nugget=0
            param$sill=1
            ##----------###
         simu=tbm2d_uni(coords, coordt, param, corrmodel, L, bivariate,parallel,ncores)$simu
         if(is.null(simu)) {print("TB simulation method fails\n");return(simu)}

         simu=sqrt(vvv)*(sqrt(1-tau2)*simu+sqrt(tau2)*rnorm(length(simu)))
        }
           simu=c(simu)
    }     
## Vecchia
#if(method=="Vecchia"){ 
 #  if(corrmodel=="Matern") model1="matern_isotropic"
  # simu=GpGp::fast_Gp_sim(covparms=c(as.numeric(param['sill']),
   #                               as.numeric(param['scale']),
    #                              as.numeric(param['smooth']),
     #                             as.numeric(param['nugget'])), 
      #                            covfun_name = model1, coords, m = M)
       #               }
 ## Circulant embeeding           for regular grid          
  if(method=="CE") simu= c(SimCE(numxgrid,numygrid,coordx,coordy,coordz,corrmodel,param,mean.val=0, max.ext)$X)        
}
##### spacetime case ######
if(spacetime)  {
 # if(method=="Vecchia") print("ciao") 
  if(method=="CE") simu=CE_Space_Time(coords=coords,time.seq=coordt, param=param,corrmodel= corrmodel,distance=distance)
    }
return(simu)
}
####################################################################
####################################################################
############# END internal functions ###############################
####################################################################
####################################################################
    M=NULL
    if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")
    if(!(method=="TB"||method=="CE")) stop("The method of simulation is not correct\n")
    corrmodel=gsub("[[:blank:]]", "",corrmodel)
    model=gsub("[[:blank:]]", "",model)
    distance=gsub("[[:blank:]]", "",distance)
    method=gsub("[[:blank:]]", "",method)
    numxgrid=numygrid=NULL
    spacetime_dyn=FALSE


### some check for parallel #######
 coremax=parallel::detectCores()
if(is.na(coremax)||coremax==1) parallel=FALSE

if(parallel) { if(is.null(ncores)){ n.cores <- coremax - 1 }
      else
        {  if(!is.numeric(ncores)) stop("number of cores not valid\n")
           if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
           n.cores=ncores
         }   
}
else n.cores=NULL
##################################

##############################################################################
###### extracting sp object informations if necessary              ###########
##############################################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
if(!is.null(spobj)) {
   if(space||bivariate){
        a=sp2Geo(spobj,NULL); coordx=a$coords 
       if(!a$pj) {if(distance!="Chor") distance="Geod"}
    }
   if(spacetime){
        a=sp2Geo(spobj,NULL); coordx=a$coords ; coordt=a$coordt 
        if(!a$pj) {if(distance!="Chor") distance="Geod"}
     }
}
###############################################################
###############################################################

    if(grid) {  xgrid=coordx;ygrid=coordy
                numxgrid=length(xgrid);numygrid=length(ygrid) 
               if(is.null(coordz)){ 
                          coords=as.matrix(expand.grid(xgrid,ygrid))
                          }
               else     { zgrid=coordz;numzgrid=length(zgrid) 
                          coords=as.matrix(expand.grid(xgrid,ygrid,zgrid))}  
             }
    else
    {   coords=coordx
        if(!is.null(coordy)) coords=cbind(coordx,coordy,coordz)
    }         
    coords_orig=coords
    if(!(is.null(anisopars))) coords=GeoAniso(coords,c(anisopars$angle,anisopars$ratio))    
    ##############################################################################
    ##############################################################################

    if(space)    {if(method=="CE"&&!grid) stop("CE method works only for regular grid\n")
                  #if(!(corrmodel=="Matern")) stop("Not implemented for this correlation model  \n")
                 } 
    
    if(spacetime)
       {if(!(corrmodel %in% c("Matern_Matern","GenWend_GenWend","GenWend_Matern_GenWend_Matern")))
          stop("Not implemented for this correlation model  \n")
        }
    if(bivariate) {
         if(!(corrmodel=="Bi_matern")) stop("Not implemented for this correlation model  \n")}
    if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
   ################################################################################
    unname(coordt);
    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy);;unname(coordz)}
  
  ################################################################################
  ################ setting parameters for each model #############################
  ################################################################################
     if(!bivariate) {  sel=substr(names(param),1,4)=="mean";  num_betas=sum(sel) }    ## number of covariates
     else
    {  sel1=substr(names(param),1,6)=="mean_1"; num_betas1=sum(sel1)
       sel2=substr(names(param),1,6)=="mean_2"; num_betas2=sum(sel2)
       num_betas=c(num_betas1,num_betas2)
    }  

 
     if(!length(param$mean)>1){
    if( !all(names(unlist(param)) %in% c(CorrParam(corrmodel), NuisParam(model,bivariate,num_betas=num_betas))) )
       stop("only nuisance and correlation parameters must be included in param\n")
    }
    
    if(!bivariate) { if(is.null(param$sill))  param$sill=1 }
if(!bivariate) { if(is.null(param$sill)) param$sill=1}
else   { if(is.null(param$sill_1)) param$sill_1=1
         if(is.null(param$sill_2)) param$sill_2=1}


#################################
    if(model %in% c("SkewGaussian","SkewGauss","Beta",'Kumaraswamy','Kumaraswamy2','LogGaussian',#"Binomial","BinomialNeg","BinomialNegZINB",
                    "StudentT","SkewStudentT","Poisson","TwoPieceTukeyh","PoissonZIP","PoissonGamma","PoissonGammaZIP","PoissonWeibull",
                     "TwoPieceBimodal", "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss","Tukeyh","Tukeyh2","Tukeygh","SinhAsinh",
                    "Gamma","Weibull","LogLogistic","Logistic","BinomialLogistic"))
       {
        if(spacetime_dyn){env <- new.env();if(is.list(X))  X=do.call(rbind,args=c(X),envir = env)}

  if(!bivariate){
           if(num_betas==1)  mm<-param$mean
           if(num_betas>1)   { BB=param[sel];
                               if(!is.null(BB$sill)) BB$sill=NULL
                               mm= X%*%as.numeric(BB) }
           param$mean=0;if(num_betas>1) {for(i in 1:(num_betas-1)) param[[paste("mean",i,sep="")]]=0}
        

        if((model %in% c("SkewGaussian","SkewGauss","TwoPieceGaussian","Logistic",
          "TwoPieceGauss","Gamma","Weibull","LogLogistic","Poisson","PoissonZIP","Tukeyh","Tukeyh2","PoissonGamma",
          "PoissonGammaZIP","PoissonWeibull",
          'LogGaussian',"TwoPieceTukeyh","TwoPieceBimodal", "Tukeygh","SinhAsinh",
                    "StudentT","SkewStudentT","TwoPieceStudentT")))  {vv<-param$sill; param$sill=1}

        if(model%in% c("SkewGaussian","SkewGauss","SkewStudentT","TwoPieceTukeyh","TwoPieceBimodal",
               "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss"))
               { sk<-param$skew
                 if(model%in% c("TwoPieceTukeyh")) tl<-param$tail
                 if(model%in% c("TwoPieceBimodal")) bimo<-param$shape
               }
        }
        else {
           if(num_betas[1]==1) { mm1<-param$mean_1;param$mean_1=0 }
           if(num_betas[1]>1)    mm1<- X%*%as.numeric((param[sel1]))
           if(num_betas[2]==1) { mm2<-param$mean_2;param$mean_2=0 }
           if(num_betas[2]>1)    mm2<- X%*%as.numeric((param[sel2]))

            mm=c(mm1,mm2)
            vv1<-param$sill_1;param$sill_1=1-param$nugget_1;
            vv2<-param$sill_2;param$sill_2=1-param$nugget_2;vv=c(vv1,vv2)
            sk1<-param$skew_1;sk2<-param$skew_2;sk=c(sk1,sk2)
        }
      }
#################################
  if(model %in% c("Tukeygh","SinhAsinh"))  {
         if(!bivariate){
          param$mean=0
          sk<-param$skew; tl<-param$tail}
         else {
            mm1<-param$mean_1;param$mean_1=0; mm2<-param$mean_2;param$mean_2=0;mm=c(mm1,mm2)
            vv1<-param$sill_1;param$sill_1=1;vv2<-param$sill_2;param$sill_2=1;vv=c(vv1,vv2)
            sk1<-param$skew_1;sk2<-param$skew_2;sk=c(sk1,sk2)
            tl1<-param$tail_1;tl2<-param$tail_2;sk=c(tl1,tl2)
        }}

    if(model %in% c("Tukeyh"))  {
         if(!bivariate){
          param$mean=0
          tl<-param$tail}
         else {
            mm1<-param$mean_1;param$mean_1=0; mm2<-param$mean_2;param$mean_2=0;mm=c(mm1,mm2)
            vv1<-param$sill_1;param$sill_1=1;vv2<-param$sill_2;param$sill_2=1;vv=c(vv1,vv2)
            tl1<-param$tail_1;tl2<-param$tail_2;sk=c(tl1,tl2)
        }}

    if(model %in% c("Tukeyh2"))  {
         if(!bivariate){
          param$mean=0
          t1l<-param$tail1
          t2l<-param$tail2
           }
     }
#################################
   if(model %in% c("Wrapped"))  {
        if(!bivariate){
            if(num_betas==1) mm<-2*atan(param$mean)+pi;
            if(num_betas>1)  mm<-2*atan(X%*%as.numeric((param[sel])))+pi;
            param$mean=0
            if(num_betas>1) {for(i in 1:(num_betas-1)) param[[paste("mean",i,sep="")]]=0}}
        else {
            mm1<-2*atan(param$mean_1)+pi;param$mean_1=0;
            mm2<-2*atan(param$mean_2)+pi;param$mean_2=0;
            mm=c(mm1,mm2)
            if(num_betas1>1) {for(i in 1:(num_betas1-1)) param[[paste("mean_1",i,sep="")]]=0}
            if(num_betas2>1) {for(i in 1:(num_betas2-1)) param[[paste("mean_2",i,sep="")]]=0}
        }}


SIM=list()

###################################################
### starting number of replicates
###################################################
if(progress&&nrep>1){
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    pb <- progressr::progressor(along = 1:nrep)
}
for( LL in 1:nrep){
    if(progress&&nrep>1){pb(sprintf("LL=%g", LL))}
    k=1;npoi=1
################################# how many random fields ################
    if(model %in% c("SkewGaussian","LogGaussian","TwoPieceGaussian","TwoPieceTukeyh")) k=1
    if(model %in% c("Weibull","Wrapped")) k=2
    if(model %in% c("LogLogistic","Logistic")) k=4
    if(model %in% c("Binomial"))   k=max(round(n))
    if(model %in% c("BinomialLogistic"))   k=2*max(round(n))
    if(model %in% c("Geometric","BinomialNeg","BinomialNegZINB")){ k=99999;
                                                 if(model %in% c("Geometric")) {model="BinomialNeg";n=1}
                                               }
    if(model %in% c("Poisson","PoissonZIP")) {k=2;npoi=999999999}
    if(model %in% c("PoissonGamma","PoissonGammaZIP")) {k=2+2*param$shape;npoi=999999999}
    if(model %in% c("PoissonWeibull")) {k=4;npoi=999999999}
    if(model %in% c("PoissonZIP","BinomialNegZINB","PoissonGammaZIP")) {param$nugget=param$nugget1}
    if(model %in% c("Gamma"))  {if(!bivariate) k=round(param$shape)
                                #if(bivariate)  k=max(param$shape_1,param$shape_2)
                               }
    if(model %in% c("Beta"))  {k=round(param$shape1)+round(param$shape2);}
    if(model %in% c("Kumaraswamy","Kumaraswamy2"))  k=4
    if(model %in% c("StudentT"))  k=round(1/param$df)+1
    if(model %in% c("TwoPieceBimodal"))  k=round(param$df)+1
    if(model %in% c("SkewStudentT","TwoPieceStudentT"))  k=round(1/param$df)+2
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
       if(ncol(coords)==2) {coordx <- coords[,1]; coordy <- coords[,2];coordz <- NULL}
         if(ncol(coords)==3) {coordx <- coords[,1]; coordy <- coords[,2];coordz <- coords[,3]}
       dime=sum(ns)
   }
   else { dime=ddim(coordx,coordy,coordz,coordt)
          if(bivariate) {ns=c(length(coordx),length(coordx))/2}
        }
   if(!bivariate) dd=array(0,dim=c(dime,1,k))
   if(bivariate)  dd=array(0,dim=c(dime,2,k))
   cumu=NULL;#s=0 # for negative binomial  case
 #########################################
numtime=1
if(!is.null(coordt)) numtime=length(coordt)
if(spacetime_dyn) numtime=1
  numcoord=nrow(coords);
  dime<-numcoord*numtime
if(bivariate) numtime=2
#########################################################
KK=1;sel=NULL;ssp=double(dime)
if(model%in% c("SkewGaussian","StudentT","SkewStudentT","TwoPieceTukeyh",
               "TwoPieceStudentT","TwoPieceGaussian"))
{
      simD=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel,n.cores)
      if(!spacetime&&!bivariate) simDD <- c(simD)
      else simDD <- matrix(simD, nrow=numtime, ncol=numcoord,byrow=TRUE)
      param$nugget=0 #ojo
}


while(KK<=npoi) {
for(i in 1:k) {

################# here the approximated simulation  ##################################################
simd=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel,n.cores)
######################################################################################################

namesnuis<-NuisParam(model, bivariate,num_betas=num_betas,copula=NULL)
nuisance<-param[namesnuis]
if(i==1&&(model=="SkewGaussian"||model=="SkewGauss")&&bivariate) param["pcol"]=0
####################################
sim<- RFfct1(numcoord,numtime,spacetime,bivariate,
                  dime,nuisance,simd,X,ns)
     
####################################
####### starting cases #############
####################################
if(model %in% c("Binomial", "BinomialNeg","BinomialNegZINB")) {

        simdim <- dim(sim)
        sim <- as.numeric(sim>0)
        dim(sim) <- simdim
         }
####################################
if(model %in% c("Weibull","SkewGaussian","SkewGauss","Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonWeibull","PoissonZIP","PoissonGammaZIP","Beta","Kumaraswamy","Kumaraswamy2",
              "LogGaussian","TwoPieceTukeyh",
                "Gamma","LogLogistic","Logistic","StudentT",
                "SkewStudentT","TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss","TwoPieceBimodal")) {
       if(!bivariate) dd[,,i]=t(sim)
       if(bivariate)  dd[,,i]=t(sim)
     }
 ####################################
if(model %in% c("BinomialNeg","BinomialNegZINB")){
                 cumu=rbind(cumu,c(t(sim)));
                 if(sum(colSums(cumu)>=n)==dime) {break;}### ## stopping rule
               }
    }
 ####################################
  if(model %in% c("poisson","Poisson","PoissonZIP"))   {
   pois1=0.5*(dd[,,1]^2+dd[,,2]^2)
   ssp=ssp+c(pois1)
   sel=rbind(sel,ssp<=c(exp(mm)))
   if(sum(apply(sel,2,prod))==0) break  ## stopping rule
 }
  ####################################
 if(model %in% c("PoissonGamma","PoissonGammaZIP"))   {
   if(KK==1){sim3=NULL;for(i in 3:k)  {sim3=cbind(sim3,dd[,,i]^2)}}
   #################################
   pois1=0.5*(dd[,,1]^2+dd[,,2]^2)
   ssp=ssp+c(pois1)
   sel=rbind(sel,ssp<=c(exp(mm)*rowSums(sim3)/(k-2)))
   if(sum(apply(sel,2,prod))==0) break  ## stopping rule
 }
if(model %in% c("PoissonWeibull"))   {
  if(KK==1){sim3=NULL;for(i in 3:k){sim3=cbind(sim3,dd[,,i]^2)}}
   #################################
   pois1=0.5*(dd[,,1]^2+dd[,,2]^2)
   ssp=ssp+c(pois1)
   sel=rbind(sel,ssp<=c(exp(mm)*(rowSums(sim3)/2)^(1/param$shape)/(gamma(1+1/param$shape))))
   if(sum(apply(sel,2,prod))==0) break  ## stopping rule
 }

 KK=KK+1
}
####### end for #########################
 ###############################################################################################
 #### simulation for discrete random field based on indipendent copies  of GRF ######
 ###############################################################################################
 if(model %in% c("Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonWeibull","PoissonGammaZIP","PoissonZIP","BinomialNeg","BinomialNegZINB"))   {
   if(model %in% c("poisson","Poisson","PoissonGamma","PoissonWeibull"))   {sim=colSums(sel);byrow=TRUE}
    if(model %in% c("PoissonZIP"))   {
      a=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel,n.cores)
     ###
      a[a<as.numeric(param$pmu)]=0;a[a!=0]=1
      sim=a*colSums(sel);
      byrow=TRUE
      }
########################################
   if(model %in% c("Binomial"))   {
                  dd1=length(dd[,,1])
                  if(length(n)==1) NN=rep(n,dd1)
                  else NN=n
                  bb=NULL; for(i in 1:k) bb=rbind(bb,dd[,,i])
                  AA=NULL; for(i in 1:dd1) AA=cbind(AA,c(rep(1,NN[i]),rep(0,k-NN[i])))
                  sim=bb*AA
                  sim=apply(sim,2,sum)
                  byrow=TRUE }
########################################
if(model %in% c("BinomialLogistic"))   {
                  dd1=length(dd[,,1])
                  if(length(n)==1) NN=rep(n,dd1)
                  else NN=n
                  bb=NULL;i=1;
                         while(i<=k) 
                          {  
                             ee=0.5*rowSums(cbind(dd[,,i]^2,dd[,,(i+1)]^2)) # exp RF
                             ss=mm+log(exp(ee)-1) # transformation
                             bb=rbind(bb,as.numeric(c(ss)>0))  
                             i=i+2
                          }
                  AA=NULL; for(i in 1:dd1) AA=cbind(AA,c(rep(1,NN[i]),rep(0,k/2-NN[i])))
                  sim=bb*AA
                  sim=apply(sim,2,sum)
                  byrow=TRUE }
#######################################
   if(model %in% c("BinomialNeg"))   {
          sim=NULL
          for(p in 1:dime) sim=c(sim,which(cumu[,p]>0,arr.ind=T)[n]-n)

          byrow=TRUE
          }
  if(model %in% c("BinomialNegZINB"))   {
           sim=NULL
          for(p in 1:dime) sim=c(sim,which(cumu[,p]>0,arr.ind=T)[n]-n)
      a=simu_approx(numxgrid,numygrid,coordx,coordy,coordz,coords,coordt,method,corrmodel,param,M,L,bivariate,spacetime,parallel,n.cores)
     ###
          a[a<as.numeric(param$pmu)]=0;a[a!=0]=1
          sim=a*sim
          byrow=TRUE
          }
#############################################
############### formatting data #############
#############################################

    if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=byrow)
        }
    else{
     
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numxgrid,numygrid))
        else                        sim <- array(sim, c(numxgrid,numygrid, numtime))
            }
}
#########################################################################################################
#### simulation for continuos random field  (on the real line) based on indipendent copies  of GRF ######
#########################################################################################################

if(model %in% c("SkewGaussian","SkewGauss","SkewStudentT","StudentT","TwoPieceGaussian","TwoPieceGauss",
  "TwoPieceTukeyh","TwoPieceBimodal","TwoPieceStudentT"))   {

if(model %in% c("SkewGaussian","SkewGauss"))   {
         if(!bivariate) aa=mm+sk*c(abs(dd[,,1]))+sqrt(vv)*c(t(simDD))
        if(bivariate)  {aa=cbind(mm[1]+sk[1]*abs(dd[,,1][,1])+sqrt(vv[1])*dd[,,2][,1],
                                  mm[2]+sk[2]*abs(dd[,,1][,2])+sqrt(vv[2])*dd[,,2][,2])}
        }
################################################
if(model %in% c("SkewStudentT"))   {
     sim=NULL
     for(i in 1:(k-2))  sim=cbind(sim,dd[,,i]^2)
        bb= sk*abs(dd[,,k-1])+sqrt(1-sk^2)*t(simDD)
        aa=mm+sqrt(vv)*(bb/sqrt(rowSums(sim)/(k-2)))
        }
################################################
if(model %in% c("StudentT"))   {
     sim=NULL
     for(i in 1:(k-1))  sim=cbind(sim,dd[,,i]^2)
        aa=mm+sqrt(vv)*(c(t(simDD))/sqrt(rowSums(sim)/(k-1)))
        }
################################################
if(model %in% c("TwoPieceGaussian"))   {
       # sim=dd[,,1]
        sim=t(simDD)
        discrete=dd[,,1]
        pp=qnorm((1-sk)/2)
        sel=(discrete<=pp);discrete[sel]=1-sk;discrete[!sel]=-1-sk;
        aa=mm+c(sqrt(vv)*(abs(sim)*discrete))
       
        }
################################################
if(model %in% c("TwoPieceTukeyh"))   {
        #sim=dd[,,1]
        sim=t(simDD)
        sim=sim*exp(tl*sim^2/2)
        discrete=dd[,,1]
        pp=qnorm((1-sk)/2)
        sel=(discrete<=pp);discrete[sel]=1-sk;discrete[!sel]=-1-sk;
        aa=mm+c(sqrt(vv)*(abs(sim)*discrete))
        }
################################################
if(model %in% c("TwoPieceBimodal"))   {
     sim=NULL
     for(i in 1:(k-1))  sim=cbind(sim,dd[,,i]^2)
        alpha=2*(bimo+1)/(k-1)
        sim=rowSums(sim)/2^(1-alpha/2);
        pp=qnorm((1-sk)/2)
        discrete=dd[,,k]
        sel=(discrete<=pp);discrete[sel]=1-sk;discrete[!sel]=-1-sk;
        aa=mm+c(sqrt(vv)*(sim)^(1/alpha)*discrete)
        #aa=mm+sqrt(vv)*(sim)^(1/bimo)*discrete
        }
################################################
if(model %in% c("TwoPieceStudentT"))   {
     sim=NULL
     for(i in 1:(k-2))  sim=cbind(sim,dd[,,i]^2)

        #aa=(c(dd[,,k-1])/sqrt(rowSums(sim)/(k-2)))
        aa=(c(t(simDD))/sqrt(rowSums(sim)/(k-2)))
        pp=qnorm((1-sk)/2)
        discrete=dd[,,k]
        sel=(discrete<=pp);discrete[sel]=1-sk;discrete[!sel]=-1-sk;
        aa=mm+c(sqrt(vv)*(abs(aa)*discrete))
        }
#############################################
############### formatting data #############
#############################################
    if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(aa)
                else                       sim <- matrix(aa, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{
        if(!spacetime&&!bivariate)  sim <- array(aa, c(numxgrid,numygrid))
        else                        sim <- array(aa, c(numxgrid,numygrid, numtime))
            }
}
#########################################################################################################
#### simulation for continuos random field  (on the positive real line) based on indipendent copies  of GRF ######
#########################################################################################################
if(model %in% c("LogLogistic","Logistic"))   {
      sim1=sim2=NULL
    for(i in 1:2)  sim1=cbind(sim1,dd[,,i]^2)
    for(i in 3:4)  sim2=cbind(sim2,dd[,,i]^2)
     sim1=rowSums(sim1)/2; sim2=rowSums(sim2)/2;
     ######################################################
      if(model %in% c("LogLogistic"))
       sim=exp(mm)*(sim1/sim2)^((1/param$shape))/(gamma(1+1/param$shape)*gamma(1-1/param$shape))
      # sim=exp(mm)*(exp(sim1)-1)^((1/param$shape))/(gamma(1+1/param$shape)*gamma(1-1/param$shape))
    if(model %in% c("Logistic"))
       {
      sim=mm+log(sim1/sim2)      *(vv)^(0.5)
      #sim=mm+log(exp(sim2)-1)*(vv)^(0.5)
     }


  if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numxgrid,numygrid))
        else                        sim <- array(sim, c(numxgrid,numygrid, numtime))
            }
}

#######################################
if(model %in% c("Gamma","Weibull"))   {
      sim=sim1=sim2=NULL;
    if(!bivariate) for(i in 1:k)  sim=cbind(sim,dd[,,i]^2)
    if(bivariate)  {for(i in 1:k)  sim1=cbind(sim1,dd[,,i][,1]^2)
                    for(i in 1:k)  sim2=cbind(sim2,dd[,,i][,2]^2)
                   }
######################################################
      if(model %in% c("Weibull"))
           {
             if(!bivariate)   sim=exp(mm)*(rowSums(sim)/2)^(1/param$shape)/(gamma(1+1/param$shape))
             if(bivariate)    sim=cbind(
                                  exp(mm[1])*(rowSums(sim1)/2)^(1/param$shape_1)/(gamma(1+1/param$shape_1)),
                                  exp(mm[2])*(rowSums(sim2)/2)^(1/param$shape_2)/(gamma(1+1/param$shape_2)))
           }
      if(model %in% c("Gamma"))
      {
      if(!bivariate) sim=exp(mm)*rowSums(sim)/k

      if(bivariate){
        if(param$shape_1==param$shape_2){
                  sim=cbind(exp(mm[1])*rowSums(sim1)/param$shape_1,
                               exp(mm[2])*rowSums(sim2)/param$shape_2)}

        if(param$shape_1>param$shape_2){
                  aa=0
                  for(cc in 1:(param$shape_2)) aa=aa+sim2[,cc]
                  sim=cbind(exp(mm[1])*rowSums(sim1)/param$shape_1,
                            exp(mm[2])* aa/param$shape_2)}
       if(param$shape_1<param$shape_2){
                  aa=0
                  for(cc in 1:(param$shape_1)) aa=aa+sim1[,cc]
                  sim=cbind( exp(mm[1])* aa/param$shape_1,
                             exp(mm[2])*rowSums(sim2)/param$shape_2)}

        }
      }
#############################################
############### formatting data #############
#############################################
         if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numxgrid,numygrid))
        else                        sim <- array(sim, c(numxgrid,numygrid, numtime))
            }
}
#########################################################################################################
#### simulation for continuos random field  based  on a compact support based on indipendent copies  of GRF ######
#########################################################################################################
if(model %in% c("Beta","Kumaraswamy","Kumaraswamy2"))   {
     sim1=NULL;sim2=NULL
      i=1
    if(model=="Beta")
    {
    while(i<=round(param$shape1))  {sim1=cbind(sim1,dd[,,i]^2);i=i+1}
    while(i<=(round(param$shape1)+round(param$shape2)))  {sim2=cbind(sim2,dd[,,i]^2);i=i+1}
    aa=rowSums(sim1)
    #sim=aa/(aa+rowSums(sim2))
   sim=param$min + (param$max-param$min)*aa/(aa+rowSums(sim2))
    }
     if(model=="Kumaraswamy"||model=="Kumaraswamy2")
    {
    while(i<=2)  {sim1=cbind(sim1,dd[,,i]^2);i=i+1}
    while(i<=4)  {sim2=cbind(sim2,dd[,,i]^2);i=i+1}
    aa=rowSums(sim1)
    sim=aa/(aa+rowSums(sim2))
   # sim=( (1-(1-sim)^(1/param$shape1))^(1/param$shape2) )
    if(model=="Kumaraswamy")
      sim=param$min + (param$max-param$min)*( (1-(sim)^(1/param$shape1))^(1/param$shape2) )
    if(model=="Kumaraswamy2")
     {
      aa=log(1-((1+exp(-mm))^(-param$shape2)))/log(0.5)
      sim=param$min + (param$max-param$min)*( (1-(sim)^aa)^(1/param$shape2) )
     }
    }
         if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
        }
         else{
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numxgrid,numygrid))
        else                        sim <- array(sim, c(numxgrid,numygrid, numtime))
            }
        }
    #######################################
    if(model %in% c("Wrapped"))   {
        if(spacetime) mm=matrix(mm,nrow=nrow(sim),ncol=ncol(sim),byrow=TRUE)
        sim=(sim+mm)%%(2*pi)
      }
 ###########################################################
 #### simulation based on a transformation of ONE standard (bivariate) GRF ######
 ###########################################################
if(model %in% c("Gaussian","LogGaussian","LogGauss","Tukeygh","Tukeyh","Tukeyh2","SinhAsinh"))
{

  if(model %in% c("Gaussian")) {sim=c(sim);byrow=FALSE}
  if(model %in% c("LogGaussian","LogGauss"))   {
        sim=c(t(sim))
        sim=exp(mm) *  (exp(sqrt(vv)*sim)/(exp( vv/2))) ## note the parametrization
        byrow=TRUE
        }
#################################################################################
 if(model %in% c("Tukeygh"))   {
     sim=c(t(sim))
     if(!sk && !tl) sim= mm+sqrt(vv)* sim
     if(!sk && tl)  sim= mm+sqrt(vv)* sim*exp(tl*sim^2/2)
     if(!tl && sk)  sim= mm+sqrt(vv)* (exp(sk*sim)-1)/sk
     if(tl&&sk)     sim= mm+sqrt(vv)* (exp(sk*sim)-1)*exp(0.5*tl*sim^2)/sk
     byrow=TRUE
    }
##############################################################################
  if(model %in% c("Tukeyh"))   {
     sim=c(t(sim))
     if(!tl) sim= mm+sqrt(vv)*sim
     if(tl)  sim= mm+sqrt(vv)*sim*exp(tl*sim^2/2)
     byrow=TRUE
   }
  if(model %in% c("Tukeyh2"))   {
       sim=c(t(sim))
       sel=sim>0
       bb=sim*exp(t1l*sim^2/2)*as.numeric(sel);  bb[bb==0]=1
       aa=sim*exp(t2l*sim^2/2)*as.numeric(!sel); aa[aa==0]=1
       sim= mm+sqrt(vv)*(aa*bb)
      byrow=TRUE
   }
#########################################
  if (model %in% c("SinhAsinh"))
  {sim=c(t(sim));sim=mm+sqrt(vv)*sinh( (1/tl)*(asinh(sim)+sk));byrow=TRUE
    }
 ### formatting data
  if(!grid)  {
                if(!spacetime&&!bivariate) sim <- c(sim)
                else                       sim <- matrix(sim, nrow=numtime,
                                                          ncol=numcoord,byrow=byrow)
        }
         else{
        if(!spacetime&&!bivariate)  sim <- array(sim, c(numxgrid,numygrid))
        else                        sim <- array(sim, c(numxgrid,numygrid, numtime))
            }

}


##################################################################
###########. formatting data for space time dynamic case. #########
    if(spacetime_dyn) {
                    sim_temp=list()
                    for(k in 1:length(coordt))
                       { if(k==1) {indx=1:(sum(ns[1:k]))}
                         if(k>1)    {indx=(sum(ns[1:(k-1)])+1):(sum(ns[1:k]))}
                         sim_temp[[k]]=c(sim)[indx] }
                      sim=sim_temp
    }
##################################################################
##################################################################
##################################################################

    SIM[[LL]]=sim
} 
######################################################################
##########  end of replicates ########################################
######################################################################

if(nrep==1) SIM=SIM[[1]] 

    #######################################
    if(bivariate)   numtime=1
    # Delete the global variables:
    # Return the objects list:
    GeoSim <- list(bivariate = bivariate,
    coordx = coords_orig[,1],
    coordy = coords_orig[,2],
    coordz = coordz,
    coordt = coordt,
    coordx_dyn =coordx_dyn,
    corrmodel = corrmodel,
    data = SIM,
    distance = distance,
    grid = grid,
    model = model,
    method=method,
    n=n,
    numcoord = numcoord,
    numtime = numtime,
    param = param,
    radius = radius,
    nrep=nrep,
    spacetime = spacetime,
    X=X)
#}
##############################################
    structure(c(GeoSim, call = call), class = c("GeoSim"))
}
