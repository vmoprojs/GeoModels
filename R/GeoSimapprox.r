####################################################
### File name: GeoSimapprox.r
####################################################


# Simulate approximate spatial and spatio-temporal random felds:
GeoSimapprox <- function(coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl",GPU=NULL, grid=FALSE,
     local=c(1,1),method="Vecchia",M=30, L=NULL,model='Gaussian', n=1, param, radius=6371,X=NULL)
{
####################################################################
############ internal function #####################################
####################################################################
    RFfct1<- function(numcoord,numtime,grid,spacetime,bivariate,
          dime,nuisance,simd,X,ns)
    {
        if(!bivariate) {if(is.null(dim(X))) {X=as.matrix(rep(1,numcoord*numtime))}}  ## in the case of no covariates
        if( bivariate) {if(is.null(dim(X))) {X=as.matrix(rep(1,ns[1]+ns[2]))}}


        if(!bivariate) {
                               sel=substr(names(nuisance),1,4)=="mean";
                               num_betas=sum(sel);mm=NULL
                               if(num_betas==1) mm=nuisance$mean
                               if(num_betas>1)  mm=c(mm,as.numeric((nuisance[sel])))
                               
                               sim <- X%*%mm+simd
                              }
                if(bivariate)  {
                  sel1=substr(names(nuisance),1,6)=="mean_1";
                  sel2=substr(names(nuisance),1,6)=="mean_2";
                  num_betas1=sum(sel1);mm1=NULL;
                  num_betas2=sum(sel2);mm2=NULL;

                   if(num_betas1==1) mm1=nuisance$mean_1
                   if(num_betas1>1)  mm1=c(mm1,as.numeric((nuisance[sel1])))
                   if(num_betas2==1) mm2=nuisance$mean_2
                   if(num_betas2>1)  mm2=c(mm2,as.numeric((nuisance[sel2])))


                  X11=as.matrix(X[1:ns[1],]);
                  X22=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]);
                   if(is.null(ns))  {sim <- c(X11%*%mm1,
                                              X22%*%mm2) + simd }
                  else            sim <- c(rep(as.numeric(nuisance['mean_1']),ns[1]),
                                         rep(as.numeric(nuisance['mean_2']),ns[2])) + simd
                  }

            if(!spacetime&&!bivariate) sim <- c(sim)
            else sim <- matrix(sim, nrow=numtime, ncol=numcoord,byrow=TRUE)
       #   }

        return(sim)
    }
###############################
ddim<-function(coordx,coordy,coordt)
{
dimt=1
if(is.null(coordy))  dims=dim(coordx)[1]
else                 dims=length(coordx)*length(coordy)
if(!is.null(coordt)) dimt=length(coordt)
return(dims*dimt)
}


###############################
tbm2d <- function(coord, a0, nu0, L , model){
  # Preparing parameters to use
  # =============================================

  if(model == "Matern") {model=0}

  parametersg <- list("C" = 1, "a" = a0, "nu1" = nu0, "nu2" = c())
  parameters=parametersg

  n <- dim(coord)[1]
  sequen <- c(seq(0,n-0.5, by = ceiling(1e6/2)),n)
  index <- cbind(1:(2*L))
  coord_n <- coord[(sequen[1]+1):sequen[2], ]
  m <- dim(coord_n)[1]
      
  # Generate random frequencies and random phases
  # =============================================
  S <- ceiling(1e7*runif(3))
  set.seed(S[1])
  G <- matrix(rgamma(2*L,nu0,scale=1),2*L,2)   ## matriz (2L x 2)
  set.seed(S[2])
  u <- matrix(rnorm(4*L),2*L,2)/sqrt(G*2)/a0/(2*pi) ## matriz (2L x 2)
  set.seed(S[3])
  phi <- 2*pi*runif(2*L) ## vector (2L x 1)
  
  # Spectral density using C code
  # ================================
  
 
  f <- .C("spectral_density", L=as.integer(2*L),model=as.integer(model),p=as.integer(length(parameters$a)),
                  matrix = as.double(u), matrix_out =as.double(rep(0, 2*length(u)*length(parameters$a))),
                  C=as.double(parameters$C), a = as.double(parameters$a), nu1 = as.double(parameters$nu1),
                  Cg=as.double(parametersg$C), ag = as.double(parametersg$a), nu1g = as.double(parametersg$nu1),
                  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)

  aux_vec <- matrix(f$matrix_out,2,0.5*length(f$matrix_out))
  
  # Eigenvalues/vectors preparation
  # ================================
  
  matrix_to_eigen <- lapply(seq(1,ncol(aux_vec),2), function(i){as.matrix(cbind(aux_vec[,i],aux_vec[,i+1]),2,2)})
  
  VD <- lapply(matrix_to_eigen ,eigen)
  aux_VD <- as.matrix(unlist(VD))+1e-100
  
  vector_to_get_eigen_vals <- rep(c(rep(1,2),rep(0,4)), times = length(VD))
  vector_to_get_eigen_vec <- rep(c(rep(0,2),rep(1,4)), times = length(VD))
  
  eigen_vals <- aux_VD[subset(aux_VD*vector_to_get_eigen_vals != 0)] 
  eigen_vec <- aux_VD[subset(aux_VD*vector_to_get_eigen_vec != 0)]
  
  neg <- which(eigen_vals<0)
  if(length(neg)>0){eigen_vals[neg] = 0}
  
  lambda2 <- as.double(sqrt(eigen_vals[seq(2,length(eigen_vals),2)]))
  lambda1 <- as.double(sqrt(eigen_vals[seq(1,length(eigen_vals)-1,2)]))
  
  Mvec <- matrix(eigen_vec,nrow = length(eigen_vec), ncol = 1)
  
  vec11 <- as.double(Mvec[seq(1,length(eigen_vec),4)])
  vec21 <- as.double(Mvec[seq(2,length(eigen_vec),4)])
  vec12 <- as.double(Mvec[seq(3,length(eigen_vec),4)])
  vec22 <- as.double(Mvec[seq(4,length(eigen_vec),4)])

  tempc <- .C("matrix_temp", N = as.integer(length(lambda1)), matrix = as.double(rep(0, 4*length(lambda1))), 
               l1 = lambda1, l2 = lambda2 ,                                           # Valores propios
               v11 = vec11 , v12 = vec12 , v21 = vec21, v22 = vec22,
               PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE )                 # Vectores propios

  A <- matrix(tempc$matrix,nrow=2, ncol=(0.5*length(tempc$matrix)))
  vts <- .C("vector_to_select", N = as.integer(2*L), matrix = as.double(rep(0, (2*L)+1)),
     PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)

  AMatrix <- A[,head( vts$matrix, -1)]
  # Simulation at target locations with c code
  # ==========================================
  
  simu <- .C("simu_on_coords", Ndim = as.integer(dim(coord)[1]) ,Mcoords = as.integer(dim(coord)[1]),
                Mu = as.integer(dim(u[index,])[1]) , coords = as.double(coord_n),
                amatrix = as.double(AMatrix[,index]),
                matrix_phi = as.double(phi[index]), matrix_u = as.double(u[index,]),
                matrix_out = as.double(rep(0, 2*dim(coord)[1])),Mcoords,
                 PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
  
  simu <- matrix(simu$matrix_out, m, 2)/sqrt(L)
  
  return(simu)
}



simu_approx=function(coords,coordt,method,corrmodel,param,M,L)
{
## Turning Bands
if(method=="TB")    { 
    simu=sqrt(as.numeric(param['sill']))*tbm2d(coords, as.numeric(param['scale']), 
                                                       as.numeric(param['smooth']), 
                                                       L=L, 
                                                       model = corrmodel)
                   }     
## Vecchia
if(method=="Vecchia"){ 
if(corrmodel=="Matern") model1="matern_isotropic"
simu=GpGp::fast_Gp_sim(covparms=c(as.numeric(param['sill']),
                                  as.numeric(param['scale']),
                                  as.numeric(param['smooth']),
                                  as.numeric(param['nugget'])), 
                                  covfun_name = model1, coords, m = M)
}
return(simu)
}

####################################################################
############# END internal functions ###############################
####################################################################


    if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")
    if(!(method=="Vecchia"||method=="TB")) stop("The method of simulation is not correct\n")
    corrmodel=gsub("[[:blank:]]", "",corrmodel)
    model=gsub("[[:blank:]]", "",model)
    distance=gsub("[[:blank:]]", "",distance)
    method=gsub("[[:blank:]]", "",method)

    if(grid) { xgrid=coordx;ygrid=coordy;
               numxgrid=length(xgrid);numygrid=length(ygrid) 
               coords=as.matrix(expand.grid(xgrid,ygrid))
             }
    else
    {   coords=coordx
        if(!is.null(coordy)) coords=cbind(coordx,coordy)
    }         
         
    spacetime_dyn=FALSE
    ##############################################################################
    ##############################################################################
    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    spacetime<-CheckST(CkCorrModel(corrmodel))
    if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
   ################################################################################
    unname(coordt);
    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy)}
  ################################################################################
  ################ setting parameters for each model #############################
  ################################################################################
     if(!bivariate)
    {  sel=substr(names(param),1,4)=="mean";
       num_betas=sum(sel)   ## number of covariates
    }
   if(bivariate)
    {  sel1=substr(names(param),1,6)=="mean_1";
       num_betas1=sum(sel1)
       sel2=substr(names(param),1,6)=="mean_2";
       num_betas2=sum(sel2)
     num_betas=c(num_betas1,num_betas2)
    }

    k=1
#################################
    if(model %in% c("SkewGaussian","SkewGauss","Beta",'Kumaraswamy','Kumaraswamy2','LogGaussian',#"Binomial","BinomialNeg","BinomialNegZINB",
                    "StudentT","SkewStudentT","Poisson","TwoPieceTukeyh","PoissonZIP","PoissonGamma","PoissonWeibull",
                     "TwoPieceBimodal", "TwoPieceStudentT","TwoPieceGaussian","TwoPieceGauss","Tukeyh","Tukeyh2","Tukeygh","SinhAsinh",
                    "Gamma","Weibull","LogLogistic","Logistic","BinomialLogistic"))
       {
        if(spacetime_dyn){
                       env <- new.env()
                       #coords=do.call(rbind,args=c(coordx_dyn),envir = env)
                       if(is.list(X))  X=do.call(rbind,args=c(X),envir = env)
                     }

  if(!bivariate){
  
           if(num_betas==1)  mm<-param$mean
           if(num_betas>1)   mm<- X%*%as.numeric((param[sel]))
           param$mean=0;if(num_betas>1) {for(i in 1:(num_betas-1)) param[[paste("mean",i,sep="")]]=0}


        if((model %in% c("SkewGaussian","SkewGauss","TwoPieceGaussian","Logistic",
          "TwoPieceGauss","Gamma","Weibull","LogLogistic","Poisson","PoissonZIP","Tukeyh","Tukeyh2","PoissonGamma","PoissonWeibull",
          'LogGaussian',"TwoPieceTukeyh","TwoPieceBimodal", "Tukeygh","SinhAsinh",
                    "StudentT","SkewStudentT","TwoPieceStudentT","Gaussian")))   ##
        {
          vv<-param$sill; param$sill=1
        }
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
            vv2<-param$sill_2;param$sill_2=1-param$nugget_2;;vv=c(vv1,vv2)
            sk1<-param$skew_1;sk2<-param$skew_2;sk=c(sk1,sk2)
        }}
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
            if(num_betas1>1) {for(i in 1:(num_betas1-1)) param[[paste("mean_1",i,sep="")]]=0}
            if(num_betas2>1) {for(i in 1:(num_betas2-1)) param[[paste("mean_2",i,sep="")]]=0}
        }}

    npoi=1
################################# how many random fields ################
    if(model %in% c("SkewGaussian","LogGaussian","TwoPieceGaussian","TwoPieceTukeyh")) k=1
    if(model %in% c("Weibull")) k=2
    if(model %in% c("LogLogistic","Logistic")) k=4
    if(model %in% c("Binomial"))   k=max(round(n))
    if(model %in% c("BinomialLogistic"))   k=2*max(round(n))
    if(model %in% c("Geometric","BinomialNeg","BinomialNegZINB")){ k=99999;
                                                 if(model %in% c("Geometric")) {model="BinomialNeg";n=1}
                                               }
    if(model %in% c("Poisson","PoissonZIP")) {k=2;npoi=999999999}
    if(model %in% c("PoissonGamma")) {k=2+2*round(param$shape);npoi=999999999}
    if(model %in% c("PoissonWeibull")) {k=4;npoi=999999999}
    if(model %in% c("PoissonZIP","BinomialNegZINB")) {param$nugget=param$nugget1}
    if(model %in% c("Gamma"))  {
                             if(!bivariate) k=round(param$shape)
                             if(bivariate)  k=max(param$shape_1,param$shape_2)
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
       coordx <- coords[,1]; coordy <- coords[,2]
       dime=sum(ns)
   }
   else { dime=ddim(coordx,coordy,coordt)
          if(bivariate) {ns=c(length(coordx),length(coordx))/2}
        }
   
   if(!bivariate) dd=array(0,dim=c(dime,1,k))
   if(bivariate)  dd=array(0,dim=c(dime,2,k))
   cumu=NULL;#s=0 # for negative binomial  case
 #########################################

#### 



numtime=1
if(!is.null(coordt)) numtime=length(coordt)
if(spacetime_dyn) numtime=1
  numcoord=nrow(coords);
  dime<-numcoord*numtime
  xx=double(dime)

#########################################################
KK=1;sel=NULL;ssp=double(dime)

if(model%in% c("SkewGaussian","StudentT","SkewStudentT","TwoPieceTukeyh",
               "TwoPieceStudentT","TwoPieceGaussian"))
{
     
      simD=simu_approx(coords,coordt,method,corrmodel,param,M,L)
      if(!spacetime&&!bivariate) simDD <- c(simD)
      else simDD <- matrix(simD, nrow=numtime, ncol=numcoord,byrow=TRUE)
      param$nugget=0 #ojo
}

  while(KK<=npoi) {
  for(i in 1:k) {

    simd=simu_approx(coords,coordt,method,corrmodel,param,M,L)
   
    ####################################################################### 
    namesnuis<-NuisParam(model, bivariate,num_betas=num_betas,copula=NULL)
    nuisance<-param[namesnuis]
    if(i==1&&(model=="SkewGaussian"||model=="SkewGauss")&&bivariate) param["pcol"]=0
    ####################################
   
    sim<- RFfct1(numcoord,numtime,grid,spacetime,bivariate,
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
    if(model %in% c("Weibull","SkewGaussian","SkewGauss","Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonWeibull","PoissonZIP","Beta","Kumaraswamy","Kumaraswamy2",
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
 if(model %in% c("PoissonGamma"))   {
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
   #sel=rbind(sel,ssp<=c(exp(mm)*rowSums(sim3)^(1/param$shape)/((k-2)*gamma(1+1/param$shape))))
   if(sum(apply(sel,2,prod))==0) break  ## stopping rule
 }

 KK=KK+1
}

 ####### end for #########################
 ###############################################################################################
 #### simulation for discrete random field based on indipendent copies  of GRF ######
 ###############################################################################################
 if(model %in% c("Binomial","BinomialLogistic","Poisson","PoissonGamma","PoissonWeibull","PoissonZIP","BinomialNeg","BinomialNegZINB"))   {

   if(model %in% c("poisson","Poisson","PoissonGamma","PoissonWeibull"))   {sim=colSums(sel);byrow=TRUE}
    
    if(model %in% c("PoissonZIP"))   {
      a=simu_approx(coords,coordt,method,corrmodel,param,M,L)
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
  
      a=simu_approx(coords,coordt,method,corrmodel,param,M,L)
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
    #######################################
    if(bivariate)   numtime=1

    # Delete the global variables:
    # Return the objects list:
    GeoSim <- list(bivariate = bivariate,
    coordx = coords[,1],
    coordy = coords[,2],
    coordt = coordt,
    coordx_dyn =coordx_dyn,
    corrmodel = corrmodel,
    data = sim,
    distance = distance,
    grid = grid,
    model = model,
    method=method,
    n=n,
    numcoord = numcoord,
    numtime = numtime,
    param = param,
    radius = radius,
    randseed=.Random.seed,
    spacetime = spacetime,
    X=X)
#}
##############################################
    structure(c(GeoSim, call = call), class = c("GeoSim"))
}
