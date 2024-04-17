
####################################################
### File name: GeoVarestboostrap.r
####################################################

   
GeoVarestbootstrap=function(fit,K=100,sparse=FALSE,GPU=NULL,  local=c(1,1),optimizer="Nelder-Mead",
  lower=NULL, upper=NULL,method="cholesky", alpha=0.95, M=30,L=3000)
{

if(length(fit$coordt)==1) fit$coordt=NULL

if(is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")
if(is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")

if(!(method=="cholesky"||method=="Vecchia"||method=="TB"||method=="CE")) stop("The method of simulation is not correct")




if(is.numeric(alpha)) if(!(alpha<1&&alpha>0) ) stop(" alpha must be  numeric between 0 and 1")

model=fit$model

cat("Parametric bootstrap can be time consuming ...\n")

if(fit$missp)  ### misspecification
 {if(fit$model=="StudentT")     model="Gaussian_misp_StudentT"
  if(fit$model=="Poisson")      model="Gaussian_misp_Poisson"
  if(fit$model=="PoissonZIP")   model="Gaussian_misp_PoissonZIP"
  if(fit$model=="SkewStudentT") model="Gaussian_misp_SkewStudenT"
  if(fit$model=="Tukeygh")      model="Gaussian_misp_Tukeygh"
 }

dimat=fit$numtime*fit$numcoord;

tempX=fit$X
if(sum(fit$X[1:dimat]==1)==dimat&&!dim(fit$X)[2]>1) fit$X=NULL 


k=1;res=NULL



  coords=cbind(fit$coordx,fit$coordy)
  N=nrow(coords)
  pp=NULL
######## simulation ##########################################
if(is.null(fit$copula)){     ### non copula models
   cat("Performing",K,"simulations....\n")
    if(method=="cholesky")
    { 

    
      data_sim = GeoSim(coordx=coords,coordt=fit$coordt,coordx_dyn=fit$coordx_dyn, anisopars=fit$anisopars,
      corrmodel=fit$corrmodel,model=fit$model,param=append(fit$param,fit$fixed),
      GPU=GPU,  local=local,sparse=sparse,grid=fit$grid, X=fit$X,n=fit$n,method=method,
      distance=fit$distance,radius=fit$radius,nrep=K)}

   if(method=="TB"||method=="Vecchia"||method=="CE")   
     { data_sim = GeoSimapprox(coordx=coords,coordt=fit$coordt,coordx_dyn=fit$coordx_dyn, anisopars=fit$anisopars,
      corrmodel=fit$corrmodel,model=fit$model,param=append(fit$param,fit$fixed),
      GPU=GPU,  local=local,grid=fit$grid,X=fit$X,n=fit$n,method=method,
      M=M, L=L,distance=fit$distance,radius=fit$radius,nrep=K)}
}
else{    ### copula models
       cat("Performing simulation....\n")
        if(method=="cholesky")
     { data_sim = GeoSimCopula(coordx=coords,coordt=fit$coordt,coordx_dyn=fit$coordx_dyn, anisopars=fit$anisopars,
       corrmodel=fit$corrmodel,model=fit$model,copula=fit$copula,param=append(fit$param,fit$fixed),
       GPU=GPU,  local=local,sparse=sparse,grid=fit$grid,X=fit$X,n=fit$n,method=method,
       distance=fit$distance,radius=fit$radius,nrep=K)}
}
###############################################################
  cat("Performing estimation...\n")
  pb <- txtProgressBar(min = 0, max = K, style = 3)

while(k<=K){
    Sys.sleep(0.1)

res_est=GeoFit( data=data_sim$data[[k]], start=fit$param,fixed=fit$fixed,
   coordx=coords, coordt=fit$coordt, coordx_dyn=fit$coordx_dyn,
   copula=fit$copula,sensitivity=FALSE,anisopars=fit$anisopars,est.aniso=fit$est.aniso,
   lower=lower,upper=upper,memdist=TRUE,neighb=fit$neighb,
   corrmodel=fit$corrmodel, model=model, sparse=FALSE,n=fit$n,
   GPU=GPU,local=local,  maxdist=fit$maxdist, maxtime=fit$maxtime, optimizer=optimizer,
   grid=fit$grid, likelihood=fit$likelihood, type=fit$type,
   X=fit$X, distance=fit$distance, radius=fit$radius)


if(res_est$convergence=='Successful'&&res_est$logCompLik<1.0e8) 
 {
 res=rbind(res,unlist(res_est$param))

setTxtProgressBar(pb, k)
}   
 k=k+1            
close(pb)
}



numparam=length(fit$param)
invG=var(res); G=try(solve(invG),silent=TRUE);if(!is.matrix(G)) print("Bootstrap estimated Godambe matrix is singular")
stderr=sqrt(diag(invG))


if((fit$likelihood=="Marginal"&&(fit$type=="Independence"))||(fit$likelihood=="Marginal"&&(fit$type=="Pairwise"))||fit$likelihood=="Conditional"&&(fit$type=="Pairwise"))
{

H=fit$sensmat

penalty <- sum(diag(H%*%invG))

claic = -2*fit$logCompLik + 2*penalty
clbic = -2*fit$logCompLik + log(dimat)*penalty
fit$varimat=H%*%invG %*%H
}    
if( fit$likelihood=="Full"&&fit$type=="Standard"){
claic <- -2 * fit$logCompLik + 2*numparam
clbic <- -2 * fit$logCompLik + log(dimat)*2*numparam
}
fit$claic=claic
fit$clbic=clbic
fit$stderr=stderr
fit$varcov=invG
fit$estimates=res
fit$X=tempX


stderr=sqrt(diag(invG))
aa=qnorm(1-(1-alpha)/2)*stderr
pp=as.numeric(fit$param)
low=pp-aa; upp=pp+aa
fit$conf.int=rbind(low,upp)
fit$pvalues=2*pnorm(-abs(pp/stderr))
return(fit)
}
