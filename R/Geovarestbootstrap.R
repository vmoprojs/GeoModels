####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate.
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Institutions: 
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: Geo_varest_bootstrap.r
### Description:
### This file contains function to compute stderr estimation
### for composite likelihoodestimation
### and model selection critera given an object of class GeoFit
### Last change: 28/03/2020.
####################################################
   
GeoVarestbootstrap=function(fit,K=100,sparse=FALSE,GPU=NULL,  local=c(1,1),optimizer="Nelder-Mead",
  lower=NULL, upper=NULL,memdist=TRUE, seed=1)
{

k=1;res=NULL;#H=list();
if(length(fit$coordt)==1) fit$coordt=NULL
print("Parametric bootstrap can be time consuming ...")
if(is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")
model=fit$model
if(fit$model=="Gaussian_misp_StudentT") fit$model="StudentT"
if(fit$model=="Gaussian_misp_Poisson") fit$model="Poisson"
set.seed(round(seed))
while(k<=K){
data_sim = GeoSim(coordx=cbind(fit$coordx,fit$coordy),coordt=fit$coordt,
     coordx_dyn=fit$coordx_dyn,
     corrmodel=fit$corrmodel,model=fit$model,
	 param=as.list(c(fit$param,fit$fixed)),
	 GPU=GPU,  local=local,grid=fit$grid, X=fit$X,n=fit$n,method="cholesky",
	 distance=fit$distance,radius=fit$radius)

estimation=GeoFit( data=data_sim$data, start=as.list(fit$param),fixed=as.list(fit$fixed),
   coordx=cbind(fit$coordx,fit$coordy), coordt=fit$coordt, coordx_dyn=fit$coordx_dyn,
   lower=lower,upper=upper,memdist=memdist,
   corrmodel=fit$corrmodel, model=model, sparse=FALSE,n=fit$n,
   GPU=GPU,local=local,  maxdist=fit$maxdist, maxtime=fit$maxtime, optimizer=optimizer,
   grid=fit$grid, likelihood=fit$likelihood, type=fit$type,
   X=fit$X, distance=fit$distance, radius=fit$radius,
   vartype='SubSamp', weighted=FALSE, 
   taper=fit$taper, tapsep=fit$tapsep, method="cholesky",onlyvar=FALSE)

if(estimation$convergence=="Successful"){
  #print((estimation$param))
res=rbind(res,estimation$param)
print(k)
k=k+1

}}

dimat=fit$numtime*fit$numcoord; numparam=length(fit$param)
invG=var(res); G=try(solve(invG),silent=TRUE);if(!is.matrix(G)) print("Bootstrap estimated Godambe matrix is singular")
stderr=sqrt(diag(invG))

if((fit$likelihood=="Marginal"&&(fit$type=="Pairwise"))||fit$likelihood=="Conditional"&&(fit$type=="Pairwise"))
{

H=fit$sensmat

penalty <- sum(diag(H%*%invG))
claic <- -2 * fit$logCompLik + 2*sum(diag(penalty))
clbic <- -2 * fit$logCompLik + log(dimat)*sum(diag(penalty))
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
return(fit)

}
