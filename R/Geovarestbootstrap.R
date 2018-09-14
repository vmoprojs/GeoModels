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
### Last change: 28/03/2018.
####################################################
   
GeoVarestbootstrap=function(fit,K=100,sparse=FALSE,GPU=NULL,  local=c(1,1))
{

k=1;res=NULL;#H=list();
if(fit$coordt==0) fit$coordt=NULL
print("Parametric bootstrap can be time consuming ...")
while(k<=K){
data_sim = GeoSim(coordx=cbind(fit$coordx,fit$coordy),coordt=fit$coordt,
     coordx_dyn=fit$coordx_dyn,
     corrmodel=fit$corrmodel,model=fit$model,
	 param=as.list(c(fit$param,fit$fixed)),
	 GPU=GPU,  local=local,grid=fit$grid, X=fit$X,n=fit$n,method="cholesky",
	 distance=fit$distance,radius=fit$radius)

estimation=GeoFit( data=data_sim$data, start=as.list(fit$param),fixed=as.list(fit$fixed),
   coordx=cbind(fit$coordx,fit$coordy), coordt=fit$coordt, coordx_dyn=fit$coordx_dyn,
   corrmodel=fit$corrmodel, model=fit$model, sparse=FALSE,n=fit$n,
   GPU=GPU,local=local,  maxdist=fit$maxdist, maxtime=fit$maxtime, 
   grid=fit$grid, likelihood=fit$likelihood, type=fit$type,
   X=fit$X, distance=fit$distance, radius=fit$radius,
   varest=FALSE, 
   vartype='SubSamp', weighted=FALSE, winconst=NULL,
   taper=fit$taper, tapsep=fit$tapsep, winstp=NULL,winconst_t=NULL, winstp_t=NULL,method="cholesky",onlyvar=FALSE)

if(estimation$convergence=="Successful"){
res=rbind(res,estimation$param)
k=k+1
}}

dimat=fit$numtime*fit$numcoord; numparam=length(fit$param)
invG=var(res); G=solve(invG)
stderr=sqrt(diag(invG))

if((fit$likelihood=="Marginal"&&(fit$type=="Pairwise"))||fit$likelihood=="Conditional"&&(fit$type=="Pairwise"))
{
#HH=
#penalty <- sum(diag(HH%*%invG))/K
penalty <- sum(diag(fit$sensmat%*%invG))
claic <- -2 * fit$logCompLik + 2*sum(diag(penalty))
clbic <- -2 * fit$logCompLik + log(dimat)*sum(diag(penalty))
}    
if( fit$likelihood=="Full"&&fit$type=="Standard"){
claic <- -2 * fit$logCompLik + 2*numparam
clbic <- -2 * fit$logCompLik + log(dimat)*2*numparam
}


return(list(claic = claic,clbic =clbic,stderr=stderr,varcov=invG))

}
