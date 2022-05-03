
####################################################
### File name: GeoVarestboostrap.r
####################################################

   
GeoVarestbootstrap=function(fit,K=100,sparse=FALSE,GPU=NULL,  local=c(1,1),optimizer="Nelder-Mead",
  lower=NULL, upper=NULL,method="cholesky",memdist=TRUE, M=30,L=500,seed=1)
{


if(fit$coordt==0)  fit$coordt=NULL

if(is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")
if(is.null(fit$sensmat)) stop("Sensitivity matrix is missing: use sensitivity=TRUE in GeoFit")

if(!(method=="cholesky"||method=="Vecchia"||method=="TB")) stop("The method of simulation is not correct")

if(!is.numeric(seed)) stop(" seed must be numeric")
model=fit$model

print("Parametric bootstrap can be time consuming ...")

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
set.seed(seed)

  pb <- txtProgressBar(min = 0, max = K, style = 3)
  coords=cbind(fit$coordx,fit$coordy)
  N=nrow(coords)
  pp=NULL
while(k<=K){
Sys.sleep(0.1)
if(method=="cholesky") {

  if(is.null(fit$copula))
data_sim = GeoSim(coordx=coords,coordt=fit$coordt,
     coordx_dyn=fit$coordx_dyn, anisopars=fit$anisopars,
     corrmodel=fit$corrmodel,model=fit$model,
	 #param=as.list(c(fit$param,fit$fixed)),
      param=append(fit$param,fit$fixed),
	 GPU=GPU,  local=local,sparse=sparse,#grid=fit$grid, 
   X=fit$X,n=fit$n,method=method,
	 distance=fit$distance,radius=fit$radius)
else
data_sim = GeoSimCopula(coordx=coords,coordt=fit$coordt,
     coordx_dyn=fit$coordx_dyn, anisopars=fit$anisopars,
     corrmodel=fit$corrmodel,model=fit$model,
   copula=fit$copula,
      param=append(fit$param,fit$fixed),
   GPU=GPU,  local=local,sparse=sparse,#grid=fit$grid, 
   X=fit$X,n=fit$n,method=method,
   distance=fit$distance,radius=fit$radius)
}



#print(append(fit$param,fit$fixed))
if(method=="Vecchia"||method=="TB") {
                data_sim = GeoSimapprox(coordx=coords,coordt=fit$coordt, coordx_dyn=fit$coordx_dyn, corrmodel=fit$corrmodel,model=fit$model, 
                param=append(fit$param,fit$fixed),method=method,M=M,L=L,GPU=GPU,  local=local,#grid=fit$grid, 
                X=fit$X,n=fit$n,distance=fit$distance,radius=fit$radius)
            }


res_est=GeoFit( data=data_sim$data, start=fit$param,fixed=fit$fixed,#start=as.list(fit$param),fixed=as.list(fit$fixed),
   coordx=coords, coordt=fit$coordt, coordx_dyn=fit$coordx_dyn,
   copula=fit$copula,sensitivity=FALSE,anisopars=fit$anisopars,est.aniso=fit$est.aniso,
   lower=lower,upper=upper,memdist=memdist,neighb=fit$neighb,
   corrmodel=fit$corrmodel, model=model, sparse=FALSE,n=fit$n,
   GPU=GPU,local=local,  maxdist=fit$maxdist, maxtime=fit$maxtime, optimizer=optimizer,
   grid=fit$grid, likelihood=fit$likelihood, type=fit$type,
   X=fit$X, distance=fit$distance, radius=fit$radius)



if(res_est$convergence=='Successful'){
 
 res=rbind(res,unlist(res_est$param))
 k=k+1
setTxtProgressBar(pb, k)
}               
close(pb)
}



numparam=length(fit$param)
invG=var(res); G=try(solve(invG),silent=TRUE);if(!is.matrix(G)) print("Bootstrap estimated Godambe matrix is singular")
stderr=sqrt(diag(invG))


if((fit$likelihood=="Marginal"&&(fit$type=="Pairwise"))||fit$likelihood=="Conditional"&&(fit$type=="Pairwise"))
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
#set.seed(sample(1:10000,1))
return(fit)

}
