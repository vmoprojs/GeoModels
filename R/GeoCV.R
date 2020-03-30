GeoCV=function(fit, K=100, n.fold=0.05, sparse=FALSE, which=1,seed=1)
{

if(n.fold>0.99||n.fold<0.01) stop("n.fold must be beween 0.01 and 0.99")
print("Cross-validation  kriking can be time consuming ...")
if(ncol(fit$X)==1) {X=Xloc=NULL}
mae=rmse=NULL
space_dyn=FALSE

if(is.list(fit$data)) space_dyn=TRUE
spacetime<-CheckST(CkCorrModel(fit$corrmodel))
bivariate<-CheckBiv(CkCorrModel(fit$corrmodel))
K=round(K)
if(K<2) stop("K must be grater or equal  to 2")
if(K>10000) stop("K is  too large")
if(bivariate)
   {if(!(which==1||which==2))
          stop("which must be 1 or 2")}

i=1
print(paste("Starting iteration from 1 to",K," ..."))
############################################################
########### spatial case ###################################
############################################################
if(!spacetime&&!bivariate)
{
N=length(fit$data)
coords=cbind(fit$coordx,fit$coordy)
data=fit$data

set.seed(round(seed))
while(i<=K){
sel_data = sample(1:N,round(N*(1-n.fold)))  
# data and coord used for prediction
datanew   = data[sel_data]
coordsnew = coords[sel_data,]
if(!is.null(X)) {
                X=fit$X[sel_data,]
                Xloc=fit$X[-sel_data,]
                }
# data and coord to predict
data_to_pred  = data[-sel_data]
loc_to_pred   = coords[-sel_data,]
pr=GeoKrig(data=datanew, coordx=coordsnew,  
	       corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=loc_to_pred, #ok
	          model=fit$model, n=fit$n, #ok
           param=as.list(c(fit$param,fit$fixed)), 
           radius=fit$radius, sparse=sparse, X=X,Xloc=Xloc) #ok
err=data_to_pred-pr$pred
N2=length(err)
rmse=c(rmse,sqrt(sum(err^2)/N2))
mae= c(mae,      sum(abs(err))/N2)
print(i)
i=i+1
}}
############################################################
########### spatial bivariate case ###################################
############################################################
if(bivariate)
{
ns=fit$ns
if(space_dyn) {data1=fit$data[[1]]
               data2=fit$data[[2]]
               coords1=fit$coordx_dyn[[1]]
               coords2=fit$coordx_dyn[[2]]
               }
if(!space_dyn) {data1=fit$data[1,];data2=fit$data[2,]
                coords=cbind(fit$coordx,fit$coordy)
                coords1=coords; coords2=coords
               }
if(!is.null(X)) {
                if(!is.list(fit$X)){X1=fit$X[1:ns[1],];X2=fit$X[(ns[1]+1):(ns[1]+ns[2]),]}
                if( is.list(fit$X)){X1=fit$X[[1]];X2=fit$X[[2]]}
                }
set.seed(round(seed))
while(i<=K){
#######################################	
if(which==1) {
	          sel_data = sample(1:ns[1],round(ns[1]*(1-n.fold))) 
              cc1= coords1[sel_data,]
              loc_to_pred   = coords1[-sel_data,]
              d1   = data1[sel_data]
              data_to_pred  = data1[-sel_data]
              if(!is.null(X)) { X=rbind(X1[sel_data,],X2); Xloc=rbind(X1[-sel_data,],X2)}
datanew=list();datanew[[1]]=d1;datanew[[2]]=data2;
coordsnew=list();coordsnew[[1]]=cc1;coordsnew[[2]]=coords2;
          }
if(which==2) {
	          sel_data = sample(1:ns[2],round(ns[2]*(1-n.fold))) 
              cc2= coords2[sel_data,]
              loc_to_pred   = coords2[-sel_data,]
              d2   = data2[sel_data]
              data_to_pred  = data2[-sel_data]
              if(!is.null(X)) {X=rbind(X1,X2[sel_data,]); Xloc=rbind(X1,X2[-sel_data,])}
datanew=list();datanew[[1]]=data1;datanew[[2]]=d2;
coordsnew=list();coordsnew[[1]]=coords1;coordsnew[[2]]=cc2;
            }
#####################################
pr=GeoKrig(data=datanew, coordx=NULL,   coordt=NULL, coordx_dyn=coordsnew,  #ok
	       corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=loc_to_pred, #ok
	          model=fit$model, n=fit$n, #ok
           param=as.list(c(fit$param,fit$fixed)), 
           radius=fit$radius, sparse=sparse,   time=NULL, 
             which=which, X=X,Xloc=Xloc) #ok    

err=data_to_pred-pr$pred  
   
N2=length(err)
rmse=c(rmse,sqrt(sum(err^2)/N2))
mae= c(mae,      sum(abs(err))/N2)
print(i)
i=i+1
}                
} ## end bivariate
############################################################
########### spatio temporal case ###########################
############################################################

if(spacetime)
{
coords=cbind(fit$coordx,fit$coordy)
ns=fit$ns
coordt=fit$coordt
T=length(coordt)
NT=sum(ns)
if(is.null(X)) rep(1,NT)

if(!space_dyn){
data=coordx_dyn=data_tot=list()
for(k in 1:T) {
              data[[k]]=fit$data[k,]
              coordx_dyn[[k]]=coords
              data_tot[[k]]=cbind(coords,rep(k,ns[k]),fit$data[k,])
              }
}
env <- new.env()
data_tot=do.call(rbind,args=c( data_tot),envir = env) 
if(is.list(X))  X=do.call(rbind,args=c(X),envir = env)

data_tot=cbind(data_tot,X)
set.seed(round(seed))
while(i<=K){
####################################### 
sel_data = sample(1:NT,round(NT*(1-n.fold))) 
data_sel=data_tot[sel_data,]
data_to_pred=data_tot[-sel_data,]

data_sel_ord=data_sel[order(data_sel[,3]),]
data_to_pred_ord=data_to_pred[order(data_to_pred[,3]),]

#pr=GeoKrig(data=datanew,    coordt=timenew, coordx_dyn=coordx_dyn,  #ok
 #        corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=loc_to_pred, #ok
  #          model=fit$model, n=fit$n, #ok
   #        param=as.list(c(fit$param,fit$fixed)), 
    #       radius=fit$radius, sparse=sparse, time=time_to_pred, X=X,Xloc=Xloc) #ok

err=c(data_to_pred)-c(pr$pred)  
   
N2=length(err)
rmse=c(rmse,sqrt(sum(err^2)/N2))
mae= c(mae,      sum(abs(err))/N2)
print(i)
i=i+1
}




} ## end spacetime

return(list(rmse=rmse,mae=mae))
}
