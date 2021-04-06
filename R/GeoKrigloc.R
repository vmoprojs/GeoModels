####################################################
### File name: GeoKrigloc.r
####################################################

GeoKrigloc= function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc,neighb=NULL,
              maxdist=NULL,maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE,  param, radius=6371, sparse=FALSE, 
               time=NULL, type="Standard",type_mse=NULL, type_krig="Simple",weigthed=TRUE, which=1,copula=NULL, X=NULL,Xloc=NULL)


{

## X and more stuuffs..
spacetime=FALSE
bivariate=FALSE
if(!is.null(coordt)) spacetime=TRUE
if(!is.null(dim(data))) if(nrow(data)==2&&is.null(coordt)) bivariate=TRUE


space=!spacetime&&!bivariate

coords=coordx
if(!is.null(coordy)){
 if(!grid)  coords=cbind(coordx,coordy) 
 if(grid)   coords=as.matrix(expand.grid(coordx,coordy))
}

Nloc=nrow(loc)
Tloc=length(time)
if(bivariate)  Tloc=1


#####################################################################
if(space){
         ### computing spatial neighborhood
         neigh=GeoNeighborhood(data, coordx=coords,distance=distance,loc=loc,neighb=neighb,maxdist=maxdist,X=X)
         res1=res2=NULL
         for(i in 1: Nloc)
          {
            pr=GeoKrig(loc=loc[i,],coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,
                X=neigh$X[[i]],Xloc= Xloc[i,],
                model=model, param=param,mse=mse, data=neigh$data[[i]],copula=copula)
                res1=c(res1,pr$pred)
                res2=c(res2,pr$mse)
          }
}
######################################################################
if(spacetime)
{  
       ### computing spatio-temporal neighborhood
         neigh=GeoNeighborhood(data, coordx=coords,coordt=coordt,distance=distance,neighb=neighb,
                  loc=loc,time=time,maxdist=maxdist,maxtime=maxtime,X=X)
         res1=res2=NULL
         k=1
         for(i in 1: Nloc){
          for(j in 1: Tloc){
            pr=GeoKrig(loc=loc[i,],time=time[j],coordx=neigh$coordx[[i]],coordt=neigh$coordt[[j]],
               X=neigh$X[[i]],Xloc= Xloc[i+(Nloc)*(j-1),],
             corrmodel=corrmodel,distance=distance, model=model, param=param,mse=mse, data=neigh$data[[k]],copula=copula)
            res1=c(res1,pr$pred)
            res2=c(res2,pr$mse)
            k=k+1
          }}
}
if(bivariate)
{ 
neigh=GeoNeighborhood(data, coordx=coords,distance=distance,loc=loc,maxdist=maxdist,neighb=neighb,bivariate=TRUE,X=X)
        res1=res2=NULL
         for(i in 1: Nloc)
          {
            pr=GeoKrig(loc=matrix(loc[i,],ncol=2),coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,
                X=neigh$X,,Xloc= Xloc[i,],which=which,
                model=model, param=param,mse=mse, data=neigh$data[[i]],copula=copula)
                res1=c(res1,pr$pred)
                res2=c(res2,pr$mse)
          }
}
varpred=NULL
  if(spacetime||bivariate) {
            pred=matrix(t(res1),nrow=Tloc,ncol=Nloc);
            varpred=matrix(c(res2),nrow=Tloc,ncol=Nloc);
    } 
  else{pred=c(res1);varpred=c(res2)}
              

if(Tloc==1)  {c(pred);c(varpred)}
    # Return the objects list:
    Kg = list(  #  bivariate=bivariate,
                    coordx = coordx,
                    coordy = coordy,
                    coordt = coordt,
                  # coordx_dyn=covmatrix$coordx_dyn,
                    corrmodel = corrmodel,
                    data=data,
                    distance = distance,
                    grid=grid,
                    loc=loc,
                    copula=copula,
              #     ns=pr$ns,
                   numcoord = nrow(coords),
                   numloc= Nloc,
                   numtime = length(coordt),
                   numt = Tloc,
                   maxdist=maxdist,
                   maxtime=maxtime,
                   model=model,
                   param = param,
                   pred=pred,
                   radius=radius,
                   spacetime = spacetime,
                   time=time,
              #     type=type,
              #     type_krig=type_krig,
                   mse=varpred)
              #     mse2=varpred2)
    structure(c(Kg, call = call), class = c("Kg"))
return(Kg)
}