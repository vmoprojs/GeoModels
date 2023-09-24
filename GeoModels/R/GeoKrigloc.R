####################################################
### File name: GeoKrigloc.r
####################################################

GeoKrigloc= function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc,neighb=NULL,
              maxdist=NULL,maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE,  param, anisopars=NULL, 
              radius=6371, sparse=FALSE, time=NULL, type="Standard",
              type_mse=NULL, type_krig="Simple",weigthed=TRUE, which=1,copula=NULL, X=NULL,Xloc=NULL,Mloc=NULL)


{

## X and more stuuffs..
M=NULL
spacetime=FALSE
bivariate=FALSE
if(!is.null(coordt)) spacetime=TRUE
if(!is.null(dim(data))) if(nrow(data)==2&&is.null(coordt)) bivariate=TRUE


space=!spacetime&&!bivariate

if(is.null(coordx_dyn)){
coords=coordx
if(!is.null(coordy)){
 if(!grid)  coords=cbind(coordx,coordy) 
 if(grid)   coords=as.matrix(expand.grid(coordx,coordy))
}
}
else{coordx=NULL;coordy=NULL;coords=NULL}

Nloc=nrow(loc)
Tloc=length(time)
if(bivariate)  Tloc=1

if(length(param$mean)>1) M=param$mean #### non constant mean



#####################################################################
if(space){
         ### computing spatial neighborhood
         neigh=GeoNeighborhood(data, coordx=coords,distance=distance,loc=loc,neighb=neighb,maxdist=maxdist,X=X,M=M)
         res1=res2=NULL
         #  pb <- txtProgressBar(min = 0, max = Nloc, style = 3)
         for(i in 1: Nloc)
          {
              #update mean
         if(!is.null(M)) param$mean=neigh$M[[i]]
        #    Sys.sleep(0.1)
            
            pr=GeoKrig(loc=loc[i,], data=neigh$data[[i]],coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,n=n,
                X=neigh$X[[i]],Xloc= Xloc[i,],Mloc=Mloc[i],
                model=model, param=param,anisopars=anisopars, mse=mse,copula=copula)
                res1=c(res1,pr$pred)
                if(mse) res2=c(res2,pr$mse)

         #   setTxtProgressBar(pb, i)
         #   close(pb)
    
          }
}
######################################################################
if(spacetime)
{  

       ### computing spatio-temporal neighborhood
         neigh=GeoNeighborhood(data, coordx=coords,coordt=coordt,distance=distance,neighb=neighb,
                  loc=loc,time=time,maxdist=maxdist,maxtime=maxtime,X=X,M=M)
         res1=res2=NULL
         k=1
        
        # pb <- txtProgressBar(min = 0, max = Nloc*Tloc, style = 3)
         for(i in 1: Nloc){
          for(j in 1: Tloc){
             if(!is.null(M)) param$mean=neigh$M[[k]]
            pr=GeoKrig(data=neigh$data[[k]],coordx=neigh$coordx[[k]],coordt=neigh$coordt[[k]],loc=loc[i,],time=time[j], #ok
               X=neigh$X[[k]],  Mloc=Mloc[i+(Nloc)*(j-1)], #ok
               Xloc= Xloc[i+(Nloc)*(j-1),],
             corrmodel=corrmodel,distance=distance, model=model, param=param,anisopars=anisopars, mse=mse,copula=copula,n=n)
            res1=c(res1,pr$pred)
            if(mse) res2=c(res2,pr$mse)
            k=k+1
         #    setTxtProgressBar(pb, k)
          #         close(pb)
          }}
}
if(bivariate)
{ 
neigh=GeoNeighborhood(data, coordx=coords,distance=distance,loc=loc,maxdist=maxdist,neighb=neighb,bivariate=TRUE,X=X)
        res1=res2=NULL
                 #  pb <- txtProgressBar(min = 0, max = Nloc, style = 3)
         for(i in 1: Nloc)
          {
           
            pr=GeoKrig(loc=matrix(loc[i,],ncol=2),coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,n=n,
                X=neigh$X,,Xloc= Xloc[i,],which=which,
                model=model, param=param,anisopars=anisopars, mse=mse, data=neigh$data[[i]],copula=copula)
                res1=c(res1,pr$pred)
               if(mse) res2=c(res2,pr$mse)
               # setTxtProgressBar(pb, i)
               # close(pb)
              
          }
}
varpred=NULL
  if(spacetime||bivariate) {
            pred=matrix(t(res1),nrow=Tloc,ncol=Nloc);
            if(mse) varpred=matrix(c(res2),nrow=Tloc,ncol=Nloc);
    } 
  else{pred=c(res1)
       if(mse)varpred=c(res2)}
              

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