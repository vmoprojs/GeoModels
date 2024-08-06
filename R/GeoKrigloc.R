####################################################
### File name: GeoKrigloc.r
####################################################

GeoKrigloc= function(estobj=NULL,data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, corrmodel, distance="Eucl", grid=FALSE, loc,neighb=NULL,
              maxdist=NULL,maxtime=NULL, method="cholesky", model="Gaussian", n=1,nloc=NULL, mse=FALSE,  param, anisopars=NULL, 
              radius=6371, sparse=FALSE, time=NULL, type="Standard",
              type_mse=NULL, type_krig="Simple",weigthed=TRUE, which=1,copula=NULL, X=NULL,Xloc=NULL,Mloc=NULL,
              spobj=NULL,spdata=NULL,parallel=FALSE,ncores=NULL)


{
call=match.call()    

###############################
## checking if there is a  GeoFit object
if(!is.null(estobj)){
   if(!inherits(estobj,"GeoFit"))
               stop("need  a 'GeoFit' object as input\n")
   data=estobj$data

if(!estobj$grid){  #not regular grid 

 if(!estobj$bivariate){  if(is.null(estobj$coordx_dyn)) coordx=cbind(estobj$coordx,estobj$coordy)
                         else cord=estobj$coordx_dyn
                      } ## spatial (temporal) non regular case
 else  {    if(is.null(estobj$coordx_dyn))  { coordx=estobj$coordx[1:estobj$ns[1]]    # bivariate not dynamic    
                                              coordy=estobj$coordy[1:estobj$ns[2]] 
                                            }  
            else {coordx_dyn=estobj$coordx_dyn}                                      # bivariate  dynamic  
       }
 }
else  { coordx=estobj$coordx; 
        coordy=estobj$coordy
      }

   if(length(estobj$coordt)==1) coordt=NULL
   else coordt=estobj$coordt
   coordx_dyn=estobj$coordx_dyn
   corrmodel=estobj$corrmodel
   model=estobj$model
   distance=estobj$distance
   grid=estobj$grid
   n=estobj$n
   param=append(estobj$param,estobj$fixed)
   radius=estobj$radius
   copula=estobj$copula
   anisopars=estobj$anisopars
   if(ncol(estobj$X)==1) X=NULL
   else X=estobj$X
}
##################################
## X and more stuuffs..
M=NULL
spacetime=FALSE
bivariate=FALSE

########################################################
####### extracting sp objects if necessary 
########################################################
bivariate<-CheckBiv(CkCorrModel(corrmodel))
spacetime<-CheckST(CkCorrModel(corrmodel))
space=!spacetime&&!bivariate
if(!is.null(spobj)) {
   if(space||bivariate){
        a=sp2Geo(spobj,spdata); coordx=a$coords 
       if(!a$pj) {if(distance!="Chor") distance="Geod"}
    }
   if(spacetime){
        a=sp2Geo(spobj,spdata); coordx=a$coords ; coordt=a$coordt 
        if(!a$pj) {if(distance!="Chor") distance="Geod"}
     }
   if(!is.null(a$Y)&&!is.null(a$X)) {data=a$Y ; X=a$X }
}
########################################################



if(is.null(coordx_dyn)){
coords=coordx
if(!is.null(coordy)){
 if(!grid)  coords=cbind(coordx,coordy) 
 if(grid)   coords=as.matrix(expand.grid(coordx,coordy))
}
}
else{coordx=NULL;coordy=NULL;coords=NULL}


if(is.null(dim(loc)))  loc= matrix(loc,1,2)


Nloc=nrow(loc)
if(is.null(Nloc)) Nloc=1
Tloc=length(time)
if(bivariate)  Tloc=1

if(length(param$mean)>1) M=param$mean #### non constant mean


coremax=parallel::detectCores()
if(is.na(coremax)||coremax==1) parallel=FALSE


#####################################################################
if(space){
         ### computing spatial neighborhood
         neigh=GeoNeighborhood(data, coordx=coords,distance=distance,loc=loc,neighb=neighb,maxdist=maxdist,X=X,M=M,parallel=FALSE,ncores=ncores)
         res1=res2=NULL

################## not parallel version #######################################
if(!parallel) {
     res1=double(Nloc); if(mse) res2=double(Nloc)
         for(i in 1: Nloc)
          {
              #update mean
         if(!is.null(M)) param$mean=neigh$M[[i]]         
            pr=GeoKrig(estobj=NULL,loc=loc[i,], data=neigh$data[[i]],coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,n=n,
                X=neigh$X[[i]],Xloc= Xloc[i,],Mloc=Mloc[i], type_krig=type_krig,sparse=sparse,
                model=model, param=param,anisopars=anisopars, mse=mse,copula=copula)
                res1[i]=pr$pred
                if(mse) res2[i]=pr$mse
          }
}
##################### parallel version ####################################
if(parallel) {
      
        if(is.null(ncores)){ n.cores <- coremax - 1 }
        else
        {  if(!is.numeric(ncores)) stop("number of cores not valid\n")
           if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
           n.cores=ncores
        }
        future::plan(multisession, workers = n.cores)
        progressr::handlers(global = TRUE) 
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:Nloc)
        cat("Performing local kriging using ",n.cores," cores...\n")

        oopts=options(future.globals.maxSize = 8000 * 1024^2)
        on.exit(options(oopts))
 xx=foreach::foreach(i= 1:Nloc,.combine = rbind,
    .options.future = list(seed = TRUE,
    globals = structure(TRUE, add = c("param")))) %dofuture% 
        { 
            pb(sprintf("i=%g", i))
        if(!is.null(M)) param$mean=neigh$M[[i]]
            pr=GeoKrig(estobj=NULL,loc=loc[i,], data=neigh$data[[i]],coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,n=n,
                X=neigh$X[[i]],Xloc= Xloc[i,],Mloc=Mloc[i], type_krig=type_krig,sparse=sparse,
                model=model, param=param,anisopars=anisopars, mse=mse,copula=copula)
            pr$data=pr$coordx=pr$coordy=pr$coordt=NULL
            c(pr$pred,pr$mse)
        }
   
xx=matrix(xx,nrow=Nloc)
res1=as.numeric(xx[,1])
if(mse) res2=as.numeric(xx[,2])
rm(xx)
future::plan(sequential)
}
 #######################################################################
} # end space
######################################################################
if(spacetime)
{  
       ### computing spatio-temporal neighborhood
         neigh=GeoNeighborhood(data, coordx=coords,coordt=coordt,distance=distance,neighb=neighb,
                  loc=loc,time=time,maxdist=maxdist,maxtime=maxtime,X=X,M=M,parallel=FALSE,ncores=ncores)
         res1=res2=double(Nloc*Tloc)
         k=1
    if(!parallel)    {
        # pb <- txtProgressBar(min = 0, max = Nloc*Tloc, style = 3)
         for(i in 1: Nloc){
          for(j in 1: Tloc){
             if(!is.null(M)) param$mean=neigh$M[[k]]
            pr=GeoKrig(estobj=NULL, data=neigh$data[[k]],coordx=neigh$coordx[[k]],coordt=neigh$coordt[[k]],loc=loc[i,],time=time[j], #ok
               X=neigh$X[[k]],  Mloc=Mloc[i+(Nloc)*(j-1)], #ok
               Xloc= Xloc[i+(Nloc)*(j-1),], type_krig=type_krig,sparse=sparse,
             corrmodel=corrmodel,distance=distance, model=model, param=param,anisopars=anisopars, mse=mse,copula=copula,n=n)
            res1[k]=pr$pred
            if(mse) res2[k]=pr$mse
            k=k+1
         #    setTxtProgressBar(pb, k)
          #         close(pb)
          }}
      # print(str(res1))
     }
   if(parallel) {
       if(is.null(ncores)){ n.cores <- coremax - 1 }
        else
        {  if(!is.numeric(ncores)) stop("number of cores not valid\n")
           if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
           n.cores=ncores
        }
     doFuture::registerDoFuture()
        future::plan(multisession, workers = n.cores)
 ############################# 
        # Define your code
        compute_krig <- function(k, neigh, loc, time, Mloc, Xloc, type_krig, corrmodel, distance, model, param, anisopars, mse, copula, n,sparse) {
          
          if (!is.null(M)) param$mean <- neigh$M[[k]]
          pr <- GeoKrig(estobj = NULL, data = neigh$data[[k]], coordx = neigh$coordx[[k]], coordt = neigh$coordt[[k]], loc = loc, time = time,
                        X = neigh$X[[k]], Mloc = Mloc, Xloc = Xloc, type_krig = type_krig,sparse=sparse,
                        corrmodel = corrmodel, distance = distance, model = model, param = param, anisopars = anisopars, mse = mse, copula = copula, n = n)
          # cat("k ",k,"\n")
          # print(pr$pred)
          return(list(pred = pr$pred, mse = ifelse(mse, pr$mse, NA)))
        }
        
        # Flatten the loops
        loop_indices <- expand.grid(i = 1:Nloc, j = 1:Tloc)
        loop_indices <- loop_indices[order(loop_indices$i),]
        rownames(loop_indices) <- NULL
        
        # Parallelize the computation using %dopar%
        progressr::handlers(global = TRUE) 
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:(Nloc * Tloc))
        cat("Performing local kriging using ",n.cores," cores...\n")

            oopts=options(future.globals.maxSize = 8000 * 1024^2)
        on.exit(options(oopts))
        
        output <- foreach::foreach(k = 1:(Nloc * Tloc), .combine = rbind,.options.future = list(seed = TRUE)) %dofuture% {
          pb(sprintf("k=%g", k))
          idx <- loop_indices[k, ]
          compute_krig(k, neigh, loc[idx$i, ], time[idx$j], Mloc[idx$i + (Nloc) * (idx$j - 1)], Xloc[idx$i + (Nloc) * (idx$j - 1), ], 
                       type_krig, corrmodel, distance, model, param, anisopars, mse, copula, n,sparse)
        }
        
        # Extract the predictions and mse values
        res1 <- as.numeric(unlist(output[, "pred"]))
        res2 <- as.numeric(unlist(output[, "mse"]))
  #################################     

        future::plan(sequential)     

   }    
} #### end spacetime


if(bivariate)
{ 
neigh=GeoNeighborhood(data, coordx=coords,distance=distance,loc=loc,maxdist=maxdist,neighb=neighb,bivariate=TRUE,X=X,parallel=FALSE,ncores=ncores)
        res1=res2=NULL
############# not parallel################################        
                 #  pb <- txtProgressBar(min = 0, max = Nloc, style = 3)
if(!parallel){           
         for(i in 1: Nloc)
          {
            pr=GeoKrig(loc=loc[i,],coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,n=n,
                X=neigh$X,,Xloc= Xloc[i,],which=which,type_krig=type_krig,sparse=sparse,
                model=model, param=param,anisopars=anisopars, mse=mse, data=neigh$data[[i]],copula=copula)
                res1=c(res1,pr$pred)
               if(mse) res2=c(res2,pr$mse)
               # setTxtProgressBar(pb, i)
               # close(pb)  
          }
}
############## parallel ###############################
if(parallel) {
        if(is.null(ncores)){ n.cores <- coremax - 1 }
    else
        {  if(!is.numeric(ncores)) stop("number of cores not valid\n")
           if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
           n.cores=ncores
        }
        future::plan(multisession, workers = n.cores)
        progressr::handlers(global = TRUE) 
        progressr::handlers("txtprogressbar")
        pb <- progressr::progressor(along = 1:Nloc)
        cat("Performing local kriging using ",n.cores," cores...\n")

        oopts=options(future.globals.maxSize = 8000 * 1024^2)
        on.exit(options(oopts))

 xx=foreach::foreach(i= 1:Nloc,.combine = rbind,
    .options.future = list(seed = TRUE,
    globals = structure(TRUE, add = c("param")))) %dofuture% {
          pb(sprintf("i=%g", i))
      pr=GeoKrig(loc=loc[i,],coordx=neigh$coordx[[i]],corrmodel=corrmodel,distance=distance,n=n,
         X=neigh$X,,Xloc= Xloc[i,],which=which,type_krig=type_krig,sparse=sparse,
         model=model, param=param,anisopars=anisopars, mse=mse, data=neigh$data[[i]],copula=copula)
    
             pr$data=pr$coordx=pr$coordy=pr$coordt=NULL
            c(pr$pred,pr$mse)
        }
res1=as.numeric(xx[,1])
if(mse) res2=as.numeric(xx[,2])
rm(xx)
future::plan(sequential)

}
}###### end bivariate #######################################
varpred=NULL
  if(spacetime||bivariate) {
            pred=matrix(t(res1),nrow=Tloc,ncol=Nloc);
            if(mse) varpred=matrix(c(res2),nrow=Tloc,ncol=Nloc);
    } 
  else{pred=c(res1)
       if(mse)varpred=c(res2)}
        
if(Tloc==1)  {c(pred);c(varpred)}
    # Return the objects list:
   GeoKrigloc = list(   
                    bivariate=bivariate,
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
                   n=n,
                   param = param,
                   pred=pred,
                   radius=radius,
                   spacetime = spacetime,
                   time=time,
              #     type=type,
                  type_krig=type_krig,
                   mse=varpred)
              #     mse2=varpred2)
    structure(c(GeoKrigloc, call = call), class = c("GeoKrigloc"))

}