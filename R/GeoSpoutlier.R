####################################################
### File name: GeoSpoutlier.r
####################################################

GeoSpoutlier <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, 
                         distance="Eucl", grid=FALSE,  neighb=10,alpha=0.001,
                         method="Z-Median",radius=6371, bivariate=FALSE,X=NULL)

{
    
    if(abs(alpha)>=1) stop("alpha must be between 0 and 1")
    if(!is.logical(bivariate))   stop("bivariate must be logical")
    if(is.null(neighb))     stop("maxdist (maxtime) or neighb must  be specified")
    if(method!="Z-Median") stop("method is not valid")
#########################
spacetime=FALSE
dyn=FALSE
  if(!is.null(coordx_dyn))  dyn=TRUE  
## handling spatial coordinates
    if(is.null(coordy)) coords=as.matrix(coordx)
    else{
    if(grid) coords=as.matrix(expand.grid(coordx,coordy))
    else     coords=cbind(coordx,coordy)  
    }

#####################################
if(is.numeric(coordt)) if(length(coordt)>1) {spacetime=TRUE}
space=!spacetime&&!bivariate 
##################################################################


   if(grid)     {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2] }

    if(is.null(coordx_dyn))
    {
      if(!is.null(coordy)){coordy <- coordx[,2]
                          coordx <- coordx[,1]
                          coords=cbind(coordx,coordy)}
      else 
      {
        if(!bivariate) coords=coordx
        if(bivariate)  coords=rbind(coordx,coordx)} 
      data=c(t(data))   
      numcoord=nrow(coords)               
      ns<-rep(numcoord,length(coordt))
      if(bivariate) ns=ns/2
    }
    else
    {
       env <- new.env()
       coords=do.call(rbind,args=c(coordx_dyn),envir = env) 
       data=unlist(data)
       ns=lengths(coordx_dyn)/2 
    }

  if(distance=="Geod"||distance=="Chor")
{
   coords_p=coords; 
   prj=mapproj::mapproject(coords[,1], coords[,2], projection="sinusoidal") 
   coords=radius*cbind(prj$x,prj$y)
}

#####################################################################
###  space
if(!bivariate&&!spacetime)
 {
       N=nrow(coords)
  ####     
  if(method=="Z-Median"){
       aa=RANN::nn2(coords,k = neighb)
       a=list()
       for(i in 1:N) a[[i]]=data[aa$nn.idx[i,]][-1]
       res=matrix(unlist(a),ncol=neighb-1,byrow=T)
       g=apply(res,1,median)
       h=data-g
       hh=(h-mean(h))/sqrt(var(h))
       sel=(hh>qnorm(1-alpha/2))|(hh<qnorm(alpha/2))
 
    if(distance=="Eucl"){
       if(!is.null(X)) a=cbind(coords,data,X)[sel,]
       if( is.null(X)) a=cbind(coords,data)[sel,]
                        }
    if(distance=="Geod"||distance=="Chor"){
       if(!is.null(X)) a=cbind(coords_p,data,X)[sel,]
       if( is.null(X)) a=cbind(coords_p,data)[sel,]
                       }

     }
   }    
#####################################################################
###  spacetime
if(spacetime)
 {     
 }
#####################################################################
## bivariate
if(bivariate)
 {
 
 }

return(a) 
}
