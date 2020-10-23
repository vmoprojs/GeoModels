GeoNeighborhood = function(data=NULL, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,bivariate=FALSE, 
                  distance="Eucl", grid=FALSE, loc, neighb=NULL,maxdist=NULL,maxtime=NULL, radius=6371, time=NULL, X=NULL)
{
  XX=NULL
  numtime=1
  sel_ss=1
  sel_tt=1
  if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
  if(!is.matrix(loc))   stop("loc parameter must be a matrix")
  if(!is.logical(bivariate))   stop("bivariate must be logical")
    #if(!(ncol(loc)==2))   stop("loc parameter must be a matrix  N X 2")
  if(is.null(neighb)&&is.null(maxdist))     stop("maxdist (maxtime) or neighb must  be specified")
  if(is.numeric(maxdist)&&is.numeric(neighb))   stop("maxdist (maxtime) or neighb must  be specified")
 
   spacetime=FALSE
  corrmodel="Exponential"
  
  if(!is.null(coordt)) {spacetime=TRUE;corrmodel="Exp_Exp"}
  if(spacetime) if(!is.vector(time))  stop("time parameter is missing")
 # checkinput = CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting",
 #                            NULL, grid, 'Marginal', maxdist, maxtime, 'Gaussian', 1,
 #                             'Nelder-Mead', NULL, radius, NULL, NULL, NULL, 
 #                         'Pairwise', FALSE, 'SubSamp', FALSE, X)
 # if(!is.null(checkinput$error)) stop(checkinput$error)
  dyn=FALSE
  if(!is.null(coordx_dyn))  dyn=TRUE  
## handling spatial coordinates
    if(is.null(coordy)) coords=as.matrix(coordx)
    else{
    if(grid) coords=as.matrix(expand.grid(coordx,coordy))
    else     coords=cbind(coordx,coordy)  
    }

Nloc=nrow(loc) # number of location sites
#####################################
sel_tt=NULL
colnames(loc)=NULL;colnames(coords)=NULL;
space=!spacetime&&!bivariate 
##################################################################

  if(distance=="Geod"||distance=="Chor")
{
   coords_p=coords; loc_p=loc
   prj=mapproj::mapproject(coords_p[,1], coords_p[,2], projection="sinusoidal") 
   coords_p=radius*cbind(prj$x,prj$y)
   prjloc=mapproj::mapproject(loc_p[,1], loc_p[,2], projection="sinusoidal")
   loc_p=radius*cbind(prjloc$x,prjloc$y)
}
##################################################################################
#####################  # case fixed radius #######################################
##################################################################################
if(is.null(neighb))             
{
  if(distance=="Geod"||distance=="Chor")
        out<- LatticeKrig::LKDist(coords_p,loc_p,delta=maxdist,distance.type="Euclidean")
  if(distance=="Eucl")
        out<- LatticeKrig::LKDist(coords,loc,delta=maxdist,distance.type="Euclidean")
##################################################################################
if(space){
    sel_ss=data_sel=numpoints=XX=list()
    ## checkinf if  there is at leas one empty neigh...
    for(i in 1:Nloc)
    {
     sel=out$ind[,1][out$ind[,2]==i]
     sel_ss[[i]]=coords[sel,]
     numpoints[[i]]=nrow(sel_ss[[i]])
     if(!is.null(data)) data_sel[[i]]=data[sel]
     if(!is.null(X)) XX[[i]]=X[sel,]
     }
 #if(dim(as.matrix(sel_ss,nrow=dimat2))[1]==0) stop("spatial distance for local kringing is too small")
}
#####################################################################################
if(spacetime)
{
   ##  it works only for fixed spatial locations
   Tloc=length(time); TT=length(coordt)
   coordt1=cbind(coordt,rep(0,TT))
   time1=cbind(time,rep(0,Tloc))
   out_t<- LatticeKrig::LKDist(coordt1,time1,delta=maxtime,distance.type="Euclidean") # temporal distance
   sel_ss=numpoints=data_sel=sel_tt=list()
   k=1
   for(i in 1:Nloc){
      sel_s=out$ind[,1][out$ind[,2]==i]
      sel_ss[[i]]=matrix(coords[sel_s,],ncol=2)
      for(j in 1:Tloc){
       sel_t=out_t$ind[,1][out_t$ind[,2]==j]
       sel_tt[[j]]=coordt[sel_t]
      if(!is.null(data)) data_sel[[k]]=data[sel_t,sel_s]
       k=k+1
    }
    if(!is.null(X)) {
                      sss=NULL;
                      for (l in 1:Tloc) sss=c(sss,sel_s+(l-1)*Nloc)
                      XX[[i]]=X[sss,]
                      }
  }
}
#####################################################################################
if(bivariate)
  {
   Nloc=nrow(loc)
   if(dyn) coords=rbind(coords,coords)
   sel_ss=numpoints=data_sel=sel_tt=XX=list()
  for(i in 1:Nloc){
      sel=out$ind[,1][out$ind[,2]==i]
      sel_ss[[i]]=coords[sel,]
      numpoints[[i]]=ncol(sel_ss[[i]])
      if(!is.null(data))data_sel[[i]]=matrix(data[,sel],nrow=2)
      if(!is.null(X))   XX[[i]]=X[c(sel,2*sel),]
                }
   }
}     
##################################################################################
#####################  # case for Neighborhood of order max.points ###############
##################################################################################
if(!is.null(neighb))             
{    
   ## computing neigh indexes
 if(distance=="Geod"||distance=="Chor")  out<- RANN::nn2(coords_p,loc_p, k=neighb)
 if(distance=="Eucl")                    out<- RANN::nn2(coords,loc,     k=neighb)
#################################
 if(space){
    sel_ss=data_sel=numpoints=XX=list()
    for(i in 1:Nloc)
    {
     ss=out$nn.idx[i,];
     sel_ss[[i]]=coords[ss,]
     numpoints[[i]]=nrow(sel_ss[[i]])
     if(!is.null(data)) data_sel[[i]]=data[ss]
     if(!is.null(X)) XX[[i]]=X[ss,]
    }
 }
 #####################################################################################
if(spacetime)
{
  Tloc=length(time); TT=length(coordt)
  coordt1=cbind(coordt,rep(0,TT))
  time1=cbind(time,rep(0,Tloc))
  out_t<- LatticeKrig::LKDist(coordt1,time1,delta=maxtime,distance.type="Euclidean") # temporal distance
  sel_ss=numpoints=data_sel=sel_tt=XX=list()
  k=1
   for(i in 1:Nloc){
      sel_s=out$nn.idx[i,];
      sel_ss[[i]]=matrix(coords[sel_s,],ncol=2)
      for(j in 1:Tloc){
       sel_t=out_t$ind[,1][out_t$ind[,2]==j]
       sel_tt[[j]]=coordt[sel_t]
      if(!is.null(data)) data_sel[[k]]=data[sel_t,sel_s]
       k=k+1
    }
     if(!is.null(X)) {
                      sss=NULL;
                      for (l in 1:Tloc) sss=c(sss,sel_s+(l-1)*Nloc)
                      XX[[i]]=X[sss,]
                      }
  }
}

#####################################################################################
if(bivariate)
  {
   Nloc=nrow(loc)
   if(dyn) coords=rbind(coords,coords)
   sel_ss=numpoints=data_sel=sel_tt=XX=list()
  for(i in 1:Nloc){
      sel=out$nn.idx[i,];
      sel_ss[[i]]=coords[sel,]
      numpoints[[i]]=ncol(sel_ss[[i]])
      if(!is.null(data))data_sel[[i]]=matrix(data[,sel],nrow=2)
      if(!is.null(X))   XX[[i]]=X[c(sel,2*sel),]
                }
   }

}

##################################################################################
##################################################################################
##################################################################################
##################################################################################

if(length(XX)==0) XX=NULL
if(length(data_sel)==0) data_sel=NULL
return(list(data=data_sel,coordx=sel_ss,coordt=sel_tt,distance=distance, 
      numpoints=numpoints,numtime=numtime,radius=radius,spacetime=spacetime,X=XX))
} 
