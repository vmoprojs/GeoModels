GeoNeighborhood = function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, distance="Eucl", grid=FALSE, 
                  loc, maxdist=NULL,maxtime=NULL, radius=6371, time=NULL, X=NULL)
{
  XX=NULL
  numtime=1
  sel_ss=1
  sel_tt=1
  if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
  if(!is.matrix(loc))   stop("loc parameter must be a matrix")
    #if(!(ncol(loc)==2))   stop("loc parameter must be a matrix  N X 2")
  spacetime=FALSE
  corrmodel="Exponential"
  
  if(!is.null(coordt)) {spacetime=TRUE;corrmodel="Exp_Exp"}
  if(spacetime) if(!is.vector(time))  stop("time parameter is missing")
  checkinput = CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting",
                             NULL, grid, 'Marginal', maxdist, maxtime, 'Gaussian', 1,
                              'Nelder-Mead', NULL, radius, NULL, NULL, NULL, 
                          'Pairwise', FALSE, 'SubSamp', FALSE, X)
  if(!is.null(checkinput$error)) stop(checkinput$error)
  spacetime_dyn=FALSE
  if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
    
## handling spatial coordinates
    if(is.null(coordy)) coords=as.matrix(coordx)
    else{
    if(grid) coords=as.matrix(expand.grid(coordx,coordy))
    else     coords=cbind(coordx,coordy)  
    }

############### total dimension #####
NN=nrow(coords);
Nloc=nrow(loc)

#####################################
sel_tt=NULL
colnames(loc)=NULL;colnames(coords)=NULL;
a_s=rbind(loc,coords)
## type of dist
if(distance=="Eucl") dd_s=as.matrix(dist(a_s))
if(distance=="Geod") dd_s=fields::rdist.earth(a_s,miles=F,R=radius)
if(distance=="Chor") dd_s=2*sin(fields::rdist.earth(a_s,miles=F,R=radius)/2)
space=!spacetime#&!bivariate

## spatial
if(space){
 sel=dd_s<maxdist
 ss=matrix((sel)[1:(Nloc),(Nloc+1):(NN+Nloc)],nrow=Nloc,ncol=NN)
 sel_ss=list()
 data_sel=list()
 XX=list()
 numpoints=double(Nloc)
 for(i in 1:Nloc)
 {
 sel_ss[[i]]=coords[ss[i,],]
 numpoints[i]=nrow(sel_ss[[i]])
 data_sel[[i]]=data[ss[i,]]
 if(!is.null(X)) XX[[i]]=X[ss[i,],]
}
 #if(dim(as.matrix(sel_ss,nrow=dimat2))[1]==0) stop("spatial distance for local kringing is too small")
}
## space-time
if(spacetime)
{
  ## spatial  part

  sel1=dd_s<maxdist
  ss=matrix((sel1)[1:(Nloc),(Nloc+1):(NN+Nloc)],nrow=Nloc,ncol=NN)
  sel_ss=list()
  data_sel1=X_sel1=list()
  
  
  for(i in 1:Nloc)  {sel_ss[[i]]=coords[ss[i,],];  
                    data_sel1[[i]]=as.matrix(data[,ss[i,]])
                     if(!is.null(X)) X_sel1[[i]]=as.matrix(X[,ss[i,]])
                    }
  ######## temporal  part
  TT=length(coordt)
  Tloc=length(time)
  
  a_t=c(time,coordt);
  dd_t=as.matrix(dist(a_t))
  sel2=dd_t<maxtime
  tt=matrix((sel2)[1:(Tloc),(Tloc+1):(TT+Tloc)],nrow=Tloc,ncol=TT)
  sel_tt=list()
  numtime=double(Tloc)
   for(i in 1:Tloc)  sel_tt[[i]]=coordt[which(tt[i,]>0)]
XX=list()
k=1
data_sel=list()
numtime=numpoints=double(Nloc*Tloc)
for(i in 1:Nloc){
for(j in 1:Tloc){
data_sel[[k]]=as.matrix(data_sel1[[i]][which(tt[j,]>0),])
 numtime[k]=nrow(data_sel[[k]])
 numpoints[k]=ncol(data_sel[[k]])
 XX[[k]]=
k=k+1
}}


}
return(list(data=data_sel,coordx=sel_ss,coordt=sel_tt,distance=distance, 
      numpoints=numpoints,numtime=numtime,radius=radius,spacetime=spacetime,X=XX))
} 
