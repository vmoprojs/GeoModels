GeoNeighborhood = function(data=NULL, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,bivariate=FALSE, 
            distance="Eucl", grid=FALSE, loc, neighb=NULL,maxdist=NULL,maxtime=NULL,
                 radius=6371, time=NULL, X=NULL,M=NULL)
{
  #####################
    nbor=function(coords,loc,distance,maxdist,neighb){
  ### to improve!
  if(distance=="Geod"||distance=="Chor")
  {
    coords_p=coords; loc_p=loc
    prj=mapproj::mapproject(coords_p[,1], coords_p[,2], projection="sinusoidal") 
    coords_p=radius*cbind(prj$x,prj$y)
    prjloc=mapproj::mapproject(loc_p[,1], loc_p[,2], projection="sinusoidal")
    loc_p=radius*cbind(prjloc$x,prjloc$y)
  }
  ##################################################################################
  searchtype = "standard"
  if(!is.null(maxdist)){searchtype = "radius"} else {maxdist=0}
  if(is.null(neighb)) neighb=min(100, NN)
  
  ## computing neigh indexes
  if(distance=="Geod"||distance=="Chor")  {
                                           #out<- RANN::nn2(coords_p,loc_p, k=neighb,searchtype=searchtype,radius=maxdist)
                                           out<- nabor::knn(coords_p,loc_p, k=neighb,radius=maxdist)
                                           }
  if(distance=="Eucl")                    {
                                          # out<- RANN::nn2(coords,loc,     k=neighb,searchtype=searchtype,radius=maxdist)
                                           out<- nabor::knn(coords,loc, k=neighb,radius=maxdist)
                                          }
    return(out)                                       
    }
  #####################  
  
  XX=NULL
  MM=NULL
  numtime=1
  sel_ss=1
  sel_tt=1
  if(is.vector(loc))    loc=t(as.matrix(loc)) ## case of 1 location sites given as vector
  if(!is.matrix(loc))   stop("loc parameter must be a matrix")
  if(!is.logical(bivariate))   stop("bivariate must be logical")
  #if(!(ncol(loc)==2))   stop("loc parameter must be a matrix  N X 2")
  if(is.null(neighb)&&is.null(maxdist))     stop("maxdist (maxtime) and/or neighb must  be specified")
  if(!is.null(neighb)) neighb=round(neighb)

  spacetime=FALSE
  if(!is.null(coordt)) {spacetime=TRUE}
  if(spacetime) if(!is.vector(time))  stop("time parameter is missing")
  
  dyn=FALSE
  if(!is.null(coordx_dyn))  dyn=TRUE  

 

  ## handling spatial coordinates
  if(is.null(coordy)) {coords=as.matrix(coordx)}else{
    if(grid) {coords=as.matrix(expand.grid(coordx,coordy))}
    else    { coords=cbind(coordx,coordy)  }
  }
  
  Nloc=nrow(loc) # number of location sites
  NN=nrow(coords)
  #####################################
  sel_tt=NULL
  colnames(loc)=NULL;colnames(coords)=NULL;
  space=!spacetime&&!bivariate 
  ##################################################################





  #################################
  if(space){

    sel_ss=data_sel=numpoints=XX=MM=list()
    out= nbor(coords,loc,distance,maxdist,neighb)
    for(i in 1:Nloc)
    {
      ss=out$nn.idx[i,];
      sel_ss[[i]]=coords[ss,]
      numpoints[[i]]=nrow(sel_ss[[i]])
      if(!is.null(data)) data_sel[[i]]=data[ss]
      if(!is.null(M)) MM[[i]]=M[ss]
      if(!is.null(X)) XX[[i]]=X[ss,]
    }
  }
  #####################################################################################
  if(spacetime) 
  {


if(!dyn){  # not dynamic

    Tloc=length(time); TT=length(coordt)   
    out= nbor(coords,loc,distance,maxdist,neighb)
    if(!is.null(M)) { Mma=matrix(M,nrow=TT,ncol=NN,byrow=T)}
    if(!is.null(X)) { Xxa=list(); for(i in 1:ncol(X)) {Xxa[[i]]=matrix(X[,i],nrow=TT,ncol=NN,byrow=T)}}
  


    out_t <- nabor::knn(coordt,time,k=length(coordt),radius=maxtime) # with distance

    #out_t <- nabor::knn(coordt,time,k=round(maxtime)+1) #with neighboord 


    qq=out_t$nn.idx[1,];ss=which(qq!=0)
    out_t <- list(ind=cbind(as.vector(out_t$nn.idx[,ss]),time))

    sel_ss=numpoints=data_sel=sel_tt=XX=MM=list()
 
    k=1
    for(i in 1:Nloc)
    {
      for(j in 1:Tloc)
       {
       sel_s=out$nn.idx[i,];
       sel_ss[[k]]=matrix(coords[sel_s,],ncol=2)
       numpoints[[k]]=nrow(sel_ss[[i]])
       
          sel_t=out_t$ind[,1][out_t$ind[,2]==time[j]]
          pp=order(coordt[sel_t])
          sel_tt[[k]]=coordt[sel_t][pp]
          if(!is.null(data)) data_sel[[k]]=(data[sel_t,sel_s])[pp,]
          if(!is.null(M))   { QQ=Mma[sel_t,sel_s]; MM[[k]]=as.numeric(t(QQ[pp,]))  }
          if(!is.null(X))   {   
                              KK=NULL
                              for(l in 1:ncol(X)) 
                                  { AA=Xxa[[l]]; QQ=AA[sel_t,sel_s]; KK=cbind(KK,as.numeric( t(QQ[pp,] ))) }
                              XX[[k]]=KK       
                            }
       k=k+1
       }
    }
}  # end not dynamic
  if(dyn) print("to fix")  #dynamic

}
  #####################################################################################
  if(bivariate)   #### to fix case dyn
  {
    Nloc=nrow(loc)
    if(dyn) coords=rbind(coords[[1]],coords[[2]])
    sel_ss=numpoints=data_sel=sel_tt=XX=MM=list()
    if(dyn)      {out1= nbor(coords[[1]],loc,distance,maxdist,neighb)
                 out2= nbor(coords[[2]],loc,distance,maxdist,neighb)}
    else  {out= nbor(coords,loc,distance,maxdist,neighb)}
     
    for(i in 1:Nloc){
      sel=out$nn.idx[i,]
      sel_ss[[i]]=coords[sel,]
      numpoints[[i]]=ncol(sel_ss[[i]])
      if(!is.null(M)) {MM[[i]]=cbind(M[ss,1],M[ss,2])}
      if(!is.null(data))data_sel[[i]]=matrix(data[,sel],nrow=2)
      if(!is.null(X))   XX[[i]]=X[c(sel,2*sel),]
    }
  }
  ##################################################################################
  ##################################################################################
  ##################################################################################
  if(length(MM)==0)  MM=NULL
  if(length(XX)==0) XX=NULL
  if(length(data_sel)==0) data_sel=NULL
 
  return(list(data=data_sel,coordx=sel_ss,coordt=sel_tt,distance=distance, 
              numpoints=numpoints,numtime=numtime,radius=radius,spacetime=spacetime,X=XX,M=MM))
} 