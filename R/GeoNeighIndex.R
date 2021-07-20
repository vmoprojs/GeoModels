
####################################################
### File name: GeoNeighIndex.r
####################################################

GeoNeighIndex<-function(coordx,coordy=NULL,coordx_dyn=NULL,coordt=NULL,
                              distance="Eucl",neighb=5,maxtime=1,radius=6371,bivariate=FALSE)
{

indices <- function(X,Y)
 {
             res = NULL;res_d = NULL
             for(i in 2:ncol(X))
             {
                sol = cbind(X[,1],X[,i])
                res = rbind(res,sol)
                sol_d = cbind(Y[,1],Y[,i])
                res_d = rbind(res_d,sol_d)
             }
            xx=as.numeric(res[,1]); yy=as.numeric(res[,2])
            sol=xx+yy+xx*yy+xx^2+yy^2
            ids <- !duplicated(sol)
            return(list(xy = res[ids,],d = res_d[ids,][,2]))
 }
##########################################
nn2Geo <- function(x, K = 1,distance,radius)  
  {
           
            nearest = RANN::nn2(x,k = K)
            #########  cases geod (2) or chordal (1) distances :  to improve this  code!!
            if(distance==2||distance==1){
                  agc=NULL;
                  nnn=nrow(x)
                  for(i in 1:nnn){   ## can we improve that?
                  a=fields::rdist.earth(matrix(x[i,],ncol=2), x[nearest$nn.idx[i,2:K],], miles = FALSE, R = 1)
                  agc=rbind(agc,a)
                 }
             if(distance==2)  agc=radius*agc   # geodesic
             if(distance==1)  agc=2*radius*sin(0.5*agc)   # chordal  
             nearest$nn.dists=cbind(rep(0,nnn),agc)
             }
            ########################################### 
            sol = indices(nearest$nn.idx,nearest$nn.dists)
            lags <- sol$d;rowidx <- sol$xy[,1];colidx <- sol$xy[,2]
         return(list (lags=lags, rowidx = rowidx, colidx = colidx))
   }
#########################
spacetime_index=function(coords,coordx_dyn,N,K,coordt,numtime,maxtime,distance,radius)
  {
##############
m_s=list();m_t=m_st=NULL;
##############         
## building marginal spatial indexes
if(is.null(coordx_dyn)) 
   {
        
         inf=nn2Geo(coords,K+1,distance,radius)
         for(i in 1:numtime) {
                  aa=cbind(inf$rowidx,inf$colidx)   ## spatial index (fixed coordinates)
                  m_s[[i]]=cbind(aa+N*(i-1),0,inf$lags)
                  }
   }
if(!is.null(coordx_dyn))
  {        ns=lengths(coordx_dyn)/2 
           for(i in 1:numtime){
                  inf=nn2Geo(coordx_dyn[[i]],K+1,distance,radius)
                  aa=cbind(inf$rowidx,inf$colidx)
                  m_s[[i]]=cbind( aa+ns[i]*(i-1),0,inf$lags)    ## spatial index (dynamic coordinates)
                  }
  }
         ## building  temporal  and spatiotemporal indexes
         ## temporal distances (not zero distance)
         nn=sort(unique(c(RANN::nn2(coordt,k=maxtime+1)$nn.dists)))[-1]  
         tnn=length(nn)   
         for(j in 1:tnn){
          for(k in 1:(numtime-tnn)){
            bb=nrow(m_s[[k]])
           m_t =rbind(m_t, cbind( m_s[[k]][,1], m_s[[k+j]][,1], rep(nn[j],bb)) )
           m_st=rbind(m_st,cbind( m_s[[k]][,1], m_s[[k+j]][,2], rep(nn[j],bb), m_s[[k]][,4]) )
         }}
        ######
        TT=cbind(m_t,rep(0,nrow(m_t)))
        SS=do.call(rbind,args=c(m_s));
        ST=m_st
        ##final space-time indexes and distances
        final=rbind(SS,TT,ST)
        return(final)
  }



##########################################

##########################
######### start ########## 
    spatial=TRUE
    spacetime=FALSE
    ### Check the parameters given in input
    # Checks if its a spatial or spatial-temporal random field:
    if(bivariate) coordt=c(0,1)
    if(!is.null(coordt))
    if(is.numeric(coordt)&&is.numeric(maxtime)) if(length(coordt)>1&&length(maxtime)>=1)  spacetime=TRUE
    distance=CheckDistance(distance)

K=neighb
## for spacetime or bivariate
      if(!is.null(coordy)){coordy <- coordx[,2]
                          coordx <- coordx[,1]
                          coords=cbind(coordx,coordy)
                          }
      else {
                       if(!bivariate) coords=coordx
                       if(bivariate)  coords=rbind(coordx,coordx)
           }  
      numcoord=nrow(coords)               
##########################
########################## 

if(!spacetime&&!bivariate)   #  spatial case
{
##########################################
  sol = nn2Geo(coords,K+1 ,distance,radius) ##### K or K+1
  gb=list(); gb$colidx=sol$colidx;
             gb$rowidx=sol$rowidx ;
             gb$lags=sol$lags
             gb$lagt=NULL
} #### end spatial case 

##############################################   
if(spacetime)   #  space time  case
{ 
  numtime=length(coordt)
  sol=spacetime_index(coords,coordx_dyn, numcoord,K,coordt,numtime,maxtime,distance,radius)
  
  gb=list(); gb$colidx=sol[,2];
             gb$rowidx=sol[,1] ;
             gb$lags=sol[,4]
             gb$lagt=sol[,3]
} #### end spacetime case


if(bivariate)  {
  #if(is.null(coordx_dyn)) sol = nn2Geo(coords,K+1 ,distance,radius) ##### K 
  #else {sol1 = nn2Geo(coordx_dyn[[1]],K+1 ,distance,radius) 
  #      sol2 = nn2Geo(coordx_dyn[[2]],K+1 ,distance,radius)
  #     }
  #gb=list();gb$colidx=c(sol$colidx

} #  spatial  bivariate case
return(gb)
}

