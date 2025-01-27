####################################################
GeoNeighIndex<-function(coordx,coordy=NULL,coordz=NULL,coordt=NULL,coordx_dyn=NULL,
                              distance="Eucl",neighb=4,maxdist=NULL,maxtime=1,radius=6371,bivariate=FALSE)
{

#########################################
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
         
            return(list(xy = res,d = res_d[,2]))
 }
##########################################
nn2Geo <- function(x,y, K = 1,distance=0,maxdist=NULL,radius=6371)  
  {
            if(is.null(maxdist)) 
               {
               #nearest = RANN::nn2(x,y,k = K,treetype = c("kd"))} ### case neighboord
               nearest = nabor::knn(x,y,k = K)} ### case neighboord
            else     {



  if(distance==2||distance==1){
                                   prj=mapproj::mapproject(x[,1], x[,2], projection="sinusoidal") 
                                   x1=y1=radius*cbind(prj$x,prj$y)
                                   K=round(nrow(x)/2 )
                                   nearest = nabor::knn(x1,y1,radius = maxdist,k=K)
                               }     
else{

                     K=min(K-1,nrow(x)) # case of  maxdist 
                       nearest = nabor::knn(x,y,radius = maxdist,k=K  )}
                     }
            #########  cases geod (2) or chordal (1) distances :  to improve this  code!!
     
            if(distance==2||distance==1){
                  nn=nrow(x); 
                  nnd=ncol(nearest$nn.dists)
                  mm=matrix(0,nrow=nn,ncol=nnd)
                  for(i in 1:nn){   ## can we improve that?
                  si=nearest$nn.idx[i,];sel1=si[si>0]
                  a=fields::rdist.earth.vec(x1=matrix(x[i,],ncol=2),
                                            x2=matrix(x[sel1,],ncol=2), miles = FALSE, R = 1)
                  mm[i,][1:length(a)]=a
                
                 }

              mm[,1]=0 # just to be sure
             if(distance==2)  mm=radius*mm              # geodesic
             if(distance==1)  mm=2*radius*sin(0.5*mm)   # chordal  
             nearest$nn.dists=mm
             }
            ########################################### 
            sol = indices(nearest$nn.idx,nearest$nn.dists)
            if(is.null(maxdist)) {lags <- sol$d;rowidx <- sol$xy[,1];colidx <- sol$xy[,2]}
            if(!is.null(maxdist)){
                                    sel = sol$xy[,2]>0
                                    lags=sol$d[sel];rowidx <- sol$xy[,1][sel];colidx <- sol$xy[,2][sel]
                                 }
         return(list (lags=lags, rowidx = rowidx, colidx = colidx))
   }
##############################################################
spacetime_index=function(coords,coordx_dyn=NULL,N,K=4,coordt=NULL
                         ,numtime,maxtime=1,maxdist=NULL,distance="Eucl",radius=6371)
{
  ##############
  m_s=list();m_t=m_st=NULL;
  ##############         
  ## building marginal spatial indexes
  if(is.null(coordx_dyn)) 
  {
    inf=nn2Geo(coords,coords,K+1,distance,maxdist,radius)
    aa=cbind(inf$rowidx,inf$colidx)   ## spatial index (fixed coordinates) 
    for(i in 1:numtime) {
      # i = 1
      # repito las coordenadas numtime veces en elementos de una lista
      m_s[[i]]=data.frame(cbind(aa+N*(i-1),0,inf$lags))
      }
  }
 else  ## spatial index (dynamic coordinates)
  {   
           if(ncol(coords)==2)  ns=lengths(coordx_dyn)/2 
           if(ncol(coords)==3)  ns=lengths(coordx_dyn)/3
  for(i in 1:numtime){
    inf=nn2Geo(coordx_dyn[[i]],coordx_dyn[[i]],K+1,distance,maxdist,radius)
    aa=cbind(inf$rowidx,inf$colidx)
    m_s[[i]]=cbind( aa+ns[i]*(i-1),0,inf$lags)    
  }
  }
  ## building  temporal  and spatiotemporal indexes
  

    ## second way
      a=sort(unique(c(nabor::knn(coordt,coordt,k=round(maxtime)+1)$nn.dists)))
   nn=a[a>0]
  tnn=length(nn)   
  # sol <- NULL
  m_t <- list()
  m_st <- list()
  contador <- 1
  
  for(j in 1:tnn){
    for(k in 1:(numtime-tnn)){
      # j = 1;k = 1
      bb=nrow(m_s[[k]])
      m_t[[contador]] =data.frame(cbind( m_s[[k]][,1], m_s[[k+j]][,1], rep(nn[j],bb),rep(0,bb)) )
      m_st[[contador]]=data.frame(cbind( m_s[[k]][,1], m_s[[k+j]][,2], rep(nn[j],bb), m_s[[k]][,4]) )
      contador  <- contador +1
    }
  }
  ######
  # SS = data.table::rbindlist(m_s)
  # TT = data.table::rbindlist(m_t)
  # ST = data.table::rbindlist(m_st)
  # ##final space-time indexes and distances
  # final=data.table::rbindlist(list(SS,TT,ST))



  SS <- do.call(rbind, m_s)
  TT <- do.call(rbind, m_t)
  ST <- do.call(rbind, m_st)

  # Final space-time indexes and distances
  final <- do.call(rbind, list(SS, TT, ST))
  
  return(as.matrix(final))
}
##############################################################
bivariate_index=function(coords,coordx_dyn,N,K,maxdist,distance,radius)
  {
 
 if(length(K)==3) {K1=K[1];K2=K[2];K3=K[3]}
 else   K1=K2=K3=K

 
 if(length(maxdist)==3) {maxdist1=maxdist[1];maxdist2=maxdist[2];maxdist3=maxdist[3]}
else  maxdist1=maxdist2=maxdist3=maxdist




if(is.null(coordx_dyn)) 
   {  
         cm=coords[1:(N/2),]
         inf1=nn2Geo(cm,cm,K1+1,distance,maxdist1,radius)
         inf2=nn2Geo(cm,cm,K3+1,distance,maxdist3,radius)
         inf3=nn2Geo(cm,cm,K2+1,distance,maxdist2,radius)
         inf4=inf3
          aa1=cbind(inf1$rowidx,  inf1$colidx,0,0,inf1$lags)
          aa2=cbind(inf2$rowidx+N/2,inf2$colidx+N/2,1,1,inf2$lags)
          aa3=cbind(inf3$rowidx  ,inf3$colidx+N/2,0,1,inf3$lags)
          aa4=cbind(inf4$rowidx+N/2  ,inf4$colidx,1,0,inf4$lags)

          #aa5=cbind(rep(1:(N/2),K1),rep(1:(N/2),K1)+N,0,1,0)
          # aa6=cbind(rep(1:(N/2),K1)+N,rep(1:(N/2),K1),1,0,0)
           a5 =nabor::knn(cm,k = K2)
           aa5=cbind(a5$nn.idx[,1],a5$nn.idx[,1]+N/2,0,1,0)
           aa6=cbind(a5$nn.idx[,1]+N/2,a5$nn.idx[,1],1,0,0)
     
          SS= cbind(rbind(aa1,aa2,aa3,aa4,aa5,aa6))  

          
      # inf=nn2Geo(A0, A0,K2+1,distance,maxdist2,radius)
      # sel=!inf$lags
      # inf$lags[sel];inf$rowidx[sel];inf$colidx[sel]
      # aa5=cbind(inf$rowidx ,inf$colidx,    1,0,inf$lags)[1:(2*N),]
     #  print(A0)
     #  print(aa5)
        

         # aa1=cbind(inf1$rowidx     ,inf1$colidx,    0,0,inf1$lags)
         # aa2=cbind(inf2$rowidx+N/2 ,inf2$colidx+N/2,1,1,inf2$lags)
         # aa3=cbind(inf3$rowidx     ,inf3$colidx+N/2,0,1,inf3$lags)
         # aa4=cbind(inf4$rowidx+N/2,inf4$colidx,    1,0,inf4$lags)
         # SS= as.matrix(rbind(aa1,aa2,aa3,aa4))     
        
   }
else #### coordx_dyn case
  {       
     if(ncol(coords)==2)  ns=lengths(coordx_dyn)/2
     if(ncol(coords)==3)  ns=lengths(coordx_dyn)/3
    inf1=nn2Geo(coordx_dyn[[1]],  coordx_dyn[[1]],  K1+1,distance,maxdist1,radius)
    inf2=nn2Geo(coordx_dyn[[2]],  coordx_dyn[[2]],  K3+1,distance,maxdist3,radius)
    inf3=nn2Geo(coordx_dyn[[1]],  coordx_dyn[[2]],  K2+1,distance,maxdist2,radius)
    inf4=nn2Geo(coordx_dyn[[2]],  coordx_dyn[[1]],  K2+1,distance,maxdist2,radius)
          aa1=cbind(inf1$rowidx,      inf1$colidx,0,0,inf1$lags)
          aa2=cbind(inf2$rowidx+ns[1],inf2$colidx+ns[1],1,1,inf2$lags)
          aa3=cbind(inf3$rowidx  ,    inf3$colidx+ns[1],0,1,inf3$lags)
          aa4=cbind(inf4$rowidx+ns[1],inf4$colidx,1,0,inf4$lags)
          SS= as.matrix(rbind(aa1,aa2,aa3,aa4))   
  }
##final bivariate  indexes and distances
return(SS)
  }
##############################################################
######################################
######### start ######################
######################################

    spatial=TRUE
    spacetime=FALSE
    
    if(!is.null(coordx_dyn))  
                if(!is.list(coordx_dyn)) stop(" coordx_dyn must be a list")
    ### Check the parameters given in input
    # Checks if its a spatial or spatial-temporal random field:
    if(!is.null(coordt))
    if(is.numeric(coordt)&&is.numeric(maxtime)) 
    if(length(coordt)>1&&length(maxtime)>=1)  spacetime=TRUE
    distance=CheckDistance(distance)
    spatial= !spacetime&&!bivariate

K=neighb
## for spacetime or bivariate
#if(!is.null(coordx_dyn)){
      if(!is.null(coordy)){
               if(is.null(coordz)){
                          coordy <- coordx[,2]
                          coordx <- coordx[,1]
                          coords=cbind(coordx,coordy)
                                  }
                else {    coordz <- coordx[,3]
                          coordy <- coordx[,2]
                          coordx <- coordx[,1]
                              coords=cbind(coordx,coordy,coordz)
                     }        
                numcoord=nrow(coords)
           }
      else {
                       if(!bivariate)
                       { coords=coordx;  numcoord=nrow(coords) }
                       if(bivariate) {
                        if(is.null(coordx_dyn)) {coords=coordx; 
                                                 numcoord=nrow(coords) }
                        else {coords=1;numcoord=1}
                      }
           }  
     
#}             
##########################
########################## 

if(spatial)   #  spatial case
{
##########################################

  sol = nn2Geo(coords,coords,K+1 ,distance,maxdist,radius) ##### 
  gb=list(); gb$colidx=sol$colidx;
             gb$rowidx=sol$rowidx ;
             gb$lags=sol$lags
            # gb$lagt=NULL
             gb$maxdist=maxdist
             gb$neighb=neighb
} #### end spatial case 

##############################################   
if(spacetime)   #  space time  case
{ 
  numtime=length(coordt)
  sol=spacetime_index(coords,coordx_dyn, numcoord,K,coordt,numtime,maxtime,maxdist,distance,radius)
  gb=list(); gb$colidx=sol[,2];
             gb$rowidx=sol[,1] ;
             gb$lags=sol[,4]
             gb$lagt=sol[,3]
} 
if(bivariate)  { #space bivariate  case
   sol=bivariate_index(coords,coordx_dyn,numcoord,K,maxdist,distance,radius)

   gb=list(); gb$colidx=sol[,2];
             gb$rowidx=sol[,1] ;
             gb$lags=sol[,5]
             #gb$lagt=NULL
             gb$first=sol[,3]
             gb$second=sol[,4]
             gb$maxdist=maxdist
             gb$neighb=neighb
} 
return(gb)
}

