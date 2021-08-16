####################################################
### File name: GeoScatterplot.r
####################################################

GeoScatterplot <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, 
                           distance="Eucl", grid=FALSE, maxdist=NULL, neighb=NULL,
                           times=NULL, numbins=4, radius=6371, bivariate=FALSE)

{
    call <- match.call()
    opar=par(no.readonly = TRUE)
    model="Gaussian"
    corrmodel <- 'exponential'
    spatial=TRUE
    spacetime=FALSE
    maxtime=ns=NULL
    ### Check the parameters given in input
    # Checks if its a spatial or spatial-temporal random field:
    if(bivariate) {coordt=c(0,1);corrmodel="Bi_Matern_sep"}
    if(!is.null(coordt))
      if(is.numeric(coordt)&&is.numeric(times)) if(length(coordt)>1&&length(times)>=1) {corrmodel <- 'gneiting'; spacetime=TRUE}
    # Checks the input:
    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting", NULL, grid,
                             'None', maxdist, maxtime, model,NULL, 'Nelder-Mead', NULL,
                             radius,  NULL, NULL,NULL, 'GeoWLS', FALSE, 'SubSamp', FALSE,NULL,NULL)

    # Checks if there are errors in the input:
  #if(is.null(maxdist))         stop('maxdist must be specified\n')
  if(!is.null(checkinput$error))
      stop(checkinput$error)
  if(spatial||bivariate){
         if(!is.null(numbins) & !is.integer(numbins))
           if(numbins < 0)
              stop('insert a positive integer value for the number of bins\n')
              }
  if(spacetime)
     {
     if(length(times)>3)   stop('insert no more than three temporal  instants\n')
     if(maxtime<min(dist((coordt))))
                     stop('maximum temporal distance is too small\n')
     }

         n=1
    if(is.null(neighb)) {initparam <- StartParam(coordx, coordy, coordt,coordx_dyn, corrmodel, data,distance, "Fitting",
                           NULL, grid, 'None', maxdist,neighb,
                           maxtime, model, n, NULL, NULL, FALSE, radius, 
                           NULL, NULL, NULL, 'GeoWLS', 'GeoWLS', FALSE,
                           'SubSamp', FALSE, 1, 1,1,1,NULL,NULL,FALSE)}

  
    ### END -- Specific checks of the Empirical Variogram    
    if(grid)     {a=as.matrix(expand.grid(coordx,coordy));coordx=a[,1];coordy=a[,2] }

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

##################################################    
##### scatterplot based on distances 
##################################################
if(is.null(neighb)) {


#########################
# spatial univariate case 
########################
if(!bivariate&&!spacetime)
 {
  
    numcoords   = nrow(coords)
    squares     = ceiling(sqrt(numbins))
    squares_2   = round(sqrt(numbins))
    if(length(maxdist)==1) bins = seq(0,maxdist,maxdist/numbins)
    else                   bins=c(0,maxdist)
    n_pairs = (numcoords*(numcoords-1))*0.5
    v0 = rep(-1,n_pairs)
    v1 = rep(exp(-99),n_pairs)
    v2 = rep(exp(-99),n_pairs)


    V  = .C("pairs",as.integer(numcoords),as.double(data),as.double(coords[,1]),as.double(coords[,2]),as.double(numbins), 
            as.double(bins),as.double(v0),as.double(v1),as.double(v2),as.double(maxdist),PACKAGE='GeoModels', DUP = TRUE,NAOK = TRUE) 
    
    v0 = as.numeric(unlist(V[7]));v0 = v0[v0 != -1]
    v1 = as.numeric(unlist(V[8]));v1 = v1[v1 != exp(-99)]
    v2 = as.numeric(unlist(V[9]));v2 = v2[v2 != exp(-99)]        
    vvv = data.frame(v0,v1,v2)

    par(mfrow = c(squares_2,squares))
    for(i in 1:numbins){
       v111 = vvv[vvv$v0 == bins[i],2]
       v222 = vvv[vvv$v0 == bins[i],3]

       if(length(v111) == 0 || length(v222) == 0  ){print(paste("There is no points between ",signif(bins[i],3),"and ", signif(bins[i+1],3)))}
       else{
          main = c(paste("(",signif(bins[i],3)," , ",signif(bins[i+1],3),"]"))
         #main = c(paste("cor =",signif(cor(v111,v222),3)),paste("(",signif(bins[i],3)," , ",signif(bins[i+1],3),"]"))
         plot(v111,v222,col = "#481567FF",xlab="",ylab="",main=main) 
         abline(0,1)
          }
    }
}
#########################
# spatial bivariate case 
########################
if(bivariate)
 {
  if(length(maxdist)==1) maxdist=c(maxdist,maxdist)
  squares     = ceiling(sqrt(numbins))
  squares_2   = round(sqrt(numbins))
  data1=data[1:ns[1]]
  data2=data[(ns[1]+1):(ns[1]+ns[2])]
  coords1=coords[1:ns[1],]
  coords2=coords[(ns[1]+1):(ns[1]+ns[2]),]
 
    numcoords1   = nrow(coords1)
    bins1 = seq(0,maxdist[1],maxdist[1]/numbins)
    n_pairs1 = (numcoords1*(numcoords1-1))*0.5
    v01 = rep(-1,n_pairs1)
    v11 = rep(exp(-99),n_pairs1)
    v21 = rep(exp(-99),n_pairs1)
 V1  = .C("pairs",as.integer(numcoords1),as.double(data1),as.double(coords1[,1]),as.double(coords1[,2]),as.double(numbins), 
            as.double(bins1),as.double(v01),as.double(v11),as.double(v21),as.double(maxdist[1]),PACKAGE='GeoModels', DUP = TRUE,NAOK = TRUE) 
 ################################################################
    numcoords2 = nrow(coords2)
    bins2 = seq(0,maxdist[2],maxdist[2]/numbins)
    n_pairs2 = (numcoords2*(numcoords2-1))*0.5
    v02 = rep(-1,n_pairs2)
    v12 = rep(exp(-99),n_pairs2)
    v22 = rep(exp(-99),n_pairs2)
 V2  = .C("pairs",as.integer(numcoords2),as.double(data2),as.double(coords2[,1]),as.double(coords2[,2]),as.double(numbins), 
            as.double(bins2),as.double(v02),as.double(v12),as.double(v22),as.double(maxdist[2]),PACKAGE='GeoModels', DUP = TRUE,NAOK = TRUE) 
   
    v01 = as.numeric(unlist(V1[7]));v01 = v01[v01 != -1]
    v11 = as.numeric(unlist(V1[8]));v11 = v11[v11 != exp(-99)]
    v21 = as.numeric(unlist(V1[9]));v21 = v21[v21 != exp(-99)]        
    vvv1 = data.frame(v01,v11,v21)

    v02 = as.numeric(unlist(V2[7]));v02 = v02[v02 != -1]
    v12 = as.numeric(unlist(V2[8]));v12 = v12[v12 != exp(-99)]
    v22 = as.numeric(unlist(V2[9]));v22 = v22[v22 != exp(-99)]        
    vvv2 = data.frame(v02,v12,v22)

##############
    par(mfrow = c(squares_2,squares*2))
    for(i in 1:numbins){
       v111_1 = vvv1[vvv1$v01 == bins1[i],2]
       v222_1 = vvv1[vvv1$v01 == bins1[i],3]
       main = c(paste("cor =",signif(cor(v111_1,v222_1),3)),paste("(",signif(bins1[i],3)," , ",signif(bins1[i+1],3),"]"))
       plot(v111_1,v222_1,col = "#481567FF",xlab="",ylab="",main=main) 
       abline(0,1)}
    for(i in 1:numbins){
       v111_2 = vvv2[vvv2$v02 == bins2[i],2]
       v222_2 = vvv2[vvv2$v02 == bins2[i],3]
       main = c(paste("cor =",signif(cor(v111_2,v222_2),3)),paste("(",signif(bins2[i],3)," , ",signif(bins2[i+1],3),"]"))
       plot(v111_2,v222_2,col = "#481567FF",xlab="",ylab="",main=main) 
       abline(0,1)
    }
 }
#########################
# spatio temporal  case 
########################
 if(spacetime){}
}

##################################################    
##### scatterplot based on neighboord
##################################################

if(!is.null(neighb)) {


#########################
# spatial univariate case 
########################
ln=length(neighb)
squares     = ceiling(sqrt(ln))
squares_2   = round(sqrt(ln))


if(!bivariate&&!spacetime)
 {
 par(mfrow = c(squares_2,squares))
for(i in 1:ln){
  sel=GeoNeighIndex(coordx=coords,neighb=neighb[i])  
  data1=data[sel$colidx]; data2=data[sel$rowidx]
   plot(data1,data2,col = "#481567FF",xlab="",ylab="",main =paste("Neighb =",neighb[i])) 
  }

 }   
#########################
# spatial bivariate case 
########################
 if(bivariate){}


#########################
# spatio temporal  case 
########################
 if(spacetime){}

}
 par(opar)
###### end scatterplot based on neighb
}
