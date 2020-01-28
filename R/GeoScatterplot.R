GeoScatterplot <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, 
                           distance="Eucl", grid=FALSE, maxdist=NULL, 
                           times=NULL, numbins=NULL, radius=6371, bivariate=FALSE)

{
    call <- match.call()
    model="Gaussian"
    corrmodel <- 'exponential'
    spatial=TRUE
    spacetime=FALSE
    maxtime=NULL

    ### Check the parameters given in input
    # Checks if its a spatial or spatial-temporal random field:
    if(bivariate) coordt=c(0,1)
    if(!is.null(coordt))
      if(is.numeric(coordt)&&is.numeric(times)) if(length(coordt)>1&&length(times)>=1) {corrmodel <- 'gneiting'; spacetime=TRUE}
    # Checks the input:
    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting", NULL, grid,
                             'None', maxdist, maxtime, model,NULL, 'Nelder-Mead', NULL,
                             radius,  NULL, NULL,NULL, 'GeoWLS', FALSE, 'SubSamp', FALSE,NULL)
    # Checks if there are errors in the input:
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
    ### END -- Specific checks of the Empirical Variogram    
    n=1
    initparam <- StartParam(coordx, coordy, coordt,coordx_dyn, corrmodel, data,distance, "Fitting",
                           NULL, grid, 'None', maxdist,
                           maxtime, model, n, NULL, NULL, FALSE, radius, 
                           NULL, NULL, NULL, 'GeoWLS', 'GeoWLS', FALSE,
                           'SubSamp', FALSE, 1, 1,1,1,NULL)
    ##################
    # 1. spatial case 
    ##################
    if(grid)     {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2];coords=cbind(coordx,coordy) }
    else coords = coordx
    numcoords   = nrow(coords)
    squares     = ceiling(sqrt(numbins))
    squares_2   = round(sqrt(numbins))
    bins = seq(0,maxdist,maxdist/numbins)
    n_pairs = (numcoords^2-length(numcoords))*0.5
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
   	   main = c(paste("cor =",signif(cor(v111,v222),3)),paste("(",signif(bins[i],3)," , ",signif(bins[i+1],3),"]"))
       plot(v111,v222,col = "#481567FF",xlab="",ylab="",main=main) 
       abline(0,1)
}




}