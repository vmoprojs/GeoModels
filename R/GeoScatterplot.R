GeoScatterplot <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL, distance="Eucl",
                       grid=FALSE, maxdist=NULL, times=NULL, numbins=NULL,
                       radius=6371,bivariate=FALSE)

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

print("ciao")

}