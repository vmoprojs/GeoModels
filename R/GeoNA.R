GeoNA=function(data,coordx,coordy=NULL, coordt=NULL, coordx_dyn=NULL, grid=FALSE, X=NULL,setting="spatial")
{

## handling spatial coordinates
   if(is.null(coordy)) coords=as.matrix(coordx)
    else{
    if(grid) coords=as.matrix(expand.grid(coordx,coordy))
    else     coords=cbind(coordx,coordy)  
    }
####    
N=length(c(unlist(data)))  #length of data


###########################
if (setting=='spatial'){
nasel=(is.nan(data)|is.infinite(data)|is.na(data))
perc=sum(nasel)/N; sel=!nasel
if(perc>0){
  data=data[sel]
  coords=coords[sel,]
  if(!is.null(X)) X=X[sel,]
}
}
############################
if (setting=='spacetime'){}


############################
if (setting=='bivariate'){
}
nonan = list( coordx = coords,
              coordy = coordy,
              coordt = coordt,
              coordx_dyn=coordx_dyn,              
              data=data,
              grid=grid,
              perc=perc,
              setting=setting,
              X=X)

return(nonan)

}