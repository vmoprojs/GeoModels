sp2Geo <- function(spobj,spdata = NULL)
{
  
  # spobj: sp  or spacetime object
  # spdata: name of dependent varible, 
  # one of the columns in @data slot must be the dependent variable
  
  X=Y=NULL
  coordt=NULL
  
  cls <- class(spobj)[1]
  
  # if(!(cls!="SpatialPoints" & cls!="SpatialPointsDataFrame")||(cls!="STFDF" & cls!="STIDF")) 
  #  stop("It is not a sp or st object ")
  
  if (!(cls %in% c("SpatialPoints","SpatialPointsDataFrame","STFDF","STIDF")) ) 
    stop("It is not a sp or st object ")
  
  
  if(cls=="SpatialPoints")
  {
    coords <- sp::coordinates(spobj)
    pj <- sp::is.projected(spobj)
    if(is.na(pj)){warning("coordinate reference system is not defined")}
  }
  if(cls=="SpatialPointsDataFrame")
  {
    if(is.null(spdata)){stop("Dependent variable name must be declared")}
    coords <- sp::coordinates(spobj)
    X <- spobj@data[,names(spobj@data)!=spdata]
    Y <- spobj@data[,spdata]
    pj <- sp::is.projected(spobj)
    if(is.na(pj)){warning("coordinate reference system is not defined")}
    
  }
  if(cls=="STFDF")
  {
    coords <- sp::coordinates(spobj@sp)
    pj <- sp::is.projected(spobj@sp)
    coordt <- seq(spobj@time)
    Y <- spobj@data[,spdata]
    X <- spobj@data[,names(spobj@data)!=spdata]
    if(is.na(pj)){warning("coordinate reference system is not defined")}
  }
  if(cls=="STIDF")
  {
    spobj <- as(spobj, "STFDF")
    coords <- sp::coordinates(spobj@sp)
    pj <- sp::is.projected(spobj@sp)
    coordt <- seq(spobj@time)
    Y <- spobj@data[,spdata]
    X <- spobj@data[,names(spobj@data)!=spdata]
    if(is.na(pj)){warning("coordinate reference system is not defined")}
  }
  
  return(list(coords = coords,coordt=coordt,X = X,Y=Y,projected = pj))
  
}
