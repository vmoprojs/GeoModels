


GeoNosymindices<- function(X,Y)
 {
#######################
fxy <- function(x,y, tol = 15){
 

    xx = atan(x)/(pi/2)
    yy = atan(y)/(pi/2)
   
    xdig = (as.numeric(strsplit(as.character(xx), "")[[1]][-(1:2)]))
    ydig = (as.numeric(strsplit(as.character(yy), "")[[1]][-(1:2)]))
    # length(xdig);length(ydig);
    if (length(xdig) < tol){
      xdig = (as.numeric(strsplit(as.character(xx-10^(-tol)), "")[[1]][-(1:2)]))
    }
    if (length(ydig) < tol){
      ydig = (as.numeric(strsplit(as.character(yy-10^(-tol)), "")[[1]][-(1:2)]))
    }
    xdig = xdig[1:tol]
    ydig = ydig[1:tol]
    if (y>=x){
     
      z = paste0(c('0.',as.vector(rbind(xdig,ydig))), collapse = "")
    }else{
      z = paste0(c('0.',as.vector(rbind(ydig,xdig))), collapse = "")
     
      }  

  return(z)
}
fxy <- Vectorize(fxy)

             
            xx=as.numeric(X[,1]); yy=as.numeric(X[,2])
            sol=fxy(xx,yy)
            ids <- !duplicated(sol)
            return(list(xy = X[ids,],d = Y[ids]))
 }
