GeoAniso <-function (coords, anisopars=c(0,1), inverse = FALSE)
{
  
    if(!is.matrix(coords))  stop(" coords is not a matrix ")
    ncc=ncol(coords)

    if(ncc!=2&&ncc!=3)  stop(" coords dimension is wrong ")
    
    if(ncc==3)
    {
    if(sum(coords[,3])==0)  {    ### anisotropy in three dimension but it is actually two dimensions 
    if (length(anisopars) != 2)
        stop(" anisopars parameters must be a vector with 2 elements (anisotropy angle and anisotropy ratio) ")
    angle <- anisopars[1]; stretch <- anisopars[2]
    if (stretch < 1) {
        stretch <- round(stretch, digits = 6)
        if (stretch < 1) stop("anisotropy ratio must be greater or equal  than 1 ")
    }
    if( !(angle>=0&&angle<=pi))    stop("the anisotropy angle  must be between 0 and pi ")
    rm <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),
        ncol = 2); rm=rbind(cbind(rm,c(0,0)),c(0,0,0))
    tm <- diag(c(1, 1/stretch)); tm=rbind(cbind(tm,c(0,0)),c(0,0,0))
    if (inverse) coordstransf <- coords %*% solve(rm %*% tm)
    else coordstransf  <- coords %*% rm %*% tm
    }
    else    {  coordstransf=1 #  ### anisotropy in three dimension  todo  ..  
              }    
   }
############################
   if(ncc==2)       ### anisotropy in  two dimensions 
   {  
    if (length(anisopars) != 2)
        stop(" anisopars parameters must be a vector with 2 elements (anisotropy angle and anisotropy ratio) ")
    angle <- anisopars[1]; stretch <- anisopars[2]
    if (stretch < 1) {
        stretch <- round(stretch, digits = 6)
        if (stretch < 1) stop("anisotropy ratio must be greater or equal  than 1 ")
    }
    if( !(angle>=0&&angle<=pi))    stop("the anisotropy angle  must be between 0 and pi ")
    rm <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)),
        ncol = 2);  tm <- diag(c(1, 1/stretch)); 
    if (inverse) coordstransf <- coords %*% solve(rm %*% tm)
    else coordstransf  <- coords %*% rm %*% tm
   }   
   ############################ 
    return(coordstransf)          
}