GeoCovariogram2 <- function(fitted,vario,variogram=TRUE, ...)

{

if(!inherits(fitted,"GeoFit"))       stop("Enter an object obtained of clas  GeoFit\n")
if(!inherits(vario,"GeoVariogram"))       stop("Enter an object of class GeoVariogram\n")
space=!fitted$bivariate&&!fitted$spacetime

maxdist=vario$maxdist
h=seq(9e-5,maxdist,0.01)



if(space){

cc= GeoCorrFct(x=h, corrmodel=fitted$corrmodel, covariance=TRUE,variogram=variogram,
               param=append(fitted$param,fitted$fixed),model=fitted$model)
plot(cc,type="l",...)
points(vario$centers,vario$variograms,...)
box()
}
}