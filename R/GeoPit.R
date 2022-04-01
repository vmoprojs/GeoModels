GeoPit=function(fit,type="Uniform")
{
if(class(fit)!="GeoFit") stop("A GeoFit object is needed as input\n")
if(!(type=="Uniform"||type=="Gaussian")) stop("The type parameter can be Uniform or Gaussian")


model=fit$model        #type of model
fit$param=unlist(fit$param)
fit$fixed=unlist(fit$fixed)
pp=c(fit$param,fit$fixed)
dd=fit$data


if(!fit$bivariate){

## spatial and spatio temporal case
if(model %in% c("Gaussian"))
{
mm=pp["mean"]
vv=pp["sill"]
data=pnorm(dd,mean=mm,sd=sqrt(vv))
}  

if(model %in% c("Weibull"))
{
mm=pp["mean"]
sh=pp["shape"]
data= pweibull(dd,shape=sh,scale=exp(mm)/(gamma(1+1/sh)))
}


if(model %in% c("Beta2"))
{
MM=pp["mean"]
mm=1/(1+exp(-MM))
sh=pp["shape"]
pmin=pp["min"];pmax=pp["max"];
data=pbeta((dd-pmin)/(pmax-pmin),shape1=mm*sh,shape2=(1-mm)*sh)
}


if(model %in% c("Kumaraswamy2")){
MM=pp["mean"]
mm=1/(1+exp(-MM))
sh=as.numeric(pp["shape"])
pmin=as.numeric(pp["min"]);pmax=as.numeric(pp["max"]);
aa=log(1-mm^(sh))/log(0.5)
shape1=log(0.5)/log1p(-mm^(sh));
data=(1-(1-((dd-pmin)/(pmax-pmin))^(sh))^(shape1))
 }



if(model %in% c("SkewGaussian"))
{
   MM=pp["mean"]
   omega=as.numeric(sqrt((pp["skew"]^2 + pp["sill"])/pp["sill"]))
   alpha=as.numeric(pp["skew"]/pp["sill"]^0.5)
   data=sn::psn((dd-MM)/sqrt(pp["sill"]),xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
}


if(model%in%c("StudentT","Gaussian_misp_StudentT"))
{
  MM=pp["mean"]
  data=pt((dd-MM)/sqrt(pp["sill"]),df=as.numeric(round(1/pp["df"])))
}


if(model%in%c("SkewStudentT","Gaussian_misp_SkewStudentT"))
{
  MM=pp["mean"]
  alpha=as.numeric(pp["skew"])
  nu=as.numeric(round(1/pp["df"]))
  data=sn::pst((dd-MM)/sqrt(pp["sill"]), xi=0, omega=1, alpha=alpha, nu=nu)
}


if(model %in% c("Gamma"))
{
   MM=pp["mean"]
   shape=pp["shape"]
   data=pgamma(dd,shape=shape/2,rate=shape/(2*exp(MM)))
}


if(model %in% c("LogGaussian"))    ##revisar
{
   MM=pp["mean"]
   VV=pp["sill"]
   data = plnorm(dd, exp(MM)-VV/2, sqrt(VV))
}


if(model %in% c("LogLogistic"))
{
MM=pp["mean"]
shape=pp["shape"]
cc=gamma(1+1/shape)*gamma(1-1/shape)
data = actuar::pllogis(dd,shape = shape,scale=exp(MM)/cc)
}


if(model %in% c("Logistic"))  
{
  MM=pp["mean"]
  VV=pp["sill"]
  data = (1+exp(-(dd-MM)/sqrt(VV)))^(-1)
}




if(model %in% c("SinhAsinh"))
{
MM=pp["mean"]
VV=pp["sill"]
tail = as.numeric(pp["tail"])
skew = as.numeric(pp["skew"])
data=pnorm(sinh(tail *asinh((dd-MM)/sqrt(VV))-skew))
}






}
else{
stop("The spatial bivariate case is not implemented yet")
}



if(type=="Gaussian") data=qnorm(data)
fit$data=data
return(fit)
}