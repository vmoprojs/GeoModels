GeoPit=function(fit,type="Uniform")
{



if(!inherits(fit,"GeoFit"))  stop("A GeoFit object is needed as input\n")
if(!(type=="Uniform"||type=="Gaussian")) stop("The type parameter can be Uniform or Gaussian")


model=fit$model        #type of model
fit$param=unlist(fit$param)
fit$fixed=unlist(fit$fixed)
pp=c(fit$param,fit$fixed)
dd=fit$data


if(!fit$bivariate){

## spatial and spatio temporal case
#######################################   OK
if(model %in% c("Gaussian","Gaussian_misp_Binomial",
              "Gaussian_misp_Poisson","Gaussian_misp_BinomialNeg"))
{
mm=pp["mean"]
vv=pp["sill"]
data=pnorm(dd,mean=mm,sd=sqrt(vv))
}  
#######################################   OK
if(model %in% c("Weibull"))
{
mm=pp["mean"]
sh=pp["shape"]
data= pweibull(dd,shape=sh,scale=exp(mm)/(gamma(1+1/sh)))
} 
#######################################   OK

if(model %in% c("Beta2"))
{
MM=pp["mean"]
mm=1/(1+exp(-MM))
sh=pp["shape"]
pmin=pp["min"];pmax=pp["max"];
data=pbeta((dd-pmin)/(pmax-pmin),shape1=mm*sh,shape2=(1-mm)*sh)
}

#######################################   OK
if(model %in% c("Kumaraswamy2")){
MM=pp["mean"]
mm=1/(1+exp(-MM))
sh=as.numeric(pp["shape"])
pmin=as.numeric(pp["min"]);pmax=as.numeric(pp["max"]);
aa=log(1-mm^(sh))/log(0.5)
shape1=log(0.5)/log1p(-mm^(sh));
data=(1-(1-((dd-pmin)/(pmax-pmin))^(sh))^(shape1))
 }

#######################################   OK

if(model %in% c("SkewGaussian"))
{
   MM=pp["mean"]
   omega=as.numeric(sqrt((pp["skew"]^2 + pp["sill"])/pp["sill"]))
   alpha=as.numeric(pp["skew"]/pp["sill"]^0.5)
   data=sn::psn((dd-MM)/sqrt(pp["sill"]),xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
}

#######################################   OK
if(model%in%c("StudentT","Gaussian_misp_StudentT"))
{
  MM=pp["mean"]
  data=pt((dd-MM)/sqrt(pp["sill"]),df=as.numeric(round(1/pp["df"])))
}

#######################################   OK
if(model%in%c("SkewStudentT","Gaussian_misp_SkewStudentT"))
{
  MM=pp["mean"]
  alpha=as.numeric(pp["skew"])
  nu=as.numeric(round(1/pp["df"]))
  data=sn::pst((dd-MM)/sqrt(pp["sill"]), xi=0, omega=1, alpha=alpha, nu=nu)
}

#######################################   OK
if(model %in% c("Gamma"))
{
   MM=pp["mean"]
   shape=pp["shape"]
   data=pgamma(dd,shape=shape/2,rate=shape/(2*exp(MM)))
}
#######################################   revisar

if(model %in% c("LogGaussian"))    
{ 
   MM=pp["mean"]
   VV=pp["sill"]
  data = pnorm((dd-exp(MM)-VV/2)/sqrt(VV))
 # data = plnorm(dd, exp(MM)-VV/2, sqrt(VV))
}

#######################################   OK
if(model %in% c("LogLogistic"))
{
MM=pp["mean"]
shape=pp["shape"]
cc=gamma(1+1/shape)*gamma(1-1/shape)
data = actuar::pllogis(dd,shape = shape,scale=exp(MM)/cc)
}
#######################################   OK

if(model %in% c("Logistic"))   
{ 
  MM=pp["mean"]
  VV=pp["sill"]
  data = (1+exp(-(dd-MM)/sqrt(VV)))^(-1)
}

#######################################   OK

if(model %in% c("SinhAsinh"))
{
MM=pp["mean"]
VV=pp["sill"]
tail = as.numeric(pp["tail"])
skew = as.numeric(pp["skew"])
data=pnorm(sinh(tail *asinh((dd-MM)/sqrt(VV))-skew))
}

#######################################   OK
if(model %in% c("Tukeyh"))
{
MM=pp["mean"]
VV=pp["sill"]
tail = as.numeric(pp["tail"])
inverse_lamb=function(x,tail)
{
  value = sqrt(VGAM::lambertW(tail*x*x)/tail);
  return(sign(x)*value);
}
x=(dd-MM)/sqrt(VV)
data=pnorm(inverse_lamb(x,tail));
}

####################################### 
if(model %in% c("Tukeyh2"))
{
MM=pp["mean"]
VV=pp["sill"]
tail1 = as.numeric(pp["tail1"]);tail2 = as.numeric(pp["tail2"])
ll=seq(min(dd),max(dd),0.1)
inverse_lamb=function(x,tail)
{
  value = sqrt(VGAM::lambertW(tail*x*x)/tail);
   return(sign(x)*value);
}

pdfTukeyh22= function(x,tail1,tail2){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1];        
  x2=x[sel2]
  ds1=pnorm(inverse_lamb(x1,tail))
  ds2=pnorm(inverse_lamb(x2,tail))
  return(c(ds2,ds1))
}
x=(dd-MM)/sqrt(VV)
data=pdfTukeyh22(x,tail1,tail2)
}

#######################################   OK
if(model %in% c("TwoPieceGaussian"))
{
MM=pp["mean"]
VV=pp["sill"]
skew = as.numeric(pp["skew"])
ptpG=function(x,eta){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1 =eta+(1-eta)*pnorm(x1/(1-eta));       
  ds2 =(1+eta)*pnorm(x2/(1+eta)); 
  return(c(ds2,ds1))
}
x=(dd-MM)/sqrt(VV)
data=ptpG(x,skew) 
}

####################################### OK
if(model %in% c("TwoPieceStudentT"))
{
MM=pp["mean"]
VV=pp["sill"]
skew = as.numeric(pp["skew"])
df   = 1/as.numeric(pp["df"])
ptpt = function(x,skew,df){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1 =skew+(1-skew)*pt(x1/(1-skew),df=df);       
  ds2 =(1+skew)*pt(x2/(1+skew),df=df); 
  return(c(ds2,ds1))
}
x=(dd-MM)/sqrt(VV)
data=ptpt(x,skew,df)
}

####################################### OK
if(model %in% c("TwoPieceTukeyh"))
{
MM=pp["mean"]
VV=pp["sill"]
skew= as.numeric(pp["skew"])
tail= as.numeric(pp["tail"])

inverse_lamb=function(x,tail)
{
  value = sqrt(VGAM::lambertW(tail*x*x)/tail);
  return(sign(x)*value);
}

pTukeyh= function(x,tail){
c=pnorm(inverse_lamb(x,tail))
return(res=c)
}

pTTukeyh= function(x,tail,skew){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1 =skew+(1-skew)*pTukeyh(x1/(1-skew),tail);       
  ds2 =(1+skew)*pTukeyh(x2/(1+skew),tail);
  return(c(ds2,ds1))
}
x=(dd-MM)/sqrt(VV)
data=pTTukeyh(x,tail,skew)
}

#######################################  OK
if(model %in% c("TwoPieceBimodal"))
{
MM=pp["mean"]
VV=pp["sill"]
skew = as.numeric(pp["skew"])
df   = as.numeric(pp["df"])
delta= as.numeric(pp["shape"])

pdfbimodal = function(x,skew,delta,df){  
  alpha=2*(delta+1)/df
  nn=2^(1-alpha/2)
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1=0.5*(1+skew)+0.5*(1-skew)*pgamma((x1/(1-skew))^(alpha),shape=df,scale=1/nn)
  ds2=0.5*(1+skew)*(1-pgamma((-x2/(1+skew))^(alpha),shape=df,scale=1/nn))
  return(c(ds2,ds1))
}
x=(dd-MM)/sqrt(VV)
data=pdfbimodal(x,skew,delta,df)
}

#######################################  discrete
if(model %in% c("Binomial")) {
   MM=pp["mean"]
   data=pbinom(dd, size=fit$n, prob=pnorm(MM))
}
if(model %in% c("BinomialNeg")) {
   MM=pp["mean"]
   data=pnbinom(dd, size=fit$n, prob=pnorm(MM))
}
if(model %in% c("Poisson")) {
   MM=pp["mean"]
   data=ppois(dd, lambda=exp(MM))
}


}
else{
stop("The spatial bivariate case is not implemented yet")
}



if(type=="Gaussian") data=qnorm(data)
fit$data=data
return(fit)
}