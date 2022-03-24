####################################################
### File name: GeoQQ.r
####################################################
GeoQQ<-function(fit,type="Q",add=FALSE,ylim=c(0,1),breaks=10,...)
{


if(class(fit)!="GeoFit") stop("A GeoFit object is needed as input\n")
if(type!="Q"&type!="D") stop("Type can be Q or D \n")

model=fit$model        #type of model
if(!is.null(fit$copula))copula=fit$copula
else copula=NULL

fit$param=unlist(fit$param)
fit$fixed=unlist(fit$fixed)
pp=c(fit$param,fit$fixed)
MM=as.numeric(pp["mean"])
VV=as.numeric(pp["sill"])

opar=par(no.readonly = TRUE)
########################################  
#### starting qq plot
########################################  
if(type=="Q") {


xlab="Theoretical Quantiles"
ylab="Sample Quantiles"






##########################################################
##########################################################
if(!fit$bivariate){

if(is.list(fit$coordx_dyn)) dd=unlist(fit$data)
else dd=c(t(fit$data))

N= length(dd)
probabilities= (1:N)/(N+1)
probabilities1= c(0.25,0.75)
q_e=quantile(dd,probabilities)
q_e1=quantile(dd,probabilities1)
#######################################
#######################################
if(!is.null(copula)){
if(copula %in% c("Gaussian","Clayton")){

#######################################  OK
if(model %in% c("Beta2")){
mm=1/(1+exp(-MM))
sh=as.numeric(pp["shape"])
pmin=as.numeric(pp["min"]);pmax=as.numeric(pp["max"]);
q_t=pmin+(pmax-pmin)*qbeta(probabilities,shape1=mm*sh,shape2=(1-mm)*sh)
q_t1=pmin+(pmax-pmin)*qbeta(probabilities1,shape1=mm*sh,shape2=(1-mm)*sh)
plot(q_t,q_e, main ="Beta qq-plot",...)
 }
####################################### 
if(model %in% c("Kumaraswamy2")){
mm=1/(1+exp(-MM))
sh=as.numeric(pp["shape"])
pmin=as.numeric(pp["min"]);pmax=as.numeric(pp["max"]);
aa=log(1-mm^(sh))/log(0.5)
shape1=log(0.5)/log1p(-mm^(sh));
q_t=pmin+(pmax-pmin)*(1-(1-(probabilities)^(sh))^(shape1))
q_t1=pmin+(pmax-pmin)*(1-(1-(probabilities1)^(sh))^(shape1))
plot(q_t,q_e, main ="Kumaraswamy qq-plot",...)
 }
####################################### 

}
}
else{    ### non copula models
if(model %in% c("Gaussian",                                           
              "Gaussian_misp_Binomial",
              "Gaussian_misp_Poisson","Gaussian_misp_BinomialNeg")) {
  q_t=qnorm(probabilities,mean=MM,sd=sqrt(VV))
  q_t1=qnorm(probabilities1,mean=MM,sd=sqrt(VV))
  plot(q_t,q_e,main="Gaussian qq-plot",...)
  #qqnorm(dd,main="Gaussian qq-plot",xlab=xlab,ylab=ylab)
}
if(model %in% c("Binomial")) {
 q_t=qbinom(probabilities, size=fit$n, prob=pnorm(pp["mean"]))
 q_t1=qbinom(probabilities1, size=fit$n, prob=pnorm(pp["mean"]))
 plot(q_t,q_e,main="Binomial qq-plot",...)
}
if(model %in% c("BinomialNeg")) {
 q_t=qnbinom(probabilities, size=fit$n, prob=pnorm(pp["mean"]))
 q_t1=qnbinom(probabilities1, size=fit$n, prob=pnorm(pp["mean"]))
 plot(q_t,q_e,main="Binomial Neg qq-plot",...)
}
if(model %in% c("Poisson")) {
 q_t=qpois(probabilities, lambda=exp(pp["mean"]))
 q_t1=qpois(probabilities1, lambda=exp(pp["mean"]))
 plot(q_t,q_e,main="Poisson qq-plot",...)
}
#######################################  OK
if(model %in% c("SkewGaussian"))
{
   omega=as.numeric(sqrt((pp["skew"]^2 + pp["sill"])/pp["sill"]))
   alpha=as.numeric(pp["skew"]/pp["sill"]^0.5)
   q_t=MM+sqrt(VV)*sn::qsn(probabilities,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
   q_t1=MM+sqrt(VV)*sn::qsn(probabilities1,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))  
  plot(q_t,q_e,main="Skew Gaussian qq-plot",...)
}
#######################################   OK
if(model%in%c("StudentT","Gaussian_misp_StudentT"))
{
  q_t=MM+sqrt(VV)*qt(probabilities,df=as.numeric(round(1/pp["df"])))
  q_t1=MM+sqrt(VV)*qt(probabilities1,df=as.numeric(round(1/pp["df"])))
  plot(q_t,q_e,main="t qq-plot",...)
      #limma::qqt(dd,df=as.numeric(round(1/pp["df"])),main="t qq-plot",xlab=xlab,ylab=ylab)
}
#######################################  OK
if(model%in%c("SkewStudentT","Gaussian_misp_SkewStudentT"))
{
  alpha=as.numeric(pp["skew"])
  nu=as.numeric(round(1/pp["df"]))
  q_t=MM+sqrt(VV)*sn::qst(probabilities, xi=0, omega=1, alpha=alpha, nu=nu)
  q_t1=MM+sqrt(VV)*sn::qst(probabilities1, xi=0, omega=1, alpha=alpha, nu=nu)
  plot(q_t,q_e,main="skewt qq-plot",...)
}
#######################################  OK
if(model %in% c("Weibull"))
{
   shape=pp["shape"]
   q_t=exp(MM)*qweibull(probabilities,shape=shape,scale=1/(gamma(1+1/shape )))
   q_t1=exp(MM)*qweibull(probabilities1,shape=shape,scale=1/(gamma(1+1/shape )))
   plot(q_t,q_e, main ="Weibull qq-plot ",...)
}
#######################################   OK
if(model %in% c("Gamma"))
{
   shape=pp["shape"]
   q_t=exp(MM)*qgamma(probabilities,shape=shape/2,rate=shape/2)
   q_t1=exp(MM)*qgamma(probabilities1,shape=shape/2,rate=shape/2)
   plot(q_t,q_e, main ="Gamma qq-plot ",...)
}
#######################################  OK
if(model %in% c("LogGaussian"))
{ 
   q_t = qlnorm(probabilities, exp(MM)-VV/2, sqrt(VV))
   q_t1 = qlnorm(probabilities1, exp(MM)-VV/2, sqrt(VV))
   plot(q_t,q_e,xlab=xlab,ylab=ylab,main = "LogGaussian qq-plot",...)
}
#######################################  OK
if(model %in% c("LogLogistic"))
{
shape=pp["shape"]
cc=gamma(1+1/shape)*gamma(1-1/shape)
q_t = actuar::qllogis(probabilities,shape = shape,scale=exp(MM)/cc)
q_t1 = actuar::qllogis(probabilities1,shape = shape,scale=exp(MM)/cc)
q_e=quantile(dd,probabilities)
plot(q_t,q_e,xlab=xlab,ylab=ylab,main = "LogLogistic qq-plot",...)
}
if(model %in% c("Logistic"))
{
q_t = qlogis(probabilities, location = MM, scale = sqrt(VV))
q_t1 = qlogis(probabilities1,location = MM, scale = sqrt(VV))
q_e=quantile(dd,probabilities)
plot(q_t,q_e,xlab=xlab,ylab=ylab,main = "Logistic qq-plot",...)
}
#######################################  OK
if(model %in% c("SinhAsinh"))
{
tail = as.numeric(pp["tail"])
skew = as.numeric(pp["skew"])
q_t=MM+sqrt(VV)*sinh(1/tail * asinh(qnorm(probabilities))+skew/tail)
q_t1=MM+sqrt(VV)*sinh(1/tail * asinh(qnorm(probabilities1))+skew/tail)
plot(q_t,q_e,main="Sas qq-plot",...)
}
#######################################   OK
if(model %in% c("Tukeyh"))
{
tail = as.numeric(pp["tail"])
uu=qnorm(probabilities);uu1=qnorm(probabilities1);
q_t=MM+sqrt(VV)*uu*exp(0.5*tail*uu^2);
q_t1=MM+sqrt(VV)*uu1*exp(0.5*tail*uu1^2)
plot(q_t,q_e,main="Tukey-h qq-plot",...)
}
#######################################  OK
if(model %in% c("Tukeyh2"))
{
qtpTukeyh221= function(x,tail1,tail2){
  ll=1:length(x)
  sel1=I(x>=0)*ll
  sel2=I(x<0)*ll
  x1=x[sel1];        
  x2=x[sel2]
  uu1<-qnorm(x1,0,1);
  uu2<-qnorm(x2,0,1)
  qq1=uu1*exp(0.5*tail1*uu1^2)
  qq2=uu2*exp(0.5*tail2*uu2^2)
  return(c(qq2,qq1))
}
tail1 = as.numeric(pp["tail1"]);tail2 = as.numeric(pp["tail2"])
q_t =MM+sqrt(VV)*qtpTukeyh221(probabilities,tail1,tail2)
q_t1 =MM+sqrt(VV)*qtpTukeyh221(probabilities1,tail1,tail2)
plot(q_t,q_e,main="Tukey-hh qq-plot",...)
}
#######################################  OK
if(model %in% c("Tukeygh","Gaussian_misp_Tukeygh"))
{
tail = as.numeric(pp["tail"]);skew = as.numeric(pp["skew"]);
uu=qnorm(probabilities);uu1=qnorm(probabilities1)
q_t=MM+sqrt(VV)*(exp(skew*uu)-1)*exp(0.5*tail*uu^2)/skew
q_t1=MM+sqrt(VV)*(exp(skew*uu1)-1)*exp(0.5*tail*uu1^2)/skew
plot(q_t,q_e,main="Tukey-gh qq-plot",...)
}
#######################################   OK
if(model %in% c("TwoPieceGaussian"))
{
qtpGaussian = function(x,skew){
  ll=1:length(x)
  sel1=I(x>0)*I(x<0.5*(1+skew))*ll
  sel2=I(x<=1)*I(x>=0.5*(1+skew))*ll
  x1=x[sel1]
  x2=x[sel2]
    qq1=(1+skew)*qnorm(x1/(1+skew))
    qq2=(1-skew)*qnorm((x2-skew)/(1-skew))
  return(c(qq1,qq2))
}
skew = as.numeric(pp["skew"])
q_t =MM+sqrt(VV)*qtpGaussian(probabilities,skew)
q_t1 =MM+sqrt(VV)*qtpGaussian(probabilities1,skew)
plot(q_t,q_e,main="Two-Piece Gaussian qq-plot",...)
}
#######################################  OK
if(model %in% c("TwoPieceBimodal"))
{
ptpbimodal = function(x,skew,delta,df){  
  alpha=2*(delta+1)/df
  nn=2^(1-alpha/2)
  ll=1:length(x)
  sel1=I(x<0)*ll
  sel2=I(x>=0)*ll
  x1=x[sel1];x2=x[sel2];
  pp1=(0.5*(1+skew)*as.numeric(zipfR::Igamma(df/2,nn*(-x1)^(alpha)/(2*((1+skew)^(alpha))),lower=FALSE))/(gamma(df/2)))
  pp2=(0.5*(skew+1)+(0.5*(1-skew)*as.numeric(zipfR::Igamma(df/2,nn*(x2)^(alpha)/(2*((1-skew)^(alpha))),lower=TRUE))/(gamma(df/2))))
  return(c(pp1,pp2))
}
skew = as.numeric(pp["skew"])
df   = as.numeric(pp["df"])
delta= as.numeric(pp["shape"])
f = function(x) ptpbimodal(x,skew = skew,delta=delta,df=df)
f.inv = GoFKernel::inverse(f,lower = -4,upper = 4)
q_t = MM+sqrt(VV)*sort(as.numeric(lapply(probabilities,f.inv)))
q_t1 = MM+sqrt(VV)*sort(as.numeric(lapply(probabilities1,f.inv)))
plot(q_t,q_e,main="Two-Piece Bimodal qq-plot",...)
}
#######################################  OK
if(model %in% c("TwoPieceStudentT"))
{
qtpt1 = function(x,skew,df){
  ll=1:length(x)
  sel1=I(x>0)*I(x<0.5*(1+skew))*ll
  sel2=I(x<=1)*I(x>=0.5*(1+skew))*ll
  x1=x[sel1]
  x2=x[sel2]
    qq1=(1+skew)*qt(x1/(1+skew),df=df)
    qq2= (1-skew)*qt((x2-skew)/(1-skew),df=df)
  return(c(qq1,qq2))
}
skew = as.numeric(pp["skew"])
df   = 1/as.numeric(pp["df"])
q_t =MM+sqrt(VV)*qtpt1(probabilities,skew,df)
q_t1 =MM+sqrt(VV)*qtpt1(probabilities1,skew,df)
plot(q_t,q_e,main="Two-Piece Student qq-plot",...)
}
#######################################   OK
if(model %in% c("TwoPieceTukeyh"))
{
qtukh=function(xx,tail)
{
  uu=qnorm(xx)
  q_t=uu*exp(0.5*tail*uu^2)
return(q_t)
}
qtptukey = function(x,skew,tail){
  ll=1:length(x)
  sel1=I(x>0)*I(x<0.5*(1+skew))*ll
  sel2=I(x<=1)*I(x>=0.5*(1+skew))*ll
  x1=x[sel1]
  x2=x[sel2]
    qq1=(1+skew)*qtukh(xx=x1/(1+skew),tail=tail)
    qq2= (1-skew)*qtukh(xx=(x2-skew)/(1-skew),tail=tail)
  return(c(qq1,qq2))
}
skew= as.numeric(pp["skew"])
tail= as.numeric(pp["tail"])
q_t =MM+sqrt(VV)*qtptukey(probabilities,skew,tail)
q_t1 =MM+sqrt(VV)*qtptukey(probabilities1,skew,tail)
plot(q_t,q_e,main="Two-Piece Tukey-h qq-plot",...)
}
#######
}
########################################
#aa=lm(q_e~1+q_t)
#abline(as.numeric(aa$coefficients[1]),as.numeric(aa$coefficients[2]))
## code from qqline of R
x <- q_t1
y<-  q_e1

slope <- diff(y)/diff(x);int <- y[1L] - slope * x[1L]

    abline(int, slope, pch=20)
}





##########################################################
##########################################################

if(fit$bivariate){
par(mfrow=c(1,2))

if(is.list(fit$coordx_dyn)){ dd1=fit$data[[1]];dd2=fit$data[[2]]}
else  {dd1=fit$data[1,];dd2=fit$data[2,];}
#N1= length(dd1);N2= length(dd2)
#probabilities1= (1:N1)/(N1+1); probabilities2= (1:N2)/(N2+1);


##########################################################

if(model %in% c("Gaussian")) { qqnorm(dd1,main="First Gaussian qq-plot");abline(0,1)
                               qqnorm(dd2,main="Second Gaussian qq-plot");abline(0,1)}

##########################################################

if(model %in% c("SkewGaussian"))
{
   omega1=sqrt((pp["skew_1"]^2 + pp["sill_1"])/pp["sill_1"])
   alpha1=pp["skew_1"]/pp["sill_1"]^0.5
   q_t1=sn::qsn(probabilities,xi=0,omega= as.numeric(omega1),alpha= as.numeric(alpha1))
   q_e1=quantile(dd1,probabilities)
   plot(q_t1,q_e1,main="First Skew Gaussian qq-plot",...)
   aa=lm(q_e1~1+q_t1)
   abline(as.numeric(aa$coefficients[1]),as.numeric(aa$coefficients[2]))

   omega2=sqrt((pp["skew_2"]^2 + pp["sill_2"])/pp["sill_2"])
   alpha2=pp["skew_2"]/pp["sill_2"]^0.5
   q_t2=sn::qsn(probabilities,xi=0,omega= as.numeric(omega2),alpha= as.numeric(alpha2))
   q_e2=quantile(dd2,probabilities)
   plot(q_t2,q_e2,main="Second Skew Gaussian qq-plot",...)
   aa=lm(q_e2~1+q_t2)
   abline(as.numeric(aa$coefficients[1]),as.numeric(aa$coefficients[2]))
}

##########################################################


  }
}

##########################################################
#### starting density plot##############################
##########################################################

if(type=="D") {
  if(!fit$bivariate){

if(is.list(fit$coordx_dyn)) dd=unlist(fit$data)
else dd=c(t(fit$data))


if(!is.null(copula)){
#################################
if(copula %in% c("Gaussian","Clayton")){

if(model %in% c("Beta2")){
mm=1/(1+exp(-MM))
sh=as.numeric(pp["shape"])
pmin=as.numeric(pp["min"]);pmax=as.numeric(pp["max"]);
ll=seq(min(dd),max(dd),  (max(dd)-min(dd))/100 )
ds=dbeta((ll-pmin)/(pmax-pmin),shape1=mm*sh,shape2=(1-mm)*sh)/(pmax-pmin)
if(!add) hist(dd,freq=F,xlim=c(pmin,pmax),xlab="",main="Beta Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...)
}

###########################################################  OK

if(model %in% c("Kumaraswamy2")){
sh=as.numeric(pp["shape"])
pmin=as.numeric(pp["min"]);pmax=as.numeric(pp["max"]);
ll=seq(min(dd),max(dd),0.1)
dkuma = function(x,MM,sh,pmin,pmax){ 
q=(x-pmin)/(pmax-pmin);k=1-q^(sh);
m1=1/(1+exp(-MM));
shapei=log(0.5)/log1p(-m1^(sh));
res=log(shapei)+log(sh)+(sh-1)*log(q)+(shapei-1)*log(k)-log(pmax-pmin);
return(exp(res))
}
ds=dkuma(ll,MM,sh,pmin,pmax)
if(!add) hist(dd,freq=F,xlim=c(pmin,pmax),xlab="",main="Kumaraswamy Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...)
}
}
}
else{

#######################################  OK
if(model%in%c("Gaussian"))
{
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Gaussian Histogram",ylim=ylim,breaks=breaks)
lines(seq(min(dd),max(dd),0.1),dnorm(seq(min(dd),max(dd),0.1),mean=MM,sd=sqrt(VV)),...)
}

#######################################  OK
if(model%in%c("StudentT","Gaussian_misp_StudentT")) 
{
df=as.numeric(round(1/pp["df"]))
ll=seq(min(dd),max(dd),0.1)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Student T Histogram",ylim=ylim,breaks=breaks)
lines(ll,dt((ll-MM)/sqrt(VV),df=df),...)/sqrt(VV)
}

#######################################   OK 
if(model%in%c("SkewStudentT","Gaussian_misp_SkewStudentT"))
{
  alpha=as.numeric(pp["skew"])
  nu=as.numeric(round(1/pp["df"]))
  ll=seq(min(dd),max(dd),0.1)
  d_st=sn::dst((ll-MM)/sqrt(VV), xi=0, omega=1, alpha=alpha, nu=nu)/sqrt(VV)
  if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Skew-T Histogram",ylim=ylim,breaks=breaks)
  lines(ll,d_st,...)
}

#######################################  OK
if(model %in% c("SkewGaussian"))       
{
   skew= as.numeric(pp["skew"])
   sill= as.numeric(pp["sill"])
   omega=sqrt((skew^2 + sill)/sill)
   alpha=skew/sill^0.5
   ll=seq(min(dd),max(dd),0.1)
   d_sn=sn::dsn((ll-MM)/sqrt(VV), xi=0, omega=omega,alpha=alpha)/sqrt(VV)
   if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Skew Gaussian Histogram",ylim=ylim,breaks=breaks)
   lines(ll,d_sn,...)
}


####################################### OK
if(model %in% c("Weibull"))
{
   shape=pp["shape"]
   ll=seq(min(dd),max(dd),0.1)
   d_w=dweibull(ll,shape=shape,scale=exp(MM)/(gamma(1+1/shape )))
   if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Weibull Histogram",ylim=ylim,breaks=breaks)
   lines(ll,d_w,...)
}

#######################################  OK
if(model %in% c("Gamma"))
{
   shape=pp["shape"]
   ll=seq(min(dd),max(dd),0.1)
   d_g=dgamma(ll,shape=shape/2,rate=shape/(2*exp(MM)))
   if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Gamma Histogram",ylim=ylim,breaks=breaks)
   lines(ll,d_g,...) 

}
#######################################  OK 
if(model %in% c("LogGaussian"))
{
   ll=seq(min(dd),max(dd),0.1)
   qtpsas=function(x,MM,VV){
   q=x*exp(VV/2);
   a=-0.5*(log(q)-MM)^2/VV-log(q)-log(sqrt(VV))-0.5*log(2*pi)+VV/2;
   return(exp(a))
   }
   d_l = dlnorm(ll,MM,VV)
   if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="LogGaussian Histogram",ylim=ylim,breaks=breaks)
   lines(ll,d_l,...) 
}
#######################################  OK
if(model %in% c("Logistic"))
{
ll=seq(min(dd),max(dd),0.1)
d_l = dlogis(ll,location = MM, scale = sqrt(VV))
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Logistic Histogram",ylim=ylim)
lines(ll,d_l,...) 
}
#######################################  OK
if(model %in% c("LogLogistic"))
{
shape=pp["shape"]
cc=gamma(1+1/shape)*gamma(1-1/shape)
ll=seq(min(dd),max(dd),0.1)
d_l = actuar::dllogis(ll,shape = shape,scale=exp(MM)/cc)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="LogLogistic Histogram",ylim=ylim,breaks=breaks)
lines(ll,d_l,...) 

}
####################################### 
if(model %in% c("SinhAsinh"))
{
tail = as.numeric(pp["tail"])
skew = as.numeric(pp["skew"])
ll=seq(min(dd),max(dd),0.1)
qtpsas1=function(x,skew,tail){
s=sinh(tail*asinh(x)-skew)
a=(2*pi*(1+x^2))^(-0.5)
c=sqrt(1+s^2)
d=exp(-0.5*s^2) 
return(tail*c*d*a)
}
ds=qtpsas1((ll-MM)/sqrt(VV),skew,tail)/sqrt(VV)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="SAS Histogram",ylim=ylim,breaks=breaks)
   lines(ll,ds,...) 
}

####################################### OK
if(model %in% c("Tukeyh"))
{
tail = as.numeric(pp["tail"])
ll=seq(min(dd),max(dd),0.1)
inverse_lamb=function(x,tail)
{
  value = sqrt(VGAM::lambertW(tail*x*x)/tail);
   return(sign(x)*value);
}

qth=function(x,tail,VV){
a= x*(1 + VGAM::lambertW(tail*x*x))
b=inverse_lamb(x,tail)
c=dnorm(inverse_lamb(x,tail),0,1)
return(b*c/(a*sqrt(VV)))
}
ds=qth((ll-MM)/sqrt(VV),tail,VV)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Tukey-h Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...) 
}


####################################### OK
if(model %in% c("Tukeyh2"))
{
tail1 = as.numeric(pp["tail1"]);tail2 = as.numeric(pp["tail2"])
ll=seq(min(dd),max(dd),0.1)
inverse_lamb=function(x,tail)
{
  value = sqrt(VGAM::lambertW(tail*x*x)/tail);
   return(sign(x)*value);
}

qtpTukeyh22= function(x,tail1,tail2,VV){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1];        
  x2=x[sel2]
  a1= x1*(1 + VGAM::lambertW(tail1*x1*x1))
  b1=inverse_lamb(x1,tail1)
  c1=dnorm(inverse_lamb(x1,tail1),0,1)
  ds1=b1*c1/(a1*sqrt(VV))
  a2= x2*(1 + VGAM::lambertW(tail2*x2*x2))
  b2=inverse_lamb(x2,tail2)
  c2=dnorm(inverse_lamb(x2,tail2),0,1)
  ds2=b2*c2/(a2*sqrt(VV))
  return(c(ds2,ds1))
}
ds=qtpTukeyh22((ll-MM)/sqrt(VV),tail1,tail2,VV)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Tukey-hh Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...) 
}


#######################################   OK
if(model %in% c("TwoPieceGaussian"))
{
skew = as.numeric(pp["skew"])
ll=seq(min(dd),max(dd),0.1)
qtpGaussian1=function(x,eta,VV){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1 =dnorm(x1/(1-eta),0,1)/sqrt(VV);       
  ds2 =dnorm(x2/(1+eta),0,1)/sqrt(VV); 
  return(c(ds2,ds1))
}
ds=qtpGaussian1((ll-MM)/sqrt(VV),skew,VV)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Two-Piece Gaussian Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...) 
}


####################################### OK
if(model %in% c("TwoPieceStudentT"))
{
skew = as.numeric(pp["skew"])
df   = 1/as.numeric(pp["df"])
ll=seq(min(dd),max(dd),0.1)
qtpt = function(x,skew,df,VV){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1 =dt(x1/(1-skew),df=df)/sqrt(VV);       
  ds2 =dt(x2/(1+skew),df=df)/sqrt(VV); 
  return(c(ds2,ds1))
}
ds=qtpt((ll-MM)/sqrt(VV),skew,df,VV)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Two-Piece Student Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...)
}


####################################### OK
if(model %in% c("TwoPieceTukeyh"))
{
skew= as.numeric(pp["skew"])
tail= as.numeric(pp["tail"])
ll=seq(min(dd),max(dd),0.1)

inverse_lamb=function(x,tail){
  value = sqrt(VGAM::lambertW(tail*x*x)/tail);
   return(sign(x)*value)
   }

dTukeyh= function(x,tail){
a= x*(1 + VGAM::lambertW(tail*x*x))
b=inverse_lamb(x,tail)
c=dnorm(inverse_lamb(x,tail),0,1)
return(res=b*c/a)
}

dTTukeyh= function(x,tail,skew,VV){
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1 =dTukeyh(x1/(1-skew),tail)/sqrt(VV);       
  ds2 =dTukeyh(x2/(1+skew),tail)/sqrt(VV);
  return(c(ds2,ds1))
}
ds=dTTukeyh((ll-MM)/sqrt(VV),tail,skew,VV)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Two-Piece Tukey-h Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...)
}

#######################################  OK
if(model %in% c("TwoPieceBimodal"))
{
skew = as.numeric(pp["skew"])
df   = as.numeric(pp["df"])
delta= as.numeric(pp["shape"])
ll=seq(min(dd),max(dd),0.1)
ptpbimodal1 = function(x,skew,delta,df,VV){  
  alpha=2*(delta+1)/df
  nn=2^(1-alpha/2)
  aa=1:length(x)
  sel1=I(x>=0)*aa
  sel2=I(x<0)*aa
  x1=x[sel1]        
  x2=x[sel2]
  ds1=0.5*alpha*x1^(alpha-1)*(1-skew)^(1-alpha)*dgamma((x1/(1-skew))^(alpha),shape=df,scale=1/nn)/sqrt(VV)
  ds2=0.5*alpha*(-x2)^(alpha-1)*(1+skew)^(1-alpha)*dgamma((-x2/(1+skew))^(alpha),shape=df,scale=1/nn)/sqrt(VV)
  return(c(ds2,ds1))
}
ds=ptpbimodal1((ll-MM)/sqrt(VV),skew,delta,df,VV)
if(!add) hist(dd,freq=F,xlim=c(min(dd),max(dd)),xlab="",main="Two-Piece Bimodal Histogram",ylim=ylim,breaks=breaks)
lines(ll,ds,...)
}




###############################################  OK
if(model %in% c("BinomialNeg")) {
ll=as.numeric(table(dd)/length(dd))
y=sort(unique(as.numeric(dd)))#min(dd):max(dd)
ds=dnbinom(y, size=fit$n, prob=pnorm(MM))
if(!add) plot(y,ll,type = "h", col = "blue",main="Binomial Negative Histogram", xlab="",  ylab="",lwd = 2)
points(y,ds,type = "p", col = "black", lwd = 3)
lines(y,ds)
}


###############################################  OK
if(model %in% c("Binomial")) {
ll=as.numeric(table(dd)/length(dd))
y=sort(unique(as.numeric(dd)))
ds=dbinom(y, size=fit$n, prob=pnorm(MM))
if(!add) plot(y,ll,type = "h", col = "blue",main="Binomial Histogram", xlab="",  ylab="",lwd = 2)
points(y,ds,type = "p", col = "black", lwd = 3)
lines(y,ds)
}


###############################################  OK
if(model %in% c("Poisson")) {
ll=as.numeric(table(dd)/length(dd))
y=sort(unique(as.numeric(dd)))
ds=dpois(y, lambda=exp(MM))
if(!add) plot(y,ll,type = "h", col = "blue",main="Poisson Histogram", xlab="",  ylab="",lwd = 2)
points(y,ds,type = "p", col = "black", lwd = 3)
lines(y,ds)
}

}







}
 
}

##########################################################
par(opar) 


}
