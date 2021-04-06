

####################################################
### File name: GeoQQ.r
####################################################
GeoQQ<-function(fit)
{

if(class(fit)!="GeoFit") stop("A GeoFit object is needed as input\n")
model=fit$model        #type of model

xlab="Theoretical Quantiles"
ylab="Sample Quantiles"


pp=c(fit$param,fit$fixed)

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
if(model %in% c("Gaussian")) {
  q_t=qnorm(probabilities)
  q_t1=qnorm(probabilities1)
  plot(q_t,q_e,main="Gaussian qq-plot",xlab=xlab,ylab=ylab,pch=20)
  #qqnorm(dd,main="Gaussian qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("SkewGaussian"))
{
   omega=as.numeric(sqrt((pp["skew"]^2 + pp["sill"])/pp["sill"]))
   alpha=as.numeric(pp["skew"]/pp["sill"]^0.5)
   q_t=sn::qsn(probabilities,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha)) 
   q_t1=sn::qsn(probabilities1,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))  
  plot(q_t,q_e,main="Skew Gaussian qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model%in%c("StudentT","Gaussian_misp_StudentT")) 
{
  q_t=qt(probabilities,df=as.numeric(round(1/pp["df"])))
  q_t1=qt(probabilities1,df=as.numeric(round(1/pp["df"])))
  plot(q_t,q_e,main="t qq-plot",xlab=xlab,ylab=ylab,pch=20)
      #limma::qqt(dd,df=as.numeric(round(1/pp["df"])),main="t qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model%in%c("SkewStudentT","Gaussian_misp_SkewStudentT")) 
{
  alpha=as.numeric(pp["skew"])
  nu=as.numeric(round(1/pp["df"]))
  q_t=sn::qst(probabilities, xi=0, omega=1, alpha=alpha, nu=nu)
  q_t1=sn::qst(probabilities1, xi=0, omega=1, alpha=alpha, nu=nu)
  plot(q_t,q_e,main="skewt qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("Weibull"))
{
   shape=pp["shape"]
   q_t=qweibull(probabilities,shape=shape,scale=1/(gamma(1+1/shape )))
   q_t1=qweibull(probabilities1,shape=shape,scale=1/(gamma(1+1/shape )))
   plot(q_t,q_e, main ="Weibull qq-plot ",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("Gamma"))
{
   shape=pp["shape"]
   q_t=qgamma(probabilities,shape=shape/2,rate=shape/2)
   q_t1=qgamma(probabilities1,shape=shape/2,rate=shape/2)
   plot(q_t,q_e, main ="Gamma qq-plot ",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("LogGaussian"))
{
   SS=pp["sill"]; mm=pp["mean"]
   q_t = qlnorm(probabilities, mm/exp(SS/2), sqrt(SS/exp(SS)))
   q_t1 = qlnorm(probabilities1, mm/exp(SS/2), sqrt(SS/exp(SS)))
   plot(q_t,q_e,xlab=xlab,ylab=ylab,main = "LogGaussian qq-plot",pch=20)
}
#######################################
if(model %in% c("LogLogistic"))
{
shape=pp["shape"]
cc=gamma(1+1/shape)*gamma(1-1/shape)
q_t = actuar::qllogis(probabilities,shape = shape,scale=1/cc)
q_t1 = actuar::qllogis(probabilities1,shape = shape,scale=1/cc)
q_e=quantile(dd,probabilities)
plot(q_t,q_e,xlab=xlab,ylab=ylab,main = "LogLogistic qq-plot",pch=20) 
}
#######################################
if(model %in% c("SinhAsinh"))
{
tail = as.numeric(pp["tail"])
skew = as.numeric(pp["skew"])
q_t=sinh(1/tail * asinh(qnorm(probabilities))+skew/tail)
q_t1=sinh(1/tail * asinh(qnorm(probabilities1))+skew/tail)
plot(q_t,q_e,main="Sas qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("Tukeyh"))
{
tail = as.numeric(pp["tail"])
uu=qnorm(probabilities);uu1=qnorm(probabilities1);
q_t=uu*exp(0.5*tail*uu^2);
q_t1=uu1*exp(0.5*tail*uu1^2)
plot(q_t,q_e,main="Tukey-h qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("Tukeyh2"))
{
qtpTukeyh22= function(x,tail1,tail2){ 
  ll=1:length(x)
  sel1=I(x>0)*ll
  sel2=I(x<=0)*ll
  x1=x[sel1];         
  x2=x[sel2] 
  uu1<-qnorm(x1,0,1); 
  uu2<-qnorm(x2,0,1)
  qq1=uu1*exp(0.5*tail1*uu1^2)
  qq2=uu2*exp(0.5*tail2*uu2^2)
  return(c(qq1,qq2))
}
tail1 = as.numeric(pp["tail1"]);tail2 = as.numeric(pp["tail2"])
q_t =qtpTukeyh22(probabilities,tail1,tail2)
q_t1 =qtpTukeyh22(probabilities1,tail1,tail2)
plot(q_t,q_e,main="Tukey-hh qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("Tukeygh","Gaussian_misp_Tukeygh")) 
{
tail = as.numeric(pp["tail"]);skew = as.numeric(pp["skew"]);
uu=qnorm(probabilities);uu1=qnorm(probabilities1)
q_t=(exp(skew*uu)-1)*exp(0.5*tail*uu^2)/skew
q_t1=(exp(skew*uu1)-1)*exp(0.5*tail*uu1^2)/skew
plot(q_t,q_e,main="Tukey-gh qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("TwoPieceGaussian"))
{
qtpGaussian = function(x,skew){
  ll=1:length(x)
  sel1=I(x>0)*I(x<0.5*(1+skew))*ll
  sel2=I(x<=1)*I(x>=0.5*(1+skew))*ll
  x1=x[sel1]
  x2=x[sel2]
    qq1=(1+skew)*qnorm(x1/(1+skew)) 
    qq2= (1-skew)*qnorm((x2-skew)/(1-skew))
  return(c(qq1,qq2))
}
skew = as.numeric(pp["skew"])
q_t =qtpGaussian(probabilities,skew)
q_t1 =qtpGaussian(probabilities1,skew)
plot(q_t,q_e,main="Two-Piece Gaussian qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
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
q_t = sort(as.numeric(lapply(probabilities,f.inv)))
q_t1 = sort(as.numeric(lapply(probabilities1,f.inv)))
plot(q_t,q_e,main="Two-Piece Bimodal qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
if(model %in% c("TwoPieceStudentT"))
{
qtpt = function(x,skew,df){
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
q_t =qtpt(probabilities,skew,df)
q_t1 =qtpt(probabilities1,skew,df)
plot(q_t,q_e,main="Two-Piece Student qq-plot",xlab=xlab,ylab=ylab,pch=20)
}
#######################################
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
q_t =qtptukey(probabilities,skew,tail)
q_t1 =qtptukey(probabilities1,skew,tail)
plot(q_t,q_e,main="Two-Piece Tukey-h qq-plot",xlab=xlab,ylab=ylab,pch=20)
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
   plot(q_t1,q_e1,main="First Skew Gaussian qq-plot",xlab=xlab,ylab=ylab,pch=20)
   aa=lm(q_e1~1+q_t1) 
   abline(as.numeric(aa$coefficients[1]),as.numeric(aa$coefficients[2]))

   omega2=sqrt((pp["skew_2"]^2 + pp["sill_2"])/pp["sill_2"])
   alpha2=pp["skew_2"]/pp["sill_2"]^0.5
   q_t2=sn::qsn(probabilities,xi=0,omega= as.numeric(omega2),alpha= as.numeric(alpha2))
   q_e2=quantile(dd2,probabilities)
   plot(q_t2,q_e2,main="Second Skew Gaussian qq-plot",xlab=xlab,ylab=ylab,pch=20)
   aa=lm(q_e2~1+q_t2) 
   abline(as.numeric(aa$coefficients[1]),as.numeric(aa$coefficients[2]))
}

##########################################################

par(mfrow=c(1,1))
  }
}
