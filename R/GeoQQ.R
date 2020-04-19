


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

#######################################
if(model %in% c("Gaussian")) 
{ pp=pnorm(dd,mean = pp['mean'], sd = sqrt(pp['sill']))
  rr=as.numeric(quantile(dd,pp))  
  plot(dd,rr,main="Gaussian qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("SkewGaussian"))
{
   omega=as.numeric(sqrt((pp["skew"]^2 + pp["sill"])/pp["sill"]))
   alpha=as.numeric(pp["skew"]/pp["sill"]^0.5)

   pp=sn::psn(dd,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))  
   rr=as.numeric(quantile(dd,pp))
   plot(dd,rr,main="Skew Gaussian qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model%in%c("StudentT","Gaussian_misp_StudentT")) 
{
      limma::qqt(dd,df=as.numeric(round(1/pp["df"])),main="t qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("Weibull"))
{
   shape=pp["shape"]
   pp=pweibull(dd,shape=shape,scale=1/(gamma(1+1/shape )))
   rr=as.numeric(quantile(dd,pp))
   plot(dd,rr, main ="Weibull qq-plot ",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("Gamma"))
{
   shape=pp["shape"]
   pp=pgamma(dd,shape=shape/2,rate=shape/2)
   rr=as.numeric(quantile(dd,pp))
   plot(dd,rr, main ="Gamma qq-plot ",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("LogGaussian"))
{
   SS=pp["sill"]; mm=pp["mean"]
   pp = plnorm(dd, mm/exp(SS/2), sqrt(SS/exp(SS)))
   rr=as.numeric(quantile(dd,pp))
   plot(dd,rr,xlab = "",ylab = "",main = "LogGaussian qq-plot")
}
#######################################
if(model %in% c("LogLogistic"))
{
shape=pp["shape"]
cc=gamma(1+1/shape)*gamma(1-1/shape)
pp = actuar::pllogis(dd,shape = shape,scale=1/cc)
rr=as.numeric(quantile(dd,pp))
plot(dd,rr,xlab = "",ylab = "",main = "LogLogistic qq-plot") 
}
#######################################
if(model %in% c("SinhAsinh"))
{
tail = as.numeric(pp["tail"])
skew = as.numeric(pp["skew"])
sas.quantiles=sinh(1/tail * asinh(qnorm(dd))+skew/tail)
plot(sas.quantiles,sort(c(dd)),main="Sas qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("Tukeyh"))
{
pTukeyh = function(x,tail){
    t = sqrt(VGAM::lambertW(tail*x*x)/tail)*sign(x)
    pnorm(t)
    }
tail = as.numeric(pp["tail"])
pp= pTukeyh(dd,tail = tail)
rr=as.numeric(quantile(dd,pp))
plot(dd,rr,main="Tukey-h qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("TwoPieceGaussian"))
{

ptpGaussian = function(x,skew){
    (1+skew)*pnorm(x/(1+skew))*I(x<0) + (skew + (1-skew)*pnorm(x/(1-skew)))*I(x>=0)
}
skew = as.numeric(pp["skew"])
pp= ptpGaussian(dd,skew = skew)
rr=as.numeric(quantile(dd,pp))
plot(dd,rr,,main="Two-Piece Gaussian qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("TwoPieceBimodal"))
{
ptpbimodal = function(x,skew,delta,df){  
  alpha=2*(delta+1)/df
  nn=2^(1-alpha/2)
  ll=length(x)
  qq=double(ll)
  for( i in 1:ll){
  if(x[i]<=0) qq[i]=(0.5*(1+skew)*as.numeric(zipfR::Igamma(df/2,nn*(-x[i])^(alpha)/(2*((1+skew)^(alpha))),lower=FALSE))/(gamma(df/2)))
  if(x[i]>0)  qq[i]=(0.5*(skew+1)+(0.5*(1-skew)*as.numeric(zipfR::Igamma(df/2,nn*(x[i])^(alpha)/(2*((1-skew)^(alpha))),lower=TRUE))/(gamma(df/2))))
  }
return(as.numeric(qq))
}
skew = as.numeric(pp["skew"])
df   = as.numeric(pp["df"])
delta= as.numeric(pp["shape"])
pp=ptpbimodal(dd,skew = skew,delta=delta,df=df)
rr=as.numeric(quantile(dd,pp))
plot(dd,rr,main="Two-Piece Bimodal qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("TwoPieceStudentT"))
{
ptpStudent = function(x,skew,df){
    (1+skew)*pt(x/(1+skew),df = df)*I(x<0) + (skew + (1-skew)*pt(x/(1-skew),df = df))*I(x>=0)
}
skew = as.numeric(pp["skew"])
df   = 1/as.numeric(pp["df"])
pp= ptpStudent(dd,skew = skew,df = df)
rr=as.numeric(quantile(dd,pp))
plot(dd,rr,main="Two-Piece Student qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("TwoPieceTukeyh"))
{
ptpTukeyh = function(x,skew,tail){
    t = sqrt(VGAM::lambertW(tail*x*x)/tail)*sign(x)
    (1+skew)*pnorm(t/(1+skew))*I(t<0) + (skew + (1-skew)*pnorm(t/(1-skew)))*I(t>=0)
}
skew= as.numeric(pp["skew"])
tail= as.numeric(pp["tail"])
pp= ptpTukeyh(dd,skew = skew,tail = tail)
rr=as.numeric(quantile(dd,pp))
plot(dd,rr,main="Two-Piece Tukey-h qq-plot",xlab=xlab,ylab=ylab)
}
########################################
abline(0,1)
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
   pp1=sn::psn(dd1,xi=0,omega= as.numeric(omega1),alpha= as.numeric(alpha1))
   rr1=as.numeric(quantile(dd1,pp1))
   plot(dd1,rr1,main="First Skew Gaussian qq-plot",xlab=xlab,ylab=ylab)
   abline(0,1)

   omega2=sqrt((pp["skew_2"]^2 + pp["sill_2"])/pp["sill_2"])
   alpha2=pp["skew_2"]/pp["sill_2"]^0.5
   pp2=sn::psn(dd2,xi=0,omega= as.numeric(omega2),alpha= as.numeric(alpha2))
   rr2=as.numeric(quantile(dd2,pp2))
   plot(dd2,rr2,main="Second Skew Gaussian qq-plot",xlab=xlab,ylab=ylab)
   abline(0,1)
}

##########################################################

par(mfrow=c(1,1))
  }
}
