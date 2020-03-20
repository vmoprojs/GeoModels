


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

#######################################
if(model %in% c("Gaussian")) qqnorm(dd,main="Gaussian qq-plot")
#######################################
if(model %in% c("SkewGaussian"))
{
   omega=sqrt((pp["skew"]^2 + pp["sill"])/pp["sill"])
   alpha=pp["skew"]/pp["sill"]^0.5
   skgauss.quantiles=sn::qsn(probabilities,xi=0,
                       omega= as.numeric(omega),alpha= as.numeric(alpha),
                       solver= "RFB")
   plot(sort(skgauss.quantiles),sort(c(dd)),main="Skew Gaussian qq-plot",xlab=xlab,ylab=ylab)
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
   weibull.quantiles=qweibull(probabilities,shape=shape,scale=1/(gamma(1+1/shape )))
   plot(sort(weibull.quantiles),sort(c(dd)), main ="Weibull qq-plot ",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("Gamma"))
{
   shape=pp["shape"]
   gamma.quantiles=qgamma(probabilities,shape=shape/2,scale=shape/2)
   plot(sort(gamma.quantiles),sort(c(dd)), main ="Gamma qq-plot ",xlab=xlab,ylab=ylab)
}

if(model %in% c("LogGaussian"))
{
   SS=pp["sill"]; mm=pp["mean"]
   lognormal.quantiles = qlnorm(probabilities, mm/exp(SS/2), sqrt(SS/exp(SS)))
   plot(sort(lognormal.quantiles),sort(c(dd)),xlab = "",ylab = "",main = "LogGaussian qq-plot")
}
if(model %in% c("LogLogistic"))
{
shape=pp["shape"]
cc=gamma(1+1/shape)*gamma(1-1/shape)
loglogistic.quantiles = actuar::qllogis(probabilities,shape = shape,scale=1/cc)
plot(sort(loglogistic.quantiles),sort(dd),xlab = "",ylab = "",main = "LogLogistic qq-plot") 
}



#######################################
#######################################
if(model %in% c("Tukeyh"))
{

tail = as.numeric(pp["tail"])
pTukeyh = function(x,tail){
    t = sqrt(VGAM::lambertW(tail*x*x)/tail)*sign(x)
    pnorm(t)
    }

f = function(x) pTukeyh(x,tail = tail)
f.inv = GoFKernel::inverse(f,lower = -Inf,upper = Inf)
tukeyh.quantiles = sort(as.numeric(lapply(probabilities,f.inv)))
plot(tukeyh.quantiles,sort(c(dd)),main="Tukey-h qq-plot",xlab=xlab,ylab=ylab)
}

#######################################
if(model %in% c("TwoPieceGaussian"))
{
skew = as.numeric(pp["skew"])

ptpGaussian = function(x,skew){
    (1+skew)*pnorm(x/(1+skew))*I(x<0) + (skew + (1-skew)*pnorm(x/(1-skew)))*I(x>=0)
}

f = function(x) ptpGaussian(x,skew = skew)
f.inv = GoFKernel::inverse(f,lower = -Inf,upper = Inf)

twopieceGaussian.quantiles = sort(as.numeric(lapply(probabilities,f.inv)))

plot(twopieceGaussian.quantiles,sort(c(dd)),main="Two-Piece Gaussian qq-plot",xlab=xlab,ylab=ylab)

}
#######################################
if(model %in% c("TwoPieceStudentT"))
{
qtwopieceStudent = function(probability,skew,df){
    qt(probability/(1+skew),df = df)*(1+skew)*I(probability < (skew+1)*0.5) + (qt((probability-skew)/(1-skew),df = df)*(1-skew))*I(probability >= (skew+1)*0.5)
}

skew = as.numeric(pp["skew"])
df   = 1/as.numeric(pp["df"])

ptpStudent = function(x,skew,df){
    (1+skew)*pt(x/(1+skew),df = df)*I(x<0) + (skew + (1-skew)*pt(x/(1-skew),df = df))*I(x>=0)
}

f = function(x) ptpStudent(x,skew = skew,df = df)
f.inv = GoFKernel::inverse(f,lower = -Inf,upper = Inf)

twopieceStudent.quantiles = sort(as.numeric(lapply(probabilities,f.inv)))
plot(twopieceStudent.quantiles,sort(c(dd)),main="Two-Piece Student qq-plot",xlab=xlab,ylab=ylab)
}
#######################################
if(model %in% c("TwoPieceTukeyh"))
{
skew             = as.numeric(pp["skew"])
tail             = as.numeric(pp["tail"])


ptpTukeyh = function(x,skew,tail){
    t = sqrt(VGAM::lambertW(tail*x*x)/tail)*sign(x)
    (1+skew)*pnorm(t/(1+skew))*I(t<0) + (skew + (1-skew)*pnorm(t/(1-skew)))*I(t>=0)
}

f = function(x) ptpTukeyh(x,skew = skew,tail = tail)
f.inv = GoFKernel::inverse(f,lower = -Inf,upper = Inf)

twopieceTukeyh.quantiles = sort(as.numeric(lapply(probabilities,f.inv)))
plot(twopieceTukeyh.quantiles,sort(c(dd)),main="Two-Piece Tukey-h qq-plot",xlab=xlab,ylab=ylab)
}



abline(0,1)
 }





##########################################################
##########################################################

if(fit$bivariate){
par(mfrow=c(1,2))

if(is.list(fit$coordx_dyn)){ dd1=fit$data[[1]];dd2=fit$data[[2]]}
else  {dd1=fit$data[1,];dd2=fit$data[2,];}
N1= length(dd1);N2= length(dd2)
probabilities1= (1:N1)/(N1+1); probabilities2= (1:N2)/(N2+1); 


##########################################################

if(model %in% c("Gaussian")) { qqnorm(dd1,main="First Gaussian qq-plot");abline(0,1)
                               qqnorm(dd2,main="Second Gaussian qq-plot");abline(0,1)}

##########################################################

if(model %in% c("SkewGaussian"))
{
   omega1=sqrt((pp["skew_1"]^2 + pp["sill_1"])/pp["sill_1"])
   alpha1=pp["skew_1"]/pp["sill_1"]^0.5
   skgauss.quantiles1=sn::qsn(probabilities1,xi=0,
                       omega= as.numeric(omega1),alpha= as.numeric(alpha1),
                       solver= "RFB")
   plot(sort(skgauss.quantiles1),sort(c(dd1)),main="First Skew Gaussian qq-plot",xlab=xlab,ylab=ylab)
   abline(0,1)

   omega2=sqrt((pp["skew_2"]^2 + pp["sill_2"])/pp["sill_2"])
   alpha2=pp["skew_2"]/pp["sill_2"]^0.5
   skgauss.quantiles2=sn::qsn(probabilities2,xi=0,
                       omega= as.numeric(omega2),alpha= as.numeric(alpha2),
                       solver= "RFB")
   plot(sort(skgauss.quantiles2),sort(c(dd2)),main="Second Skew Gaussian qq-plot",xlab=xlab,ylab=ylab)
   abline(0,1)
}

##########################################################



par(mfrow=c(1,1))
  }
}
