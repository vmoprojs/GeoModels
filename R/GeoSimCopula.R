####################################################
### File name: GeoSimCopula.r
####################################################


# Simulate spatial and spatio-temporal random felds:
GeoSimCopula <- function(coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl",GPU=NULL, grid=FALSE,
     local=c(1,1),method="cholesky",model='Gaussian', n=1, param,anisopars=NULL, radius=6371, sparse=FALSE,copula="Gaussian",seed=NULL,X=NULL)
{

if(!is.null(seed))  set.seed(seed)
if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")

if(is.null(CkModel(model))) stop("The name of the  model  is not correct\n")
    corrmodel=gsub("[[:blank:]]", "",corrmodel)
    model=gsub("[[:blank:]]", "",model)
    distance=gsub("[[:blank:]]", "",distance)
    method=gsub("[[:blank:]]", "",method)

if((copula!="Clayton")&&(copula!="Gaussian")) stop("the type of copula is wrong")



#### corr parameters
paramcorr=param[CorrParam(corrmodel)]
####Gaussian copula #############################################
if(copula=="Gaussian")
{
param1=c(list(mean=0,sill=1,nugget=param$nugget),paramcorr)

sim=GeoSim(coordx=coordx, coordy=coordy,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance,GPU=GPU, grid=grid,
     local=local,method=method,model='Gaussian', n=1, param=param1,anisopars=anisopars, radius=radius, sparse=sparse)
unif=pnorm(sim$data,mean=0,sd=1);
}
####beta copula #############################################
if(copula=="Clayton")
{
pp=round(as.numeric(param['nu']))
param1=c(list(shape1=pp,shape2=2,sill=1,mean=0,min=0,max=1,nugget=param$nugget),paramcorr)
sim=GeoSim(coordx=coordx, coordy=coordy,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance,GPU=GPU, grid=grid,
     local=local,method=method,model='Beta', n=1, param=param1,anisopars=anisopars, radius=radius, sparse=sparse)
unif=(sim$data)^(pp/2)
}
####################################################################
####################################################################
if(sim$spacetime||sim$bivariate) DD=dim(simcop)
  if(!sim$bivariate){
           if(is.null(dim(X))) {X=as.matrix(rep(1,sim$numcoord*sim$numtime))}
           sel=substr(names(param),1,4)=="mean";
           num_betas=sum(sel) 
           if(num_betas==1)  mm<-as.numeric(param$mean)
           if(num_betas>1)   mm<- X%*%as.numeric((param[sel]))
           param$mean=0;if(num_betas>1) {for(i in 1:(num_betas-1)) param[[paste("mean",i,sep="")]]=0}
    }       
##############################
##############################
#######     models     #######
##############################
if(!sim$bivariate) {}



####################################
############ discrete  RF ##########
####################################
if(model=="Binomial") {
 simcop=qbinom(unif, size=n, prob=pnorm(mm))
}
if(model=="BinomialNeg") {
 simcop=qnbinom(unif, size=n, prob=pnorm(mm))
}
if(model=="BinomialNegZINB") {
 simcop=qzinegbin(unif, size=n, #prob = pnorm(mm), 
        munb=   n*(1-pnorm(mm))/pnorm(mm),#,n/pnorm(mm)-n,
        pstr0 = as.numeric(param$pmu))
}
if(model=="Poisson") {
 simcop=qpois(unif, lambda=exp(mm))
}
if(model=="PoissonZIP") {
 simcop=qzipois(unif, lambda=exp(mm), pstr0 = as.numeric(param$pmu) ) # require package VGAM
}
####################################
############ positive real  RF #####
####################################
if(model=="Gamma") 
         {
p2=as.numeric(param$shape)
simcop=exp(mm)*qgamma(unif,shape=p2/2,scale=p2/2)
         }
############
if(model=="Weibull") 
         {
p2=as.numeric(param$shape)
simcop= exp(mm)*qweibull(unif,shape=p2,scale=1/(gamma(1+1/p2 )))
         }
########################
if(model=="LogLogistic") 
         {
p2=as.numeric(param$shape)
cc=gamma(1+1/p2)*gamma(1-1/p2)
simcop=actuar::qllogis(unif,shape = p2,scale=exp(mm)/cc)            
         }
##############################

if(model=="LogGaussian")
{ 
    vv=as.numeric(param$sill)
   simcop = qlnorm(unif, exp(mm)-vv/2, sqrt(vv))
}
##############################
##############################
############ Real RF #########
##############################

if(model=="Gaussian") {
         simcop=qnorm(unif,mean=mm,sd=sqrt(as.numeric(param$sill)))
         }
if(model=="Logistic") 
         {
         simcop=qlogis(unif,location=mm,scale=sqrt(as.numeric(param$sill)))
         }
if(model=="StudentT") {
vv=as.numeric(param$sill)
dd=as.numeric(param$df)
simcop=mm+sqrt(vv)*qt(unif,df=round(1/dd))
         }
############
if(model=="SkewGaussian") 
         {
vv=as.numeric(param$sill)  
sk=as.numeric(param$skew)  
omega=as.numeric(sqrt((sk^2 + vv)/vv)) 
alpha=as.numeric(sk/(vv^0.5))         
simcop=mm+sqrt(vv)*sn::qsn(unif,xi=0,omega= as.numeric(omega),alpha= as.numeric(alpha))
         }
#######################################  
if(model=="SkewStudentT")
{
  vv=as.numeric(param$sill) 
  sk=as.numeric(param$skew)
  dd=as.numeric(round(1/pp["df"]))
 simcop=mm+sqrt(vv)*sn::qst(unif, xi=0, omega=1, alpha=sk, nu=dd)
}
#######################################
if(model=="SinhAsinh")
{
vv=as.numeric(param$sill) 
tail = as.numeric(param$tail) 
skew = as.numeric(param$skew)
simcop=mm+sqrt(vv)*sinh(1/tail * asinh(qnorm(unif))+skew/tail)
}
#######################################   OK
if(model=="Tukeyh")
{
vv=as.numeric(param$sill) 
tail = as.numeric(param$tail) 
uu=qnorm(unif)
simcop=mm+sqrt(vv)*uu*exp(0.5*tail*uu^2);
}
#######################################  OK
if(model=="Tukeyh2")
{
qtpTukeyh22= function(x,tail1,tail2){
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
tail1 =as.numeric(param$tail1);tail2 = as.numeric(param$tail2)
vv=as.numeric(param$sill) 
simcop =mm+sqrt(vv)*qtpTukeyh22(unif,tail1,tail2)
}
#######################################  OK
if(model=="Tukeygh")
{
tail = as.numeric(param$tail);skew =  as.numeric(param$skew) 
vv=as.numeric(param$sill) 
uu=qnorm(unif)
simcop=mm+sqrt(vv)*(exp(skew*uu)-1)*exp(0.5*tail*uu^2)/skew
}
#######################################   OK
if(model=="TwoPieceGaussian")
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
skew = as.numeric(param$skew) 
vv=as.numeric(param$sill) 
simcop =mm+sqrt(vv)*qtpGaussian(unif,skew)
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
vv=as.numeric(param$sill) 
skew = as.numeric(param$skew)
df   = as.numeric(param$df)
delta= as.numeric(param$shape)
f = function(x) ptpbimodal(x,skew = skew,delta=delta,df=df)
f.inv = GoFKernel::inverse(f,lower = -4,upper = 4)
simcop = mm+sqrt(vv)*sort(as.numeric(lapply(unif,f.inv)))
}
#######################################  
if(model=="TwoPieceStudentT")
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
vv=as.numeric(param$sill) 
skew = as.numeric(param$skew)
df   = 1/as.numeric(param$df)
simcop=mm+sqrt(vv)*qtpt(unif,skew,df)
}
#######################################   OK
if(model=="TwoPieceTukeyh")
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
vv=as.numeric(param$sill) 
skew = as.numeric(param$skew)
tail= as.numeric(param$tail)
simcop =mm+sqrt(vv)*qtptukey(unif,skew,tail)
}




####################################
############ bounded  real RF    ###
####################################
if(model=="Kumaraswamy") 
{
p1=param$shape1;p2=param$shape2
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
simcop=pmin + (pmax-pmin)*((1-unif^(1/p1))^(1/p2))
}
############
if(model=="Kumaraswamy2") 
{ # parametrization using beta median  regression
mm=1/(1+exp(-mm))
p2=as.numeric(param$shape)
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
aa=log(1-mm^(p2))/log(0.5)
simcop=pmin + (pmax-pmin)*((1-unif^(aa))^(1/p2))
}
############
if(model=="Beta") 
{ # parametrization using beta regression
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
simcop=pmin + (pmax-pmin)*qbeta(unif,shape1=param$shape1,shape2=param$shape2)
}
############
if(model=="Beta2") 
{ # parametrization using beta mean  regression

mm=1/(1+exp(-mm))
p2=as.numeric(param$shape)
pmin=as.numeric(param$min);pmax=as.numeric(param$max);
simcop=pmin + (pmax-pmin)*qbeta(unif,shape1=mm*p2,shape2=(1-mm)*p2)
}





############
if(sim$spacetime||sim$bivariate) {dim(simcop)=DD}
else {if (!grid) simcop=c(simcop)
     }
##############################
##############################
##############################

    GeoSim_Copula <- list(bivariate = sim$bivariate,
    coordx = sim$coordx,
    coordy = sim$coordy,
    coordt = sim$coordt,
    coordx_dyn =sim$coordx_dyn,
    corrmodel = corrmodel,
    data = simcop,
    distance = sim$distance,
    grid = sim$grid,
    model = sim$model,
    method=method,
    n=sim$n,
    numcoord = sim$numcoord,
    numtime = sim$numtime,
    param = param,
    radius = radius,
    randseed=.Random.seed,
    spacetime = sim$spacetime,
    sparse=sim$sparse,
    copula=copula,
    X=X)
#}
##############################################
  structure(c(GeoSim_Copula, call = call), class = c("GeoSim_Copula"))
}
