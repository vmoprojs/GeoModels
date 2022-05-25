
####################################################
### File name: GeoResiduals.r
####################################################

GeoResiduals<-function(fit)

{
if(!inherits(fit,"GeoFit"))  stop("A GeoFit object is needed as input\n")
######

fit$param=unlist(fit$param)
fit$fixed=unlist(fit$fixed)
model=fit$model        #type of model
num_betas=fit$numbetas  #number of mean parameters

if(!fit$bivariate)
{
## extracting mean parameters
namescorr <- CorrParam(fit$corrmodel) 
namesnuis <- NuisParam(fit$model,fit$bivariate,num_betas)
param <- c(fit$param, fit$fixed)
namesparam<- names(param)
paramcorr <- param[namescorr]
nuisance <- param[namesnuis]
sel=substr(names(nuisance),1,4)=="mean"
beta2=as.numeric(nuisance[sel])
## names of estimated and fixed parameters
nm=names(fit$param)
nf=names(fit$fixed)
copula=fit$copula

#################################
#### computing mean ########
#################################

mu=fit$X%*%beta2  
#################################
#### computing residuals ########
#################################
if(is.list(fit$coordx_dyn)) dd=unlist(fit$data)
else dd=c(t(fit$data))


############################################################
if(!is.null(copula)){    #### copula models
if(copula=="Clayton"||copula=="Gaussian")
         {
if(model=="Beta2")
             {
              mm=c(1/(1+exp(-mu)))
              sh=as.numeric(param['shape']);  
              res1=pbeta(dd,mm*sh,(1-mm)*sh)
             }
if(model=="Kumaraswamy2") {
             mm=c(1/(1+exp(-mu))); 
             sh=as.numeric(param['shape']);
             ga= (log(1-mm^sh)/log(0.5))^{-1}
             res1=(1-dd^ga)^(sh) }
          }
}
else {           #### non-copula models
###  positive multiplicative models
if(model %in% c("Gamma","Weibull","LogLogistic","LogGaussian"))
res1=dd/exp(c(mu))
### additive  models  on the real line
if(model %in% c("Gaussian","SkewGaussian","Logistic", 
               "Tukeyh","Tukeyh2","SinhAsinh","Tukeygh","Gaussian_misp_Tukeygh",
               "StudentT",  "Gaussian_misp_StudentT","Gaussian_misp_SkewStudentT","SkewStudentT",
               "TwoPieceGaussian","TwoPieceTukeyh","TwoPieceGauss","TwoPieceStudentT","TwoPieceBimodal"))
{res1=(dd-c(mu))/sqrt(as.numeric(param['sill']))}

if(model=="Gaussian_misp_Binomial")
{
    aa=pnorm(c(mu))
    hh=fit$n*aa
    res1=(dd-hh)/sqrt(hh*(1-aa))
}
if(model=="Gaussian_misp_BinomialNeg")
{
    aa=pnorm(c(mu))
    hh=fit$n*(1-aa)/aa
    res1=(dd-hh)/sqrt(hh/aa)
}
if(model=="Gaussian_misp_Poisson")
{
    aa=exp(c(mu))
    res1=(dd-aa)/sqrt(aa)
}
#########
}



fit$X=as.matrix(rep(1,length(dd)))

#### updating  object
mm=0
names(mm)="mean"
nuis_update=c(mm,nuisance[!sel])
fit$param=c(nuis_update,paramcorr)
fit$numbetas=1
fit$X=as.matrix(rep(1,length(c(fit$data))))

if(!is.null(copula)){
#############################################
if(copula=="Clayton"||copula=="Gaussian")
{
if(model %in% c("Beta2")) {fit$param['shape']=2;fit$param['mean']=0; fit$param['max']=1; fit$param['min']=0}
if(model %in% c("Kumaraswamy2")) {fit$param['shape']=1;fit$param['mean']=0; fit$param['max']=1; fit$param['min']=0}
}
}
else
{
#####################################################
if(model %in% c("Gaussian","Logistic","Tukeyh","Tukeyh2","Tukeygh","SinhAsinh", "Gaussian_misp_StudentT","Gaussian_misp_Tukeygh",
         "StudentT","TwoPieceGauss","TwoPieceStudentT","TwoPieceGaussian","TwoPieceTukeyh","TwoPieceBimodal"))
{fit$param['sill']=1;fit$param['mean']=0}


if(model %in% c("SkewGaussian")) 
{param['mean']=0;
 fit$param['skew']=as.numeric(param['skew'])/sqrt(as.numeric(param['sill']))
 fit$param['sill']=1
}

if(model %in% c("Gaussian_misp_SkewStudentT","SkewStudentT")) 
{param['mean']=0;
 fit$param['skew']=as.numeric(param['skew'])
 fit$param['df']= fit$param['df']
 fit$param['sill']=1
}
##
}

fit$param=fit$param[nm]
fit$fixed=fit$fixed[nf]

### deleting NA
fit$param=fit$param[!is.na(fit$param)]
fit$fixed=fit$fixed[!is.na(fit$fixed)]


### formatting data
if(fit$spacetime) 
{if(!is.list(fit$coordx_dyn)) 
          data_res=matrix(res1,nrow=nrow(fit$data),ncol=ncol(fit$data),byrow=TRUE)
 else{ 
     ns=unlist(lapply(fit$coordx_dyn,nrow))
     sim_temp=list()
                    for(k in 1:length(fit$coordt))
                       { if(k==1) {indx=1:(sum(ns[1:k]))}
                         if(k>1)    {indx=(sum(ns[1:(k-1)])+1):(sum(ns[1:k]))}
                         sim_temp[[k]]=c(res1)[indx] }
    data_res=sim_temp  
 }
}   
else   {data_res=as.vector(res1)}

}
#####################################################################
#####################################################################
if(fit$bivariate)
{
 X=fit$X
 ns=fit$ns

namescorr <- CorrParam(fit$corrmodel) 
namesnuis <- NuisParam(fit$model,fit$bivariate,num_betas)
param <- c(fit$param, fit$fixed)
namesparam<- names(param)
paramcov <- param[namescorr]
nuisance <- param[namesnuis]
###
if(is.list(fit$coordx_dyn)) dd=unlist(fit$data)
else dd=c(t(fit$data))
dd1=dd[1:ns[1]]
dd2=dd[(ns[1]+1):(ns[2]+ns[1])]
###

###
 sel1=substr(names(nuisance),1,6)=="mean_1"
 beta1=as.numeric(nuisance[sel1])
 sel2=substr(names(nuisance),1,6)=="mean_2"
 beta2=as.numeric(nuisance[sel2])

 X1=as.matrix(X[1:ns[1],]);
 X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 

mu1=X1%*%beta1
mu2=X2%*%beta2

if(model %in% c("Gaussian")){
       res1=(dd1-c(mu1))/sqrt(as.numeric(paramcov['sill_1']))
       res2=(dd2-c(mu2))/sqrt(as.numeric(paramcov['sill_2']))
        }



if(!is.list(fit$coordx_dyn)) 
           {data_res=rbind(res1,res2)}
else {data_res=list(); data_res[[1]]=res1;data_res[[2]]=res2}
}



#####################################################################
#####################################################################

fit$data=data_res
### geofit object
GeoFit <- list(bivariate=fit$bivariate,
                         claic = fit$claic,
                         clbic = fit$clbic,
                         coordx = fit$coordx,
                         coordy = fit$coordy,
                         coordt = fit$coordt,
                         coordx_dyn=fit$coordx_dyn,
                         copula=fit$copula,
                         convergence = fit$convergence,
                         corrmodel = fit$corrmodel,
                         data= fit$data,
                         distance = fit$distance,
                         fixed = as.list(fit$fixed),
                         grid = fit$grid,
                         iterations = fit$counts,
                         likelihood = fit$likelihood,
                         logCompLik = fit$logCompLik,
                         message = fit$message,
                         model = fit$model,
                         n=fit$n,
                         ns=fit$ns,
                         numbetas=fit$numbetas,
                         numcoord=fit$numcoord,
                         numtime=fit$numtime,
                         optimizer=fit$optimizer,
                         lower=fit$lower,
                         upper=fit$upper,
                         param = as.list(fit$param),
                         nozero = fit$setup$nozero,
                         score = fit$score,
                         maxdist =fit$maxdist,
                         maxtime = fit$maxtime,
                         missp=fit$missp,
                         radius = fit$radius,
                         spacetime = fit$spacetime,
                         stderr = fit$stderr,
                         sensmat = fit$sensmat,
                         varcov = fit$varcov,
                         varimat = fit$varimat,
                         vartype = fit$vartype,
                         type = fit$type,
                         weighted=fit$weighted,
                         winconst = fit$winconst,
                         winstp = fit$winstp,
                         winconst_t = fit$winconst_t,
                         winstp_t = fit$winstp_t,
                         X = fit$X)
    structure(c(GeoFit, call = call), class = c("GeoFit"))
#########################
}
