


GeoResiduals<-function(fit)

{
if(class(fit)!="GeoFit") stop("A GeoFit object is needed as input\n")
######
num_betas=fit$numbetas
model=fit$model

## extracting mean parameters
namescorr <- CorrParam(fit$corrmodel) 
namesnuis <- NuisParam(fit$model,fit$bivariate,num_betas)
param <- c(fit$param, fit$fixed)
namesparam<- names(param)
paramcorr <- param[namescorr]
nuisance <- param[namesnuis]
sel=substr(names(nuisance),1,4)=="mean"
beta2=as.numeric(nuisance[sel])   

#################################
#### computing mean ########
#################################
mu=fit$X%*%beta2  
#################################
#### computing residuals ########
#################################

dd=c(t(fit$data))
### multiplicative models
if(model %in% c("Gamma","Weibull","LogLogistic","LogGaussian"))
res1=dd/exp(mu)  
### additive models
if(model %in% c("Gaussian","Logistic","SkewGaussian","TwoPieceGaussian",
         "StudentT","TwoPieceGauss","TwoPieceStudentT"))
res1=dd-mu  

#if(binomial or binomialneg or geom or bernoulli)
#
#............

fit$X=as.matrix(rep(1,length(dd)))

#### updating  object
mm=0
names(mm)="mean"
nuis_update=c(mm,nuisance[!sel])
fit$param=c(nuis_update,paramcorr)
fit$numbetas=1
fit$X=as.matrix(rep(1,length(c(fit$data))))

### formatting data
if(fit$spacetime) data_res=matrix(res1,nrow=nrow(fit$data),ncol=ncol(fit$data),byrow=TRUE)
else   data_res=as.vector(res1)

fit$data=data_res
#if(fit$bivariate.....)


### geofit object
GeoFit <- list(bivariate=fit$bivariate,
                         claic = fit$claic,
                         clbic = fit$clbic,
                         coordx = fit$coordx,
                         coordy = fit$coordy,
                         coordt = fit$coordt,
                         coordx_dyn=fit$coordx_dyn,
                         convergence = fit$convergence,
                         corrmodel = fit$corrmodel,
                         data= fit$data,
                         distance = fit$distance,
                         fixed = fit$fixed,
                         grid = fit$grid,
                         iterations = fit$counts,
                         likelihood = fit$likelihood,
                         logCompLik = fit$logCompLik,
                         message = fit$message,
                         model = fit$model,
                         n=fit$n,
                         numbetas=fit$numbetas,
                         numcoord=fit$numcoord,
                         numtime=fit$numtime,
                         param = fit$param,
                         nozero = fit$setup$nozero,
                         score = fit$score,
                         maxdist =fit$maxdist,
                         maxtime = fit$maxtime,
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
