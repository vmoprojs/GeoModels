####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate.
### Email: moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: GeoCorrFct.r
### Description:
### This file contains a set of procedures
### for computing a correlation model
### given a set of spatial(temporal) distances
### Last change: 28/05/2021.
####################################################


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

GeoCorrFct<- function(x,t=NULL,corrmodel, model="Gaussian",distance="Eucl",  
                                  param, radius=6371,n=1)

{
  

CorrelationFct <- function(bivariate,corrmodel, lags, lagt, numlags, numlagt, mu,model, nuisance,param,N)
    {
       if(!bivariate) { 
                             p=.C('VectCorrelation', corr=double(numlags*numlagt), as.integer(corrmodel), as.double(lags),
                             as.integer(numlags), as.integer(numlagt), as.double(mu),as.integer(model),as.double(nuisance),as.double(param),
                             as.double(lagt),as.integer(N), PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
                             cc=p$corr
                    }
        else    {
                             p=.C('VectCorrelation_biv', corr=double(numlags*4),vario=double(numlags*4), as.integer(corrmodel), as.double(lags),
                             as.integer(numlags), as.integer(numlagt),  as.double(mu),as.integer(model),as.double(nuisance), as.double(param),
                             as.double(lagt), as.integer(N),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
                             cc=c(p$corr,p$vario)   

                    }
        return(cc)
    }
  #############################################################################################
  #################### end internal function ##################################################
  #############################################################################################
    # Check the user input
    if(is.null(CkCorrModel(corrmodel)))   stop("The correlation model is not valid\n")
    if(is.null(CkModel(model)))   stop("The  model is not valid\n")
    if(!is.numeric(x)) stop("Distances must be numeric\n")
    if(sum(x<0)>=1) stop("Distances must be positive\n")
    spacetime<-CheckST(CkCorrModel(corrmodel))
    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    if(is.null(CkCorrModel (corrmodel))) stop("The name of the coorelation model  is not correct\n")
  
mu=0;nuisance=0
mm=0
num_beta=c(1,1)

nx=length(x)
if(spacetime) nt=length(t)
else nt=1

num_betas=c(1,1)
if(sum((names(param)=='mean'))==0) param$mean=0 # adding mean if missing
print(param)
  ## selecting nuisance mean annd corr parameters
      if(!bivariate){
        parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
        nuisance <- c(param)[NuisParam(model,FALSE,num_betas)]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=nuisance[sel]
        nuisance=nuisance[!sel]
        }
      if(bivariate){
        parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
        nuisance <- c(param)[NuisParam(model,FALSE,num_betas)]
    }
correlation <- CorrelationFct(bivariate,CkCorrModel(corrmodel), x, t, nx, nt,mu,
                                     CkModel(model), nuisance,parcorr,n)

#### non-Gaussian case
return(correlation)
}

