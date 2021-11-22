####################################################
### File name: GeoCorrFct_Cop.r
####################################################

GeoCorrFct_Cop<- function(x,t=NULL,corrmodel, model="Gaussian",copula="Gaussian",distance="Eucl",  
                                  param, radius=6371,n=1,covariance=FALSE)

{
############################################################################
############################################################################
## C functions
biv_unif_CopulaClayton<- function(a,b,c,d)
{
  sol = .C("biv_unif_CopulaClayton_call", as.double(a), as.double(b),
           as.double(c),as.double(d), ress = as.double(0),
           PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
  return(exp(sol$ress))
}
biv_unif_CopulaGauss<- function(a,b,c)
{
  sol = .C("biv_unif_CopulaGauss_call", as.double(a), as.double(b),
           as.double(c),ress = as.double(0),
           PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
  return(sol$ress)
}
###################
biv_unif_CopulaClayton<-Vectorize(biv_unif_CopulaClayton,vectorize.args=c("a","b"))
biv_unif_CopulaGauss<-Vectorize(biv_unif_CopulaGauss,vectorize.args=c("a","b"))
###################
arg_clayton<-function(x,y,rho,nu) return(q1(x)*q2(y)*biv_unif_CopulaClayton(x,y,rho,nu))
arg_gaussian<-function(x,y,rho)   return(q1(x)*q2(y)*biv_unif_CopulaGauss(x,y,rho))

###########################################    
corr_copula<-function(rho,copula,q1,q2,e1,e2,v1,v2,nu){
  if (copula=="Clayton"){
    res<-vector()
    for (i in 1:length(rho)){
          if (rho[i]>0.9999 & rho[i]<=1) res[i]=1
          else                           res[i]=(pracma::integral2(arg_clayton,0,1,0,1,rho=rho[i],nu=nu)$Q-e1*e2)/sqrt(v1*v2)
      
    }
    return(res)
  }
  if (copula=="Gaussian"){
    res<-vector()
    for (i in 1:length(rho)){
      if (rho[i]>0.9999 & rho[i]<=1) res[i]=1
      else                           res[i]=(pracma::integral2(arg_gaussian,0,1,0,1,rho=rho[i])$Q-e1*e2)/sqrt(v1*v2)  
    }
    return(res)
  }
}
###########################################
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
   
mu=0;nuisance=0
mm=0
num_beta=c(1,1)

nx=length(x)
if(spacetime) nt=length(t)
else nt=1

num_betas=c(1,1)
if(sum((names(param)=='mean'))==0) param$mean=0 # adding mean if missing


mu=as.numeric(param$mean)
  ## selecting nuisance mean annd corr parameters
      if(!bivariate){
       
        parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
        nuisance <- c(param)[NuisParam(model,FALSE,num_betas)]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        nuisance=nuisance[!sel]
        }
      if(bivariate){
        parcorr <- c(param)[CorrelationPar(CkCorrModel(corrmodel))]
        nuisance <- c(param)[NuisParam(model,FALSE,num_betas)]
    }

correlation <- CorrelationFct(bivariate,CkCorrModel(corrmodel), x, t, nx, nt,mu,
                                     CkModel(model), nuisance,parcorr,n)
cc=correlation*(1-as.numeric(nuisance['nugget'] )  )
######################################
######################################
######################################
######################################
if(model=="Beta2")        { if(bivariate) {} 
                        else {
                        delta=as.numeric(nuisance["shape"])
                       # min=as.numeric(param['min']);max=as.numeric(param['max'])
                       # x=(x-min)/(max-min)
                       # y=(y-min)/(max-min)
                        q1<-function(x) qbeta(x,mu1*delta,(1-mu1)*delta)
                        q2<-function(x) qbeta(x,mu2*delta,(1-mu2)*delta)
                        mu1=1/(1+exp(-mm));mu2=1/(1+exp(-mm))
                        a1=mu1*delta;a2=mu2*delta
                        b1=(1-mu1)*delta;b2=(1-mu2)*delta
                        e1=a1/(a1+b1); e2=a2/(a2+b2)
                        v1=a1*b1/((a1+b1)^2*(a1+b1+1));v2=a2*b2/((a2+b2)^2*(a2+b2+1))
                        vs=sqrt(v1*v2)
                     }
}        
            
################################################################
################################################################                     
if(copula=="Clayton")    cc=corr_copula(cc,"Clayton" ,q1,q2,e1,e2,v1,v2,as.numeric(param['nu']))
if(copula=="Gaussian")   cc=corr_copula(cc,"Gaussian",q1,q2,e1,e2,v1,v2,0)

if(covariance) cc=cc*vs
return(cc)
}

