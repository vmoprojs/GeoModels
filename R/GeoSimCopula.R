####################################################
### File name: GeoSimCopula.r
####################################################


# Simulate spatial and spatio-temporal random felds:
GeoSimCopula <- function(coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl",GPU=NULL, grid=FALSE,
     local=c(1,1),method="cholesky",model='Gaussian', n=1, param, radius=6371, sparse=FALSE,copula="Gaussian",X=NULL)
{

if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")
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
     local=local,method=method,model='Gaussian', n=1, param=param1, radius=radius, sparse=sparse)
unif=pnorm(sim$data,mean=0,sd=1);
}
####beta copula #############################################
if(copula=="Clayton")
{
pp=round(as.numeric(param['nu']))
param1=c(list(shape1=pp,shape2=2,sill=1,mean=0,min=0,max=1,nugget=param$nugget),paramcorr)
sim=GeoSim(coordx=coordx, coordy=coordy,coordt=coordt, coordx_dyn=coordx_dyn,corrmodel=corrmodel, 
    distance=distance,GPU=GPU, grid=grid,
     local=local,method=method,model='Beta', n=1, param=param1, radius=radius, sparse=sparse)
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


if(model=="Gaussian") {
         simcop=qnorm(unif,mean=mm,sd=sqrt(as.numeric(param$sill)))
         }
if(model=="Logistic") 
         {
         simcop=qlogis(unif,location=mm,scale=sqrt(as.numeric(param$sill)))
         }
if(model=="Gamma") 
         {
        # simcop=qgamma(unif,shape=as.numeric(param$shape),scale=1/mm)
         }
#######
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
