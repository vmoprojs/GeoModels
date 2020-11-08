####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate.
### Email: moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: GeoCovmatrix.r
### Description:
### This file contains a set of procedures
### for computing a covariance (tapered) matrix for a given
### space(time) covariance model.
### Last change: 28/05/2020.
####################################################

### decomposition of a square  matrix
MatDecomp<-function(mtx,method)    {
        if(method=="cholesky")  {
            mat.decomp <- try(chol(mtx), silent=TRUE)
            #mat.decomp <- try(Rfast::cholesky(mtx,parallel=TRUE))
            if (inherits(mat.decomp , "try-error")) return (FALSE)
        }
        if(method=="svd")      {
            mat.decomp <- svd(mtx)
            cov.logdeth <- try(sum(log(sqrt(mat.decomp$d))), silent=TRUE)
            if (inherits(cov.logdeth, "try-error"))  return (FALSE)
        }
        return(mat.decomp)
    }
### square root of a square matrix
MatSqrt<-function(mat.decomp,method)    {  
        if(method=="cholesky")  varcov.sqrt <- mat.decomp
        if(method=="svd")       varcov.sqrt <- t(mat.decomp$v %*% sqrt(diag(mat.decomp$d)))  #sqrt(diag(mat.decomp$d))%*%t(mat.decomp$u)
        return(varcov.sqrt)
    }  
### inverse a square matrix given a decomposition
MatInv<-function(mat.decomp,method)    {

        if(method=="cholesky")  varcov.inv <- chol2inv(mat.decomp)
        if(method=="svd")       
        { 
              tol = sqrt(.Machine$double.eps)
              e <- mat.decomp$d;e[e > tol] <- 1/e[e > tol] 
              varcov.inv<-mat.decomp$v %*% diag(e,nrow=length(e)) %*% t(mat.decomp$u) 
        } 
        return(varcov.inv)
    } 
### determinant a square matrix given a decomposition    
MatLogDet<-function(mat.decomp,method)    {
        if(method=="cholesky")  det.mat <- 2*sum(log(diag(mat.decomp)))
        if(method=="svd")       det.mat <- sum(log(mat.decomp$d))
        return(det.mat)
    }      


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

GeoCovmatrix <- function(coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
                       maxdist=NULL, maxtime=NULL, model="Gaussian", n=1, param, radius=6371, 
                       sparse=FALSE,taper=NULL, tapsep=NULL, type="Standard",X=NULL)

{
  ########################################################################################################
  ##########  Internal function: computing covariance matrix #############################################
  ########################################################################################################

    Cmatrix <- function(bivariate, coordx, coordy, coordt,corrmodel, dime, n, ns, NS, nuisance, numpairs,
                           numpairstot, model, paramcorr, setup, radius, spacetime, spacetime_dyn,type,X)
    {
   


###################################################################################
############### computing correlation #############################################
###################################################################################

if(model %in% c(1,9,34,12,18,39,27,38,29,21,26,24,10,22,40,28,33,42))
{
  if(type=="Standard") {
      fname <-"CorrelationMat2"
      if(spacetime) fname <- "CorrelationMat_st_dyn2"
        if(bivariate) {
            if(model==1) fname <- "CorrelationMat_biv_dyn2"
            if(model==10)fname <- "CorrelationMat_biv_skew_dyn2" }
corr=double(numpairstot)
  #      cr=.C(fname, corr=corr,  as.double(coordx),as.double(coordy),as.double(coordt),
   #       as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius), 
   #       as.integer(ns),as.integer(NS),
   #       PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",  "integer","double","double","double","integer","integer"),  
    corr=corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,radius,ns,NS,
 INTENT = c("w","r","r","r","r","r","r","r", "r", "r"),
            PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE) 
      }
###############################################################
   if(type=="Tapering")  {
        fname <- "CorrelationMat_tap"
        if(spacetime) fname <- "CorrelationMat_st_tap"
       if(bivariate) fname <- "CorrelationMat_biv_tap"
       corr=double(numpairs)
       #print(paramcorr)
        #cr=.C(fname,  corr=corr, as.double(coordx),as.double(coordy),as.double(coordt),
        #  as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius),as.integer(ns),
        #   as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",  "integer","double","double","double","integer","integer"),  
     corr=corr, coordx,coordy,coordt,corrmodel,nuisance, paramcorr,radius,ns,NS,
  INTENT = c("w","r","r","r","r","r","r","r", "r", "r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE) 
     ## deleting correlation equual  to 1 because there are problems  with hipergeometric function
        sel=(abs(cr$corr-1)<.Machine$double.eps);cr$corr[sel]=0  
      }
    }
###################################################################################
###################################################################################
###################################################################################


if(model==1)   ## gaussian case  
{  
  if(!bivariate)
     {
      corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
      vv=nuisance['sill']
     }
if(bivariate){}
}
###############################################################
if(model==9)  {  ## TukeyGH

if(!bivariate)   {             
          corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
          h=as.numeric(nuisance['tail'])
          g=as.numeric(nuisance['skew'])
          ss=as.numeric(nuisance['sill']) 

  if(!g&&!h) { vv=ss } 

  if(g&&!h){  #  
              varcov <-  diag(dime)
              aa=( -exp(g^2)+exp(g^2*2))*g^(-2)
              corr <- (( -exp(g^2)+exp(g^2*(1+corr)))*g^(-2))/aa
              vv=  aa*ss
              } 
   if(!g&&h){ ##
              varcov <-  diag(dime)
              aa=(1-2*h)^(-1.5) # variance
              corr <- corr/(aa*( (1-h)^2-h^2*corr^2 )^(1.5))
              vv=aa*ss
               } 

  if(h&&g){ # ok
              varcov <-  diag(dime)
                  rho=corr; rho2=rho*rho;
                  tail2=tail*tail
                  tail=h; eta=g
                  eta2=eta*eta; u=1-tail;a=1+rho;
                  A1=exp(a*eta2/(1-tail*a));
                  A2=2*exp(0.5*eta2*  (1-tail*(1-rho2))  / (u*u- tail2*rho2)  );
                  A3=eta2*sqrt(u*u- rho2*tail2)
                  mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                  cova=ss*((A1-A2+1)/A3-mu*mu)
                  vari=ss*((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*
                           sqrt(1-2*tail))-mu*mu);
                  corr <-cova/vari
                  vv=vari 
            }      
}
if(bivariate){}
}
###############################################################
if(model==40)  {  ## TukeyH2  
          
if(!bivariate)    {          
    corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
    hr=as.numeric(nuisance['tail1'])
    hl=as.numeric(nuisance['tail2'])

x1=1-(1-corr^2)*hr
x2=(1-hr)^2-(corr*hr)^2
y1=1-(1-corr^2)*hl
y2=(1-hl)^2-(corr*hl)^2
g=1-hl-hr+(1-corr^2)*hl*hr
h1=sqrt(1-corr^2/(x1^2))+(corr/x1)*asin(corr/x1)
h2=sqrt(1-corr^2/(y1^2))+(corr/y1)*asin(corr/y1)

h3=sqrt(1-corr^2/(x1*y1))+sqrt(corr^2/(x1*y1))*asin(sqrt(corr^2/(x1*y1)))
p1=x1*h1/(2*pi*(x2)^(3/2))+corr/(4*(x2)^(3/2))
p2=y1*h2/(2*pi*(y2)^(3/2))+corr/(4*(y2)^(3/2))
p3=-(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2))+corr/(4*(g)^(3/2))

   mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
   vv1=0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)-(mm)^2
   corr=(p1+p2+2*p3-mm^2)/vv1
  vv=as.numeric(nuisance['sill'])* vv1
       }
if(bivariate){}
} 

###############################################################
if(model==34)  {  ## TukeyH  
          
if(!bivariate)    {             
       corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
    h=nuisance['tail']
   if(!h)     vv=as.numeric(nuisance['sill']) 
   if(h){     aa=(1-2*h)^(-1.5) # variance
              #varcov <-  diag(dime)
              corr <- (-corr/((1+h*(corr-1))*(-1+h+h*corr)*(1+h*(-2+h-h*corr^2))^0.5))/aa
              vv=aa*as.numeric(nuisance['sill'])
              }      
       }
if(bivariate){}
} 
######################################################################        
if(model==12)   ##  student case 
    {
if(!bivariate){
   corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
nu=as.numeric(1/nuisance['df']) 
#corr=exp(log(nu-2)+2*lgamma(0.5*(nu-1))-log(2)-2*lgamma(nu/2)+log(Re(hypergeo::hypergeo(0.5,0.5, nu/2,corr^2)))+log(corr))
if(nu<170) corr=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,corr^2))*cr$corr)/(2*gamma(nu/2)^2)
vv=as.numeric(nuisance['sill'])*(nu)/(nu-2)
}
if(bivariate){}
}

############################################################### 
  if(model==18)   ##  skew student case 
    { 
     if(!bivariate)
{
    corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 

    nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])
    skew2=sk*sk;l=nu/2; f=(nu-1)/2; w=sqrt(1-skew2);y=corr;
    CorSkew=(2*skew2/(pi*w*w+skew2*(pi-2)))*(sqrt(1-y*y)+y*asin(y)-1)+w*w*cr$corr/(w*w+skew2*(1-2/pi)) ;
    mm=sqrt(nu)*gamma(f)*sk/(sqrt(pi)*gamma(l));
    corr=(pi*(nu-2)*gamma(f)^2/(2*(pi*gamma(l)^2-skew2*(nu-2)*gamma(f)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,l,y*y))
                                *((1-2*skew2/pi)*CorSkew+2*skew2/pi)-2*skew2/pi);
    vv=as.numeric(nuisance['sill'])*(nu/(nu-2)-mm*mm);
}
if(bivariate){}
}
############################################################### 
 if(model==39)   ##  two piece bimodal case
    { 
if(!bivariate)
{
 corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
 nu=as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew'])
 delta=as.numeric(nuisance['shape'])
 alpha=2*(delta+1)/nu
 nn=2^(1-alpha/2)

 ll=qnorm((1-sk)/2)
 p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
 corr2=corr^2;sk2=sk^2
 a1=Re(hypergeo::hypergeo(-1/alpha ,-1/alpha,nu/2,corr2))
 a3=3*sk2 + 2*sk + 4*p11 - 1

 MM=(2^(2/alpha)*(gamma(nu/2 + 1/alpha))^2) 
 vari=2^(2/alpha)*(gamma(nu/2 + 2/alpha))*gamma(nu/2)* (1+3*sk2) - sk2*2^(2/alpha+2)*gamma(nu/2+1/alpha)^2 
 corr= MM*(a1*a3-4*sk2)/vari
 vv=as.numeric(nuisance['sill'])*vari/(nn^(2/alpha)*gamma(nu/2)^2)
}

if(bivariate){}
}
###############################################################           
if(model==27)   ##  two piece student case case
    {  
  if(!bivariate)
{
          corr=cr$corr*(1-as.numeric(nuisance['nugget']))
          nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])
          corr2=corr^2;sk2=sk^2
          a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
          a2=corr*asin(corr) + (1-corr2)^(0.5)
          ll=qnorm((1-sk)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          a3=3*sk2 + 2*sk + 4*p11 - 1
          KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
          corr= KK*(a1*a2*a3-4*sk2);
          ttemp=gamma(0.5*(nu-1))/gamma(0.5*nu)
          vv=as.numeric(nuisance['sill'])*((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*ttemp^2)
       #   corr=cr$corr*(1-nuisance['nugget'] )
       #   nu=1/as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew']);sill=as.numeric(nuisance['sill'])
       #   sk2=sk^2
       #  corr2=corr^2;sk2=sk^2
       #   a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
       #   a2=corr *asin(corr) + (1-corr2)^(0.5)
       #   ll=qnorm((1-sk)/2)
       #   p11=pbivnorm::pbivnorm(ll,ll, rho = corr, recycle = TRUE)
       #  a3=3*sk2 + 2*sk + 4*p11 - 1
       #   KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
       #   corr= KK*(a1*a2*a3-4*sk2);
       #   vv=sill* ((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*(gamma(0.5*(nu-1))/gamma(0.5*nu))^2)
}

if(bivariate){}
}
############################################################### 
if(model==38)   ##  two piece tukey h  case
    {
      if(!bivariate)
{
          corr=cr$corr*(1-as.numeric(nuisance['nugget']))
          tail=as.numeric(nuisance['tail']); sk=as.numeric(nuisance['skew'])
          corr2=corr^2;sk2=sk^2;
          gg2=(1-(1-corr2)*tail)^2
          xx=corr2/gg2
          A=(asin(sqrt(xx))*sqrt(xx)+sqrt(1-xx))/(1-xx)^(1.5)
          ll=qnorm((1-sk)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          a3=3*sk2 + 2*sk + 4*p11 - 1
          mm=8*sk2/(pi*(1-tail)^2); 
          ff=(1+3*sk2)/(1-2*tail)^(1.5)
          M=(2*(1-corr2)^(3/2))/(pi*gg2)
          corr=  (M*A*a3-mm)/( ff- mm)
          vv= as.numeric(nuisance['sill'])*((1-2*tail)^(-1.5)* (1+3*(sk2)) - 4*(sk2)*2/(pi*(1-tail)^2))
  }
if(bivariate){}
}
###############################################################  
if(model==29)   ##  two piece gaussian case
    {
      if(!bivariate)
{

          corr1=cr$corr*(1-as.numeric(nuisance['nugget']))
          sk=as.numeric(nuisance['skew']);
          corr2=sqrt(1-corr1^2); sk2=sk^2
          ll=qnorm((1-sk)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          KK=3*sk2+2*sk+ 4*p11 - 1
          corr=(2*((corr2 + corr1*asin(corr1))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )
          vv=as.numeric(nuisance['sill'])*(1+3*sk2-8*sk2/pi)  
  }
if(bivariate){}
}
###############################################################
if(model==10)  {  ##  skew Gaussian case

          if(!bivariate){ 
              #corr=cr$corr*(1-as.numeric(nuisance['nugget']))
              #sk=as.numeric(nuisance['skew'])
              #corr2=corr^2; ; sk2=sk^2; vv=as.numeric(nuisance['sill'])
              #corr=(2*sk2)*(sqrt(1-corr2) + corr*asin(corr)-1)/(pi*vv+sk2*(pi-2)) + (cr$corr*vv)/(vv+sk2*(1-2/pi))
              #vv=vv+as.numeric(nuisance['skew'])^2*(1-2/pi)
              corr=cr$corr*(1-as.numeric(nuisance['nugget']))
              sk=as.numeric(nuisance['skew'])
              corr2=cr$corr^2; ; sk2=sk^2; vv=as.numeric(nuisance['sill'])
              corr=(2*sk2)*(sqrt(1-corr2) + cr$corr*asin(cr$corr)-1)/(pi*vv+sk2*(pi-2)) + (corr*vv)/(vv+sk2*(1-2/pi))
              vv=vv+as.numeric(nuisance['skew'])^2*(1-2/pi)
     }
      if(bivariate){}
}
#################################################################################
################ covariance matrix for models defined on the real line ##########
#################################################################################
if(model %in% c(1,9,34,12,18,39,27,38,29,10,40)){
if(!bivariate)
{
 if(type=="Standard"){
     # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        varcov=varcov*vv
      }
    if(type=="Tapering")  {
          vcov <- vv*corr; 
          varcov <- new("spam",entries=vcov,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
          diag(varcov)=vv

        }
}     
#####
    if(bivariate)      {
       if(type=="Standard"){
          corr <- cr$corr
          varcov<-diag(dime)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
          varcov <- t(varcov)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
        }
        if(type=="Tapering")  {
          varcov <-new("spam",entries=cr$corr,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
        }
      }

}

#################################################################################
############# covariance models defined on a bounded support of the real line ###
#################################################################################
if(model==28)   ##  beta case
    {
      if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget']))
         corr2=corr^2
         shape1=as.numeric(nuisance['shape1']);     
         shape2=as.numeric(nuisance['shape2']);  
         ssup=as.numeric(nuisance['max'])-as.numeric(nuisance['min'])
         cc=0.5*(shape1+shape2)
         vv=ssup^2*shape1*shape2/((cc+1)*(shape1+shape2)^2)        ## variance remember min and max!
idx=which(abs(corr2)>1e-10);corr22=corr2[idx]
######################
  #nu=shape1;alpha=shape2
  nu2=shape1/2;alpha2=shape2/2
  res=0;ss=0;k=0
  while(k<=100){
    p1=2*(lgamma(cc+k)-lgamma(cc)+lgamma(nu2+1+k)-lgamma(nu2+1))
    p2=lgamma(k+1)+(lgamma(nu2+k)-lgamma(nu2))+2*(lgamma(cc+1+k)-lgamma(cc+1))
    b1=p1-p2
    b2=log(hypergeo::genhypergeo(U=c(cc+k,cc+k,alpha2), L=c(cc+k+1,cc+k+1), polynomial=TRUE,maxiter=1000, z=corr22))
    b3=k*log(corr22)
    sum=exp(b1+b2+b3)
    res=res+sum
    if (all(sum<1e-6)){
      break
    } else{
      A=res
    }
    k=k+1
  }
######################
         corr[idx]=A
         corr=shape1*(cc + 1 ) * ((1-corr2)^(cc) *corr -1)/shape2 ## correlation
         corr[-idx]=0
        }
     if(bivariate){}     
}

if(model==33)   ##  Kumaraswamy case
    { 
}

if(model==42)   ##  Kumaraswamy 2 case
    {     
}
#################################################################################
# covariance matrix for models defined on a bounded support of the  real line  ##
#################################################################################
if(model %in% c(28,33,42)){
if(!bivariate)
{
 if(type=="Standard"){
     # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
        varcov=varcov*vv
      }
    if(type=="Tapering")  {
          vcov <- vv*corr; 
          varcov <- new("spam",entries=vcov,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
          diag(varcov)=vv
        }
}     
#####
    if(bivariate)      {
       if(type=="Standard"){
          corr <- cr$corr
          varcov<-diag(dime)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
          varcov <- t(varcov)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
        }
        if(type=="Tapering")  {
          varcov <-new("spam",entries=cr$corr,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
        }
      }

}

#################################################################################
################ models defined on the positive real line #######################
#################################################################################
        
if(model==21)   ##  gamma case
    {
      if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget']))
         corr=corr^2   ### gamma correlation
         sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel]) 
         mu = X%*%mm   ## mean function
         vv=exp(mu)^2 * 2/as.numeric(nuisance['shape'])
      
        }
     if(bivariate){}     
}
###############################################################   
if(model==26)   ##  weibull case
    {

        if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget'])) 
         sh=as.numeric(nuisance['shape'])
         sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel]) 
         mu = X%*%mm
         vv=exp(mu)^2 * (gamma(1+2/sh)/gamma(1+1/sh)^2-1)
         # weibull correlations   
        bcorr=    gamma(1+1/sh)^2/(gamma(1+2/sh)-gamma(1+1/sh)^2)                
        corr=bcorr*((1-corr^2)^(1+2/sh)*Re(hypergeo::hypergeo(1+1/sh,1+1/sh ,1 ,corr^2))-1) 
        }
  if(bivariate) {}
}
############################################################### 
    if(model==24)   ## log-logistic case
    {

       if(!bivariate) {
         corr=cr$corr*(1-as.numeric(nuisance['nugget']))
         sh=as.numeric(nuisance['shape'])
         sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel]) 
         mu = X%*%mm   ## mean  function
         vv=(exp(mu))^2* (2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)  
    corr= ((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
             (Re(hypergeo::hypergeo(-1/sh,-1/sh ,1 ,corr^2))* Re(hypergeo::hypergeo(1/sh,1/sh ,1 ,corr^2)) -1)
        }
   if(bivariate){}     
}


############################################################### 
if(model==22)  {  ## Log Gaussian
       if(!bivariate) {
      corr=cr$corr*(1-as.numeric(nuisance['nugget']))
      vvar=as.numeric(nuisance['sill'])
      corr=(exp(vvar*corr)-1)/(exp(vvar)-1)
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel]) 
            mu = X%*%mm     
            vv=(exp(mu))^2*(exp(vvar)-1) /(exp(vvar*0.5))^2        
             }  
    if(bivariate){}    
  } 
#################################################################################
################ covariance matrix for models defined on the positive real line #
#################################################################################
if(model %in% c(24,26,21,22)){
       if(type=="Standard"){
     # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr
      }
    if(type=="Tapering")  {
          varcov <- new("spam",entries=corr,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
          diag(varcov)=1
        }
    V=vv%*%t(vv)
    varcov=varcov*sqrt(V)
}

###############################################################
################################ end continous models #########
###############################################################

###############################################################
################################ start discrete #models ########
###############################################################
if(model %in% c(2,11,30,16,14,43,45)){ #  binomial (negative)Gaussian type , Poisson (inflated)


if(!bivariate){
        
            fname <-"CorrelationMat_dis2"
            if(spacetime) fname <- "CorrelationMat_st_dyn_dis2"
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel]) 
            mu = X%*%mm
            other_nuis=as.numeric(nuisance[!sel])   
if(type=="Standard")  {
  corr=double(numpairstot)
             

              cr=.C(fname, corr=corr,  as.double(coordx),as.double(coordy),as.double(coordt),
              as.integer(corrmodel), as.double(c(mu)),as.integer(min(n)), as.double(other_nuis), as.double(paramcorr),as.double(radius),
              as.integer(ns), as.integer(NS),as.integer(model),
              PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

  #cr=dotCall64::.C64(fname,SIGNATURE = 
   #    c("double","double","double","double",  "integer","double","integer","double","double","double","integer","integer","integer"),  
   #     corr=corr, coordx,coordy,coordt,corrmodel,c(mu), min(n),other_nuis,paramcorr,radius,ns,NS,model,
  #INTENT = c("w","r","r","r","r","r","r","r", "r", "r","r", "r", "r"),
   #          PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE) 


     corr=cr$corr # ojo que corr en este caso es una covarianza y ya va con el nugget
     varcov <-  diag(dime) 
     varcov[lower.tri(varcov)] <- corr   
     varcov <- t(varcov)
     varcov[lower.tri(varcov)] <- corr     
}
############################
############################
 if(type=="Tapering")  {
        tap <-new("spam",entries=setup$taps,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
        idx=spam::triplet(tap)$indices

        fname <- "CorrelationMat_dis_tap"
        if(spacetime) fname <- "CorrelationMat_st_dis_tap"
         corr=double(numpairs)
      
        cr=.C(fname,  corr=corr, as.double(coordx),as.double(coordy),as.double(coordt),
        as.integer(corrmodel), as.double(other_nuis), as.double(paramcorr),as.double(radius),as.integer(ns),
           as.integer(NS),as.integer(min(n)),as.double(mu[idx[,1]]),as.double(mu[idx[,2]]),as.integer(model),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
    
#      cr=dotCall64::.C64(fname,SIGNATURE = 
#         c("double","double","double","double","integer","double","double","double","integer","integer","integer","double","double","integer"),  
#        corr=corr, coordx,coordy,coordt,corrmodel,other_nuis,paramcorr,radius,ns,NS,min(n),mu[idx[,1]],mu[idx[,2]],model,
#  INTENT = c("w","r","r","r","r","r","r","r", "r", "r","r", "r", "r", "r"),
#             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE) 


        varcov <-new("spam",entries=cr$corr,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))
        }

            ## updating the diagonal with variance  
  if(model %in% c(2,11)) { pg=pnorm(mu); vv=pg*(1-pg)*n; diag(varcov)=vv } 
  if(model %in% c(14))   { pg=pnorm(mu); vv=  (1-pg)/pg^2; diag(varcov)=vv }       
  if(model %in% c(16))   { pg=pnorm(mu); vv=n*(1-pg)/pg^2; diag(varcov)=vv }
  if(model %in% c(30))   { vv=exp(mu); diag(varcov)=vv }

  if(model %in% c(43))   { mm=exp(mu); pg=pnorm(param$pmu)
                           vv=(1-pg)*mm*(1+pg*mm)
                           diag(varcov)=vv }

  if(model %in% c(45))   {  
                            p=pnorm(param$pmu)
                            MM=pnorm(mu)
                            vv=m*(1-MM)*(1-p)*(1+m*p*(1-MM)) /MM^2
                            diag(varcov)=vv }
}
if(bivariate) {  fname <- "CorrelationMat_biv_dyn_dis2"}      
}        
###############################################################
################################ end discrete #models #########
###############################################################

return(varcov)
}


  #############################################################################################
  #################### end internal function ##################################################
  #############################################################################################
    # Check the user input
    spacetime<-CheckST(CkCorrModel(corrmodel))
    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    ## setting zero mean and nugget if no mean or nugget is fixed
    if(!bivariate){
    if(is.null(param$mean)) param$mean<-0
    if(is.null(param$nugget)) param$nugget<-0  }
    else{
    if(is.null(param$mean_1)) param$mean_1<-0
    if(is.null(param$mean_2)) param$mean_2<-0
    if(is.null(param$nugget_1)) param$nugget_1<-0 
    if(is.null(param$nugget_2)) param$nugget_2<-0 }
    unname(coordt)
    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy)}
    #if the covariance is compact supported  and option sparse is used
    #then set the code as a tapering and an object spam is returned
if(sparse) {
    covmod=CkCorrModel(corrmodel)
    if(covmod %in% c(10,11,13,15,19,6,7,
                     63,64,65,66,67,68,
                     69,70,71,72,73,74,75,76,77,
                     111,112,129,113,114,131,132,130,134,
                     115,116,120))
    {
      type="Tapering"
      if(bivariate){
      #taper="unit_matrix_biv"
      if(covmod %in% c(111,113,115,130)) maxdist=c(param$scale,param$scale,param$scale) 
      if(covmod %in% c(112,114,116,134)) maxdist=c(param$scale_1,param$scale_12,param$scale_2)
      if(covmod %in% c(120,129,131,132)) maxdist=c(param$scale_1,0.5*(param$scale_1+param$scale_2),param$scale_2)  
      }
      if(spacetime)
      {
      maxdist=param$scale_s;maxtime=param$scale_t
       if(covmod==63||covmod==65||covmod==67) {  tapsep=c(param$power2_s,param$power_t,param$scale_s,param$scale_t,param$sep) }
       if(covmod==64||covmod==66||covmod==68) {  tapsep=c(param$power_s,param$power2_t,param$scale_s,param$scale_t,param$sep) }
    }
      if(!(spacetime||bivariate)){
        maxdist=param$scale
        #print(param)
        if(covmod==6)  maxdist=as.numeric(param$scale*exp((lgamma(2*param$smooth+1/param$power2+1)-lgamma(1/param$power2))/ (1+2*param$smooth) ))
        if(covmod==7)  maxdist=as.numeric(param$scale*exp((lgamma(2*param$smooth+param$power2+1)-lgamma(param$power2))/ (1+2*param$smooth) ))
      }
  }
  taper=corrmodel
}

    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, NULL, distance, "Simulation",
                             NULL, grid, NULL, maxdist, maxtime,  model=model, n,  NULL,
                              param, radius, NULL, taper, tapsep,  "Standard", NULL, NULL, NULL,X)
  
    if(!is.null(checkinput$error)) stop(checkinput$error)
    spacetime_dyn=FALSE
    if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
    # Initialising the parameters:
    initparam <- StartParam(coordx, coordy, coordt,coordx_dyn, corrmodel, NULL, distance, "Simulation",
                           NULL, grid, NULL, maxdist, NULL,maxtime, model, n, 
                           param, NULL, NULL, radius, NULL, taper, tapsep,  type, type,
                           NULL, NULL, FALSE, NULL, NULL,NULL,NULL,X,FALSE)
   
    cc=cbind(initparam$coordx,initparam$coordy)  

    if(!spacetime_dyn) dime=initparam$numcoord*initparam$numtime 
    else               dime=sum(initparam$ns)

    if(!initparam$bivariate) numpairstot=dime*(dime-1)*0.5
    if(initparam$bivariate)  numpairstot=dime*(dime-1)*0.5+dime
    if(!is.null(initparam$error)) stop(initparam$error)
    setup<-initparam$setup
    if(initparam$type=="Tapering")
    {
      if(initparam$spacetime) fname= "CorrelationMat_st_tap"
      if(!initparam$spacetime) fname= "CorrelationMat_tap"
      if(initparam$bivariate) fname= "CorrelationMat_biv_tap"
      corr <- double(initparam$numpairs)
      #tapmod <- setup$tapmodel
      ### unit taperssss ####
      if(sparse){
           if(spacetime) tapmod=230
           if(bivariate) tapmod=147
           if(!(spacetime||bivariate)) tapmod=36
           }
      else(tapmod=CkCorrModel(taper))
    #######################
    tp=.C(fname,tapcorr=double(initparam$numpairs),as.double(cc[,1]),as.double(cc[,2]),as.double(initparam$coordt),as.integer(tapmod),
      as.double(1),as.double(tapsep),as.double(1),
      PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
        setup$taps<-tp$tapcorr
    }
    if(is.null(X))  initparam$X=as.matrix(rep(1,dime))

    if(bivariate) {if(is.null(X))  initparam$X=as.matrix(rep(1,initparam$ns[1]+initparam$ns[2])) }
  
    if(spacetime||bivariate){
          initparam$NS=cumsum(initparam$ns);
            if(spacetime_dyn){  initparam$NS=c(0,initparam$NS)[-(length(initparam$ns)+1)]}
            else{               initparam$NS=rep(0,initparam$numtime)}
    }
    if(is.null(initparam$NS)) initparam$NS=0
 
    covmatrix<- Cmatrix(initparam$bivariate,cc[,1],cc[,2],initparam$coordt,initparam$corrmodel,dime,n,initparam$ns,
                        initparam$NS,
                        initparam$param[initparam$namesnuis],
                        initparam$numpairs,numpairstot,initparam$model,
                        initparam$param[initparam$namescorr],setup,initparam$radius,initparam$spacetime,spacetime_dyn,initparam$type,initparam$X)
   
    
    if(type=="Tapering") sparse=TRUE

    # Delete the global variables:
               .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)



    if(initparam$bivariate)   initparam$numtime=2
    # Return the objects list:
    CovMat <- list(bivariate =  initparam$bivariate,
                   coordx = initparam$coordx,
                   coordy = initparam$coordy,
                   coordt = initparam$coordt,
                   coordx_dyn = coordx_dyn,
                   covmatrix=covmatrix,
                   corrmodel = corrmodel,
                   distance = distance,
                   grid=   grid,
                   nozero=initparam$setup$nozero,
                   maxdist = maxdist,
                   maxtime = maxtime,
                   n=n,
                   ns=initparam$ns,
                   NS=initparam$NS,
                   model=initparam$model,
                   namescorr = initparam$namescorr,
                   namesnuis = initparam$namesnuis,
                   namessim = initparam$namessim,
                   numblock = initparam$numblock,
                   numcoord = initparam$numcoord,
                   numcoordx = initparam$numcoordx,
                   numcoordy = initparam$numcoordy,
                   numtime = initparam$numtime,
                   param = initparam$param,
                   radius = initparam$radius,
                   setup=setup,
                   spacetime = initparam$spacetime,
                   sparse=sparse,
                   tapmod=taper,
                   tapsep=tapsep,
                   X=initparam$X)
    structure(c(CovMat, call = call), class = c("CovMat"))
}

