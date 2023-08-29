####################################################
### File name: GeoCorrFct.r
####################################################

GeoCorrFct<- function(x,t=NULL,corrmodel, model="Gaussian",distance="Eucl",  
                                  param, radius=6371,n=1,covariance=FALSE,variogram=FALSE)

{
  
  #############################################################################################
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
 
    if(is.null(CkCorrModel (corrmodel))) stop("The name of the correlation model  is not correct\n")
    if(is.null(CkModel(model)))   stop("The name of the  model  is not correct\n")
    if(!is.numeric(x)) stop("Distances must be numeric\n")
    if(sum(x<0)>=1) stop("Distances must be positive\n")
    if(!is.list(param)) stop("param  must be a list\n")

    spacetime<-CheckST(CkCorrModel(corrmodel))
    bivariate<-CheckBiv(CkCorrModel(corrmodel))

    if(!bivariate) {if(is.null(param$sill)) param$sill=1}




   
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
        #print(nuisance)
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



####### Gaussian
   if(model=="Gaussian"){
        vs=as.numeric(nuisance["sill"])
        
        if(bivariate){  cova11 <- correlation[(0*length(nx)+1):(1*length(nx))]
                        cova12 <- correlation[(1*length(nx)+1):(2*length(nx))]
                        cova22 <- correlation[(3*length(nx)+1):(4*length(nx))]
                               }
        else { 
        cova <- correlation*(1-as.numeric(nuisance["nugget"]))

        }
}

##########################
###### non Gaussian cases
##########################

   if(model=="SkewGausssian") {    
            if(bivariate) {}
              else {
              correlation1=(1-as.numeric(nuisance['nugget']) )*correlation  
              vv=as.numeric(nuisance['sill']);sk=as.numeric(nuisance['skew']);sk2=sk^2;corr2=correlation^2;  
              cc=(2*sk2)*(sqrt(1-corr2) + correlation*asin(correlation)-1)/(pi*vv+sk2*(pi-2)) + (correlation1*vv)/(vv+sk2*(1-2/pi))
              vs=(vv+sk2*(1-2/pi))
              cova=cc; 
               }
                   }
  ##########################################
   if(model=="StudentT")        { if(bivariate) {}
                        else {
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              nu=1/as.numeric(nuisance['df']);sill=as.numeric(nuisance['sill'])
                              vs=sill*(nu)/(nu-2)
                              cc=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,correlation^2))*correlation1)/(2*gamma(nu/2)^2)
                              cova=cc;
                               }
                  }
##########################################

 if(model=="Tukeyh")        { if(bivariate) {}
                        else {
                              correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                              h=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              vs=  sill*(1-2*h)^(-1.5)    
                              cc=(correlation*(1-2*h)^(1.5))/((1-h)^2-(h*correlation)^2)^(1.5)
                              cova=cc;
                             }
                  }
   if(model=="Tukeyh2")        { if(bivariate) {}
                        else {
                              correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                              hr=as.numeric(nuisance['tail1']); hl=as.numeric(nuisance['tail2'])
                              sill=as.numeric(nuisance['sill'])
                              corr=correlation
                              corr[corr>=0.9999999999]=0.9999999999
                             x1=1-(1-corr^2)*hr; y1=1-(1-corr^2)*hl
                             x2=(1-hr)^2-(corr*hr)^2;y2=(1-hl)^2-(corr*hl)^2
                             g=1-hl-hr+(1-corr^2)*hl*hr
                             h1=sqrt(1-corr^2/(x1^2))+(corr/x1)*asin(corr/x1);
                             h2=sqrt(1-corr^2/(y1^2))+(corr/y1)*asin(corr/y1)
                             h3=sqrt(1-corr^2/(x1*y1))+sqrt(corr^2/(x1*y1))*asin(sqrt(corr^2/(x1*y1)))
                             p1=x1*h1/(2*pi*(x2)^(3/2))+corr/(4*(x2)^(3/2))
                             p2=y1*h2/(2*pi*(y2)^(3/2))+corr/(4*(y2)^(3/2))
                             p3=-(x1*y1)^(1/2)*h3/(2*pi*(g)^(3/2))+corr/(4*(g)^(3/2))

                             mm=(hr-hl)/(sqrt(2*pi)*(1-hl)*(1-hr))
                             vv1=0.5*(1-2*hl)^(-3/2)+0.5*(1-2*hr)^(-3/2)-(mm)^2
                             cc=(p1+p2+2*p3-mm^2)/vv1  # correlation
                             vs=as.numeric(nuisance['sill'])*vv1
                              cova=cc; 
                             } 
                  }   
##########################################
  if(model=="Tukeygh")       { if(bivariate) {}
                        else {
                              correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                              rho=correlation
                              tail=as.numeric(nuisance['tail'])
                              sill=as.numeric(nuisance['sill'])
                              eta=as.numeric(nuisance['skew'])
                              rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;

                        if(tail>1e-05&&abs(eta)>1e-05){
                               rho2=rho*rho; eta2=eta*eta; tail2=tail*tail;
                               u=1-tail; a=1+rho;
                               mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                               vs=sill*(exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*sqrt(1-2*tail))-mu*mu;
                               A1=exp(a*eta2/(1-tail*a));
                               A2=2*exp(0.5*eta2*  (1-tail*(1-rho2))  / (u*u- tail2*rho2)  );
                               A3=eta2*sqrt(u*u- rho2*tail*tail);
                               cc=((A1-A2+1)/A3-mu*mu)/vs;
                               cova=cc;
                            } 
                        if(tail<=1e-05&&abs(eta)>1e-05){
                              vs=sill*( -exp(eta^2)+exp(eta^2*2))*eta^(-2)
                              cc= (( -exp(eta^2)+exp(eta^2*(1+rho)))*eta^(-2))/vs
                              cova=cc;
                            } 
                        if(tail>1e-05&&abs(eta)<=1e-05){
                              vs=  sill*(1-2*tail)^(-1.5)     
                              cc=(-rho/((1+h*(tail-1))*(-1+tail+tail*rho)*(1+tail*(-2+tail-tail*rho^2))^0.5))/vs
                              cova=cc;
                            } 
                        if(tail<=1e-05&&abs(eta)<=1e-05){
                            cova=rho;
                            vs=sill
                            }  
                }}


if(model=="Kumaraswamy"||model=="Kumaraswamy2")  { if(bivariate) {}
                            else {
correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
if(model=="Kumaraswamy"){
ga=as.numeric(nuisance['shape2'])
eta=as.numeric(nuisance['shape1'])
}
if(model=="Kumaraswamy2"){
ga=as.numeric(nuisance['shape2'])
eta=log1p(-(1+exp(mu))^(-ga))/log(0.5)
}
mm=eta*beta(1+1/ga,eta)
sill=eta*beta(1+2/ga,eta)-mm^2
###
NN=length(correlation);
res=double(NN)
bb=.C("corr_kuma_vec",as.double(correlation),as.double(eta),as.double(ga), 
     res=as.double(res),
    as.integer(NN),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
rho=bb$res
cova=rho;
vs=sill
 }
}

##########################
 if(model=="SinhAsinh") { if(bivariate) {}
                        else {
                          correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                          corr=correlation
                          d=as.numeric(nuisance['tail']); 
                          e=as.numeric(nuisance['skew']); 
                          sill=as.numeric(nuisance['sill'])
                  
    mm=sinh(e/d)*exp(0.25)*(besselK(.25,(d+1)/(2*d))+besselK(.25,(1-d)/(2*d)))/(sqrt(8*pi))
    vs=sill*(cosh(2*e/d)*exp(0.25)*(besselK(.25,(d+2)/(2*d))+besselK(0.25,(2-d)/(2*d)))/(sqrt(32*pi))-0.5-mm^2)
####### starting extra functions
c1<-function(e,d,n,r){
   U=c(0.5-0.5/d,-0.5/d);L=c(1-1/d,0.5-0.5/d-n/2+r)
   res=exp(-e/d)*2^(-0.5+1.5/d+n/2-r)*gamma(0.5+0.5/d+n/2-r)*hypergeo::genhypergeo(U=U, L=L,0.5)
   return(res)
  }
c2<-function(e,d,n,r){
   U=c(0.5+0.5/d,0.5/d);L=c(1+1/d,0.5+0.5/d-n/2+r)
   res=exp(-e/d)*2^(-0.5-1.5/d+n/2-r)*gamma(0.5-0.5/d+n/2-r)*hypergeo::genhypergeo(U=U, L=L,0.5)
   return(res)
  }
c3<-function(e,d,n,r){
   U=c(0.5+n/2-r,1+n/2-r);L=c(1.5-0.5/d+n/2-r,1.5+0.5/d+n/2-r)
   r1=exp(-e/d)*pi*gamma(1+n-2*r)*hypergeo::genhypergeo(U=U, L=L,0.5)
   r2=d*gamma(1.5-0.5/d+n/2-r)*gamma(1.5+0.5/d+n/2-r)
   return(r1/r2)
   }
I1<-function(e,d,n,r){
   a1=cosh(2*e/d)+sinh(2*e/d);  a2=pracma::sec(0.5*pi/d-0.5*n*pi+pi*r)+pracma::sec(0.5*pi/d+0.5*n*pi-pi*r)*a1; a3=pracma::sec(0.5*pi/d+0.5*n*pi-pi*r)+pracma::sec(0.5*pi/d-0.5*n*pi+pi*r)*a1
   r1=(-1)^(3+n-2*r)*c1(e,d,n,r)-c2(e,d,n,r)+c1(e,d,n,r)*a1; r2=(-1)^(2+n-2*r)*c2(e,d,n,r)*a1+2^(-2-n+2*r)*c3(e,d,n,r)*a2; r3=-(-0.5)^(2+n-2*r)*c3(e,d,n,r)*a3
  return(r1+r2+r3)
     }
SS<-Vectorize(I1, c("r"))
coef<-function(e,d,N){
  mat=NULL;n=1
  while(n<=N){
   r=as.vector(seq(0,trunc(n/2),1))
   res=factorial(n)*SS(e,d,n,r)*(-0.5)^r/(2*sqrt(2*pi)*factorial(n-2*r)*factorial(r))
   mat=c(mat,sum(res));n=n+1}
return(mat)}
CC<-Vectorize(coef, c("N"))
corrsas<-function(e,d,N,vv,rho1){
  mat=NULL;  j=1
  while(j<=N){
     A=CC(e,d,j)^2*rho1^(seq(1:j))/factorial(1:j)
     mat=sum(A)/vv;j=j+1}
    return(mat)}
CorrSAS<-Vectorize(corrsas, c("rho1"))
##########
corr=CorrSAS(e,d,20,vs,corr)
vs=vs*sill
cova=corr; 
                        }
          }
##########################################
 if(model=="TwopieceT")        { if(bivariate) {}
                        else {
                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              nu=1/as.numeric(nuisance['df']); sk=as.numeric(nuisance['skew']);sill=as.numeric(nuisance['sill'])
                              sk2=sk^2
                              vs=sill* ((nu/(nu-2))*(1+3*sk2) - 4*sk2*(nu/pi)*(gamma(0.5*(nu-1))/gamma(0.5*nu))^2)
                              corr2=correlation1^2;sk2=sk^2

                                  a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
                                  a2=correlation1 *asin(correlation1) + (1-corr2)^(0.5)
                                  ll=qnorm((1-sk)/2)
                                  p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                                  a3=3*sk2 + 2*sk + 4*p11 - 1
                                  KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
                                  cc= KK*(a1*a2*a3-4*sk2);
                              ##
                              cova=cc;
                               }
                  } 
  ##########################################
  if(model=="TwopieceTukeyh")        { if(bivariate) {}
                        else {

                              correlation1=correlation*(1-as.numeric(nuisance['nugget'] ))
                              tail=as.numeric(nuisance['tail']); sk=nuisance['skew'];sill=nuisance['sill']
                              corr2=correlation1^2;sk2=sk^2
                              gg2=(1-(1-corr2)*tail)^2
                              xx=corr2/gg2
                              A=(asin(sqrt(xx))*sqrt(xx)+sqrt(1-xx))/(1-xx)^(1.5)
                              ll=qnorm((1-sk)/2)
                              p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                              a3=3*sk2 + 2*sk + 4*p11 - 1
                              mm=8*sk2/(pi*(1-tail)^2); 
                              ff=(1+3*sk2)/(1-2*tail)^(1.5)
                              M=(2*(1-corr2)^(3/2))/(pi*gg2)
                              cc=  (M*A*a3-mm)/( ff- mm)
                              vs=sill*(ff- mm) 
                              cova=cc;
                               }
                  }  
##########################################
 if(model=="TwopieceGauss")        { 
                        if(bivariate) {}
                        else {        
                        correlation1=correlation*(1-as.numeric(nuisance['nugget']))
                        corr2=sqrt(1-correlation1^(2))
                        sk=as.numeric(nuisance['skew']); sk2=sk^2
                        ll=qnorm((1-sk)/2)
                        p11=pbivnorm::pbivnorm(ll,ll, rho = correlation, recycle = TRUE)
                        KK=3*sk2+2*sk+ 4*p11 - 1
                        cc=(2*((corr2 + correlation1*asin(correlation1))*KK)- 8*sk2)/(3*pi*sk2  -  8*sk2   +pi   )
                        vs= as.numeric(nuisance['sill'])*(1+3*sk2-8*sk2/pi)
                        cova=cc;
                        }
                  } 
    ##########################################
 if(model=="Gamma")        { if(bivariate) {}
                        else {
                              correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                              vs=2*exp(mm)^2/as.numeric(nuisance['shape'])
                              cc=correlation^2
                              cova=cc;
                            }
                  }
##########################################
if(model=="Weibull")        { if(bivariate) {} 
                        else {
                         
                        correlation=correlation*(1-as.numeric(nuisance['nugget'] )  )
                        sh=as.numeric(nuisance["shape"])
                        vs=exp(mm)^2*(gamma(1+2/sh)/gamma(1+1/sh)^2-1)
                        auxcorr= (gamma(1+1/sh))^2/((gamma(1+2/sh))-(gamma(1+1/sh))^2)
                        cc=auxcorr*(Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,correlation^2)) -1)
                        cova=cc;
                     }
                    }
##########################################
if(model=="Loglogistic")    { if(bivariate) {}  
                      else { 
                     correlation=correlation*(1-as.numeric(nuisance['nugget'] ))
                     sh=as.numeric(nuisance["shape"])
                     vs=exp(mm)^2*(2*sh*sin(pi/sh)^2/(pi*sin(2*pi/sh))-1)
                     cc=((pi*sin(2*pi/sh))/(2*sh*(sin(pi/sh))^2-pi*sin(2*pi/sh)))*
                                    (Re(hypergeo::hypergeo(-1/sh, -1/sh, 1,correlation^2))*
                                     Re(hypergeo::hypergeo( 1/sh,  1/sh, 1,correlation^2)) -1)
                      cova=cc;
                    }
                    }
if(model=="LogGaussian")    { if(bivariate) {}  
                      else {    
                    vvar=as.numeric(nuisance["sill"])
                    yy=vvar*correlation*(1-as.numeric(nuisance["nugget"]))
                    cc<-(exp(yy)-1)/(exp(vvar)-1)
                    vs<-(exp(vvar)-1)
                    cova=cc; 
                     }
                  }

if(model=="Binomial"||model=="BinomialNeg"||model=="Geometric"||model=="BinomialnegZINB")
   {
                    if(bivariate) {}
                    if(!bivariate) {          
                           pp=pnorm(mu)
                           if(model=="Binomial") vs=min(n)*pp*(1-pp)
                           if(model=="BinomialNeg")      vs=(n)*(1-pp)/pp^2;
                           if(model=="BinomialnegZINB") { 
                                     pg=pnorm(as.numeric(nuisance['pmu']))
                                     vs=n*(1-pp)*(1-pg)*(1+n*pg*(1-pp)) /pp^2
                                               }
                           cova=correlation
                           }
   }
if(model=="Poisson") {
                    if(bivariate) {}
                    if(!bivariate) {   
                           correlation=(1-as.numeric(nuisance['nugget']))*correlation   
                           corr2=correlation^2    
                           vs=exp(mu);
                           z=2*vs/(1-corr2)
                           cc=corr2*(1-(besselI(z,0,expon.scaled = TRUE)+besselI(z,1,expon.scaled = TRUE)))
                           cova=cc
                         }
                      }
                     
#################################################################

if(!covariance) vs=1

 res=cova*vs; 
if(variogram) res=vs*(1-cova)

return(res)

}

