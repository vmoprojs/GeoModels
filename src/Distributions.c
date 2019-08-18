#include "header.h"

// Integrand for hypergeometric computation:



/* integrand  in  hypergeometric function*/
double int_gen_hyp(double x,double a, double b,double z,double c)
{
    double res=0.0;
    res=R_pow(x,b-1)*R_pow(1-x,c-b-1)* R_pow(1-z*x,-a);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}
void integr_gen_hyp(double *x, int n, void *ex){
    int i;double a,b,c,y;
    a =    ((double*)ex)[0];  //a
    b = ((double*)ex)[1];  //b
    c=     ((double*)ex)[2];  //c
    y =     ((double*)ex)[3];  //y
    for (i=0;i<n;i++) {x[i]=int_gen_hyp(x[i],a,b,y,c);}
    return;
}
// function computing generalized wendland
double HyperG_integral(double x, double *param) {
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;         /* as instructed in WRE */
    iwork =   (int *) Calloc(subdiv, int);  /* idem */
    work = (double *) Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=0;
    upper=1;
    // Compute the integral
    Rdqags(integr_gen_hyp, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);

    Free(iwork);Free(work);
    return(result);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/* Wendland covariance */

/* integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp)
{
    double res=0.0,y;
    y=lag/supp;
    res=R_pow(1-x,mu-1)*R_pow(x*x-y*y,alpha)/beta(2*alpha+1,mu);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}
void integr_gen(double *x, int n, void *ex){
    int i;double mu,alpha,beta,y;
    mu =    ((double*)ex)[0];  //mu
    alpha = ((double*)ex)[1];  //alpha
    beta =     ((double*)ex)[2];  //csupp
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_gen(x[i],mu,alpha,y,beta);}
    return;
}
// function computing generalized wendland
double wendintegral(double x, double *param) {
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;		     /* as instructed in WRE */
    iwork =   (int *) Calloc(subdiv, int);  /* idem */
    work = (double *) Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=x/param[2];
    upper=1;
    // Compute the integral
    if(x<=param[2]) {
    Rdqags(integr_gen, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);

    }else   {result=0;}
    Free(iwork);Free(work);
    return(result);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
// Integrand function (derivatives of the Student-t cdf):
double int_pt(double x, double df)
  {
    double res=0.0, x2=0.0, y=0.0;
    x2=R_pow(x,2);
    y=1+x2/df;
    res=0.5*dt(x,df,0)*((df+1)*x2/R_pow(df,2)/y-log(y));
    return res;
  }
// Vectorised integrand function:
void integr_pt(double *x, int n, void *ex)
{
  int i=0;
  double d=0.0;
  d=*((double *)ex);
  for(i=0;i<n;i++)
    x[i]=int_pt(x[i],d);
  return;
}




// compute  bivariate normal standard pdf:
double d2norm(double x, double y, double rho)
{
  double res=0.0, omr=1-R_pow(rho,2);
  res=(1/(2*M_PI))*exp(-0.5*(R_pow(x,2)-2*rho*x*y+R_pow(y,2))/omr)/sqrt(omr);
  return(res);
}


double d2lognorm(double x, double y, double sill,double nugget, double mux,double muy,double rho)
{
  rho=(1-nugget)*rho;
  double KK=exp(sill/2);
  x=x*KK; y=y*KK;
  double res=0.0, q=0.0, omr=R_pow(sill,2)-R_pow(rho*sill,2);
  q=(sill*R_pow((log(x)-mux),2) + sill*R_pow((log(y)-muy),2)
    -2*rho*sill*(log(x)-mux)*(log(y)-muy))/omr;
  res=exp(-q/2)/(2*x*y*M_PI*sqrt(omr));
  return(res*R_pow(KK,2));
}

double biv_sinh(double corr,double zi,double zj,double mi,double mj,double skew,double tail,double vari)
{
double b1=0.0,b2=0.0,A=0.0,B=0.0,k=0.0,res=0.0,Z1,Z2;
double xi=(zi-mi)/(sqrt(vari));
double xj=(zj-mj)/(sqrt(vari));
  b1=tail * asinh(xi)-skew;
  b2=tail * asinh(xj)-skew;
  k=1-R_pow(corr,2);
  A=R_pow(2 * M_PI * R_pow(k,0.5) * (vari),-1) * cosh(b1) * cosh(b2) * R_pow(tail,2)/sqrt((R_pow(xi,2)+1) * (R_pow(xj,2)+1));
  Z1=sinh(b1);Z2=sinh(b2);
  B=exp(- (Z1*Z1 + Z2*Z2 - 2*corr*Z1*Z2)/(2*k)  );
  res=A*B;                     
  return(res);
 } 


// pochammer factorial
double ff(double a,int k) 
{
return(gammafn(k+a)/gammafn(a));
}





double  biv_Weibull2(double rho12,double zi,double zj,double mi,double mj, double shape1,double shape2)

{
double a1=0.0,a2=0.0,a=0.0,b=0.0,c=0.0,dens=.0,k=.0,mui=.0,muj=.0;
a1=gammafn(1+1/shape1);
a2=gammafn(1+1/shape2);
mui=exp(mi);muj=exp(mj);
k=1-rho12*rho12;
a=(shape1*shape2*R_pow(a1,shape1)*R_pow(a2,shape2)*R_pow(zi,shape1-1)*R_pow(zj,shape2-1))/(R_pow(mui,shape1)*R_pow(muj,shape2)*k);
b=exp(-(R_pow(a1*zi/mui,shape1) + R_pow(a2*zj/muj,shape2))/k);
c=bessel_i(2*fabs(rho12)*R_pow(a1*zi/mui,shape1/2)*R_pow(a2*zj/muj,shape2/2)/k,0,1);
dens=a*b*c;
return(dens);
}

double asy_log_besselI(double z,double nu)
{
     double val;
     double K=4*pow(nu,2);
     val =  (z-0.5*log(2*M_PI*z))+
           log((1-(K-1)/(8*z)*(1-(K-9)/(2*8*z)*(1-(K-25)/(3*8*z)))));
     return(val);
}


double biv_Weibull(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double res1=0.0,ui=0.0, uj=0.0,z=0.0,a=0.0,A=0.0,k=0.0,res=0.0,B=0.0;double ci=exp(mui);double cj=exp(muj);;
    k=pow(gammafn(1+1/shape),-1);
    ui=zi/ci;uj=zj/cj;
        a=1-R_pow(corr,2);
        z=2*fabs(corr)*pow(ui*uj,shape/2)*pow(k,-shape)/a;
        A=pow(shape,2)*pow(k,-2*shape)*pow(ui*uj,shape-1)/a;
        B= exp(-pow(k,-shape)*(pow(ui,shape)+pow(uj,shape))/a);
        if(z<700) 
               res=A*B*bessel_i(z,0,1)/(ci*cj);
        else{
               B=-pow(k,-shape)*(pow(ui,shape)+pow(uj,shape))/a;
               res1=(log(A)+B-(log(ci)+log(cj))) + asy_log_besselI(z,0);
               res=exp(res1);
             }
    return(res);
}
/*
double biv_Weibull(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double ui=0.0, uj=0.0,z=0.0,a=0.0,A=0.0,k=0.0,res=0.0,B=0.0;double ci=exp(mui);double cj=exp(muj);;
    k=pow(gammafn(1+1/shape),-1);
    ui=zi/ci;uj=zj/cj;
   // if(corr)   {
        a=1-R_pow(corr,2);
        z=2*fabs(corr)*pow(ui*uj,shape/2)/(pow(k,shape)*a);
        A=pow(shape,2)*pow(ui*uj,shape-1)/(a*pow(k,2*shape));
        B= exp(-(pow(ui,shape)+pow(uj,shape))/(a*pow(k,shape)));
        res=A*B*bessel_i(z,0,1);
    return(res/(ci*cj));
    
}*/


double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double res1=0.0,a=0.0,A=0.0,D=0.0,z=0.0,res=0.0,B=0.0,C=0.0;
    double ci=zi/exp(mui);double cj=zj/exp(muj);
    double gam = gammafn(shape/2);
        a=1-R_pow(corr,2);  

        z=shape*fabs(corr)*sqrt((ci)*(cj))/a;

        A=R_pow((ci)*(cj),shape/2-1) * R_pow(z/2,1-shape/2) ; ///ok
        C= exp(-shape*((ci)+(cj))/(2*a));//ok
        B=gam*R_pow(a,shape/2)*R_pow(2,shape)*R_pow(shape,-shape);  
        D=bessel_i(z,shape/2-1,1); //ok
        if(z<700) 
             res=(A*C*D)/(exp(muj)*exp(muj)*B);
        else{
            C=-shape*((ci)+(cj))/(2*a);
            res1=log(A)+C-(mui-muj-log(B))+asy_log_besselI(z,shape/2-1);
            res=exp(res1);
        }
        return(res);
}
/*******************************
double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double a=0.0,A=0.0,D=0.0,z=0.0,res=0.0,B=0.0,C=0.0;
    double ci=zi/exp(mui);double cj=zj/exp(muj);
    double gam = gammafn(shape/2);
    //if(corr)   {
        a=1-R_pow(corr,2);  
        z=shape*fabs(corr)*sqrt((ci)*(cj))/a;
        A=R_pow((ci)*(cj),shape/2-1) * exp(-shape*((ci)+(cj))/(2*a)); ///ok
        C=R_pow(z/2,1-shape/2); //ok
        B=gam*R_pow(a,shape/2)*R_pow(2,shape)*R_pow(shape,-shape);  
        D=bessel_i(z,shape/2-1,1); //ok
        res=(A*C*D)/(exp(muj)*exp(muj)*B);
        return(res);
}*/


void biv_gamma_call(double *corr,double *zi,double *zj,double *mui, double *muj, double *shape, double *res)
{
    *res = biv_gamma(*corr,*zi,*zj,*mui,*muj,*shape);
}

/***********/
double Hyp_conf_laplace_approx(double a, double b, double z)
{
double y=0.0,r=0.0,A=0.0;
  y=2*a/(b-z+sqrt(R_pow(z-b,2)+ 4*a*z ));

  r=(R_pow(y,2)/a)+(R_pow(1-y,2)/(b-a));

  A=R_pow(b,b-0.5)*R_pow(r,-0.5)*R_pow((y/a),a)*R_pow((1-y)/(b-a),b-a)*exp(z*y);
  return(A);
}
/*****************************************************/
double asym_aprox_1F1(double a, double b, double z)
{
double k1,k2,k3,res;
double alpha=(b-a-1)/z;
if(b<z+a+1)
{
  k1=(gammafn(b)*exp(z)*R_pow(1-alpha,a-1))/(gammafn(a)*gammafn(b-a));
  k2=sqrt(2*M_PI*alpha/z)*R_pow(alpha/exp(1),alpha*z);
  k3=(2-a*alpha)*(1-a)/(2*R_pow(1-alpha,2));
  res=k1*k2*(1+ (k3+ 1/(12*alpha))/z);  
}
if(b>z+a+1)
{
k1=gammafn(b)/(gammafn(b-a)*R_pow(b-a-z-1,a));
k2=a*(a+1)*(a+1-b)/(2*R_pow(b-a-z-1,2));
res=k1*(1-k2);
}
if(b==z+a+1)
{k1=gammafn(b)/(2*gammafn(a)*gammafn(b-a));
 k2=R_pow(2/(b-a-1),a/2);
 k3=2*gammafn(0.5*(a+3))*sqrt(2/(b-a-1))/3;
res=k1*k2*(gammafn(a/2)-k3);}
return(res);
}



double biv_gamma2(double rho,double x,double y,double shape1, double shape2, double rate)
{
double  ss,yr,xr,xi;  
double  ar,ai=0.0;
double  br,bi=0.0;
double  r1,r2,ri=0.0;
int len=1;int lnchf=0;int ip=0;

double res,beta,b,c,f,rho2,a;
int k;

rho2=R_pow(rho,2);

beta=shape1+shape2;

if(rho){
c=rate/(2*(1-rho2));
   
ar=shape2/2;
a=0.0;
for(k=0;k<=20;k++)
{
br=beta/2+k;
xr=c*rho2*x;
yr=c*rho2*y;
//r1=Hyp_conf_laplace_approx(ar,br,c*rho2*x);
//r2=Hyp_conf_laplace_approx(ar,br,c*rho2*y);
F77_NAME(chfm)(&xr,&xi,&ar,&ai,&br,&bi,&r1,&ri,&len,&lnchf,&ip);
F77_NAME(chfm)(&yr,&xi,&ar,&ai,&br,&bi,&r2,&ri,&len,&lnchf,&ip);
if(!R_FINITE(r1*r2)||ISNAN(r1*r2)){
 r1=asym_aprox_1F1(ar,br,c*rho2*x);
 r2=asym_aprox_1F1(ar,br,c*rho2*y);}

ss=ff(shape1/2,k)*r1*r2*R_pow((R_pow(c,2)*rho2*x*y),k)/(gammafn(k+1)*R_pow(ff(beta/2,k),2));
a=a+ss;
}
b=R_pow(rate/2,beta)*R_pow(x*y,beta/2-1)*exp(-c*(x+y));
f=R_pow(gammafn(beta/2),2)*R_pow(1-rho2,shape1/2);

res=b*a/f;}
else
{
  b=(R_pow(rate,beta/2)*R_pow(x,beta/2-1)*exp(-(rate*x/2))) /(R_pow(2,beta/2)*gammafn(beta/2));
  f=(R_pow(rate,beta/2)*R_pow(y,beta/2-1)*exp(-(rate*y/2)) )/(R_pow(2,beta/2)*gammafn(beta/2));
  res=b*f;
}
  return(res);
}


//bivariate skew gaussian distribution
double biv_skew(double corr,double z1,double z2,double mi,double mj,double vari,double skew)
{
   double aux1=0.0,aux11=0.0, aux2=0.0,aux21 = 0.0, aux22=0.0, pdf1=0,pdf2=0,cdf1=0,cdf2=0,quadr,zi,zj;
    zi=z1-mi;
    zj=z2-mj;
    double det,dens,det1,det2,lim1,lim2,a11,a22;
    double nu2  = R_pow(skew,2);
    double om2  = R_pow(vari,1);
    double cor = corr;
                                      // pdf 1
                                       aux1  =  om2 + nu2 ;
                                       aux2  =  cor * aux1;
                                       det1  =  R_pow(aux1,2) - R_pow(aux2,2) ; 
                                       det2  =  aux1 - aux1 * R_pow(cor,2)  ;
                                       quadr =  (1.0/det2)*( (R_pow(zj,2.0)+R_pow(zi,2.0)) - 2.0*cor*zi*zj  ) ;  
                                       pdf1  =  (0.5/M_PI) * (1/sqrt(det1)) * exp(- 0.5 * quadr  ) ;
                                       lim1  = (zi * skew)/(aux1);
                                       lim2  = (zj * skew)/(aux1);
                                       a11   = (om2)/(aux1);
                                       a22   = (cor * om2)/(aux1);
                                       cdf1  =  cdf_norm(lim1,lim2,a11,a22) ;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                       // pdf 2
                                       aux1  =  om2 + nu2 ;
                                       aux2  =  om2 * cor - nu2 * cor ;
                                       det   =  R_pow(aux1,2) - R_pow(aux2,2) ;   
                                       quadr =  (1.0/det)*( (aux1)*(R_pow(zj,2.0)+R_pow(zi,2.0)) - 2.0*aux2*zi*zj  ) ;  
                                       pdf2  =  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
                                       aux1  =  om2 + nu2 ;
                                       aux2  =  om2 * cor - nu2 * cor ;
                                       aux11 = (skew)/( R_pow(aux1,2.0) - R_pow(aux2,2));
                                       aux21 = nu2 * (1 - R_pow(cor,2.0));
                                       aux22 = om2 * (1 + R_pow(cor,2.0));
                                      
                                       lim1  = aux11 * (zi*aux21 + zi*aux22 - 2*zj*om2*cor);
                                       lim2  = aux11 * (zj*aux21 + zj*aux22 - 2*zi*om2*cor); // bien hasta aquí
                                       a11   =  1 - skew*aux11 * (nu2 + om2 - nu2 * R_pow(cor,2.0) + 3*om2 * R_pow(cor,2.0));
                                       a22   =  - cor + (skew) * aux11 * ( nu2  *(cor - R_pow(cor,3.0)) + om2*(cor + R_pow(cor,3.0)) + 2 * cor * om2  ); 
                                       cdf2  =  cdf_norm(lim1,lim2,a11,a22) ; 
  
dens = 2*(pdf1 * cdf1 + pdf2 * cdf2);
return(dens);
}


double biv_wrapped (double alfa,double u, double v,double mi,double mj,double nugget,double sill,double corr)
{
    double x,y,s1=0.0,s12=0.0,quadr=0.0,det=0.0,wrap_gauss=0.0; // 5???
    //2*atan(mean[i])-M_PI
    x=u-2*atan(mi)-M_PI; y=v-2*atan(mj)-M_PI;
    s1=nugget+sill;
    s12=sill*corr; //sill * corr
    det=R_pow(s1,2)-R_pow(s12,2);
    double k1=-alfa,k2=-alfa;
     while(k1<=alfa){
     while(k2<=alfa){
     quadr = -0.5*(1.0/det)*(s1*R_pow((y+2*k1*M_PI),2.0)+s1*R_pow((x+2*k2*M_PI),2.0)
                                                    -2.0*s12*(x+2*k2*M_PI)*(y+2*k1*M_PI));
     wrap_gauss +=  (1/2.0*M_PI)*(1/sqrt(det)*exp(quadr)) ;
     k2 = k2+1;}
     k1 = k1+1;k2 = -alfa;}
return(wrap_gauss);
}



//bivariate skew gaussian bivariate  distribution
double biv_skew2(double corr,double zi,double zj,double vari1,double vari2,double rho ,double skew1,double skew2)
{
   
    double aux1,aux2,aux3,pdf1,pdf2,quadr;
    double dens,det,fact1,fact2, fact3,fact4,factor,c1,lim1,lim2,cdf1,cdf2,a11,a12,a22;
       double taul=R_pow(vari1,0.5); double taur=R_pow(vari2,0.5);
                                       c1= 1.0 /  (1-R_pow(corr,2)) ;
                                       aux1 = R_pow(taul,2) + R_pow(skew1,2) ;
                                       aux2 =  taul*taur*corr + skew1*skew2*corr ;
                                       aux3 = R_pow(taur,2) + R_pow(skew2,2) ;
                                       det =   aux1*aux3 - R_pow(aux2,2) ;   

                                       quadr= (1.0/det)*( aux1*R_pow(zj,2.0)+aux3*R_pow(zi,2.0) - 2.0*aux2*zi*zj  ) ;  
                                       pdf1=  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
                                       //////
                                       aux1 = (R_pow(skew1/taul,2)+1)*c1;
                                       aux2 = ((skew1*skew2)/(taul*taur)+1)*c1*corr;
                                       aux3 = (R_pow(skew2/taur,2)+1)*c1;

                                       det =   aux1*aux3 - R_pow(aux2,2) ;  
                                       a11  =  (1/det)* aux3   ;
                                       a12  =  (1/det)* aux2   ;
                                       a22  =  (1/det)* aux1   ;
                                       factor = c1 /  det ;
                                       fact1 = factor * (aux3*(skew1/R_pow(taul,2)) - aux2*(skew2/(taul*taur))*corr) ;
                                       fact2 = factor * (-aux3*(skew1/(taul*taur))*corr + aux2*(skew2/R_pow(taur,2)) ) ;
                                       fact3 = factor * (aux2*(skew1/R_pow(taul,2)) - aux1*(skew2/(taul*taur))*corr) ;
                                       fact4 = factor * (-aux2*(skew1/(taul*taur))*corr + aux1*(skew2/R_pow(taur,2)) ) ;

                                       lim1 =   fact1*zi +  fact2*zj   ;
                                       lim2 =   fact3*zi +  fact4*zj   ;
                                       cdf1 = cdf_norm2(lim1,lim2,a11,a12,a22) ;
                                       //////////////////////////////////////////////////////////////////////////////////
                                       aux1 = R_pow(taul,2) + R_pow(skew1,2) ;
                                       aux2 =  taul*taur*corr - skew1*skew2*corr ;
                                       aux3 = R_pow(taur,2) + R_pow(skew2,2) ;
                                       det =   aux1*aux3 - R_pow(aux2,2) ;   
                                       quadr= (1.0/det)*( aux1*R_pow(zj,2.0)+aux3*R_pow(zi,2.0) - 2.0*aux2*zi*zj  ) ;  
                                       pdf2=  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
                                       //////
                                       aux1 = R_pow(skew1/taul,2)*c1 + c1 ;
                                       aux2 = (skew1*skew2)/(taul*taur)*c1*corr - c1*corr;
                                       aux3 = R_pow(skew2/taur,2)*c1 + c1 ;
                                       det =   aux1*aux3 - R_pow(aux2,2) ;  
                                       a11  =  (1/det)* aux3   ;
                                       a12  =  (1/det)* aux2   ;
                                       a22  =  (1/det)* aux1   ;
                                       factor = c1 /  det ;
                                       fact1 = factor * (aux3*(skew1/R_pow(taul,2)) - aux2*(skew2/(taul*taur))*corr) ;
                                       fact2 = factor * (-aux3*(skew1/(taul*taur))*corr + aux2*(skew2/R_pow(taur,2)) ) ;
                                       fact3 = factor * (aux2*(skew1/R_pow(taul,2)) - aux1*(skew2/(taul*taur))*corr) ;
                                       fact4 = factor * (-aux2*(skew1/(taul*taur))*corr + aux1*(skew2/R_pow(taur,2) )) ;
                                       lim1 =   fact1*zi +  fact2*zj   ;
                                       lim2 =   fact3*zi +  fact4*zj   ;
                                       cdf2 = cdf_norm2(lim1,lim2,a11,a12,a22) ; 
                            dens=2*(pdf1*cdf1+pdf2*cdf2);

return(dens);
}


// trivariate skew gaussian distribution/
double triv_skew(double x,double c_0i,double c_0j, double rho,double data_i,double data_j,double *nuis)
{

    int N=3,m,n,i, maxpts=3*2000,fail=100;
    double esterror=10.0,abseps=1.0e-6, releps=0.0;
    double *K;K=(double *) Calloc(N,double);
    double *lower;lower=(double *) Calloc(N,double);
    double *upper;upper=(double *) Calloc(N,double);
    double *corr;corr=(double *) Calloc(N,double);
    int *infin;infin=(int *) Calloc(N,int);
    lower[0]=0;lower[1]=0;lower[2]=0;
    infin[0]=0;infin[1]=0;infin[2]=0;

    double **Om;
    Om= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){Om[i]=(double *) Calloc(N,double);}

    double **C;
    C= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){C[i]=(double *) Calloc(N,double);}

    double **M;
    M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}

    double **cov;
    cov= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){cov[i]=(double *) Calloc(N,double);}
    
    double **inverse;
    inverse= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){inverse[i]=(double *) Calloc(N,double);}  

    double **inverse1;
    inverse1= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){inverse1[i]=(double *) Calloc(N,double);} 

    int *indx; indx=(int *) Calloc(N,int);
    double *col; col=(double *) Calloc(N,double);
    double *data;  data=(double *) Calloc(N,double);

    double det,pdf1,pdf2,pdf3,pdf4,dd;
    double dens,cdf1,cdf2,cdf3,cdf4;
    double nu2=R_pow(nuis[3],2.0),tau2 =nuis[2];
    double kk=tau2+nu2,co=tau2/nu2;
    double co1=co+1,com1=co-1;
    double ratio=nuis[3]/tau2;
  

  /**********    omega and its inverse*/
    Om[0][0]=1;   Om[0][1]=rho ;Om[0][2]=c_0i;
    Om[1][0]=Om[0][1];   Om[1][1]=1;     Om[1][2]=c_0j;
    Om[2][0]=Om[0][2];   Om[2][1]=Om[1][2];  Om[2][2]=1;
     ludcmp(Om,N,indx,&dd);    //Lu decomposition
     //for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
     for(n=0;n<N;n++) {
     for(m=0;m<N;m++) col[m]=0.0;
     col[n]=1.0;
     lubksb(Om,N,indx,col);
    for(m=0;m<N;m++) inverse1[m][n]=col[m];  // inverse using LU decomposition
      }
     /**********/

    data[0]=data_i-nuis[0] ;data[1]=data_j-nuis[0] ;data[2]=x-nuis[0] ;
    det=co*(R_pow(c_0j,2)-2*c_0i*rho*c_0j+R_pow(rho,2)+R_pow(c_0i,2)-1);
/********************************************************************************************/
/********************************************************************************************/
                    M[0][0]=kk     ;   M[0][1]=kk*rho ;M[0][2]=kk*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=kk*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf1=  dNnorm(N,M,data);
                    cov[0][0]= (co1)*(c_0j-1)*(c_0j+1)/det;   
                    cov[0][1]=-(co1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=-(co1)*(rho*c_0j-c_0i)/det;
                    cov[1][1]= (co1)*(c_0i-1)*(c_0i+1)/det;
                    cov[1][2]=(co1)*(c_0j-c_0i*rho)/det;
                    cov[2][2]=(co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  

                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                   // for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }

                    Matrix_prod(inverse,inverse1,C,N);
                    Matrix_prod_vec(C,data,K,N);

                    upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                    upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                    upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                    corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];
                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf1, &fail );
/********************************************************************************************/
/********************************************************************************************/
                    M[0][0]=kk     ;   M[0][1]=(tau2-nu2)*rho ;M[0][2]=(tau2-nu2)*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=kk*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf2= dNnorm(N,M,data);
                    cov[0][0]=(co1)*(c_0j-1)*(c_0j+1)/det;
                    cov[0][1]=(com1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=(com1)*(rho*c_0j-c_0i)/det;
                    cov[1][2]=(co1)*(c_0j-c_0i*rho)/det;
                    cov[1][1]=(co1)*(c_0i-1)*(c_0i+1)/det; 
                    cov[2][2]= (co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  
                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                    //for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }
                     Matrix_prod(inverse,inverse1,C,N);
                     Matrix_prod_vec(C,data,K,N);
                   
                   upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                   upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                   upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                   corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];

                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf2, &fail );
/********************************************************************************************/
/********************************************************************************************/
                    M[0][0]=kk     ;   M[0][1]=(tau2-nu2)*rho ;M[0][2]=kk*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=(tau2-nu2)*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf3= dNnorm(N,M,data);

                    cov[0][0]=(co1)*(c_0j-1)*(c_0j+1)/det;   
                    cov[0][1]=(com1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=-(co1)*(rho*c_0j-c_0i)/det;
                    cov[1][1]=(co1)*(c_0i-1)*(c_0i+1)/det;
                    cov[1][2]=-(com1)*(c_0j-c_0i*rho)/det;
                    cov[2][2]=(co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  

                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                   // for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }
                     Matrix_prod(inverse,inverse1,C,N);
                     Matrix_prod_vec(C,data,K,N);
                   upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                   upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                   upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                   corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];

                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf3, &fail );
/********************************************************************************************/
/********************************************************************************************/

                    M[0][0]=kk;   M[0][1]=kk*rho; M[0][2]=(tau2-nu2)*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=(tau2-nu2)*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf4= dNnorm(N,M,data);

                    cov[0][0]=(co1)*(c_0j-1)*(c_0j+1)/det;   
                    cov[0][1]=-(co1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=(com1)*(rho*c_0j-c_0i)/det ;
                    cov[1][1]=(co1)*(c_0i-1)*(c_0i+1)/det; 
                    cov[1][2]=-(com1)*(c_0j-c_0i*rho)/det;
                    cov[2][2]=(co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  
                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                    //for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }
                    
                     Matrix_prod(inverse,inverse1,C,N);
                     Matrix_prod_vec(C,data,K,N);
                 upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                 upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                 upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                 corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];
                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf4, &fail );
/********************************************************************************************/
/********************************************************************************************/
dens=2*(pdf1*cdf1 + pdf2*cdf2 + pdf3*cdf3 + pdf4*cdf4);
 Free(lower);Free(upper);Free(corr);Free(infin);Free(K);
 Free(indx);Free(col);
 for(i=0;i<N;i++)  {Free (inverse[i]);Free(inverse1[i]);Free(Om[i]);Free (M[i]);Free (cov[i]);Free (C[i]);}
Free(inverse);Free(M);Free(cov);Free(C);Free(Om);Free(inverse1);
 return(dens);
}


double log_biv_Norm(double corr,double zi,double zj,double mi,double mj,double vari, double nugget)
{
    double u,v,u2,v2,det,s1,s12,dens;
    u=zi-mi;
    v=zj-mj;
    u2=R_pow(u,2);v2=R_pow(v,2);
    s1=vari+nugget;
    s12=vari*corr;
    det=R_pow(s1,2)-R_pow(s12,2);
    dens=(-0.5*(2*log(2*M_PI)+log(det)+(s1*(u2+v2)-2*s12*u*v)/det));
return(dens);
}

/*multivariate gaussian PDF*/
double dNnorm(int N,double **M, double *dat)
{
 double dd,pdf; int m=0,n=0,i=0;
 int *indx; indx=(int *) Calloc(N,int);
 double *col; col=(double *) Calloc(N,double);
 double **inverse;
 inverse= (double **) Calloc(N,double *);
 for(i=0;i<N;i++){inverse[i]=(double *) Calloc(N,double);} 

 ludcmp(M,N,indx,&dd);    //Lu decomposition
 for(m=0;m<N;m++){ dd *= M[m][m];}  //determinant using LU
 for(n=0;n<N;n++) {
 for(m=0;m<N;m++) col[m]=0.0;
 col[n]=1.0;
 lubksb(M,N,indx,col);
 for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
 }
 pdf=R_pow(2*M_PI,-N/2)*R_pow(dd,-0.5)*exp(-0.5*QFORM(inverse,dat,dat,N)) ;   
 Free(indx);Free(col);
 for(i=0;i<N;i++)  {Free (inverse[i]);}Free(inverse);
return(pdf);
}

void mult_pmnorm( int *nvar , double *lower , double *upper , int *infin , double *corr , int *maxpts , double *abseps , double *releps , double *esterror , double *result , int *fail )
{ 
  F77_CALL(sadmvn)( nvar, lower, upper, infin, corr, maxpts, abseps, releps, esterror, result, fail ) ;
} 

// compute the trivariate normal cdf  for bernoulli  RF :
double ptnorm(int which,int *cormod, double h0,double h1,double h2, double u0, double u1,double u2, 
                 double *nuis, double *par, double thr)
{  
  int N=3; int maxpts=3*2000;
  double esterror=10.0,abseps=1.0e-6, releps=0.0,res2=0.0;int fail=100;
  double a=(nuis[0]-thr)/(sqrt(nuis[2]+nuis[1]));
  double *lower;double *upper; int *infin;
  lower=(double *) Calloc(3,double);
  upper=(double *) Calloc(3,double);
  infin=(int *) Calloc(3,int);
  switch(which) //
              {
         case 1: /// 111
         lower[0]=0;lower[1]=0; lower[2]=0; 
         upper[0]=a;upper[1]=a; upper[2]=a;
         infin[0]=0;infin[1]=0;infin[2]=0;
            break;
            case 2: /// 110
         lower[0]=0;lower[1]=0; lower[2]=a; 
         upper[0]=a;upper[1]=a; upper[2]=0 ;
         infin[0]=0;infin[1]=0;infin[2]=1;     
            break;
            case 3: /// 101
         lower[0]=0;lower[1]=a; lower[2]=0; 
         upper[0]=a;upper[1]=0 ; upper[2]=a;
         infin[0]=0;infin[1]=1;infin[2]=0;      
            break;
            case 4: /// 011
         lower[0]=a;lower[1]=0; lower[2]=0; 
         upper[0]=0;upper[1]=a; upper[2]=a;
         infin[0]=1;infin[1]=0;infin[2]=0;        
            break;
 }
  double corr[3]={nuis[2]*CorFct(cormod,h0,u0,par,0,0),
                  nuis[2]*CorFct(cormod,h1,u1,par,0,0),
                  nuis[2]*CorFct(cormod,h2,u2,par,0,0)};


  mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &res2, &fail );
  //Free(lim_inf);Free(lim_sup);Free(infin);
  Free(lower);Free(upper);Free(infin);
  return(res2);
}


// trivariate CDF gaussian 
double p3norm(int *cormod, double h0,double h1,double h2, double u0, double u1,double u2, 
                 double *nuis, double *par)
{  
  int N=3; int maxpts=3*2000;
  double esterror=10.0,abseps=1.0e-6, releps=0.0,res2=0.0;int fail=100;
  double a=1;
  double *lower;double *upper; int *infin;
  lower=(double *) Calloc(3,double);
  upper=(double *) Calloc(3,double);
  infin=(int *) Calloc(3,int);

 lower[0]=0;lower[1]=0; lower[2]=0; 
 upper[0]=a;upper[1]=a; upper[2]=a;
 infin[0]=0;infin[1]=0;infin[2]=0;

 double corr[3]={nuis[2]*CorFct(cormod,h0,u0,par,0,0),
                  nuis[2]*CorFct(cormod,h1,u1,par,0,0),
                  nuis[2]*CorFct(cormod,h2,u2,par,0,0)};

  mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &res2, &fail );
  Free(lower);Free(upper);Free(infin);
  return(res2);
}



// cdf of  a bivariate Gausssian distribution
double cdf_norm(double lim1,double lim2,double a11,double a12)
{
    double res=0;
    double  lowe[2]={0,0}, uppe[2]={lim1/sqrt(a11),lim2/sqrt(a11)}, corre[1] ={a12/a11};
    int infin[2]={0,0};
    double auxil=1-R_pow(corre[0],2);
    double det=R_pow(a11,2)-R_pow(a12,2) ;
    res= a11* sqrt(auxil/det) *  F77_CALL(bvnmvn)(lowe,uppe,infin,corre);
    return(res);
}
// cdf of  a bivariate Gausssian distribution
double cdf_norm2(double lim1,double lim2,double a11,double a12, double a22)
{
  double res=0;

  double  lowe[2]={0,0}, uppe[2]={lim1/sqrt(a11),lim2/sqrt(a22)}, corre[1] ={a12/sqrt(a11*a22)};
  int infin[2]={0,0};
  double auxil=1-R_pow(corre[0],2);
  double det=a11*a22-R_pow(a12,2) ;
  res= sqrt(a11*a22)* sqrt(auxil/det) *  F77_CALL(bvnmvn)(lowe,uppe,infin,corre);
  return(res);
}


// compute the bivariate normal cdf for the bernoulli RF:
// compute the bivariate normal cdf for the bernoulli RF:
double pbnorm(int *cormod, double h, double u, double mean1, double mean2, 
  double nugget, double var,double *par, double thr)
{
  double res=0;
  double lim_inf[2]={0,0};//lower bound for the integration
  double lim_sup[2]={mean1,mean2};
  int infin[2]={0,0};//set the bounds for the integration
  double corr[1]={(1-nugget)*CorFct(cormod,h,u,par,0,0)};
  res=F77_CALL(bvnmvn)(lim_inf,lim_sup,infin,corr);
  return(res);
}



// compute a sequence (depending on space and time) of bivariate normal cdf:
void vpbnorm(int *cormod, double *h, double *u, int *nlags, int *nlagt,
	     double *nuis, double *par, double *rho, double *thr)
{
  int i,j,t=0;
  for(j=0;j<*nlagt;j++)
    for(i=0;i<*nlags;i++){
      rho[t]=pbnorm(cormod,h[i],u[j],nuis[0],nuis[0],nuis[1],nuis[2],par,thr[0]);
      t++;}
  return;
}

/*

double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11)
{
  int a=0,i=0;
  double kk1=0.0,kk2=0.0,dens1=0.0,dens2=0.0;

    for(a=fmax_int(0,u-v+NN-1);a<=NN-2;a++){
   for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=fac(NN-1+u,1)/(fac(i,1)*fac(NN-1-i,1)*fac(a-i,1)*fac(u-a+i,1));
    kk2=fac(v-u-1,1)/(fac(v+a-NN-u+1,1)*fac(NN-a-2,1));
    dens1+=kk1*kk2*R_pow(p11,i+1)*R_pow(1+p11-(x+y),u-a+i)*
             R_pow(x-p11,NN-i-1)*R_pow(y-p11,a-i)*R_pow(1-y,v-u-NN+a+1)*R_pow(y,NN-a-1);
  }}

    for(a=fmax_int(0,u-v+NN);a<=NN-1;a++){
    for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=fac(NN-1+u,1)/(fac(i,1)*fac(NN-1-i,1)*fac(a-i,1)*fac(u-a+i,1));
    kk2=fac(v-u-1,1)/(fac(v+a-NN-u,1)*fac(NN-a-1,1));
    dens2+=kk1*kk2*R_pow(p11,i)*R_pow(1+p11-(x+y),u-a+i)*
                 R_pow(x-p11,NN-i)*R_pow(y-p11,a-i)*R_pow(1-y,v-u-NN+a)*R_pow(y,NN-a);
  }}
  return(dens1+dens2);
}

// bivariate negative  binomial
double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
double kk1=0.0,dens=0.0;int i=0;
if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);

if(u==v)            {
    for(i=fmax_int(0,NN-u-1);i<=NN-1;i++){
      kk1=fac(NN-1+u,1)/(fac(i,1)*R_pow(fac(NN-1-i,1),2)*fac(u-NN+1+i,1));
      dens+=kk1*R_pow(p11,i+1)*R_pow(1+p11-(p01+p10),u-NN+1+i)*R_pow(p01-p11,NN-1-i)*R_pow(p10-p11,NN-1-i); }
}

if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
return(dens);
}
*/


double aux_biv_binomneg_simple(int NN, int u, double p01,double p10,double p11)
{
          int i=0;
  double kk1=0.0,dens=0.0;

    for(i=fmax_int(0,NN-u-1);i<=NN-1;i++){
      kk1=exp(lgammafn(NN-1+u+1)-(lgammafn(i+1)+lgammafn(NN-i)+lgammafn(NN-i)+lgammafn(u-NN+2+i)));
      dens+=kk1*pow(p11,i+1)*pow(1+p11-(p01+p10),u-NN+1+i)*pow(p01-p11,NN-1-i)*pow(p10-p11,NN-1-i);
       }
         return(dens);
}



double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11)
{
  int a=0,i=0;
  double kk1=0.0,kk2=0.0,dens1=0.0,dens2=0.0;

    for(a=fmax_int(0,u-v+NN-1);a<=NN-2;a++){
   for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=exp(lgammafn(NN+u)-(lgammafn(i+1)+lgammafn(NN-i)+lgammafn(a-i+1)+lgammafn(u-a+i+1)));         
    kk2=exp(lgammafn(v-u)-(lgammafn(v+a-NN-u+2)+lgammafn(NN-a-1)));                        
    dens1+=kk1*kk2*pow(p11,i+1)*pow(1+p11-(x+y),u-a+i)*
             pow(x-p11,NN-i-1)*pow(y-p11,a-i)*pow(1-y,v-u-NN+a+1)*pow(y,NN-a-1);
  }}

    for(a=fmax_int(0,u-v+NN);a<=NN-1;a++){
    for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=exp(lgammafn(NN+u)-(lgammafn(i+1)+lgammafn(NN-1-i+1)+lgammafn(a-i+1)+lgammafn(u-a+i+1)));
    kk2=exp(lgammafn(v-u)-(lgammafn(v+a-NN-u+1)+lgammafn(NN-a)));                      
    dens2+=kk1*kk2*pow(p11,i)*pow(1+p11-(x+y),u-a+i)*
                 pow(x-p11,NN-i)*pow(y-p11,a-i)*pow(1-y,v-u-NN+a)*pow(y,NN-a);
  }}
  return(dens1+dens2);
}

/*
double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
double kk1=0.0,dens=0.0;int i=0;
if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);

if(u==v)            {
    for(i=fmax_int(0,NN-u-1);i<=NN-1;i++){
      kk1=exp(lgammafn(NN-1+u+1)-(lgammafn(i+1)+lgammafn(NN-i)+lgammafn(NN-i)+lgammafn(u-NN+2+i)));
      dens+=kk1*pow(p11,i+1)*pow(1+p11-(p01+p10),u-NN+1+i)*pow(p01-p11,NN-1-i)*pow(p10-p11,NN-1-i); }
}

if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
return(dens);
}*/


double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
double dens=0.0;
if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);

if(u==v)      dens=aux_biv_binomneg_simple(NN,v,p10,p01,p11);

if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
return(dens);
}

double bin_aux(int a,int NN,int u,int v,double p1, double p2,double p11)
{
  double kk,dens;
kk=exp(lgammafn(NN+1)-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)+lgammafn(NN-u-v+a+1)));
dens=kk*(R_pow(p11,a)*R_pow(p1-p11,u-a)*R_pow(p2-p11,v-a)*R_pow(1+p11-(p1+p2),NN-u-v+a));
return(dens);
}

double biv_binom(int NN, int u, int v, double p01,double p10,double p11)
{
    
int a;
double kk=0.0,dens=0.0;
for(a=fmax_int(0,u+v-NN);a<=fmin_int(u,v);a++){
//kk=exp(logfac(NN)-(logfac(a)+logfac(u-a)+logfac(v-a)+logfac(NN-u-v+a)));
kk=exp(lgammafn(NN+1)-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)+lgammafn(NN-u-v+a+1)));
dens+=kk*(R_pow(p11,a)*R_pow(p01-p11,u-a)*R_pow(p10-p11,v-a)*R_pow(1+p11-(p01+p10),NN-u-v+a));
 }
    return(dens);
}


/*
double biv_binom (int NN, int u, int v, double p01,double p10,double p11)
{
    
int a;
double kk=0,dens=0.0;
for(a=fmax_int(0,u+v-NN);a<=fmin_int(u,v);a++){
kk=fac(NN,1)/(fac(a,1)*fac(u-a,1)*fac(v-a,1)*fac(NN-u-v+a,1));
dens+=kk*(R_pow(p11,a)*R_pow(p01-p11,u-a)*R_pow(p10-p11,v-a)*R_pow(1+p11-(p01+p10),NN-u-v+a));
 }
    return(dens);
}*/



/// biv binomial type II
/*
double  biv_binom2 (int NN_i,int NN_j, int k, int u, int v, double p01,double p10,double p11)
{

int a,i,j;
double const1=0.0,const2=0.0,const3=0.0,dens=0.0,dens1;
double P10=p01-p11;
double P01=p10-p11;
double P00=1+p11-(p01+p10);
double P11=p11;

for(i=0;i<=fmin_int(NN_i-k,u);i++){
   for(j=0;j<=fmin_int(NN_j-k,v);j++){
     
for(a=fmax_int(0,u+v-k-i-j);a<=fmin_int(u-i,v-j);a++){         
       const1=fac(k,1)/(fac(a,1)*fac(u-i-a,1)*fac(v-j-a,1)*fac(k-u-v+i+j+a,1));
        dens1=const1*R_pow(P11,a)*R_pow(P00,k-u-v+i+j+a)*
             R_pow(P10,u-i-a)*R_pow(P01,v-j-a);
   
       const2=fac(NN_i-k,1)/(fac(NN_i-k-i,1)*fac(i,1));
       const3=fac(NN_j-k,1)/(fac(NN_j-k-j,1)*fac(j,1));
       //Rprintf("%f %f %f\n",const2,const3,dens1);
       dens+=dens1*const2*const3*
             R_pow(P11+P10,i) * R_pow(P11+P01,j) *
             R_pow(P00+P01,NN_i-k-i) * R_pow(P00+P10,NN_j-k-j);
  }}}
    return(dens);
}
*/


double  biv_binom2 (int NN_i,int NN_j, int k, int u, int v, double p01,double p10,double p11)
{

int a,i,j;
double const1=0.0,const2=0.0,const3=0.0,dens=0.0,dens1;
double P10=p01-p11;
double P01=p10-p11;
double P00=1+p11-(p01+p10);
double P11=p11;

for(i=0;i<=fmin_int(NN_i-k,u);i++){
   for(j=0;j<=fmin_int(NN_j-k,v);j++){
     
for(a=fmax_int(0,u+v-k-i-j);a<=fmin_int(u-i,v-j);a++){         
       const1=exp(lgammafn(k+1)-(lgammafn(a+1)+lgammafn(u-i-a+1)+lgammafn(v-j-a+1)+lgammafn(k-u-v+i+j+a+1)));
        dens1=const1*R_pow(P11,a)*R_pow(P00,k-u-v+i+j+a)*
             R_pow(P10,u-i-a)*R_pow(P01,v-j-a);
   
       const2=exp(lgammafn(NN_i-k+1)-(lgammafn(NN_i-k-i+1)+lgammafn(i+1)));
       const3=exp(lgammafn(NN_j-k+1)-(lgammafn(NN_j-k-j+1)+lgammafn(j+1)));

       dens+=dens1*const2*const3*
             R_pow(P11+P10,i) * R_pow(P11+P01,j) *
             R_pow(P00+P01,NN_i-k-i) * R_pow(P00+P10,NN_j-k-j);
  }}}
    return(dens);
}

/*
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0;
    a_u=fmin(u,v);
    dens=0;
    for(a=0;a<=a_u;a++){
        dens+=(R_pow(NN*p11,a)*R_pow(NN*(p01-p11),u-a)*R_pow(NN*(p10-p11),v-a))/(fac(a,1)*fac(u-a,1)*fac(v-a,1));
    }
    return(exp(-NN*(p01+p10-p11))*dens);
}



double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0;
    a_u=fmin_int(u,v);
    dens=0;
    pp=p11*p01*p10;
    bb=(1-p11)/p11;
    for(a=0;a<=a_u;a++){
       dens+=(R_pow(NN*((p11*(p01+p10)-p01*p10*(p11+1))/pp),a)*R_pow(NN*((p10-p11)/(p10*p11)),u-a)*R_pow(NN*bb,v-a))/(fac(a,1)*fac(u-a,1)*fac(v-a,1));
    }
    return(exp(-NN*bb)*dens);
}
*/

// bivariate pois-binomial 
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0, kk=0.0;
    a_u=fmin(u,v);
    dens=0;
    for(a=0;a<=a_u;a++){

    kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
        dens+=kk*(R_pow(NN*p11,a)*R_pow(NN*(p01-p11),u-a)*R_pow(NN*(p10-p11),v-a));
    }
    return(exp(-NN*(p01+p10-p11))*dens);
}



// bivariate pois-bineg 
double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0, kk=0.0;
    a_u=fmin_int(u,v);
    dens=0;
    pp=p11*p01*p10;
    bb=(1-p11)/p11;
    for(a=0;a<=a_u;a++){

       kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
       dens+=kk*(R_pow(NN*((p11*(p01+p10)-p01*p10*(p11+1))/pp),a)*R_pow(NN*((p10-p11)/(p10*p11)),u-a)*R_pow(NN*bb,v-a));
    }
    return(exp(-NN*bb)*dens);
}



// marginal binomial
double marg_binom(int n,double x,double p)
{
 double pr=0.0;double kk;
 kk=fac(n,n-x+1)/fac(x,1);
 pr=kk*R_pow(p,x)*R_pow(1-p,n-x);
return(pr);
}
//marginal geom
double marg_geom(int x,double p)
{
 double pr=0.0;
 pr=p*R_pow(1-p,x);
    return(pr);
}
//marginal binom neg
double marg_binomneg(int n,int x,double p)
{
 double pr=0.0;
 pr=(fac(x+n-1,1)/(fac(n-1,1)*fac(x,1)))*R_pow(p,n)*R_pow(1-p,x);
    return(pr);
}
//marginal poisson
double marg_pois(int n,double x,double p)
{
 double pr=0.0,lambda=0.0;
 lambda=n*p;
 pr=R_pow(lambda,x)*exp(-lambda)/fac(x,1);
    return(pr);
}
                    

double marg_p(double categ_0,double psm,int *model,int n)
{
    double res=0.0;
    if(model[0]==2||model[0]==11) res=marg_binom (n-1,categ_0,psm);
    if(model[0]==14)           res=marg_geom (categ_0,psm);
    if(model[0]==15)           res=marg_pois (n-1,categ_0,psm);
    return(res);
}

/************ pochammer symbols*******************/
double Poch(double q,double n)
{
    return(gammafn(q+n)/gammafn(q));
}
/*********** ratio of gamma function type 1********************/
double aprox(int k,double a, double b) 
{
return(R_pow(k,a-b)*(1+ (a-b)*(a+b-1)/(2*k)));
}
/************ ratio of gamma function type 2 *******************/
double aprox1(double a) {
  double res;
    res=R_pow(R_pow(a,5)+(5*R_pow(a,4)/4)+(25*R_pow(a,3)/32)+(35*R_pow(a,2)/128)+(75*a/2048)+3/8192,1/10);
    return(res);
 }


/**************** for bivaraite T distribution */////////////////////////
double hyp2f1( double a,double b,double c,double x)
{
    double d, d1, d2, e;
    double p, q, r, s, y, ax;
    double ia, ib, ic, id, err;
    int flag, i, aid;
    
    err = 0.0;
    ax = fabs(x);
    s = 1.0 - x;
    flag = 0;
    ia = round(a); /* nearest integer to a */
    ib = round(b);

    if (x == 0.0) {return 1.0; }
    
    if( a <= 0 )
    {
        if( fabs(a-ia) < EPS1 )    /* a is a negative integer */
            flag |= 1;
    }
    
    if( b <= 0 )
    {
        if( fabs(b-ib) < EPS1 )    /* b is a negative integer */
            flag |= 2;
    }
    
    if( ax < 1.0 )
    {
        if( fabs(b-c) < EPS1 )   /* b = c */
        {
            y = pow( s, -a ); /* s to the -a power */
            goto hypdon;
        }
        if( fabs(a-c) < EPS1 )   /* a = c */
        {
            y = pow( s, -b ); /* s to the -b power */
            goto hypdon;
        }
    }
    if( c <= 0.0 )
    {
        ic = round(c);  /* nearest integer to c */
        if( fabs(c-ic) < EPS1 )    /* c is a negative integer */
        {
            /* check if termination before explosion */
            if( (flag & 1) && (ia > ic) )
                goto hypok;
            if( (flag & 2) && (ib > ic) )
                goto hypok;
            goto hypdiv;
        }
    }
    
    if( flag )      /* function is a polynomial */
        goto hypok;
    
    if( ax > 1.0 )      /* series diverges  */
        goto hypdiv;
    
    p = c - a;
    ia = round(p); /* nearest integer to c-a */
    if( (ia <= 0.0) && (fabs(p-ia) < EPS1) ) /* negative int c - a */
        flag |= 4;
    
    r = c - b;
    ib = round(r); /* nearest integer to c-b */
    if( (ib <= 0.0) && (fabs(r-ib) < EPS1) ) /* negative int c - b */
        flag |= 8;
    
    d = c - a - b;
    id = round(d); /* nearest integer to d */
    q = fabs(d-id);
    
    /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
     * for reporting a bug here.  */
    if( fabs(ax-1.0) < EPS1 )      /* |x| == 1.0 */
    {
        if( x > 0.0 )
        {
            if( flag & 12 ) /* negative int c-a or c-b */
            {
                //Rprintf("c-a or c-b negative");
                if( d >= 0.0 )
                    goto hypf;
                else
                    goto hypdiv;
            }
            if( d <= 0.0 )
            {
                goto hypdiv;}
            //y = exp(lgammafn(c)+lgammafn(d) -(lgammafn(p) + lgammafn(r)));
            y = gammafn(c)*gammafn(d)/(gammafn(p)*gammafn(r));
            goto hypdon;
        }
        
        if( d <= -1.0 )
            goto hypdiv;
    }
    /* Conditionally make d > 0 by recurrence on c
     * AMS55 #15.2.27
     */
    if( d < 0.0 )
    {
        /* Try the power series first */
        y = hyt2f1( a, b, c, x, &err );
        if( err < ETHRESH )
            goto hypdon;
        /* Apply the recurrence if power series fails */
        err = 0.0;
        aid = 2 - id;
        e = c + aid;
        d2 = hyp2f1(a,b,e,x);
        d1 = hyp2f1(a,b,e+1.0,x);
        q = a + b + 1.0;
        for( i=0; i<aid; i++ )
        {
            r = e - 1.0;
            y = (e*(r-(2.0*e-q)*x)*d2 + (e-a)*(e-b)*x*d1)/(e*r*s);
            e = r;
            d1 = d2;
            d2 = y;
        }
        goto hypdon;
    }
    
    
    if( flag & 12 )
        goto hypf; /* negative integer c-a or c-b */
    
hypok:
    y = hyt2f1( a, b, c, x, &err ); 
hypdon:
    if( err > ETHRESH )
    {
      //Rprintf( "hyp2f1, PLOSS\n");
       }
    return(y);
    /* The transformation for c-a or c-b negative integer
     * AMS55 #15.3.3
     */
hypf:
    y = pow( s, d ) * hys2f1( c-a, c-b, c, x, &err );
    goto hypdon;
    
    /* The alarm exit */
hypdiv:
    //Rprintf( "hyp2f1, OVERFLOW\n");
    return( MAXNUM );
}


/* Apply transformations for |x| near 1
 * then call the power series
 */
double hyt2f1( double a,double b,double c,double x,double *loss )
{
    double p, q, r, s, t, y, d, err, err1;
    double ax, id, d1, d2, e, y1;
    int i, aid;
    err = 0.0;
    s = 1.0 - x;
    if( x < -0.5 )
    {
        if( b > a )
            y = pow( s, -a ) * hys2f1( a, c-b, c, -x/s, &err );
        else
            y = pow( s, -b ) * hys2f1( c-a, b, c, -x/s, &err );
        goto done;
    }
    d = c - a - b;
    id = round(d);  /* nearest integer to d */
  
    if( x > 0.9 )
    {
        if( fabs(d-id) > EPS1 ) /* test for integer c-a-b */
        {
            //  Rprintf("integer case\n");
            /* Try the power series first */
            y = hys2f1( a, b, c, x, &err );
            if( err < ETHRESH ) {
                //Rprintf("Power series failed");
                goto done;
            }
            /* If power series fails, then apply AMS55 #15.3.6 */
            q = hys2f1( a, b, 1.0-d, s, &err );
            if (d < 0)
                q *= gammafn(d) /(gammafn(c-a) * gammafn(c-b));
            else q *= exp(lgammafn(d)  - (lgammafn(c-a) + lgammafn(c-b)));
            r = pow(s,d) * hys2f1( c-a, c-b, d+1.0, s, &err1 );
            if ( d > 0)
                r *= gammafn(-d) /gammafn(a) * gammafn(b);
            else  r *= exp(lgammafn(-d) - (lgammafn(a) + lgammafn(b)));
            y = q + r;
            q = fabs(q); /* estimate cancellation error */
            r = fabs(r);
            if( q > r )
                r = q;
            err += err1 + (MACHEP*r)/y;
            
            y *= gammafn(c);
            goto done;
        }
        else
        {
            /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12 */
            if( id >= 0.0 )
            {
                //    Rprintf("id >= 0\n");
                e = d;
                d1 = d;
                d2 = 0.0;
                aid = id;
            }
            else
            {
                e = -d;
                d1 = 0.0;
                d2 = d;
                aid = -id;
            }
            
            ax = log(s);
            
            /* sum for t = 0 */
            y = digamma(1.0) + digamma(1.0+e) - digamma(a+d1) - digamma(b+d1) - ax;
            y /= gammafn(e+1.0);
            
            p = (a+d1) * (b+d1) * s / gammafn(e+2.0); /* Poch for t=1 */
            t = 1.0;
            do
            {
                r = digamma(1.0+t) + digamma(1.0+t+e) - digamma(a+t+d1)
                - digamma(b+t+d1) - ax;
                q = p * r;
                y += q;
                p *= s * (a+t+d1) / (t+1.0);
                p *= (b+t+d1) / (t+1.0+e);
                t += 1.0;
            }
            while( fabs(q/y) > EPS1 );
            
            
            if( id == 0.0 )
            {
                y *= gammafn(c)/(gammafn(a)*gammafn(b));
                goto psidon;
            }
            
            y1 = 1.0;
            
            if( aid == 1 )
                goto nosum;
            
            //  Rprintf("sum case b=%lf\n", b);
            t = 0.0;
            p = 1.0;
            for( i=1; i<aid; i++ )
            {
                r = 1.0-e+t;
                p *= s * (a+t+d2) * (b+t+d2) / r;
                t += 1.0;
                p /= t;
                y1 += p;
            }
        nosum:
            p = gammafn(c);
            //  Rprintf("e = %lf, a=%lf, b=%lf, d1= %lf, d2=%lf\n ", e, a, b, d1, d2);
            //  y1 *= gammafn(e) * p / (gammafn(a+d1) * gammafn(b+d1));
            y1 *= exp(lgammafn(e) + log(p) - (lgammafn(a+d1) - lgammafn(b+d1)));
            //  y *= p / (gammafn(a+d2) * gammafn(b+d2));
            y *= exp(log(p) - lgammafn(a+d2+ .00001))/gammafn(b+d2 + .00001);
            if( (aid & 1) != 0 )
                y = -y;
            
            q = pow( s, id ); /* s to the id power */
            if( id > 0.0 )
                y *= q;
            else
                y1 *= q;
            
            y += y1;
        psidon:
            goto done;
        }   
    }
    /* Use defining power series if no special cases */
    y = hys2f1( a, b, c, x, &err );
done:
    *loss = err;
    return(y);
}


/* Defining power series expansion of Gauss hypergeometric function */

double hys2f1( double a,double b,double c,double x,double *loss )
//double *loss; /* estimates loss of significance */
{
    double f, g, h, k, m, s, u, umax;
    int i;
    
    i = 0;
    umax = 0.0;
    f = a;
    g = b;
    h = c;
    s = 1.0;
    u = 1.0;
    k = 0.0;
    do
    {
        if( fabs(h) < EPS1 )
        {
            *loss = 1.0;
            return( MAXNUM );
        }
        m = k + 1.0;
        u = u * ((f+k) * (g+k) * x / ((h+k) * m));
        s += u;
        k = fabs(u);  /* remember largest term summed */
        if( k > umax )
            umax = k;
        k = m;
        if( ++i > 10000 ) /* should never happen */
        {
            *loss = 1.0;
            return(s);
        }
    }
    while( fabs(u/s) > MACHEP );
    
    /* return estimated relative error */
    *loss = (MACHEP*umax)/fabs(s) + (MACHEP*i);
    
    return(s);
}

double hypergeo(double a,double b,double c,double x)
{
    double sol;
    sol =(hyp2f1( a, b, c, x ));
    return(sol);
}

void hypergeo_call(double *a,double *b,double *c,double *x, double *res)
{
    *res = hypergeo(*a,*b,*c,*x);
}
/***********************************************************/
/*********** bivariate T distribution********************/ 
double biv_T(double rho,double zi,double zj,double nuu)
{
  int k=0; 
  double nu=1/nuu;
  double res0=0.0,RR=0.0,pp1=0.0,pp2=0.0;
  double bb1,bb2;
  double x=zi;double y=zj;        
  double cc=(nu+1)/2; double nu2=nu/2;
  double x1=(x*x+nu); double y1=(y*y+nu);
  double rho2=R_pow(1-rho*rho,-cc);

  double b1 = (R_pow(nu,nu))*R_pow(x1*y1,-cc)*R_pow(gammafn(cc),2);
  double c1 = M_PI*R_pow(gammafn(nu/2),2)*rho2;
  double b2 = rho*x*y*R_pow(nu,nu+2)*R_pow(x1*y1,-nu2-1);
  double c2 = 2*M_PI*rho2;

  double a1 = 0; double a2 = 0;
  double aux  = R_pow(rho*x*y,2)/(x1*y1);
  double aux1 = R_pow(rho*nu,2)/(x1*y1);
 //if(fabs(rho)<=EPS1)
  /*if(!fabs(rho))
  {
    C = lgammafn(cc)+log(R_pow((1+x*x/nu),-cc))-log(sqrt(M_PI*nu))-lgammafn(nu/2);
    B = lgammafn(cc)+log(R_pow((1+y*y/nu),-cc))-log(sqrt(M_PI*nu))-lgammafn(nu/2);
    return(exp(B)*exp(C));
  }*/
  while( k<=6000 )
    {
   // pp1=hypergeo(cc+k,cc+k,0.5,aux);
    pp1=(0.5-2*(cc+k))*log(1-aux)+log(hypergeo(0.5-(cc+k),0.5-(cc+k),0.5,aux)); //euler
    bb1=pp1+k*log(aux1)+2*(lgammafn(cc+k)-lgammafn(cc))-lgammafn(k+1)-lgammafn(nu2+k)+lgammafn(nu2);
    a1 = a1 + exp(bb1);
   // pp2=hypergeo(nu2+1+k,nu2+1+k,1.5,aux);
    pp2=(1.5-2*(nu2+1+k))*log(1-aux)+log(hypergeo(1.5-(nu2+1+k),1.5-(nu2+1+k),1.5,aux));//euler
    bb2=pp2+k*log(aux1)+2*log((1+k/nu2))+lgammafn(nu2+k)-lgammafn(k+1)-lgammafn(nu2);
    a2 = a2 + exp(bb2);
    RR=(b1/c1)*a1+(b2/c2)*a2;
   if(!R_FINITE(RR)) return(res0);
    if((fabs(RR-res0)<1e-90)  ) {break;}
    else {res0=RR;}
        k++;
    }
return(RR);
}
/*********** Appell F4 function ********/
double appellF4(double a,double b,double c,double d,double x,double y)
{
double RR=0.0,bb=0.0,res0=0.0;int k=0;
if((int)c%2==0) c=c+0.0000001;
 for (k=0;k<=15000;k=k+1)
  {
    bb=k*log(y)+(lgammafn(a+k)+lgammafn(b+k)+lgammafn(d))
               -(lgammafn(a)+lgammafn(b)+lgammafn(d+k)+lgammafn(k+1))
               +(c-(a+k)-(b+k))*log(1-x)+log(hypergeo(c-a-k,c-b-k,c,x)); //euler
              // +log(hypergeo(a+k,b+k,c,x));
    RR=RR+exp(bb);
    if(!R_FINITE(RR)) return(res0);
 if(fabs(RR-res0)<=EPS) {break;}
    else {res0=RR;}
}
return(RR);
}
/****************************************/
double appellF4_mod(double nu,double rho2,double x,double y)
{
  double xx,yy,x2,y2,arg,arg1,pp1,pp2,app;
xx=x*x;yy=y*y;
x2=xx+nu;
y2=yy+nu;
arg=(nu+1)/2; 
arg1=nu/2;
pp1=R_pow(nu,nu)*R_pow(x2*y2,-arg)*R_pow(gammafn(arg),2);
pp2=M_PI*R_pow(gammafn(arg1),2)*R_pow(1-rho2,-arg);
app=appellF4(arg,arg,0.5,arg1,rho2*xx*yy/(x2*y2), nu*nu*rho2/(x2*y2));
return(4*pp1*app/pp2);
}
/*********** bivariate two piece-T distribution********************/ 
double biv_two_pieceT(double rho,double zi,double zj,double sill,double nuu,double eta,
             double p11,double mui,double muj)
{
double res;  
double nu=1/nuu;
double etamas=1+eta;
double etamos=1-eta;
double rho2=rho*rho;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);
if(zi>=mui&&zj>=muj)
{res=          (p11/R_pow(etamos,2))*appellF4_mod(nu,rho2,zistd/etamos,zjstd/etamos);}
if(zi>=mui&&zj<muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho2,zistd/etamos,zjstd/etamas);}
if(zi<mui&&zj>=muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho2,zistd/etamas,zjstd/etamos);}
if(zi<mui&&zj<muj)
{res=    ((p11+eta)/R_pow(etamas,2))*appellF4_mod(nu,rho2,zistd/etamas,zjstd/etamas);}
return(res/sill);
}
/***** bivariate half gaussian ****/     
double biv_half_Gauss(double rho,double zi,double zj)
{
double kk=0, dens=0,a=0,b=0,rho2=1-rho*rho;
  kk=(M_PI)*sqrt(rho2);
  a=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)-2*rho*zi*zj));
  b=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)+2*rho*zi*zj));
  dens=(a+b)/kk;
  return(dens);
}
/***** bivariate two piece gaussian ****/ 
double biv_two_pieceGaussian(double rho,double zi,double zj,double sill,double eta,
             double p11,double mui,double muj)
{
double res;  
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);
if(zi>=mui&&zj>=muj)
{res=          (p11/R_pow(etamos,2))*biv_half_Gauss(rho,zistd/etamos,zjstd/etamos);}
if(zi>=mui&&zj<muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Gauss(rho,zistd/etamos,zjstd/etamas);}
if(zi<mui&&zj>=muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Gauss(rho,zistd/etamas,zjstd/etamos);}
if(zi<mui&&zj<muj)
{res=    ((p11+eta)/R_pow(etamas,2))*biv_half_Gauss(rho,zistd/etamas,zjstd/etamas);}
return(res/sill);
}



double biv_Kumara(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2)
{
  double xx=0.0,yy=0.0,ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0;
 ki=1-pow(zi,shape2); kj=1-pow(zj,shape2);
//if(rho) {
  rho2=rho*rho;
  xx=rho2*pow(ki*kj,shape1);
  yy=rho2*(1-pow(ki,shape1))*(1-pow(kj,shape1));
  p1=pow(shape1*shape2,2)*pow(zi*zj,shape2-1)*pow(ki*kj,shape1-1)*pow(1-rho2,2);
    p2= appellF4(2,2,1,1,xx,yy);
  res=p1*p2;
/*} else  {p1=shape1*shape2*pow(zi,shape2-1)* pow(ki,shape1-1);
         p2=shape1*shape2*pow(zj,shape2-1)* pow(kj,shape1-1);
         res=p1*p2;}*/
return(res);
}

/***********************************************************************************/
/************ functions for binomial or negative binomial  two piece *********/
/***********************************************************************************/
/***********************************************************************************/
/************ functions for binomial or negative binomial  two piece *********/
/***********************************************************************************/

double pbnorm22(double lim1,double lim2,double corr,double nugget)
{
    double  lowe[2]  = {0,0}, uppe[2] = {lim1,lim2}, corre[1] = {(1-nugget)*corr};
    double  value;
    int     infin[2] = {0,0};
    value            = F77_CALL(bvnmvn)(lowe,uppe,infin,corre); 
    return(value);
}

// cdf bivariate half-normal distribution
double pbhalf_gauss(double zi,double zj,double rho,double nugget)
{
  double dens = 0;
  dens = pbnorm22(zi,zj,rho,nugget) + pbnorm22(-zi,-zj,rho,nugget) -
              pbnorm22(-zi,zj,rho,nugget) - pbnorm22(zi,-zj,rho,nugget);
  return(dens); 
}
/****** cdf univariate two-piece gaussian distribution *****/
double pnorm_two_piece(double x, double eta)
{
    double cdf = 0;
    if (x <=  0) cdf = (1 + eta)*pnorm(x/(1 + eta),0,1,1,0);
    if (x > 0)   cdf = eta + (1 - eta)*pnorm(x/(1 - eta),0,1,1,0);
  return(cdf);
}
//***********************************************//
// cdf univariate half-gaussian distribution
double phalf_gauss (double z){
  double dens = 0;
  dens = 2 * pnorm(z,0,1,1,0) - 1;
  return(dens);
}

// cdf bivariate two_piece gaussian distribution
double pbnorm_two_piece(int *cormod, double h, double u, 
    double xi, double xj, double nugget, double var,double eta,double *par)
{
    double p11,g_eta=0.0,q_g_eta=0.0,dens=0.0,corr=0.0;
    double etamos,etamas;
    g_eta   = (1 - eta)/2;
    q_g_eta = qnorm(g_eta,0,1,1,0);
    corr    = CorFct(cormod,h,u,par,0,0);
    p11     = pbnorm22(q_g_eta,q_g_eta,corr,nugget);
    etamas = 1+eta;
    etamos = 1-eta;
    if(xi  <= 0 & xj <= 0) dens = (1 + pbhalf_gauss(-xi/etamas,-xj/etamas,corr,nugget) - phalf_gauss(-xi/etamas) - phalf_gauss(-xj/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) );
    if(xi  > 0  & xj <= 0) dens = (1 + pbhalf_gauss(0,-xj/etamas,corr,nugget) - phalf_gauss(-xj/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + 
      (phalf_gauss(xi/etamos) - pbhalf_gauss(xi/etamos,-xj/etamas,corr,nugget)) * (pnorm(q_g_eta,0,1,1,0) - p11);
    if(xi  <= 0 & xj >  0) dens =  (1 + pbhalf_gauss(-xi/etamas,0,corr,nugget) - phalf_gauss(-xi/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + 
      (phalf_gauss(xj/etamos) - pbhalf_gauss(-xi/etamas,xj/etamos,corr,nugget)) * (pnorm(q_g_eta,0,1,1,0) - p11);
    if(xi  > 0  & xj >  0) dens =  ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + (phalf_gauss(xi/etamos) - pbhalf_gauss(xi/etamos,0,corr,nugget))*(pnorm(q_g_eta,0,1,1,0)-p11) + 
     (phalf_gauss(xj/etamos) - pbhalf_gauss(0,xj/etamos,corr,nugget))* (pnorm(q_g_eta,0,1,1,0)-p11) + pbhalf_gauss(xi/etamos,xj/etamos,corr,nugget)* p11;
    return(dens);
}
/***********************************************************************************/
/***********************************************************************************/
/********** bivariate log-logistic **********/
double biv_LogLogistic(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double c=gammafn(1+1/shape)*gammafn(1-1/shape);
    double A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=exp(mui);double cj=exp(muj);
    double ki=R_pow(c*zi/ci,shape)+1;
    double kj=R_pow(c*zj/cj,shape)+1;
    double rho2=R_pow(corr,2);
    double kij=ki*kj;
   // if(corr)   {
        A=(R_pow(c*shape,2)/(ci*cj))*R_pow((c*c*zi*zj)/(ci*cj),shape-1)*R_pow(1-rho2,2);
        B=R_pow(kij,-2);
        C=appellF4(2,2,1,1,
        (rho2*R_pow(c*c*zi*zj,shape))/(R_pow(ci*cj,shape)*kij),
        rho2/(kij));
        res=A*B*C;
    //}else{    B=(c*shape/ci)*R_pow((c*zi/ci),shape-1)*R_pow(ki,-2);  C=(c*shape/cj)*R_pow((c*zj/cj),shape-1)*R_pow(kj,-2);
    //    res=B*C;}
    return(res);
}

/********** bivariate logistic **********/
double biv_Logistic(double corr,double zi,double zj,double mui, double muj, double sill)
{
    double a=0.0,A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=mui;double cj=muj;
    double ki=exp((zi-ci)/sqrt(sill));
    double kj=exp((zj-cj)/sqrt(sill));
    double rho2=R_pow(corr,2);
    //if(corr)   {
        a=1-rho2;
        A=(ki*kj)/(R_pow(a,-2)*sill);
        B=R_pow((ki+1)*(kj+1),-2);
        C=appellF4(2,2,1,1,
        (rho2*ki*kj)/((ki+1)*(kj+1)),
        rho2/((ki+1)*(kj+1)));
        res=A*B*C;
    //} else{ B=ki*R_pow((ki+1),-2)/sqrt(sill);C=kj*R_pow((kj+1),-2)/sqrt(sill);res=B*C;}
    return(res);
}




//********** functions for bivariate tukey h *******************//

// compute lambert w function
double LambertW(double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
  //if (dbgW) fprintf(stderr,"LambertW: z=%g\n",z);
  if (z<-em1 || isinf(z) || isnan(z)) { 
    //fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); exit(1);
      return 0.0;
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return 
     -1.0
     +2.331643981597124203363536062168*r
     -1.812187885639363490240191647568*q
     +1.936631114492359755363277457668*r*q
     -2.353551201881614516821543561516*q2
     +3.066858901050631912893148922704*r*q2
     -4.175335600258177138854984177460*q3
     +5.858023729874774148815053846119*r*q3
     -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
    p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else 
    w=log(z); /* asymptotic */
  if (z>3.0) w-=log(w); /* useful? */
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  //fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z);
    return 0.0;
  //exit(1);
}

// pdf bivariate gaussian distribution
double dbnorm(double x_i,double x_j,double mean_i,double mean_j,double sill,double corr)
{
    double  fraq = 1.0,dens = 0.0,aux1 = 1.0,z1 = 1.0,z2 = 1.0,z3 = 1.0,z = 1.0;
    fraq = 2*M_PI*sill*sqrt(1-corr*corr);
    z1   = (x_i - mean_i)*(x_i - mean_i);
    z2   = (x_j - mean_j)*(x_j - mean_j);
    z3   = 2*corr*(x_i - mean_i)*(x_j - mean_j);
    z    = (z1 + z2 - z3)/sill;
    aux1 = 2*(1-corr*corr); 
    dens = (1/fraq)*exp(-z/aux1);
    return(dens);
}


// compute the inverse lambert w transformation
double inverse_lamb(double x,double tail)
{
  double sign,value;
  value = sqrt(LambertW(tail*x*x)/tail);
  if (x > 0) sign= 1;
  if (x < 0) sign= -1;
   return(sign*value);
}


// pdf bivariate tukey h random field 
double biv_tukey_h(double corr,double data_i, double data_j, double mean_i, double mean_j, double tail, double sill)
{
  double dens = 0.0,x_i = 0.0,x_j = 0.0,est_mean_i = 0.0,est_mean_j = 0.0;
  double est_mean_ij = 1.0,extra = 1.0;

  est_mean_i = (data_i - mean_i)/sqrt(sill);
  est_mean_j = (data_j - mean_j)/sqrt(sill);
  
  x_i = inverse_lamb(est_mean_i,tail);
  x_j = inverse_lamb(est_mean_j,tail);

  est_mean_ij = 1/(est_mean_i*est_mean_j);
    extra       = 1/( (1 + LambertW(tail*est_mean_i*est_mean_i))*(1 + LambertW(tail*est_mean_j*est_mean_j)));
  dens = dbnorm(x_i,x_j,0,0,1,corr)*
              x_i*x_j * est_mean_ij * extra/sill;
  return(dens);
}








/***** bivariate half tukey h ****/     
double biv_half_Tukeyh(double rho,double ti,double tj,double tail)
{
  double dens = 0.0;
  dens = biv_tukey_h(rho,ti,tj,0,0,tail,1) + biv_tukey_h(rho,-ti,-tj,0,0,tail,1) + biv_tukey_h(rho,-ti,tj,0,0,tail,1) + biv_tukey_h(rho,ti,-tj,0,0,tail,1);
  return(dens);
}
 

/***** bivariate two piece tukey h ****/ 
double biv_two_pieceTukeyh(double rho,double zi,double zj,double sill,double eta,double tail,
             double p11,double mui,double muj)
{
double res;  
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);

if(zi>=mui&&zj>=muj)
{res=          (p11/R_pow(etamos,2))*biv_half_Tukeyh(rho,zistd/etamos,zjstd/etamos,tail);}
if(zi>=mui&&zj<muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,zistd/etamos,zjstd/etamas,tail);}
if(zi<mui&&zj>=muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,zistd/etamas,zjstd/etamos,tail);}
if(zi<mui&&zj<muj)
{res=    ((p11+eta)/R_pow(etamas,2))*biv_half_Tukeyh(rho,zistd/etamas,zjstd/etamas,tail);}

return(res/sill);
}


