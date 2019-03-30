/*###################################################
### Authors: Moreno Bevilacqua.
### Email: moreno.bevilacqua@unibg.it
### File name: Gradient.c
### Description:
### This file contains a set of procedures
### for the computation of the composite likelihood
### gradients.
### Last change: 03/02/2017.
##################################################*/

#include "header.h"

// Compute the gradient vector of the conditional pairwise log likelihood for a
void Grad_Cond_Gauss(double rho, int *flag, double *gradcor, double *grad,
                     int *npar, double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0],nugget=par[1],sill=par[2];
  double a=nugget+sill,b=sill*rho,pa=a*a,pb=b*b;
  double c=-pa+pb,d=pa+pb,da=2*a,k=1/(c*c);
  double pn=nugget*nugget,ps=sill*sill;
  double C=0.0,L=0.0,R=0.0;
  double pu=0.0,pv=0.0,su=0.0,sv=0.0,suv=0.0;
  int h=0,i=0,j=0;
  //defines useful quantities:
  u=u-mean;
  v=v-mean;
  pu=u*u;
  pv=v*v;
  R=pu+pv;L=u*v;
  su=(-1+pu/a)/da;
  sv=(-1+pv/a)/da;
  suv=su+sv;
  // Derivatives of the conditional respect with the mean
  if(flag[0]==1){grad[i]=2*(u+v)/(a+b)-(u/a+v/a);i++;}
  // Derivative of the conditional respect with the nugget
   if(flag[1]==1){grad[i]=k*(R*d-L*4*b*a-2*a*(pa-pb))-suv;i++;}
  // Derivative of the conditional respect with the sill
   if(flag[2]==1){grad[i]=-k*(2*(pa*a-pb*(2*sill+3*nugget)+
				 rho*b*(pb-pn))+R*(c+2*nugget*b*rho)+
			      2*L*rho*(ps-pn-pb))-suv;i++;}
  // Derivatives with respect to the correlation parameters
  h=0;
  C=-2*k*sill*(R*a*b-L*d+b*c);
  for(j=i;j<*npar;j++){grad[j]=C*gradcor[h];h++;}
  return;
}
// Compute the gradient vector of the conditional pairwise log likelihood for a Gaussian model:
/*void Grad_Cond_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v)
{
  // Initialization variables:
  double mean=par[0], nugget=par[1], sill=par[2];
  double a=nugget+sill,b=sill*rho;
  double pa=a*a,pb=b*b,ppa=pa*a,ps=sill*sill,pr=rho*rho;
  double c=pa+pb,d=pa+pb,da=2*a,k=-c*-c;
  double pu=0.0, pv=0.0, su=0.0, sv=0.0, suv=0;
  double C=0.0,L=0.0,R=0.0;
  int h=0,i=0,j=0;
  //defines useful quantities:
  //defines useful quantities:
  u=u-mean;
  v=v-mean;
  pu=u*u;
  pv=v*v;
  R=pu+pv;
  L=u*v;
  su=(-1+pu/a)*pow(da,-1); //first statistics: first component
  sv=(-1+pv/a)*pow(da,-1); //second statistics: second component
  suv=su+sv;
  // Derivatives of the conditional respect with the mean
  if(flag[0]==1){grad[i]=2*((u+v)/(a+b))-(u/a+v/a);i++;}

  // Derivative of the conditional respect with the nugget
  if(flag[1]==1){grad[i]=k*(R*d-L*4*b*a-2*(ppa-pb*a))-suv;i++;}

  // Derivative of the conditional respect with the sill
  if(flag[2]==1){ grad[i]=-k*(R*(-pa+b*rho*(a+nugget))+
                           2*rho*L*(ps*(1-pr)-pow(nugget,2))+
                           2*(ppa-b*rho*(ps*(2-pr+nugget*
                           (nugget+3*sill)))))-suv;i++;}


  // Derivatives with respect to the correlation parameters
  h=0;
  C=-2*k*(R*a*b-L*d-b*c);
  for(j=i;j<*npar;j++) {grad[j]=C*gradcor[h];h++;}
  return;
  }*/


// Compute the gradient vector of the difference log likelihood for a Gaussian model :
void Grad_Diff_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v)
{
  // Initialization variables:
  double nugget=par[1], sill=par[2];
  double vario=0.0, sh=0.0;
  int h=0, i=0, j=0;
  //variogram:
  vario=nugget+sill*(1-rho);
  sh=0.5*(0.5*pow(u-v,2)/vario-1)/vario;
  // Derivative of the conditional respect with the nugget
  if(flag[1]==1) { grad[i]=sh; i++; }
  // Derivative of the conditional respect with the sill
  if(flag[2]==1) { grad[i]=(1-rho)*sh; i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++) { grad[j]=-sill*gradcor[h]*sh; h++; }
  return;
}




// Compute the gradient vector of the variogram for a Gaussian model :
void Grad_Diff_Vario(double rho, int *flag, double *gradcor,
		     double *grad, int *npar, double *par)
{
  // Initialization variables:
  double nugget=par[1], sill=par[2];
  double vario=0.0, sh=0.0;
  int h=0, i=0, j=0;
  //variogram:
  vario=nugget+sill*(1-rho);
  sh=1/vario;
  // Derivative of the conditional respect with the nugget
  if(flag[1]==1){ grad[i]=sh; i++; }
  // Derivative of the conditional respect with the sill
  if(flag[2]==1){ grad[i]=(1-rho)*sh; i++; }
  // Derivatives with respect to the correlation parameters
  for(j=i;j<*npar;j++){ grad[j]=-sill*gradcor[h]*sh; h++;}
  return;
}

/*************************************************************************************/
void Grad_Pair_Twopiecegauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,qqd=0.0,p11d=0.0,p11b=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 

  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double nugget=nuis[nbetas];
  double sill=nuis[nbetas+1];
  double skew=nuis[nbetas+2];

  double rho1=(1-nugget)*rho;
  double qq=qnorm((1-skew)/2,0,1,1,0);
  double p11=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,par,0);
  double ff=log(biv_two_pieceGaussian(rho1,u,v,sill,skew,p11,ai,aj));

   
/*
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*betas[kk];
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_two_pieceGaussian(rho1,u,v,sill,skew,p11,ai_d,aj_d)) - ff)/delta; 
   i++; }
}  
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
    p11b=pbnorm(cormod,lag,lagt,qq,qq,nugget+delta,1,par,0);
  grad[i]=(log(biv_two_pieceGaussian((1-(nugget+delta))*rho,u,v,sill,skew,p11b,ai,aj)) - 
           ff)/delta; 
    i++;}
  /* Derivvativve of the difference respect with the sill */
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_two_pieceGaussian(rho1,u,v,sill+delta,skew,p11,ai,aj)) - 
             ff)/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    qqd=qnorm((1-(skew+delta))/2,0,1,1,0);
    p11d=pbnorm(cormod,lag,lagt,qqd,qqd,nugget,1,par,0);
    grad[i]=(log(biv_two_pieceGaussian(rho1,u,v,sill,skew+delta,p11d,ai,aj)) - 
             ff)/delta; 
    i++; 
  }
  /* Derivvativves with respect to the correlation parameters*/
    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       p11d=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,parC,0);
      grad[kk+i]=(log(biv_two_pieceGaussian((1-nugget)*rhod,u,v,sill,skew,p11d,ai,aj)) - 
                  ff)/delta;
    kk++;}
      h++;
    }
  return;
}

/*************************************************************************************/
void Grad_Pair_Binomneg(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0,uu,vv,N;
  uu=(int) u; vv=(int) v; N=(int) NN;
  double p1=0.0,p2=0.0,p11=0.0,p1_d=0.0,p2_d=0.0,p11_d=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  double p11_dsill,p11_dcorr;
  //double sill=nuis[nbetas+1];
  double nugget=nuis[nbetas];

  /* probabilitities*/
  p11=pbnorm(cormod,lag,lagt,ai,aj,nugget,1,par,0);
  p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
for(o=0;o<nbetas;o++) {b1[o]=betas[o];}
for(o=0;o<nbetas;o++) { delta=sqrt(EPS)*betas[o];
                        if(flag[o]==1)  b1[o]=betas[o]+delta;}
double ff=log(biv_binomneg (N,uu,vv,p1,p2,p11));
/* Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
     p1_d=pnorm(ai_d,0,1,1,0); p2_d=pnorm(aj_d,0,1,1,0);
     p11_d=pbnorm(cormod,lag,lagt,ai_d,aj_d,nugget,1,par,0);
    
   grad[i]=(log(biv_binomneg (N,uu,vv,p1_d,p2_d,p11_d)) - ff)/delta; 
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { delta=sqrt(EPS)*(nugget);
    p11_dsill=pbnorm(cormod,lag,lagt,ai,aj,nugget+delta,1,par,0);
     grad[i]=(log(biv_binomneg(N,uu,vv,p1,p2,p11_dsill)) - ff)/delta;
    i++; }
  // Derivvativve of the difference respect with the sill
  //if(flag[nbetas+1]==1) { 
  //  delta=sqrt(EPS)*(sill);
  //  p11_dsill=pbnorm(cormod,lag,lagt,ai,aj,nugget,sill+delta,par,0);
   //  grad[i]=(log(biv_binomneg(N,uu,vv,p1,p2,p11_dsill)) - log(biv_binomneg(N,uu,vv,p1,p2,p11)))/delta;
   // i++; }
  /* Derivvativves with respect to the correlation parameters*/

                 h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*(par[h]);
       parC[h]=par[h]+delta;
       p11_dcorr=pbnorm(cormod,lag,lagt,ai,aj,nugget,1,parC,0);
     grad[kk+i]=(log(biv_binomneg (N,uu,vv,p1,p2,p11_dcorr)) - ff)/delta;
    kk++;}
      h++;
    }


  return;
}
/*************************************************************************************/
void Grad_Pair_Binom(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  int h=0, i=0, j=0,kk=0,o=0,k=0,uu,vv,N;
  uu=(int) u; vv=(int) v; N=(int) NN;
  double p1=0.0,p2=0.0,p11=0.0,p1_d=0.0,p2_d=0.0,p11_d=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  double p11_dsill,p11_dcorr;
  //double sill=nuis[nbetas+1];
  double nugget=nuis[nbetas];

  /* probabilitities*/
  p11=pbnorm(cormod,lag,lagt,ai,aj,nugget,1,par,0);
  p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
for(o=0;o<nbetas;o++) {b1[o]=betas[o];}
for(o=0;o<nbetas;o++) { delta=sqrt(EPS)*betas[o];
                        if(flag[o]==1)  b1[o]=betas[o]+delta;}
double ff=log(biv_binom (N,uu,vv,p1,p2,p11));
/* Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
     p1_d=pnorm(ai_d,0,1,1,0); p2_d=pnorm(aj_d,0,1,1,0);
     p11_d=pbnorm(cormod,lag,lagt,ai_d,aj_d,nugget,1,par,0);
    
   grad[i]=(log(biv_binom (N,uu,vv,p1_d,p2_d,p11_d)) - ff)/delta; 
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { delta=sqrt(EPS)*(nugget);
    p11_dsill=pbnorm(cormod,lag,lagt,ai,aj,nugget+delta,1,par,0);
     grad[i]=(log(biv_binom(N,uu,vv,p1,p2,p11_dsill)) - ff)/delta;
    i++; }
  // Derivvativve of the difference respect with the sill
  //if(flag[nbetas+1]==1) { 
  //  delta=sqrt(EPS)*(sill);
  //  p11_dsill=pbnorm(cormod,lag,lagt,ai,aj,nugget,sill+delta,par,0);
   //  grad[i]=(log(biv_binomneg(N,uu,vv,p1,p2,p11_dsill)) - log(biv_binomneg(N,uu,vv,p1,p2,p11)))/delta;
   // i++; }
  /* Derivvativves with respect to the correlation parameters*/

                 h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*(par[h]);
       parC[h]=par[h]+delta;
       p11_dcorr=pbnorm(cormod,lag,lagt,ai,aj,nugget,1,parC,0);
     grad[kk+i]=(log(biv_binom (N,uu,vv,p1,p2,p11_dcorr)) - ff)/delta;
    kk++;}
      h++;
    }


  return;
}

void Grad_Pair_Weibull(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];
  double rhosill=(1-nugget)*rho;
  double ff=log(biv_Weibull(rhosill,u,v,ai,aj,shape));
  //Rprintf("----ff---%f\n",ff);
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_Weibull(rhosill,u,v,ai_d,aj_d,shape)) - ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=(log(biv_Weibull((1-(nugget+delta))*rho,u,v,ai,aj,shape)) - 
              ff)/delta; 
     /// Rprintf("----grad2---%f\n",grad[i]);
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/
 // if(flag[nbetas+1]==1) { 
 //   delta=sqrt(EPS)*sill;
 //     grad[i]=(log(biv_Weibull((sill+delta)*rho,u,v,ai,aj,shape)) - 
 //             ff)/delta; 
 //   i++; 
 // }
   /* Derivvativve of the difference respect with the shape*/   
   if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*shape;
     //delta=R_pow(EPS,1/3);
    grad[i]=(log(biv_Weibull(rhosill,u,v,ai,aj,shape+delta)) -
                 ff)/(delta); 
    //Rprintf("%f %f ----grad3---%f\n",u,v,grad[i]);
    i++; 
  }
    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
 grad[kk+i]=(log(biv_Weibull((1-nugget)*rhod,u,v,ai,aj,shape)) - ff)/delta;
    kk++;}
      h++;
    }
//Rprintf("%f %f %f %f %f %f %f \n",grad[0],grad[1],grad[2],grad[3],grad[4],grad[5],grad[6]);
  return;
}


void Grad_Pair_Gamma(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];
  double rhosill=(1-nugget)*rho;
  double ff=log(biv_gamma(rhosill,u,v,ai,aj,shape));
  //Rprintf("----ff---%f\n",ff);
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_gamma(rhosill,u,v,ai_d,aj_d,shape)) - ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=(log(biv_gamma((1-(nugget+delta))*rho,u,v,ai,aj,shape)) - 
              ff)/delta; 
     /// Rprintf("----grad2---%f\n",grad[i]);
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/
 // if(flag[nbetas+1]==1) { 
 //   delta=sqrt(EPS)*sill;
 //     grad[i]=(log(biv_Weibull((sill+delta)*rho,u,v,ai,aj,shape)) - 
 //             ff)/delta; 
 //   i++; 
 // }
   /* Derivvativve of the difference respect with the shape*/   
   if(flag[nbetas+2]==1) { 
    delta=R_pow(EPS,1/2)*shape;
     //delta=R_pow(EPS,1/3);
    grad[i]=(log(biv_gamma(rhosill,u,v,ai,aj,shape+delta)) -
                 ff)/(delta); 
    ///Rprintf("----grad3---%f\n",grad[i]);
    i++; 
  }
    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
 grad[kk+i]=(log(biv_gamma((1-nugget)*rhod,u,v,ai,aj,shape)) - ff)/delta;
    kk++;}
      h++;
    }
//Rprintf("%f %f %f %f %f %f %f \n",grad[0],grad[1],grad[2],grad[3],grad[4],grad[5],grad[6]);
  return;
}


void Grad_Pair_LogGauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
       b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  double nugget=nuis[nbetas];
 double sill=1-nugget;
 double ff=log(d2lognorm(u,v,sill,nugget, ai, aj,rho));

/*
  // Derivativve of the difference respect with the mean*/
   
     for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(d2lognorm(u,v,sill,nugget, ai_d, aj_d,rho)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=(log(d2lognorm(u,v,sill,nugget+delta, ai, aj,rho)) - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(d2lognorm(u,v,sill+delta,nugget, ai, aj,rho)) - ff)/delta; 
    i++; 
  }
  /* Derivvativves with respect to the correlation parameters*/

                   h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
        grad[kk+i]=(log(d2lognorm(u,v,sill,nugget, ai, aj,rhod)) - ff)/delta;
    kk++;}
      h++;
    }
  return;
}

void Grad_Pair_Logistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
 

b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 

 double sill=nuis[nbetas+1];

double ff=log(biv_Logistic(rho,u,v,ai,aj,sill));
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
                          
   grad[i]=(log(biv_Logistic(rho,u,v,ai_d,aj_d,sill)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { grad[i]=1; i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_Logistic(rho,u,v,ai,aj,sill+delta)) - ff)/delta; 
    i++; 
  }
  /* Derivvativves with respect to the correlation parameters*/

                  h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
grad[kk+i]=(log(biv_Logistic(rhod,u,v,ai,aj,sill)) - ff)/delta;
    kk++;}
      h++;
    }

  return;
}


/*************************************************************************************/
void Grad_Pair_TwopieceT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,qqd=0.0,p11d=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
    
    double df=nuis[nbetas];
    double nugget=nuis[nbetas+1];
    double sill=nuis[nbetas+2];
    double skew=nuis[nbetas+3];

    double rho1=(1-nugget)*rho;
    double qq=qnorm((1-skew)/2,0,1,1,0);
    double p11=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,par,0);
    double ff=log(biv_two_pieceT(rho1,u,v,sill,df,skew,p11,ai,aj));
 Rprintf("%d %d %d %d %d %d %f %f %f %f\n",nbetas,flag[0],flag[1],flag[2],flag[3],flag[4],df,nugget,sill,skew);
/*
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*betas[kk];
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_two_pieceT(rho1,u,v,sill,df,skew,p11,ai_d,aj_d)) - ff)/delta; 
   i++; }
}  
// Derivvativve of the difference respect with df*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*df; 
    grad[i]=(log(biv_two_pieceT(rho1,u,v,sill,df+delta,skew,p11,ai,aj)) - ff)/delta;
    i++; 
  }
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas+1]==1) {   
    delta=sqrt(EPS)*nugget; 
    p11d=pbnorm(cormod,lag,lagt,qq,qq,nugget+delta,1,par,0);
    grad[i]=(log(biv_two_pieceT((1-(nugget+delta))*rho,u,v,sill,df,skew,p11d,ai,aj)) - 
             ff)/delta;
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_two_pieceT(rho1,u,v,sill+delta,df,skew,p11,ai,aj)) - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with skew*/  
    if(flag[nbetas+3]==1) { 
    delta=sqrt(EPS)*skew;
    qqd=qnorm((1-(skew+delta))/2,0,1,1,0);
    p11d=pbnorm(cormod,lag,lagt,qqd,qqd,nugget,1,par,0);
    grad[i]=(log(biv_two_pieceT(rho,u,v,sill,df,skew+delta,p11d,ai,aj)) - ff)/delta; 
    i++; 
  }
  /* Derivvativves with respect to the correlation parameters*/
        h=0;kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       p11d=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,parC,0);
      grad[kk+i]=(log(biv_two_pieceT((1-nugget)*rhod,u,v,sill,df,skew,p11d,ai,aj)) - 
                  ff)/delta;
    kk++;}
      h++;
    }
  return;
}


void Grad_Pair_StudenT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
 b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 

 double df=nuis[nbetas];
 double nugget=nuis[nbetas+1];
 double sill=nuis[nbetas+2];
 double rho1=(1-nugget)*rho;
 double ff=log(biv_T(rho1,u,v,ai,aj,df,sill));
/*
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*betas[kk];
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_T(rho1,u,v,ai_d,aj_d,df,sill)) - ff)/delta; 
   i++; }
}
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*df; 
    grad[i]=(log(biv_T(rho1,u,v,ai,aj,df+delta,sill))- ff)/delta;
    i++; 
  }
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas+1]==1) { delta=sqrt(EPS)*nugget; 
      grad[i]=(log(biv_T((1-(nugget+delta))*rho,u,v,ai,aj,df,sill)) - ff)/delta;
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*sill; 
    grad[i]=(log(biv_T(rho1,u,v,ai,aj,df,sill+delta)) -ff )/delta;
    i++; 
  }
  /* Derivvativves with respect to the correlation parameters*/
              h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
           grad[kk+i]=(log(biv_T((1-nugget)*rhod,u,v,ai,aj,df,sill)) 
            - ff)/delta;
    kk++;}
      h++;
    }

  return;
}


void Grad_Pair_Sinh(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 



  double sill=nuis[nbetas+1];//double nugget=nuis[nbetas]; nuis[2],nuis[3],nuis[1]
  double skew=nuis[nbetas+2];
  double tail=nuis[nbetas+3];
  double ff=log(biv_sinh(rho,u,v,ai,aj,skew,tail,sill));
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*betas[kk];
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_sinh(rho,u,v,ai_d,aj_d,skew,tail,sill)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { grad[i]=1; i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_sinh(rho,u,v,ai,aj,skew,tail,sill+delta)) - ff)/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    grad[i]=(log(biv_sinh(rho,u,v,ai,aj,skew+delta,tail,sill)) - ff)/delta; 
    i++; 
  }
  if(flag[nbetas+3]==1) { 
    delta=sqrt(EPS)*tail;
    grad[i]=(log(biv_sinh(rho,u,v,ai,aj,skew,tail+delta,sill)) - ff)/delta; 
    i++; 
  }

              h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
           grad[kk+i]=(log(biv_sinh(rhod,u,v,ai,aj,skew,tail,sill)) 
            - ff)/delta;
    kk++;}
      h++;
    }




  return;
}
void Grad_Pair_LogLogistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];
  double rhosill=(1-nugget)*rho;
  double ff=log(biv_LogLogistic(rhosill,u,v,ai,aj,shape));
  //Rprintf("----ff---%f\n",ff);
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_LogLogistic(rhosill,u,v,ai_d,aj_d,shape)) - ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=(log(biv_LogLogistic((1-(nugget+delta))*rho,u,v,ai,aj,shape)) - 
              ff)/delta; 
     /// Rprintf("----grad2---%f\n",grad[i]);
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/
 // if(flag[nbetas+1]==1) { 
 //   delta=sqrt(EPS)*sill;
 //     grad[i]=(log(biv_Weibull((sill+delta)*rho,u,v,ai,aj,shape)) - 
 //             ff)/delta; 
 //   i++; 
 // }
   /* Derivvativve of the difference respect with the shape*/   
   if(flag[nbetas+2]==1) { 
    delta=R_pow(EPS,1/2)*shape;
     //delta=R_pow(EPS,1/3);
    grad[i]=(log(biv_LogLogistic(rhosill,u,v,ai,aj,shape+delta)) -
                 ff)/(delta); 
    ///Rprintf("----grad3---%f\n",grad[i]);
    i++; 
  }
    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
 grad[kk+i]=(log(biv_LogLogistic((1-nugget)*rhod,u,v,ai,aj,shape)) - ff)/delta;
    kk++;}
      h++;
    }
//Rprintf("%f %f %f %f %f %f %f \n",grad[0],grad[1],grad[2],grad[3],grad[4],grad[5],grad[6]);
  return;
}


void Grad_Pair_Skewgauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 

  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double sill=nuis[nbetas+1];double nugget=nuis[nbetas];
  double skew=nuis[nbetas+2];
/*
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*betas[kk];
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}
   grad[i]=(log(biv_skew((1-nugget)*rho,u,v,ai_d,aj_d,sill,skew)) - log(biv_skew((1-nugget)*rho,u,v,ai,aj,sill,skew)))/delta; 
   i++; }
}  
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=(log(biv_skew((1-(nugget+delta))*rho,u,v,ai,aj,sill+delta,skew)) - log(biv_skew((1-nugget)*rho,u,v,ai,aj,sill,skew)))/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_skew((1-nugget)*rho,u,v,ai,aj,sill+delta,skew)) - log(biv_skew((1-nugget)*rho,u,v,ai,aj,sill,skew)))/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    grad[i]=(log(biv_skew((1-nugget)*rho,u,v,ai,aj,sill,skew+delta)) - log(biv_skew((1-nugget)*rho,u,v,ai,aj,sill,skew)))/delta; 
    i++; 
  }
  /* Derivvativves with respect to the correlation parameters*/
        h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
      grad[kk+i]=(log(biv_skew(rhod,u,v,ai,aj,sill,skew)) - log(biv_skew(rho,u,v,ai,aj,sill,skew)))/delta;
    kk++;}
      h++;
    }
  return;
}











void Grad_Pair_Wrapped(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,alfa=2.0;
  double delta=0,*b1,*parC;
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 

  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double sill=nuis[nbetas+1];double nugget=nuis[nbetas];
 double ff=log(biv_wrapped(alfa,u,v,ai,aj,nugget,sill,rho));

/*
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
     delta=sqrt(EPS)*betas[kk];
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){ai_d=ai_d+sX[l][o]*(b1[o]);
                           aj_d=aj_d+sX[m][o]*(b1[o]);}

   grad[i]=(log(biv_wrapped(alfa,u,v,ai_d,aj_d,nugget,sill,rho)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { grad[i]=1; i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_wrapped(alfa,u,v,ai,aj,nugget,sill+delta,rho)) - ff)/delta; 
    i++; 
  }
  /* Derivvativves with respect to the correlation parameters*/

      h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
       grad[kk+i]=(log(biv_wrapped(alfa,u,v,ai,aj,nugget,sill,rhod)) - ff)/delta;
    kk++;}
      h++;
    }


  return;
}

/******************************/
void Grad_Pair_Gauss(double rho, int *flag,int *flagcor,  double *gradcor, double *grad,
                     int *npar,int *nparc, int nbetas, double *par, double u, double v, double *Xl,double *Xm)
{
  // Initialization variables:
  double nugget=par[nbetas],sill=par[nbetas+1];
  double a=nugget+sill,b=sill*rho,pa=a*a,pb=b*b;
  double c=-pa+pb,d=pa+pb,k=1/(c*c);
  double C=0.0,L=0.0,R=0.0;
  double pn=nugget*nugget,ps=sill*sill,pu=0.0, pv=0.0;
  int kk=0,h=0, i=0, j=0;
  //defines useful quantities:
  pu=pow(u,2); pv=pow(v,2);R=pu+pv;L=u*v;
  // Derivatives  respect with the mean
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
   grad[i]=(((u*a - v*b))*Xl[kk]    + ((v*a - u*b))*Xm[kk])/(c);
     i++;}
  }
  // Derivative  respect with the nugget
   if(flag[nbetas]==1){grad[i]=0.5*k*(R*d-L*4*b*a-2*a*(pa-pb));i++;}
  // Derivative respect with the sill
   if(flag[nbetas+1]==1){grad[i]=-0.5*k*(2*(pa*a-pb*(2*sill+3*nugget)+
             rho*b*(pb-pn))+R*(c+2*nugget*b*rho)+
          2*L*rho*(ps-pn-pb));i++;}
  // Derivatives with respect to the correlation parameters
  h=0;
  C=-k*sill*(R*a*b-L*d+b*c);
  for(j=i;j<*npar;j++){;
    grad[j]=C*gradcor[h];h++;}
  
  return;
}

void Grad_Pair_Gauss2(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas)
{
  

    double rhod,delta=0,*parC,zi,zj;int j=0;
  parC=(double *) Calloc(nparcT[0],double);
  for(j=0;j<nparcT[0];j++) parC[j]=par[j];
   zi=u-ai;zj=v-aj;
double nugget=nuis[nbetas],sill=nuis[nbetas+1];

  double a=nugget+sill,b=sill*rho,pa=a*a,pb=b*b;
  double c=-pa+pb,d=pa+pb,k=1/(c*c);
  double L=0.0,R=0.0;
  double pn=nugget*nugget,ps=sill*sill,pu=0.0, pv=0.0;
  int kk=0,h=0, i=0;
  //defines useful quantities:
  pu=pow(zi,2); pv=pow(zj,2);R=pu+pv;L=zi*zj;
  // Derivatives  respect with the mean
  for(kk=0;kk<nbetas;kk++){
  if(flag[kk]==1){
   grad[i]=(((zi*a - zj*b))*Xl[kk] + ((zj*a - zi*b))*Xm[kk])/(c);
     i++;}
  }
  // Derivative  respect with the nugget
   if(flag[nbetas]==1){grad[i]=0.5*k*(R*d-L*4*b*a-2*a*(pa-pb));i++;}
  // Derivative respect with the sill
   if(flag[nbetas+1]==1){grad[i]=-0.5*k*(2*(pa*a-pb*(2*sill+3*nugget)+
             rho*b*(pb-pn))+R*(c+2*nugget*b*rho)+
          2*L*rho*(ps-pn-pb));i++;}
  // Derivatives with respect to the correlation parameters
  h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
    grad[kk+i]=(log_biv_Norm(rhod,u,v, ai, aj,sill,nugget) - log_biv_Norm(rho,u,v, ai, aj,sill,nugget))/delta;
    kk++;}
      h++;
    }


//Rprintf("%f %f %f %f %f %f %f %f %f %f %f  \n",grad[0],grad[1],grad[2],grad[3],zi,zj,a,b,c,nugget,sill);

  return;
}

void Grad_Pair_Gauss_biv(double rhott, double rhotv,double rhovt,double rhovv,int *flag,
                         double *gradcortt, double *gradcortv, double *gradcorvt, double *gradcorvv,
                         double *grad,int *npar, double *par, double u, double w)
{
    // Initialization variables:
    double a,b,c,def,mean1,mean2,det;
    mean1=par[0];mean2=par[1];
    int h=0, i=0, j=0;
    //defines useful quantities:
    u=u-mean1;w=w-mean2;
    det=rhott*rhovv-rhotv*rhovt;
    // Derivatives  respect with the mean1
    if(flag[0]==1)
    {  grad[i]=(-2*rhovv*u+rhotv*u+rhovt*w)/det;i++;}
    // Derivatives  respect with the mean2
    if(flag[1]==1){
      grad[i]=(-2*rhott*w+rhovt*u+rhotv*u)/det;i++;}
    // Derivatives with respect to the correlation parameters
    h=0;
    for(j=i;j<*npar;j++){
    def=rhott*gradcorvv[h]-rhotv*gradcorvt[h]-rhovt*gradcortv[h]+rhovv*gradcortt[h];
    a=(u*(u*gradcorvv[h]-w*gradcorvt[h]) + w*(w*gradcortt[h]-u*gradcortv[h]));
    b=(def)*(u*(rhovv*u-rhovt*w))/det;
    c=(def)*(w*(rhott*w-rhotv*u))/det;
    grad[j]=-0.5*(a-b-c+def)/det;h++;}
    return;
}



void Grad_Cond_Gauss_biv(double *gradcorttii,double *gradcorvvii,double *gradcorvtii,double *gradcortvii ,
                        double *gradcorttij ,double *gradcorvvij, double *gradcorvtij,double *gradcortvij,
                        double **inverse,int *flag,double *grad,int *npar, double *par,double *dat,int N)
{
    // Initialization variables:
    int i=0,j=0,h,m=0,n=0;
    double mean1,mean2;    
    double **D;
    D= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){D[i]=(double *) Calloc(N,double);}
    
    double **P1;
    P1= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){P1[i]=(double *) Calloc(N,double);}
    
    
    double **P2;
    P2= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){P2[i]=(double *) Calloc(N,double);}
    for(m=0;m<N;m++){
        for(n=0;n<N;n++) { P2[m][n]=0; P1[m][n]=0;}}
    mean1=par[0];mean2=par[1];
    //defines useful quantities:
    dat[0]=dat[0]-mean1;dat[1]=dat[1]-mean1;
    dat[2]=dat[2]-mean2;dat[3]=dat[3]-mean2;
    // Derivatives  respect with the mean1
    if(flag[0]==1)
    {  grad[i]=2   ;i++;}
    // Derivatives  respect with the mean2
    if(flag[1]==1){
        grad[i]=3;i++;}
    
    // Derivatives with respect to the correlation parameters
    h=0;
    for(j=i;j<*npar;j++){
        D[0][0]=gradcorttii[h];    D[0][1]=gradcorttij[h];       D[0][2]=gradcortvii[h];        D[0][3]=gradcortvij[h];
        D[1][0]=D[0][1];           D[1][1]=gradcorttii[h];       D[1][2]=gradcorvtij[h];        D[1][3]=gradcorvtii[h];
        D[2][0]=D[0][2];           D[2][1]= D[1][2];             D[2][2]=gradcorvvii[h];        D[2][3]=gradcorvvij[h];
        D[3][0]=D[0][3];           D[3][1]= D[1][3];             D[3][2]=D[2][3];               D[3][3]=gradcorvvii[h];
        Matrix_prod(inverse,D,P1,N);
        Matrix_prod(P1,inverse,P2,N);
        grad[j]=-0.5*(Trace(P1,N)-QFORM(P2,dat,dat,N));
        h++;}
    return;
    Free(D);Free(P1);Free(P2);
}





double mij(double qij, double w, double pij, double p)
{
  double val=0.0;
  val=(2*pij*w*(p-1)+qij*(p+pij-2*pow(p,2)))/(pij*(p-pij)*(qij-2*w));
  return(val);
}


double nij(double qij, double w, double pij, double p)
{
  double val=0.0;
  val=(w*(1-pij)-qij*(1-p))/((p-pij)*(qij-2*w));
  return(val);
}
