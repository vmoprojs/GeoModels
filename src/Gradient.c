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






/**************************************************************************************/
/**************************************************************************************/
/****** gradient for pairwise and conditional likelihood  *****************************/
/**************************************************************************************/
/**************************************************************************************/


void Grad_Pair_Gauss2(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
    double rhod,delta=0,*parC,*b1;int j=0,o=0,k=0;


        b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  


double nugget=nuis[nbetas],sill=nuis[nbetas+1];
  double ai_d=0.0,aj_d=0.0;
  int kk=0,h=0, i=0;
 

  double ff=log_biv_Norm((1-nugget)*rho,u,v, ai, aj,sill,0);


   for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log_biv_Norm((1-nugget)*rho,u,v, ai_d, aj_d,sill,0)-ff)/delta; 
      i++; 
   }

} 

  // Derivative  respect with the nugget
   if(flag[nbetas]==1){
      delta=sqrt(EPS)*nugget;
      grad[i]=(log_biv_Norm((1-(nugget+delta))*rho,u,v, ai, aj,sill,0) -ff)/(delta);
      i++;}
  // Derivative respect with the sill
   if(flag[nbetas+1]==1){
        delta=sqrt(EPS)*sill;
         grad[i]=(log_biv_Norm((1-nugget)*rho,u,v, ai, aj,sill+delta,0) -  ff)/(delta);i++;}
  // Derivatives with respect to the correlation parameters
  h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
    grad[kk+i]=(log_biv_Norm((1-nugget)*rhod,u,v, ai, aj,sill,0) - ff)/delta;
    kk++;}
      h++;
    }

  return;
}

void Grad_Cond_Gauss2(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
    double rhod,delta=0,*parC,*b1;int j=0,o=0,k=0;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

double nugget=nuis[nbetas],sill=nuis[nbetas+1];
  double ai_d=0.0,aj_d=0.0;
  int kk=0,h=0, i=0;
 
  double l1=dnorm(u,ai,sqrt(sill),1);
  double l2=dnorm(v,aj,sqrt(sill),1);
  double ff=2*log_biv_Norm((1-nugget)*rho,u,v, ai, aj,sill,0)-(l1+l2);

   for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(2*log_biv_Norm((1-nugget)*rho,u,v, ai_d, aj_d,sill,0)-
              (dnorm(u,ai_d,sqrt(sill),1)+dnorm(v,aj_d,sqrt(sill),1))
                - ff)/delta; 
   i++; }
} 

  // Derivative  respect with the nugget
   if(flag[nbetas]==1){
      delta=sqrt(EPS)*nugget;
      grad[i]=(2*log_biv_Norm((1-(nugget+delta))*rho,u,v, ai, aj,sill,0)-(l1+l2) -ff)/(delta);
      i++;}
  // Derivative respect with the sill
   if(flag[nbetas+1]==1){
        delta=sqrt(EPS)*sill;
         grad[i]=(2*log_biv_Norm((1-nugget)*rho,u,v, ai, aj,sill+delta,0) - 
                (dnorm(u,ai,sqrt(sill+delta),1)+dnorm(v,aj,sqrt(sill+delta),1))
                  -  ff)/(delta);i++;}
  // Derivatives with respect to the correlation parameters
  h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
    grad[kk+i]=(2*log_biv_Norm((1-nugget)*rhod,u,v, ai, aj,sill,0)-(l1+l2) - ff)/delta;
    kk++;}
      h++;
    }

  return;
}

/*************************************************************************************/
void Grad_Cond_Twopiecegauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,qqd=0.0,p11d=0.0,p11b=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];
  double sill=nuis[nbetas+1];
  double skew=nuis[nbetas+2];

  double rho1=(1-nugget)*rho;
  double qq=qnorm((1-skew)/2,0,1,1,0);
  double p11=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,par,0);
  double l1=one_log_two_pieceGauss(u,ai,sill,skew);
  double l2=one_log_two_pieceGauss(v,aj,sill,skew);
  double ff=2*log(biv_two_pieceGaussian(rho1,u,v,sill,skew,p11,ai,aj))-(l1+l2);

  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log(biv_two_pieceGaussian(rho1,u,v,sill,skew,p11,ai_d,aj_d)) - 
   (one_log_two_pieceGauss(u,ai_d,sill,skew)+one_log_two_pieceGauss(v,aj_d,sill,skew)))-ff)/delta; 
   i++; }
}  
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
    p11b=pbnorm(cormod,lag,lagt,qq,qq,nugget+delta,1,par,0);
  grad[i]=((2*log(biv_two_pieceGaussian((1-(nugget+delta))*rho,u,v,sill,skew,p11b,ai,aj)) - (l1+l2))-ff)/delta; 
    i++;}
  /* Derivvativve of the difference respect with the sill */
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=((2*log(biv_two_pieceGaussian(rho1,u,v,sill+delta,skew,p11,ai,aj)) - 
      one_log_two_pieceGauss(u,ai,sill+delta,skew)+one_log_two_pieceGauss(v,aj,sill+delta,skew))-ff)/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    qqd=qnorm((1-(skew+delta))/2,0,1,1,0);
    p11d=pbnorm(cormod,lag,lagt,qqd,qqd,nugget,1,par,0);
    grad[i]=((2*log(biv_two_pieceGaussian(rho1,u,v,sill,skew+delta,p11d,ai,aj)) - 
              one_log_two_pieceGauss(u,ai,sill,skew+delta)+one_log_two_pieceGauss(v,aj,sill,skew+delta))-ff)/delta; 
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
      grad[kk+i]=((2*log(biv_two_pieceGaussian((1-nugget)*rhod,u,v,sill,skew,p11d,ai,aj)) - (l1+l2))- ff)/delta;
    kk++;}
      h++;
    }
  return;
}

/*************************************************************************************/
void Grad_Pair_Twopiecegauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,qqd=0.0,p11d=0.0,p11b=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];


  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double nugget=nuis[nbetas];
  double sill=nuis[nbetas+1];
  double skew=nuis[nbetas+2];

  double rho1=(1-nugget)*rho;
  double qq=qnorm((1-skew)/2,0,1,1,0);
  double p11=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,par,0);
  double ff=log(biv_two_pieceGaussian(rho1,u,v,sill,skew,p11,ai,aj));

  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
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
void Grad_Cond_Binomneg(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
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

  double  l1=dbinom(uu,NN,p1,1);
  double  l2=dbinom(vv,NN,p2,1);
  double ff=2*log(biv_binomneg (N,uu,vv,p1,p2,p11))-(l1+l2);
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
     p1_d=pnorm(ai_d,0,1,1,0); p2_d=pnorm(aj_d,0,1,1,0);
     p11_d=pbnorm(cormod,lag,lagt,ai_d,aj_d,nugget,1,par,0);
    
   grad[i]=((2*log(biv_binomneg (N,uu,vv,p1_d,p2_d,p11_d)) -(dbinom(uu,NN,p1_d,1)+dbinom(vv,NN,p2_d,1)))-ff)/delta; 
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { delta=sqrt(EPS)*(nugget);
    p11_dsill=pbnorm(cormod,lag,lagt,ai,aj,nugget+delta,1,par,0);
     grad[i]=((2*log(biv_binomneg(N,uu,vv,p1,p2,p11_dsill))-(l1+l2)) - ff)/delta;
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
     grad[kk+i]=((2*log(biv_binomneg (N,uu,vv,p1,p2,p11_dcorr))-(l1+l2)) - ff)/delta;
    kk++;}
      h++;
    }
  return;
}

/*************************************************************************************/
void Grad_Pair_Binomneg(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
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
double ff=log(biv_binomneg (N,uu,vv,p1,p2,p11));
/* Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
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
void Grad_Cond_Binom(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
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


    double  l1=dbinom(uu,NN,p1,1);
    double  l2=dbinom(vv,NN,p2,1);
double ff=2*log(biv_binom (N,uu,vv,p1,p2,p11))-(l1+l2);
/* Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
     p1_d=pnorm(ai_d,0,1,1,0); p2_d=pnorm(aj_d,0,1,1,0);
     p11_d=pbnorm(cormod,lag,lagt,ai_d,aj_d,nugget,1,par,0);
    
   grad[i]=((2*log(biv_binom (N,uu,vv,p1_d,p2_d,p11_d)) - (dbinom(uu,NN,p1_d,1)+dbinom(vv,NN,p2_d,1)))-ff)/delta; 
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { delta=sqrt(EPS)*(nugget);
    p11_dsill=pbnorm(cormod,lag,lagt,ai,aj,nugget+delta,1,par,0);
     grad[i]=((2*log(biv_binom(N,uu,vv,p1,p2,p11_dsill)) -(l1+l2))- ff)/delta;
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
     grad[kk+i]=((2*log(biv_binom (N,uu,vv,p1,p2,p11_dcorr)) -(l1+l2))- ff)/delta;
    kk++;}
      h++;
    }
  return;
}

/*************************************************************************************/
void Grad_Pair_Binom(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
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

double ff=log(biv_binom (N,uu,vv,p1,p2,p11));
/* Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
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

void Grad_Cond_Weibull(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];

  double l1=one_log_weibull(u,ai,shape);
  double l2=one_log_weibull(v,aj,shape);
  double ff=2*log(biv_Weibull((1-nugget)*rho,u,v,ai,aj,shape))-(l1+l2);

for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=( (2*log(biv_Weibull((1-nugget)*rho,u,v,ai_d,aj_d,shape)) -(one_log_weibull(u,ai_d,shape)+one_log_weibull(v,aj_d,shape)))- ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=((2*log(biv_Weibull((1-(nugget+delta))*rho,u,v,ai,aj,shape))-(l1+l2)) - ff)/delta; 
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
grad[i]=(2*log(biv_Weibull((1-nugget)*rho,u,v,ai,aj,shape+delta))-(one_log_weibull(u,ai_d,shape+delta)+one_log_weibull(v,aj_d,shape+delta))-ff)/(delta);  
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
 grad[kk+i]=((2*log(biv_Weibull((1-nugget)*rhod,u,v,ai,aj,shape)) -(l1+l2))- ff)/delta;
    kk++;}
      h++;
    }
  return;
}


void Grad_Pair_Weibull(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];
  double ff=log(biv_Weibull((1-nugget)*rho,u,v,ai,aj,shape));

for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_Weibull((1-nugget)*rho,u,v,ai_d,aj_d,shape)) - ff)/delta; 
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
    delta=R_pow(EPS,1/2)*shape;
     //delta=R_pow(EPS,1/3);
    grad[i]=(log(biv_Weibull((1-nugget)*rho,u,v,ai,aj,shape+delta)) -
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
 grad[kk+i]=(log(biv_Weibull((1-nugget)*rhod,u,v,ai,aj,shape)) - ff)/delta;
    kk++;}
      h++;
    }
  return;
}


void Grad_Cond_Poisson(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  double nugget=nuis[nbetas];

  double l1=dpois(u,exp(ai),1);
  double l2=dpois(v,exp(aj),1);
  double ff=2*log( biv_Poisson((1-nugget)*rho,u,v,ai,aj))-(l1+l2);

for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log( biv_Poisson((1-nugget)*rho,u,v,ai_d,aj_d)) - (dpois(u,exp(ai_d),1)+ dpois(v,exp(aj_d),1)))-ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=((2*log( biv_Poisson((1-(nugget+delta))*rho,u,v,ai,aj)) - (l1+l2))-ff)/delta; 
     /// Rprintf("----grad2---%f\n",grad[i]);
    i++; 
  }

    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
 grad[kk+i]=((2*log( biv_Poisson((1-nugget)*rhod,u,v,ai,aj)) -(l1+l2))- ff)/delta;
    kk++;}
      h++;
    }
  return;
}



void Grad_Pair_Poisson(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
  double nugget=nuis[nbetas];

  double rhosill=(1-nugget)*rho;
  double ff=log( biv_Poisson(rhosill,u,v,ai,aj));

for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log( biv_Poisson(rhosill,u,v,ai_d,aj_d)) - ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=(log( biv_Poisson((1-(nugget+delta))*rho,u,v,ai,aj)) - 
              ff)/delta; 
     /// Rprintf("----grad2---%f\n",grad[i]);
    i++; 
  }

    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
      
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
 grad[kk+i]=(log( biv_Poisson((1-nugget)*rhod,u,v,ai,aj)) - ff)/delta;
    kk++;}
      h++;
    }
  return;
}



void Grad_Cond_Gamma(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];
  double  l1=one_log_gamma(u,ai,shape);
  double  l2=one_log_gamma(v,aj,shape);
  double ff=2*log(biv_gamma((1-nugget)*rho,u,v,ai,aj,shape))-(l1+l2);

for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(  (2*log(biv_gamma((1-nugget)*rho,u,v,ai_d,aj_d,shape)) -(one_log_gamma(u,ai_d,shape)+one_log_gamma(v,aj_d,shape)) )-     ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=( (2*log(biv_gamma((1-(nugget+delta))*rho,u,v,ai,aj,shape)) -(l1+l2))- ff)/delta; 
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
    grad[i]=((2*log(biv_gamma((1-nugget)*rho,u,v,ai,aj,shape+delta)) - (one_log_gamma(u,ai,shape+delta)+one_log_gamma(v,aj,shape+delta)))-ff)/(delta); 
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
 grad[kk+i]=((2*log(biv_gamma((1-nugget)*rhod,u,v,ai,aj,shape))-(l1+l2)) - ff)/delta;
    kk++;}
      h++;
    }
 return;
}


void Grad_Pair_Gamma(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];

  double ff=log(biv_gamma((1-nugget)*rho,u,v,ai,aj,shape));

for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_gamma((1-nugget)*rho,u,v,ai_d,aj_d,shape)) - ff)/delta; 
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
    delta=sqrt(EPS)*shape;
     //delta=R_pow(EPS,1/3);
    grad[i]=(log(biv_gamma((1-nugget)*rho,u,v,ai,aj,shape+delta)) -
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
 return;
}


void Grad_Cond_LogGauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
 
  double nugget=nuis[nbetas];
double sill=nuis[nbetas+1];
 double l1=one_log_loggaussian(u,ai,sill);
 double l2=one_log_loggaussian(v,aj,sill);
 double ff=2*log(d2lognorm(u,v,sill,0, ai, aj,rho*(1-nugget)))-(l1+l2);

  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(2*log(d2lognorm(u,v,sill,0, ai_d, aj_d,rho*(1-nugget))) - (one_log_loggaussian(u,ai_d,sill)+ one_log_loggaussian(v,aj_d,sill))-ff)/delta; 
   i++; }
}

  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=((2*log(d2lognorm(u,v,sill,0, ai, aj,rho*(1-(nugget+delta))))-(l1+l2)) - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=((2*log(d2lognorm(u,v,sill+delta,0, ai, aj,rho*(1-nugget))) -(one_log_loggaussian(u,ai,sill+delta)
                             +one_log_loggaussian(v,aj,sill+delta))) - ff)/delta; 
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
 grad[kk+i]=((2*log(d2lognorm(u,v,sill,0, ai, aj,rhod*(1-nugget))) -(l1+l2))- ff)/delta;
    kk++;}
      h++;
    }

  return;
}

void Grad_Pair_LogGauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
       b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
 
  double nugget=nuis[nbetas];
double sill=nuis[nbetas+1];
 //double sill=1-nugget;
 double ff=log(d2lognorm(u,v,sill,0, ai, aj,rho*(1-nugget)));

/*
  // Derivativve of the difference respect with the mean*/
   
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(d2lognorm(u,v,sill,0, ai_d, aj_d,rho*(1-nugget))) - ff)/delta; 
   i++; }
}

  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=(log(d2lognorm(u,v,sill,0, ai, aj,rho*(1-(nugget+delta)))) - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(d2lognorm(u,v,sill+delta,0, ai, aj,rho*(1-nugget))) - ff)/delta; 
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
       
 grad[kk+i]=(log(d2lognorm(u,v,sill,0, ai, aj,rhod*(1-nugget))) - ff)/delta;
    kk++;}
      h++;
    }

  return;
}

void Grad_Cond_Logistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
 

b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];
 double sill=nuis[nbetas+1];
  double l1=one_log_logistic(u,ai,sill) ;
  double l2=one_log_logistic(v,aj,sill)  ; 

double ff=2*log(biv_Logistic((1-nugget)*rho,u,v,ai,aj,sill))-(l1+l2);
  // Derivativve of the difference respect with the mean*/
  for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
                          
   grad[i]=((2*log(biv_Logistic((1-nugget)*rho,u,v,ai_d,aj_d,sill)) - (one_log_logistic(u,ai_d,sill) +one_log_logistic(v,aj_d,sill)))-ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=((2*log(biv_Logistic((1-(nugget+delta))*rho,u,v,ai,aj,sill)) -(l1+l2)) - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=((2*log(biv_Logistic((1-nugget)*rho,u,v,ai,aj,sill+delta)) - (one_log_logistic(u,ai,sill+delta)+one_log_logistic(u,ai,sill+delta)))-ff)/delta; 
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
       
grad[kk+i]=((2*log(biv_Logistic((1-nugget)*rhod,u,v,ai,aj,sill)) - (l1+l2))-ff)/delta;
    kk++;}
      h++;
    }

  return;
}


void Grad_Pair_Logistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
 

b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

   double nugget=nuis[nbetas];
 double sill=nuis[nbetas+1];

double ff=log(biv_Logistic((1-nugget)*rho,u,v,ai,aj,sill));
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_Logistic((1-nugget)*rho,u,v,ai_d,aj_d,sill)) - ff)/delta; 
   i++; }
}
 // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=(log(biv_Logistic((1-(nugget+delta))*rho,u,v,ai,aj,sill))  - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_Logistic((1-nugget)*rho,u,v,ai,aj,sill+delta)) - ff)/delta; 
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
       
grad[kk+i]=(log(biv_Logistic((1-nugget)*rhod,u,v,ai,aj,sill)) - ff)/delta;
    kk++;}
      h++;
    }

  return;
}



/*************************************************************************************/
void Grad_Cond_TwopieceT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,qqd=0.0,p11d=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

    
    double df=nuis[nbetas];
    double nugget=nuis[nbetas+1];
    double sill=nuis[nbetas+2];
    double skew=nuis[nbetas+3];

    double rho1=rho;
    double qq=qnorm((1-skew)/2,0,1,1,0);
    double p11=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,par,0);
    double l1=one_log_two_pieceT(u,ai,sill,df,skew);
    double l2=one_log_two_pieceT(v,aj,sill,df,skew);
    double ff=2*log(biv_two_pieceT(rho1,u,v,sill,df,skew,p11,ai,aj,nugget))-(l1+l2);
/*
  // Derivativve of the difference respect with the mean*/
 for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log(biv_two_pieceT(rho1,u,v,sill,df,skew,p11,ai_d,aj_d,nugget)) -
      (one_log_two_pieceT(u,ai_d,sill,df,skew)+one_log_two_pieceT(v,aj_d,sill,df,skew)))-ff)/delta; 
   i++; }
}  
// Derivvativve of the difference respect with df*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*df; 
    grad[i]=((2*log(biv_two_pieceT(rho1,u,v,sill,df+delta,skew,p11,ai,aj,nugget)) -
        (one_log_two_pieceT(u,ai,sill,df+delta,skew)+one_log_two_pieceT(v,aj,sill,df+delta,skew)))-ff)/delta; 
    i++; 
  }
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas+1]==1) {   
    delta=sqrt(EPS)*nugget; 
    p11d=pbnorm(cormod,lag,lagt,qq,qq,nugget+delta,1,par,0);
    grad[i]=((2*log(biv_two_pieceT(rho,u,v,sill,df,skew,p11d,ai,aj,nugget+delta)) - (l1+l2))-ff)/delta;
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=((2*log(biv_two_pieceT(rho1,u,v,sill+delta,df,skew,p11,ai,aj,nugget)) - 
  (one_log_two_pieceT(u,ai,sill+delta,df,skew)+one_log_two_pieceT(v,aj,sill+delta,df,skew)))-ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with skew*/  
    if(flag[nbetas+3]==1) { 
    delta=sqrt(EPS)*skew;
    qqd=qnorm((1-(skew+delta))/2,0,1,1,0);
    p11d=pbnorm(cormod,lag,lagt,qqd,qqd,nugget,1,par,0);
    grad[i]=((2*log(biv_two_pieceT(rho,u,v,sill,df,skew+delta,p11d,ai,aj,nugget)) -
            (one_log_two_pieceT(u,ai,sill,df,skew+delta)+one_log_two_pieceT(v,aj,sill,df,skew+delta)))-ff)/delta;  
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
      grad[kk+i]=((2*log(biv_two_pieceT(rhod,u,v,sill,df,skew,p11d,ai,aj,nugget)) - (l1+l2))-ff)/delta;
    kk++;}
      h++;
    }
  return;
}


/*************************************************************************************/
void Grad_Pair_TwopieceT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,qqd=0.0,p11d=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

    
    double df=nuis[nbetas];
    double nugget=nuis[nbetas+1];
    double sill=nuis[nbetas+2];
    double skew=nuis[nbetas+3];

    double rho1=rho;
    double qq=qnorm((1-skew)/2,0,1,1,0);
    double p11=pbnorm(cormod,lag,lagt,qq,qq,nugget,1,par,0);
    double ff=log(biv_two_pieceT(rho1,u,v,sill,df,skew,p11,ai,aj,nugget));
 //Rprintf("%d %d %d %d %d %d %f %f %f %f\n",nbetas,flag[0],flag[1],flag[2],flag[3],flag[4],df,nugget,sill,skew);
/*
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_two_pieceT(rho1,u,v,sill,df,skew,p11,ai_d,aj_d,nugget)) - ff)/delta; 
   i++; }
}  
// Derivvativve of the difference respect with df*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*df; 
    grad[i]=(log(biv_two_pieceT(rho1,u,v,sill,df+delta,skew,p11,ai,aj,nugget)) - ff)/delta;
    i++; 
  }
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas+1]==1) {   
    delta=sqrt(EPS)*nugget; 
    p11d=pbnorm(cormod,lag,lagt,qq,qq,nugget+delta,1,par,0);
    grad[i]=(log(biv_two_pieceT(rho,u,v,sill,df,skew,p11d,ai,aj,nugget+delta)) - 
             ff)/delta;
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_two_pieceT(rho1,u,v,sill+delta,df,skew,p11,ai,aj,nugget)) - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with skew*/  
    if(flag[nbetas+3]==1) { 
    delta=sqrt(EPS)*skew;
    qqd=qnorm((1-(skew+delta))/2,0,1,1,0);
    p11d=pbnorm(cormod,lag,lagt,qqd,qqd,nugget,1,par,0);
    grad[i]=(log(biv_two_pieceT(rho,u,v,sill,df,skew+delta,p11d,ai,aj,nugget)) - ff)/delta; 
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
      grad[kk+i]=(log(biv_two_pieceT(rhod,u,v,sill,df,skew,p11d,ai,aj,nugget)) - 
                  ff)/delta;
    kk++;}
      h++;
    }
  return;
}

void Grad_Cond_StudenT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj
       ,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
 b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];


 double df=nuis[nbetas];
 double nugget=nuis[nbetas+1];
 double sill=nuis[nbetas+2];
 double rho1=rho;
 double l1=one_log_T(u,ai,sill,df);
 double l2=one_log_T(v,aj,sill,df);
 double ff=2*log(biv_T(rho1,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df,nugget)/sill)-(l1+l2);
/*
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log(biv_T(rho1,(u-ai_d)/sill,(v-aj_d)/sill,df,nugget)/sill) -(one_log_T(u,ai_d,sill,df)+one_log_T(v,aj_d,sill,df)))-ff)/delta; 
   i++; }
}
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*df; 
    grad[i]=((2*log(biv_T(rho1,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df+delta,nugget)/sill)- 
      (one_log_T(u,ai,sill,df+delta)+one_log_T(v,aj,sill,df+delta)) -ff))/delta;
    i++; 
  }
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas+1]==1) { delta=sqrt(EPS)*nugget; 
      grad[i]=((2*log(biv_T(rho,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df,nugget+delta)/sill)- (l1+l2)) - ff)/delta;
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*sill; 
    grad[i]=(2*log(biv_T(rho1,(u-ai)/sqrt(sill+delta),(v-aj)/sqrt(sill+delta),df,nugget)/(sill+delta)) -
   -(one_log_T(u,ai,sill+delta,df)+one_log_T(v,aj,sill+delta,df))-ff )/delta;
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
       
           grad[kk+i]=((2*log(biv_T(rhod,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df,nugget)/sill) -(l1+l2)) - ff)/delta;
    kk++;}
      h++;
    }

  return;
}

void Grad_Pair_StudenT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj
       ,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
 b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];


 double df=nuis[nbetas];
 double nugget=nuis[nbetas+1];
 double sill=nuis[nbetas+2];
 double rho1=rho;
 double ff=log(biv_T(rho1,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df,nugget)/sill);
/*
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_T(rho1,(u-ai_d)/sill,(v-aj_d)/sill,df,nugget)/sill) - ff)/delta; 
   i++; }
}
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*df; 
    grad[i]=(log(biv_T(rho1,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df+delta,nugget)/sill)- ff)/delta;
    i++; 
  }
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas+1]==1) { delta=sqrt(EPS)*nugget; 
      grad[i]=(log(biv_T(rho,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df,nugget+delta)/sill) - ff)/delta;
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*sill; 
    grad[i]=(log(biv_T(rho1,(u-ai)/sqrt(sill+delta),(v-aj)/sqrt(sill+delta),df,nugget)/(sill+delta)) -ff )/delta;
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
       
           grad[kk+i]=(log(biv_T(rhod,(u-ai)/sqrt(sill),(v-aj)/sqrt(sill),df,nugget)/sill) 
            - ff)/delta;
    kk++;}
      h++;
    }

  return;
}

void Grad_Cond_Tukeyh(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];
 
  double nugget=nuis[nbetas];
  double sill=nuis[nbetas+1];//double nugget=nuis[nbetas]; nuis[2],nuis[3],nuis[1]
  double tail=nuis[nbetas+2];
    double rho1=(1-nugget)*rho;

  double  l1=one_log_tukeyh(u,ai,sill,tail);
  double  l2=one_log_tukeyh(v,aj,sill,tail);
  double ff=2*log(biv_tukey_h(rho1,u,v,ai,aj,tail,sill))-(l1+l2);

  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log(biv_tukey_h(rho1,u,v,ai_d,aj_d,tail,sill)) -
     -(one_log_tukeyh(u,ai_d,sill,tail)+one_log_tukeyh(v,aj_d,sill,tail)))-ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {   delta=sqrt(EPS)*nugget;
    grad[i]=((2*log(biv_tukey_h((1-(nugget+delta))*rho,u,v,ai,aj,tail,sill)) -(l1+l2)) - ff)/delta; 
    i++;  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=((2*log(biv_tukey_h(rho1,u,v,ai,aj,tail,sill+delta)) 
      -(one_log_tukeyh(u,ai,sill+delta,tail)+one_log_tukeyh(v,aj,sill+delta,tail)))- ff)/delta; 
    i++; 
  }

  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*tail;
    grad[i]=((2*log(biv_tukey_h(rho1,u,v,ai,aj,tail+delta,sill)) 
   -(one_log_tukeyh(u,ai,sill,tail+delta)+one_log_tukeyh(v,aj,sill,tail+delta))) - ff)/delta; 
    i++; 
  }

              h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
           grad[kk+i]=((2*log(biv_tukey_h(rhod*(1-nugget),u,v,ai,aj,tail,sill))-(l1+l2))  - ff)/delta;
    kk++;}
      h++;
    }
   return;
}   

void Grad_Pair_Tukeyh(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];
  double sill=nuis[nbetas+1];//double nugget=nuis[nbetas]; nuis[2],nuis[3],nuis[1]
  double tail=nuis[nbetas+2];
    double rho1=(1-nugget)*rho;
  double ff=log(biv_tukey_h(rho1,u,v,ai,aj,tail,sill));

  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_tukey_h(rho1,u,v,ai_d,aj_d,tail,sill)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {   delta=sqrt(EPS)*nugget;
    grad[i]=(log(biv_tukey_h((1-(nugget+delta))*rho,u,v,ai,aj,tail,sill)) - ff)/delta; 
    i++;  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_tukey_h(rho1,u,v,ai,aj,tail,sill+delta)) - ff)/delta; 
    i++; 
  }

  if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*tail;
    grad[i]=(log(biv_tukey_h(rho1,u,v,ai,aj,tail+delta,sill)) - ff)/delta; 
    i++; 
  }

              h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
           grad[kk+i]=(log(biv_tukey_h(rhod*(1-nugget),u,v,ai,aj,tail,sill)) 
            - ff)/delta;
    kk++;}
      h++;
    }
   return;
}   


void Grad_Cond_Sinh(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];



  double nugget=nuis[nbetas];
  double sill=nuis[nbetas+1];
  double skew=nuis[nbetas+2];
  double tail=nuis[nbetas+3];
  double rho1=(1-nugget)*rho;

  double l1=one_log_sas(u,ai,skew,tail,sill);
  double l2=one_log_sas(v,aj,skew,tail,sill);
  double ff=2*log(biv_sinh(rho1,u,v,ai,aj,skew,tail,sill))-(l1+l2);


  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log(biv_sinh(rho1,u,v,ai_d,aj_d,skew,tail,sill)) -
    (one_log_sas(u,ai_d,skew,tail,sill)+one_log_sas(v,aj_d,skew,tail,sill)))- ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { delta=sqrt(EPS)*nugget;
    grad[i]=((2*log(biv_sinh(1-(nugget+delta)*rho,u,v,ai,aj,skew,tail,sill+delta))-(l1+l2)) - ff)/delta; 
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=((2*log(biv_sinh(rho1,u,v,ai,aj,skew,tail,sill+delta)) -
     (one_log_sas(u,ai,skew,tail,sill+delta)+one_log_sas(v,aj,skew,tail,sill+delta)))-ff)/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    grad[i]=((2*log(biv_sinh(rho1,u,v,ai,aj,skew+delta,tail,sill)) - 
      (one_log_sas(u,ai,skew+delta,tail,sill)+one_log_sas(v,aj,skew+delta,tail,sill)))- ff)/delta; 
    i++; 
  }
  if(flag[nbetas+3]==1) { 
    delta=sqrt(EPS)*tail;
    grad[i]=((2*log(biv_sinh(rho1,u,v,ai,aj,skew,tail+delta,sill)) -
    (one_log_sas(u,ai,skew,tail+delta,sill)+one_log_sas(v,aj,skew,tail+delta,sill)))-ff)/delta; 
    i++; 
  }

              h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
           grad[kk+i]=((2*log(biv_sinh(rhod*(1-nugget),u,v,ai,aj,skew,tail,sill)) -(l1+l2)) - ff)/delta;
    kk++;}
      h++;
    }
  return;
}

void Grad_Pair_Sinh(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];



  double nugget=nuis[nbetas];
  double sill=nuis[nbetas+1];//double nugget=nuis[nbetas]; nuis[2],nuis[3],nuis[1]
  double skew=nuis[nbetas+2];
  double tail=nuis[nbetas+3];
    double rho1=(1-nugget)*rho;
  double ff=log(biv_sinh(rho1,u,v,ai,aj,skew,tail,sill));

  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_sinh(rho1,u,v,ai_d,aj_d,skew,tail,sill)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { delta=sqrt(EPS)*nugget;
    grad[i]=(log(biv_sinh(1-(nugget+delta)*rho,u,v,ai,aj,skew,tail,sill+delta)) - ff)/delta; 
    i++; }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_sinh(rho1,u,v,ai,aj,skew,tail,sill+delta)) - ff)/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    grad[i]=(log(biv_sinh(rho1,u,v,ai,aj,skew+delta,tail,sill)) - ff)/delta; 
    i++; 
  }
  if(flag[nbetas+3]==1) { 
    delta=sqrt(EPS)*tail;
    grad[i]=(log(biv_sinh(rho1,u,v,ai,aj,skew,tail+delta,sill)) - ff)/delta; 
    i++; 
  }

              h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) { 
  if(flagcor[h]==1){
       delta=sqrt(EPS)*par[h];
       parC[h]=par[h]+delta;
       rhod=CorFct(cormod,lag,lagt,parC,0,0);
       
           grad[kk+i]=(log(biv_sinh(rhod*(1-nugget),u,v,ai,aj,skew,tail,sill)) 
            - ff)/delta;
    kk++;}
      h++;
    }
  return;
}



void Grad_Cond_LogLogistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];
  double shape=nuis[nbetas+2];
  double l1=one_log_loglogistic(u,exp(ai),shape);
  double l2=one_log_loglogistic(v,exp(aj),shape);
  double ff=2*log(biv_LogLogistic((1-nugget)*rho,u,v,ai,aj,shape))-(l1+l2);
  //Rprintf("----ff---%f\n",ff);
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log(biv_LogLogistic((1-nugget)*rho,u,v,ai_d,aj_d,shape)) - 
    (one_log_loglogistic(u,exp(ai_d),shape)+one_log_loglogistic(v,exp(aj_d),shape)))-ff)/delta; 
  //   Rprintf("----grad1---%f\n",grad[i]);
   i++; 
 }}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
    delta=sqrt(EPS)*nugget;
      grad[i]=((2*log(biv_LogLogistic((1-(nugget+delta))*rho,u,v,ai,aj,shape)) - (l1+l2))-ff)/delta; 
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
    grad[i]=((2*log(biv_LogLogistic((1-nugget)*rho,u,v,ai,aj,shape+delta)) -
      (one_log_loglogistic(u,exp(ai),shape+delta)+one_log_loglogistic(v,exp(aj),shape+delta)))-ff)/(delta); 
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
 grad[kk+i]=((2*log(biv_LogLogistic((1-nugget)*rhod,u,v,ai,aj,shape))-(l1+l2)) - ff)/delta;
    kk++;}
      h++;
    }
//Rprintf("%f %f %f %f %f %f %f \n",grad[0],grad[1],grad[2],grad[3],grad[4],grad[5],grad[6]);
  return;
}


void Grad_Pair_LogLogistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
 // double sum0=0.0,sum1=0.0,sum2=0.0,
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC;
    b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];

  double nugget=nuis[nbetas];

  double shape=nuis[nbetas+2];
  double rhosill=(1-nugget)*rho;
  double ff=log(biv_LogLogistic(rhosill,u,v,ai,aj,shape));
  //Rprintf("----ff---%f\n",ff);
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
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



void Grad_Cond_Skewgauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];


  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double sill=nuis[nbetas+1];double nugget=nuis[nbetas];
  double skew=nuis[nbetas+2];

   double l1=one_log_SkewGauss(u,ai,sill,skew);
   double l2=one_log_SkewGauss(v,aj,sill,skew);

  double  ff=2*log(biv_skew(rho,u,v,ai,aj,sill,skew,nugget))-(l1+l2);
/*
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=((2*log(biv_skew(rho,u,v,ai_d,aj_d,sill,skew,nugget))-(one_log_SkewGauss(u,ai_d,sill,skew)+one_log_SkewGauss(v,aj_d,sill,skew)))- ff)/delta; 
   i++; }
}  
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=((2*log(biv_skew(rho,u,v,ai,aj,sill+delta,skew,nugget+delta)) -(l1+l2))- ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=((2*log(biv_skew(rho,u,v,ai,aj,sill+delta,skew,nugget)) -
      (one_log_SkewGauss(u,ai,sill+delta,skew)+one_log_SkewGauss(v,aj,sill+delta,skew)))-ff)/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    grad[i]=((2*log(biv_skew(rho,u,v,ai,aj,sill,skew+delta,nugget)) -
      (one_log_SkewGauss(u,ai,sill,skew+delta)+one_log_SkewGauss(v,aj,sill,skew+delta)))-ff)/delta; 
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
       
      grad[kk+i]=((2*log(biv_skew(rhod,u,v,ai,aj,sill,skew,nugget)) -(l1+l2))- ff)/delta;
    kk++;}
      h++;
    }
  return;
}

void Grad_Pair_Skewgauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{ 
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0;
  double delta=0,*b1,*parC; 
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];


  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double sill=nuis[nbetas+1];double nugget=nuis[nbetas];
  double skew=nuis[nbetas+2];
    double  ff=log(biv_skew(rho,u,v,ai,aj,sill,skew,nugget));
/*
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }
   grad[i]=(log(biv_skew(rho,u,v,ai_d,aj_d,sill,skew,nugget)) - ff)/delta; 
   i++; }
}  
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) {  delta=sqrt(EPS)*nugget;
    grad[i]=(log(biv_skew(rho,u,v,ai,aj,sill+delta,skew,nugget+delta)) - ff)/delta; 
    i++; 
  }
  /* Derivvativve of the difference respect with the sill*/  
  if(flag[nbetas+1]==1) { 
    delta=sqrt(EPS)*sill;
    grad[i]=(log(biv_skew(rho,u,v,ai,aj,sill+delta,skew,nugget)) - ff)/delta; 
    i++; 
  }
    if(flag[nbetas+2]==1) { 
    delta=sqrt(EPS)*skew;
    grad[i]=(log(biv_skew(rho,u,v,ai,aj,sill,skew+delta,nugget)) - ff)/delta; 
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
       
      grad[kk+i]=(log(biv_skew(rhod,u,v,ai,aj,sill,skew,nugget)) - ff)/delta;
    kk++;}
      h++;
    }
  return;
}


void Grad_Cond_Wrapped(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,alfa=2.0;
  double delta=0,*b1,*parC;
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];


  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double sill=nuis[nbetas+1];double nugget=nuis[nbetas];
 double ff=log(biv_wrapped(alfa,u,v,ai,aj,nugget,sill,rho));

/*
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }

   grad[i]=(log(biv_wrapped(alfa,u,v,ai_d,aj_d,nugget,sill,rho)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
      delta=sqrt(EPS)*nugget;
      grad[i]=(log(biv_wrapped(alfa,u,v,ai,aj,nugget+delta,sill,rho)) -ff)/(delta);
      i++;}
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


void Grad_Pair_Wrapped(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double **sX,int l,int m,double *betas)
{
  // Initialization variables:
  int h=0, i=0, j=0,kk=0,o=0,k=0;
  double rhod=0.0,ai_d=0.0,aj_d=0.0,alfa=2.0;
  double delta=0,*b1,*parC;
  b1=(double *) Calloc(nbetas,double);
  parC=(double *) Calloc(nparcT[0],double);
  for(k=0;k<nparcT[0];k++) parC[k]=par[k];


  //q3,dp1dBeta,dp2dBeta,dp11dBeta,dp11dsill,p11_dcorr,p11_dsill,C;
  double sill=nuis[nbetas+1];double nugget=nuis[nbetas];
 double ff=log(biv_wrapped(alfa,u,v,ai,aj,nugget,sill,rho));

/*
  // Derivativve of the difference respect with the mean*/
for(kk=0;kk<nbetas;kk++){

    for(o=0;o<nbetas;o++) {b1[o]=betas[o];} 
  if(flag[kk]==1){
     delta=sqrt(EPS)*(betas[kk]);
     b1[kk]=betas[kk]+delta;
     ai_d=0.0;aj_d=0.0;
     for(o=0;o<nbetas;o++){
                           
                           ai_d=ai_d+sX[l][o]*(  b1[o] );
                           aj_d=aj_d+sX[m][o]*(  b1[o] );
                           }

   grad[i]=(log(biv_wrapped(alfa,u,v,ai_d,aj_d,nugget,sill,rho)) - ff)/delta; 
   i++; }
}
  // Derivvativve of the difference respect with the nugget*/
  if(flag[nbetas]==1) { 
      delta=sqrt(EPS)*nugget;
      grad[i]=(log(biv_wrapped(alfa,u,v,ai,aj,nugget+delta,sill,rho)) -ff)/(delta);
      i++;}
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
