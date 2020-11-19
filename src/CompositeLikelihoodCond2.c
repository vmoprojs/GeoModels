#include "header.h"

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

/******************************************************************************************/
void Comp_Cond_Gauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
       int i=0;
    double  weights=1.0,sill,nugget,corr,bl,l1,l2;
     sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      bl=log_biv_Norm((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,0);
                      l1= dnorm(data1[i], mean1[i],sqrt(sill),1);
                      l2= dnorm(data2[i], mean2[i],sqrt(sill),1);
                      *res+= (2*bl-l1-l2)*weights;
                    }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/******************************************************************************************/
void Comp_Cond_SinhGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
   
    int i=0;double corr,zi,zj,qi,qj,Z1,Z2,bb=0.0,weights=1.0,l1=0.0,l2=0.0,b1=0.0,b2=0.0;
           if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}

   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    qi=(zi-mean1[i])/(sqrt(nuis[1]));qj=(zi-mean2[i])/(sqrt(nuis[1]));
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    b1=nuis[3]*asinh(qi)-nuis[2];
                    b2=nuis[3]*asinh(qj)-nuis[2];
                    Z1=sinh(b1);Z2=sinh(b2);
                    l1=-0.5*log(R_pow(qi,2)+1)-0.5*log(2*M_PI*nuis[1])+log(cosh(b1))+log(nuis[3])-Z1*Z1/2;
                    l2=-0.5*log(R_pow(qj,2)+1)-0.5*log(2*M_PI*nuis[1])+log(cosh(b2))+log(nuis[3])-Z2*Z2/2;;
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bb=2*log(biv_sinh((1-nuis[0])*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]))-(l1+l2);
                    *res+= weights*bb;
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gamma2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
   
    int i;double corr,zi,zj,weights=1.0,bl=1.0,l1=0.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    l1=(nuis[2]/2)*log(nuis[2]/(2*exp(mean1[i])))+(nuis[2]/2-1)*log(zi)-(nuis[2]/(2*exp(mean1[i])))*zi-lgammafn(nuis[2]/2);
                    l2=(nuis[2]/2)*log(nuis[2]/(2*exp(mean2[i])))+(nuis[2]/2-1)*log(zj)-(nuis[2]/(2*exp(mean2[i])))*zj-lgammafn(nuis[2]/2);
                    bl=2*log(biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))
                     -(l1+l2);
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);      
                           
  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Weibull2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double corr,zi,zj,weights=1.0,bl=0.0,l1=0.0,l2=0.0;
      double nugget=nuis[0];
      double c1,c2;
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    c1=exp(mean1[i])/(gammafn(1+1/nuis[2]));
                    c2=exp(mean2[i])/(gammafn(1+1/nuis[2]));
                    l1=log(nuis[2])-nuis[2]*log(c1)+(nuis[2]-1)*log(zi)-R_pow(zi/c1,nuis[2]);
                    l2=log(nuis[2])-nuis[2]*log(c2)+(nuis[2]-1)*log(zj)-R_pow(zj/c2,nuis[2]);

                    bl=2*log(biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))
                    - (l1+l2);
                     *res+= weights*bl;
               
                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/*********************************************************/
void Comp_Cond_LogGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
   
    int i=0;double corr,zi,zj,qi,qj,weights=1.0,bl=0.0,l1=0.0,l2=0.0;
    double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){

                    zi=data1[i];zj=data2[i];
                    qi=zi*exp(sill/2);qj=zj*exp(sill/2);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    l1=-0.5*R_pow((log(qi)-mean1[i]),2)/sill-log(qi)-log(sqrt(sill))-0.5*log(2*M_PI)+sill/2;
                    l2=-0.5*R_pow((log(qj)-mean2[i]),2)/sill-log(qj)-log(sqrt(sill))-0.5*log(2*M_PI)+sill/2;
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl=2*log(d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],(1-nugget)*corr))
                     -(l1+l2);
                    *res+= weights*bl;
                    }}
   
    if(!R_FINITE(*res))*res = LOW;
    return;
}






/*********************************************************/
void Comp_Cond_Beta2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
   
    int i;double corr,zi,zj,weights=1.0,bl,l1=0.0,l2=0.0;
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]-min)/(max-min); zj=(data2[i]-min)/(max-min);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);      
   l1=(nuis[2]/2-1)*log(zi)+(nuis[3]/2-1)*log(1-zi)+lgammafn(0.5*(nuis[2]+nuis[3]))-lgammafn(nuis[2]/2)-lgammafn(nuis[3]/2)-2*log(max-min);
   l2=(nuis[2]/2-1)*log(zj)+(nuis[3]/2-1)*log(1-zj)+lgammafn(0.5*(nuis[2]+nuis[3]))-lgammafn(nuis[2]/2)-lgammafn(nuis[3]/2)-2*log(max-min);    
                  bl=2*log(biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
       
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_Kumaraswamy2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
   
    int i;double corr,zi,zj,qi,qj,weights=1.0,bl,l1=0.0,l2=0.0,k1=0.0,k2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    qi=(zi-min)/(max-min); qj=(zj-min)/(max-min);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    k1=1-pow(qi,nuis[3]); k2=1-pow(qj,nuis[3]);
                    l1=log(nuis[2])+log(nuis[3])+(nuis[3]-1)*log(qi)+(nuis[2]-1)*log(k1)-2*log(max-min);
                    l2=log(nuis[2])+log(nuis[3])+(nuis[3]-1)*log(qj)+(nuis[2]-1)*log(k2)-2*log(max-min);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);            
                  bl=2*log(biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
       
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/*********************************************************/
void Comp_Cond_Kumaraswamy22mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
   
    int i=0;double corr,zi,zj,qi,qj,weights=1.0,bl,ki=0.0,kj=0.0,
                      l1=0.0,l2=0.0,mi,mj,shapei,shapej;
    double nugget=nuis[0];
    double min=nuis[4];double max=nuis[5];
    
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
   mi=1/(1+exp(-mean1[i])); mj=1/(1+exp(-mean2[i]));
     shapei=log(0.5)/log(1-pow(mi,nuis[3]));shapej=log(0.5)/log(1-pow(mj,nuis[3]));
                    zi=(data1[i]); zj=(data2[i]);
                    qi=(zi-min)/(max-min); qj=(zj-min)/(max-min);
                    ki=1-pow(qi,nuis[3]); kj=1-pow(qj,nuis[3]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    l1=log(shapei)+log(nuis[3])+(nuis[3]-1)*log(qi)+(shapei-1)*log(ki)-log(max-min);
                    l2=log(shapej)+log(nuis[3])+(nuis[3]-1)*log(qj)+(shapej-1)*log(kj)-log(max-min);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);            
                  bl=2*log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
       
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
