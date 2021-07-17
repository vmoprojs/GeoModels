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
//Rprintf("ggh %f  %d  \n", maxdist[0],*weigthed);
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
  //Rprintf("%f %f %f  %d\n",lags[i],data1[i],data2[i],*npairs);
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                             //weights=CorFunW_gen(lags[i],6, 3, maxdist[0]);
                      bl=log_biv_Norm((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,0);
                      l1= dnorm(data1[i], mean1[i],sqrt(sill),1);
                      l2= dnorm(data2[i], mean2[i],sqrt(sill),1);
                      *res+= (2*bl-l1-l2)*weights;
                    }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Tukeyhh2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0 ,l1=0.0,l2=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];

                corr=CorFct(cormod,lags[i],0,par,0,0);
                l1=one_log_tukeyhh(zi,mean1[i],sill,h1,h2);
                l2=one_log_tukeyhh(zj,mean2[i],sill,h1,h2);
               if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
               bl=2*log(biv_tukey_hh((1-nugget)*corr,zi,zj,mean1[i],mean2[i],sill,h1,h2))-(l1+l2);
                             *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Tukeyh2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,l1=0.0,l2=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                corr=CorFct(cormod,lags[i],0,par,0,0);
                l1=one_log_tukeyh(zi,mean1[i],sill,tail);
                l2=one_log_tukeyh(zj,mean2[i],sill,tail);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                bl=2*log(biv_tukey_h((1-nugget)*corr,zi,zj,mean1[i],mean2[i],tail,sill))-(l1+l2);
                             *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_SkewGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
      double sill=nuis[1];double nugget=nuis[0],l1=0.0,l2=0.0,bb=0.0;
     // double skew2  = R_pow(nuis[2],2);
     // double vari2  = R_pow(nuis[1],1);
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}
    int i=0;double corr,zi,zj,weights=1.0;
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
    l1=one_log_SkewGauss(zi,mean1[i],nuis[1],nuis[2]);
    l2=one_log_SkewGauss(zj,mean2[i],nuis[1],nuis[2]);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
    bb=2*log(biv_skew(corr,zi,zj,mean1[i],mean2[i],nuis[1],nuis[2],nuis[0]))-(l1+l2);

    //bb=2*0.5-(l1+l2);
                  *res+= weights*bb;
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_T2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,qi,qj,weights=1.0,l1=0.0,l2=0.0;
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
    double df1=1/nuis[0];
      if( sill<0||nugget<0||nugget>=1||df<0||df>0.5){*res=LOW; return;}

   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                qi=(zi-mean1[i])/sqrt(sill);qj=(zj-mean2[i])/sqrt(sill);
                corr=CorFct(cormod,lags[i],0,par,0,0);
                  l1=one_log_T(zi,mean1[i],sill,df1);
                  l2=one_log_T(zj,mean2[i],sill,df1);
              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                bl=2*log(biv_T(corr,qi,qj,df,nugget)/sill)-(l1+l2);
                *res+= weights*bl;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gauss_misp_T2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double zi,zj,weights=1.0,corr,df=0.0,bl,l1=0.0,l2=0.0,var=0.0;

     double sill=nuis[2];
    double nugget=nuis[1];


    if( sill<0||nugget<0||nugget>=1||nuis[0]<0||nuis[0]>0.5){*res=LOW; return;}
    df=1/nuis[0];
    var=sill*df/(df-2);
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){

           corr=(1-nugget)*CorFct(cormod,lags[i],0,par,0,0);
           corr=exp(log(df-2)+2*lgammafn(0.5*(df-1))-(log(2)+2*lgammafn(df/2))+log(hypergeo(0.5,0.5, df/2,corr*corr))+log(corr*(1-nugget)));
           zi=data1[i];zj=data2[i];

         l1=dnorm(zi,mean1[i],sqrt(sill*df/(df-2)),1); l2=dnorm(zj,mean2[i],sqrt(sill*df/(df-2)),1);
         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                     bl=2*log_biv_Norm(corr,data1[i],data2[i],mean1[i],mean2[i],sill*df/(df-2),0)-(l1+l2);
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_Gauss_misp_SkewT2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double weights=1.0,sill,nugget,skew,corr,corr2,df,bl,l1,l2;


    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];

    if(df<2||fabs(skew)>1||sill<0||nugget<0||nugget>=1){*res=LOW; return;}
    //auxuliary variables
    double D1=(df-1)/2;
    double D2=df/2;
    //double delta=skew/sqrt(1-skew*skew);
    double MM=(sqrt(df)*skew)/(sqrt(M_PI))*exp(lgammafn(D1)-lgammafn(D2));
    double FF=(df/(df-2)-MM*MM);

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=CorFct(cormod,lags[i],0,par,0,0)*(1-nugget);
                     corr2=corr_skewt(corr,df,skew);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          bl=log_biv_Norm(corr2,data1[i],data2[i],mean1[i]+sqrt(sill)*MM,
                                                                  mean2[i]+sqrt(sill)*MM,
                                                                  sill*FF,0);
                          l1=dnorm(data1[i],mean1[i]+sqrt(sill)*MM,sqrt(sill*FF),1);
                          l2=dnorm(data2[i],mean2[i]+sqrt(sill)*MM,sqrt(sill*FF),1);
                        *res+= (2*bl-l1-l2)*weights;


                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
/*********************************************************/
void Comp_Cond_Gauss_misp_Tukeygh2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,corr2,zi,zj,weights=1.0,eta,tail,sill,nugget,u,eta2,mu,vv,l1,l2;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

    eta2=eta*eta;
    u=1-tail;
    mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
    vv=((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*
                           sqrt(1-2*tail))-mu*mu);
    if(fabs(eta)<1e-5)
           {
           mu=0.0;
           vv=R_pow(1-2*tail,-3/2);
           }
         if(sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
          zi=data1[i];zj=data2[i];
if(!ISNAN(zi)&&!ISNAN(zj) ){

                    corr=(1-nugget)*CorFct(cormod,lags[i],0,par,0,0);
                    corr2=corr_tukeygh(corr,eta,tail);
                 //   if(corr2<0) Rprintf("%f %f %f \n",corr2,par[0],par[1]);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                bl=log_biv_Norm(corr2,zi,zj,mean1[i]+sqrt(sill)*mu,
                                            mean2[i]+sqrt(sill)*mu, sill*vv,0);

                      l1= dnorm(zi, mean1[i]+sqrt(sill)*mu,sqrt(sill*vv),1);
                      l2= dnorm(zj, mean2[i]+sqrt(sill)*mu,sqrt(sill*vv),1);
                      *res+= (2*bl-l1-l2)*weights;

                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_SinhGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{

    int i=0;double corr,zi,zj,bb=0.0,weights=1.0,l1=0.0,l2=0.0;
           if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    l1=one_log_sas(zi,mean1[i],nuis[2],nuis[3],nuis[1]);
                    l2=one_log_sas(zj,mean2[i],nuis[2],nuis[3],nuis[1]);
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
                    l1=one_log_gamma(zi,mean1[i],nuis[2]);
                    l2=one_log_gamma(zj,mean2[i],nuis[2]);
                    bl=2*log(biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))-(l1+l2);
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
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    l1=one_log_weibull(zi,mean1[i],nuis[2]);
                    l2=one_log_weibull(zj,mean2[i],nuis[2]);

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
    int i=0;double corr,zi,zj,weights=1.0,bl=0.0,l1=0.0,l2=0.0;
    double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     zi=(data1[i]);zj=(data2[i]);
                    l1=one_log_loggaussian(zi,mean1[i],sill);
                    l2=one_log_loggaussian(zj,mean2[i],sill);
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
                    zi=data1[i]; zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    l1=one_log_beta(zi,nuis[2],nuis[3],min,max);
                    l2=one_log_beta(zj,nuis[2],nuis[3],min,max);
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
    int i;double corr,zi,zj,weights=1.0,bl,l1=0.0,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);

                    l1=one_log_kumma(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma(zj,mean2[i],nuis[2],nuis[3],min,max);

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
    int i;double corr,zi,zj,weights=1.0,bl,l1=0.0,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    l1=one_log_kumma2(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma2(zj,mean2[i],nuis[2],nuis[3],min,max);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                  bl=2*log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gauss_misp_Pois2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,l1,l2;
    double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}
double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}
    double *dat;
    dat=(double *) Calloc(N,double);
    for(i=0;i<npairs[0];i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0)*(1-nugget);
                      corr1=corr_pois(corr,mui, muj);
                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                        M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;
                      l1=dnorm(data1[i],mui,sqrt(mui),1);
                      l2=dnorm(data2[i],muj,sqrt(muj),1);;
                      bl=2*log(dNnorm(N,M,dat))-(l1+l2);

                      *res+= bl*weights;
                    }}
   for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);

    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_BinomNNGauss_misp2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, N=2,n1,n2;
    double u,v,m1,m2,l1,l2,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success

    double **M;
    M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}
    double *dat;
    dat=(double *) Calloc(N,double);

    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*corr);
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 n1=NN[i];n2=NN[i+npairs[0]];
                 m1=n1*p1;m2=n2*p2;
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 M[0][0]=m1*(1-p1);   M[1][1]=m2*(1-p2);  // var1 var2
                 M[0][1]= fmin_int(n1,n2)*(p11-p1*p2) ;       // covariance
                 M[1][0]= M[0][1];
                 dat[0]=u-m1;dat[1]=v-m2; 
                 //Rprintf("%d %f %f %f \n",fmin_int(n1,n2),p1,p2,p11 );
                 l1=dnorm(u,m1,sqrt(m1*(1-p1)),1);
                 l2=dnorm(v,m2,sqrt(m2*(1-p2)),1);;
                 bl= 2*log(dNnorm(N,M,dat)) -(l1+l2); 
                 *res+= bl*weights;       
                }}
    for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisGamma2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,l1,l2,bi,bj,vvi,vvj;
    double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}
double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}
    double *dat;
    dat=(double *) Calloc(N,double);
    for(i=0;i<npairs[0];i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    bi= nuis[2]/mui; bj= nuis[2]/muj;
                    vvi= mui*(1+1/bi); vvj= muj*(1+1/bj);
                     corr=CorFct(cormod,lags[i],0,par,0,0)*(1-nugget);
                      corr1=corr_pois_gen(corr,mui, muj, nuis[2]);
                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                        M[0][0]=vvi; M[1][1]=vvj;M[0][1]=sqrt(vvi*vvj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;
                      l1=dnorm(data1[i],mui,sqrt(vvi),1);
                      l2=dnorm(data2[i],muj,sqrt(vvj),1);;
                      bl=2*log(dNnorm(N,M,dat))-(l1+l2);

                      *res+= bl*weights;
                    }}
   for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_PoisGamma2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l1,l2;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      l1=one_log_dpoisgamma(uu,mui,nuis[2]);
                      l2=one_log_dpoisgamma(ww,muj,nuis[2]);
                      bl=2*log(biv_PoissonGamma((1-nugget)*corr,uu,ww,mui, muj,nuis[2]))
                        - (l1+l2);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Pois2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l1,l2;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      l1=dpois(uu,mui,1);
                      l2=dpois(ww,muj,1);
                      bl=2*log(biv_Poisson((1-nugget)*corr,uu,ww,mui, muj))
                        - (l1+l2);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
    p11=pbnorm22(ai,aj,(1-nugget)*corr);
    //Rprintf("p11: %f\n",p11);
                p1=pnorm(ai,0,1,1,0);
                p2=pnorm(aj,0,1,1,0);
                u=data1[i];v=data2[i];
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u; vv=(int) v;
                          l1=dbinom(uu,NN[0],p1,1);
                          l2=dbinom(vv,NN[0],p2,1);
                        bl=2*log(biv_binom (NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomLogi2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
    //Rprintf("p11: %f\n",p11);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                u=data1[i];v=data2[i];
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u; vv=(int) v;
                          l1=dbinom(uu,NN[0],p1,1);
                          l2=dbinom(vv,NN[0],p2,1);
                        bl=2*log(biv_binom (NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
void Comp_Cond_BinomNNGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                 p11=pbnorm22(ai,aj,(1-nugget)*corr);
                 p1=pnorm(ai,0,1,1,0);
                 p2=pnorm(aj,0,1,1,0);
                 u=data1[i];v=data2[i];
                 n1=NN[i];n2=NN[i+npairs[0]];
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 l1=dbinom(uu,n1,p1,1);
                 l2=dbinom(vv,n2,p2,1);
                 bl=2*log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-(l1+l2);
                 *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_BinomNNLogi2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                 p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                 p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                 u=data1[i];v=data2[i];
                 n1=NN[i];n2=NN[i+npairs[0]];
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 uu=(int) u; vv=(int) v;
                 l1=dbinom(uu,n1,p1,1);
                 l2=dbinom(vv,n2,p2,1);
                 bl=2*log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-(l1+l2);
                 *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_BinomnegGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;
                         vv=(int) v;
                         l1=one_log_negbinom_marg(uu,NN[0],p1);
                         l2=one_log_negbinom_marg(vv,NN[0],p2);
                        bl=2*log(biv_binomneg(NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Cond_BinomnegLogi2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;
                         vv=(int) v;
                         l1=one_log_negbinom_marg(uu,NN[0],p1);
                         l2=one_log_negbinom_marg(vv,NN[0],p2);
                        bl=2*log(biv_binomneg(NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Cond_TWOPIECETukeyh2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget,l1=0.0,l2=0.0;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
           zi=data1[i];zj=data2[i];


           corr=CorFct(cormod,lags[i],0,par,0,0);

            l1=one_log_two_pieceTukey(zi,mean1[i],sill,tail,eta);
            l2=one_log_two_pieceTukey(zj,mean2[i],sill,tail,eta);

           p11=pbnorm22(qq,qq,corr);
           if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=2*log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-(l1+l2);
               *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_TWOPIECET2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,qq,l1=0.0,l2=0.0;
    double eta=nuis[3];  //skewness parameter
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}
       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                corr=CorFct(cormod,lags[i],0,par,0,0);
                l1=one_log_two_pieceT(zi,mean1[i],sill,df,eta);
                l2=one_log_two_pieceT(zj,mean2[i],sill,df,eta);
                 p11=pbnorm22(qq,qq,corr);
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                 /********************************************************/
                 bl=2*log(biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget))
                    -(l1+l2);
                 /********************************************************/
                         *res+= weights*bl;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_TWOPIECEGauss2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,nugget,l1=0.0,l2=0.0;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                l1=one_log_two_pieceGauss(zi,mean1[i],sill,eta);
                l2=one_log_two_pieceGauss(zj,mean2[i],sill,eta);

                p11=pbnorm22(qq,qq,corr);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl=2*log(biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]))-l1-l2;
                    *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/********************************************************/
void Comp_Cond_TWOPIECEBIMODAL2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,alpha,weights=1.0,p11,eta,qq,sill,df,nugget,delta ,l1=0.0,l2=0.0;

    eta=nuis[4];  //skewness parameter
    delta=nuis[3];
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
    alpha=2*(delta+1)/df;
 if( fabs(eta)>1||df<0||nugget<0||nugget>=1||delta<0||sill<0) {*res=LOW;  return;}
    qq=qnorm((1-eta)/2,0,1,1,0);
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                l1=one_log_bomidal(zi,mean1[i],sill,df,delta,eta);
                l2=one_log_bomidal(zj,mean2[i],sill,df,delta,eta);
                p11=pbnorm22(qq,qq,corr);
                 if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    /********************************************************/
                  bl=2*log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-(l1+l2);
                           *res+= weights*bl;
                }}}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/********************************************************/
void Comp_Cond_BinomnegGaussZINB2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],0,par,0,0);
                    u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;
                         vv=(int) v;
                    l1=one_log_BinomnegZIP(uu,NN[0],ai,mup);
                    l2=one_log_BinomnegZIP(vv,NN[0],aj,mup);

                    bl=2*log(biv_binomnegZINB(NN[0],corr,uu,vv,ai,aj,nugget1,nugget2,mup))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;

}
/************************************************/
void Comp_Cond_PoisZIP2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{

    int i=0, uu,vv;
    double weights=1.0,corr,mui,muj,bl,l1=0.0,l2=0.0,u,v;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                        u=data1[i];v=data2[i];
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;vv=(int) v;

                    l1=one_log_PoisZIP(uu,mui,mup);
                    l2=one_log_PoisZIP(vv,muj,mup);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                     bl=2*log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-(l1+l2);
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisZIP2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double weights=1.0,corr,mui,muj,bl ,l1=0.0,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
    double p=pnorm(mup,0,1,1,0);

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}

  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);


                      l1=dnorm(data1[i],(1-p)*mui,sqrt(mui*(1-p)*(1+p*mui)),1);
                      l2=dnorm(data2[i],(1-p)*muj,sqrt(muj*(1-p)*(1+p*muj)),1);


                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      bl=2*log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-(l1+l2);

                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_LogLogistic2mem(int *cormod, double *data1,double *data2,int *NN,
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

                    l1=one_log_loglogistic(zi,exp(mean1[i]),nuis[2]);
                    l2=one_log_loglogistic(zj,exp(mean2[i]),nuis[2]);

                    bl=2*log(biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))
                     -(l1+l2);
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



void Comp_Cond_Logistic2mem(int *cormod, double *data1,double *data2,int *NN,
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

                    l1=one_log_logistic(zi,mean1[i],nuis[1]) ;
                    l2=one_log_logistic(zj,mean2[i],nuis[1])  ;

                    bl=2*log(biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1]))

                     -(l1+l2);
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPACE TIME CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
void Comp_Cond_Gauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double weights=1.0,l1=0.0,l2=0.0;
    double bl=0.0, corr=0.0;

      double  nugget=nuis[0];
    double sill=nuis[1];

      if(sill<0 || nugget<0||nugget>=1){*res=LOW; return;}
  for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                          //  s12=nuis[1]*CorFct(cormod,lags,0,par,t,v);
                           corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
              bl=log_biv_Norm((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,0);
                      l1= dnorm(data1[i], mean1[i],sqrt(sill),1);
                      l2= dnorm(data2[i], mean2[i],sqrt(sill),1);
                      *res+= (2*bl-l1-l2)*weights;
                            }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_T_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double corr,zi,zj,weights=1.0,bl,l1=0.0,l2=0.0;
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];

    if( sill<0||nugget<0||nugget>=1||df<0||df>0.5){*res=LOW; return;}
    //if( sill<0||nugget<0||nugget>=1){*res=LOW; return;}
         for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
   bl=log(biv_T(corr,(zi-mean1[i])/sqrt(sill),(zj-mean2[i])/sqrt(sill),df,nugget)/sill);

       l1=one_log_T(zi,mean1[i],sill,1/df);
       l2=one_log_T(zj,mean2[i],sill,1/df);
                             *res+= (2*bl-(l1+l2))*weights;

                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_Gauss_misp_T_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double  corr,df=0.0,u=0.0, w=0.0,weights=1.0,bl,l1=0.0,l2=0.0,var;
     double sill=nuis[2];
    double nugget=nuis[1];


    if( sill<0||nugget<0||nugget>=1||nuis[0]<0||nuis[0]>0.5){*res=LOW; return;}

      df=1/nuis[0];
        var=sill*df/(df-2);

    for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            corr=(1-nugget)*CorFct(cormod,lags[i],lagt[i],par,0,0);
                 corr=exp(log(df-2)+2*lgammafn(0.5*(df-1))-(log(2)+2*lgammafn(df/2))+log(hypergeo(0.5,0.5, df/2,corr*corr))+log(corr*(1-nugget)));
                               u=data1[i];
                                w=data2[i];
                  if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
bl=log_biv_Norm(corr,u,w,mean1[i],mean2[i],var,0);

                      l1= dnorm(u, mean1[i],sqrt(var),1);
                      l2= dnorm(w, mean2[i],sqrt(var),1);
                      *res+= (2*bl-l1-l2)*weights;

                *res+= bl*weights;
                                    }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Tukeyh_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,l1=0.0,l2=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

                l1=one_log_tukeyh(zi,mean1[i],sill,tail);
                l2=one_log_tukeyh(zj,mean2[i],sill,tail);

                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                bl=2*log(biv_tukey_h((1-nugget)*corr,zi,zj,mean1[i],mean2[i],tail,sill))-(l1+l2);
                             *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Tukeyhh_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0 ,l1=0.0,l2=0.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];

                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                l1=one_log_tukeyhh(zi,mean1[i],sill,h1,h2);
                l2=one_log_tukeyhh(zj,mean2[i],sill,h1,h2);

              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
               bl=2*log(biv_tukey_hh((1-nugget)*corr,zi,zj,mean1[i],mean2[i],sill,h1,h2))-(l1+l2);
                             *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_SkewGauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
      double sill=nuis[1];double nugget=nuis[0],l1=0.0,l2=0.0,bb=0.0;
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}

    int i=0;double corr,zi,zj,weights=1.0;
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    l1=one_log_SkewGauss(zi,mean1[i],nuis[1],nuis[2]);
                    l2=one_log_SkewGauss(zj,mean2[i],nuis[1],nuis[2]);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    bb=2*log(biv_skew(corr,zi,zj,mean1[i],mean2[i],nuis[1],nuis[2],nuis[0]))-(l1+l2);
                  *res+= weights*bb;
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_SinhGauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;double corr,zi,zj,bb=0.0,weights=1.0,l1=0.0,l2=0.0;
           if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    l1=one_log_sas(zi,mean1[i],nuis[2],nuis[3],nuis[1]);
                    l2=one_log_sas(zj,mean2[i],nuis[2],nuis[3],nuis[1]);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    bb=2*log(biv_sinh((1-nuis[0])*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]))-(l1+l2);
                    *res+= weights*bb;
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Gamma_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{

    int i;double corr,zi,zj,weights=1.0,bl=1.0,l1=0.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    l1=one_log_gamma(zi,mean1[i],nuis[2]);
                    l2=one_log_gamma(zj,mean2[i],nuis[2]);
                    bl=2*log(biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))-(l1+l2);
                   if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Weibull_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double corr,zi,zj,weights=1.0,bl=0.0,l1=0.0,l2=0.0;
      double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    l1=one_log_weibull(zi,mean1[i],nuis[2]);
                    l2=one_log_weibull(zj,mean2[i],nuis[2]);
                    bl=2*log(biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))
                    - (l1+l2);
                     *res+= weights*bl;

                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_LogGauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;double corr,zi,zj,weights=1.0,bl=0.0,l1=0.0,l2=0.0;
    double sill=nuis[1];double nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                     zi=(data1[i]);zj=(data2[i]);
                    l1=one_log_loggaussian(zi,mean1[i],sill);
                    l2=one_log_loggaussian(zj,mean2[i],sill);
                   if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    bl=2*log(d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],(1-nugget)*corr))
                     -(l1+l2);
                    *res+= weights*bl;
                    }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Beta_st2mem(int *cormod, double *data1,double *data2,int *NN,
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
                    zi=data1[i]; zj=data2[i];
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    l1=one_log_beta(zi,nuis[2],nuis[3],min,max);
                    l2=one_log_beta(zj,nuis[2],nuis[3],min,max);
                  bl=2*log(biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);

        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/

void Comp_Cond_Kumaraswamy_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{

    int i;double corr,zi,zj,weights=1.0,bl,l1=0.0,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
               zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

                    l1=one_log_kumma(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma(zj,mean2[i],nuis[2],nuis[3],min,max);

                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                  bl=2*log(biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);

        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Kumaraswamy2_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double corr,zi,zj,weights=1.0,bl,l1=0.0,l2=0.0;
    double nugget=nuis[0];
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    l1=one_log_kumma2(zi,mean1[i],nuis[2],nuis[3],min,max);
                    l2=one_log_kumma2(zj,mean2[i],nuis[2],nuis[3],min,max);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                  bl=2*log(biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max))-(l1+l2);
        *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



/*********************************************************/
void Comp_Cond_Gauss_misp_Pois_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,l1,l2;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}

    double *dat;
    dat=(double *) Calloc(N,double);
    for(i=0;i<npairs[0];i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nugget);
                      corr1=corr_pois(corr,mui, muj);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                        M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;
                      l1=dnorm(data1[i],mui,sqrt(mui),1);
                      l2=dnorm(data2[i],muj,sqrt(muj),1);;
                      bl=2*log(dNnorm(N,M,dat))-(l1+l2);
                      *res+= bl*weights;
                    }}
   for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisGamma_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl,l1,l2,bi,bj,vvi,vvj;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}

    double *dat;
    dat=(double *) Calloc(N,double);
    for(i=0;i<npairs[0];i++){

                  if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             //***********/
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                    bi= nuis[2]/mui; bj= nuis[2]/muj;
                    vvi= mui*(1+1/bi); vvj= muj*(1+1/bj);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nugget);
                      corr1=corr_pois_gen(corr,mui, muj, nuis[2]);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                        M[0][0]=vvi; M[1][1]=vvj;M[0][1]=sqrt(vvi*vvj)*corr1;;M[1][0]= M[0][1];
                        dat[0]=data1[i]-mui;dat[1]=data2[i]-muj;
                      l1=dnorm(data1[i],mui,sqrt(vvi),1);
                      l2=dnorm(data2[i],muj,sqrt(vvj),1);;
                      bl=2*log(dNnorm(N,M,dat))-(l1+l2);
                      *res+= bl*weights;
                    }}
   for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);

    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/*********************************************************/
void Comp_Cond_PoisGamma_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l1,l2;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      l1=one_log_dpoisgamma(uu,mui,nuis[2]);
                      l2=one_log_dpoisgamma(ww,muj,nuis[2]);
                      bl=2*log(biv_PoissonGamma((1-nugget)*corr,uu,ww,mui, muj,nuis[2]))
                        - (l1+l2);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_Pois_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl,l1,l2;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      l1=dpois(uu,mui,1);
                      l2=dpois(ww,muj,1);
                      bl=2*log(biv_Poisson((1-nugget)*corr,uu,ww,mui, muj))
                        - (l1+l2);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomGauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);
                p2=pnorm(aj,0,1,1,0);
                u=data1[i];v=data2[i];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u; vv=(int) v;

                           l1=dbinom(uu,NN[0],p1,1);
                           l2=dbinom(vv,NN[0],p2,1);
                        bl=2*log(biv_binom (NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomLogi_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                u=data1[i];v=data2[i];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u; vv=(int) v;

                           l1=dbinom(uu,NN[0],p1,1);
                           l2=dbinom(vv,NN[0],p2,1);
                        bl=2*log(biv_binom (NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Cond_BinomNNGauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);
                p2=pnorm(aj,0,1,1,0);
                u=data1[i];v=data2[i];
                n1=NN[i];n2=NN[i+npairs[0]];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                            uu=(int) u; vv=(int) v;
                           l1=dbinom(uu,n1,p1,1);
                           l2=dbinom(vv,n2,p2,1);
                        bl=2*log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Cond_BinomNNLogi_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0,n1,n2;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
    if( nugget>=1 || nugget<0){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                u=data1[i];v=data2[i];
                n1=NN[i];n2=NN[i+npairs[0]];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                            uu=(int) u; vv=(int) v;
                           l1=dbinom(uu,n1,p1,1);
                           l2=dbinom(vv,n2,p2,1);
                        bl=2*log(biv_binom222(n1,n2,uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_BinomnegGauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u;  vv=(int) v;
                         l1=one_log_negbinom_marg(uu,NN[0],p1);
                         l2=one_log_negbinom_marg(vv,NN[0],p2);
                        bl=2*log(biv_binomneg(NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}



/*********************************************************/
void Comp_Cond_BinomnegLogi_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
    double nugget=nuis[0];
       if( nugget>=1 || nugget<0){*res=LOW; return;}
    //compute the composite log-likelihood:

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    p11=pblogi22(log(exp(ai)+1),log(exp(aj)+1),(1-nugget)*corr);
                p1=1/(1+exp(-ai));p2=1/(1+exp(-aj));
                    u=data1[i];v=data2[i];
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u;  vv=(int) v;
                         l1=one_log_negbinom_marg(uu,NN[0],p1);
                         l2=one_log_negbinom_marg(vv,NN[0],p2);
                        bl=2*log(biv_binomneg(NN[0],uu,vv,p1,p2,p11))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}



void Comp_Cond_TWOPIECETukeyh_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget,l1=0.0,l2=0.0;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
           zi=data1[i];zj=data2[i];


           corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

            l1=one_log_two_pieceTukey(zi,mean1[i],sill,tail,eta);
            l2=one_log_two_pieceTukey(zj,mean2[i],sill,tail,eta);

           p11=pbnorm22(qq,qq,corr);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
           bl=2*log(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]))-(l1+l2);
               *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Cond_TWOPIECET_st2memm(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,qq,l1=0.0,l2=0.0;
    double eta=nuis[3];  //skewness parameter
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                l1=one_log_two_pieceT(zi,mean1[i],sill,df,eta);
                l2=one_log_two_pieceT(zj,mean2[i],sill,df,eta);
                 p11=pbnorm22(qq,qq,corr);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                 /********************************************************/
                 bl=2*log(biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget))
                    -(l1+l2);
                 /********************************************************/
                         *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/

void Comp_Cond_TWOPIECEGauss_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,nugget,l1=0.0,l2=0.0;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

                l1=one_log_two_pieceGauss(zi,mean1[i],sill,eta);
                l2=one_log_two_pieceGauss(zj,mean2[i],sill,eta);

                p11=pbnorm22(qq,qq,corr);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    bl=2*log(biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]))-(l1+l2);
                    *res+= weights*bl;
                }}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}


void Comp_Cond_TWOPIECEBIMODAL_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,alpha,weights=1.0,p11,eta,qq,sill,df,nugget,delta ,l1=0.0,l2=0.0;

    eta=nuis[4];  //skewness parameter
    delta=nuis[3];
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
    alpha=2*(delta+1)/df;
 if( fabs(eta)>1||df<0||nugget<0||nugget>=1||delta<0||sill<0) {*res=LOW;  return;}
    qq=qnorm((1-eta)/2,0,1,1,0);
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                l1=one_log_bomidal(zi,mean1[i],sill,df,delta,eta);
                l2=one_log_bomidal(zj,mean2[i],sill,df,delta,eta);
                p11=pbnorm22(qq,qq,corr);
                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                    /********************************************************/
                  bl=2*log(biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]))-(l1+l2);

                    /********************************************************/
                           *res+= weights*bl;
                }}}

    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/********************************************************/
void Comp_Cond_BinomnegGaussZINB_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0,l1=0.0,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                  ai=mean1[i];aj=mean2[i];
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                    u=data1[i];v=data2[i];
                                  if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u;
                         vv=(int) v;
                    l1=one_log_BinomnegZIP(uu,NN[0],ai,mup);
                    l2=one_log_BinomnegZIP(vv,NN[0],aj,mup);

                    bl=2*log(biv_binomnegZINB(NN[0],corr,uu,vv,ai,aj,nugget1,nugget2,mup))-(l1+l2);
                    *res+= weights*bl;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/************************************************/
void Comp_Cond_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu,vv;
    double weights=1.0,corr,mui,muj,bl,l1=0.0,l2=0.0,u,v;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];


      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                        u=data1[i];v=data2[i];
                                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                          uu=(int) u;
                         vv=(int) v;

                    l1=one_log_PoisZIP(uu,mui,mup);
                    l2=one_log_PoisZIP(vv,muj,mup);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                     bl=2*log(biv_PoissonZIP(corr,uu,vv,mui, muj,mup,nugget1,nugget2))-(l1+l2);
                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Gauss_misp_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double weights=1.0,corr,mui,muj,bl ,l1=0.0,l2=0.0;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
    double p=pnorm(mup,0,1,1,0);

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);


                      l1=dnorm(data1[i],(1-p)*mui,sqrt(mui*(1-p)*(1+p*mui)),1);
                      l2=dnorm(data2[i],(1-p)*muj,sqrt(muj*(1-p)*(1+p*muj)),1);


                                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                      bl=2*log(biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2))-(l1+l2);

                      *res+= bl*weights;
                    }}

    if(!R_FINITE(*res))  *res = LOW;
    return;
}


void Comp_Cond_LogLogistic_st2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{

    int i;double corr,zi,zj,weights=1.0,bl=1.0,l1=0.0,l2=0.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);

                    l1=one_log_loglogistic(zi,exp(mean1[i]),nuis[2]);
                    l2=one_log_loglogistic(zj,exp(mean2[i]),nuis[2]);

                    bl=2*log(biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]))
                     -(l1+l2);
             if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Logistic_st2mem(int *cormod, double *data1,double *data2,int *NN,
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
                    l1=one_log_logistic(zi,mean1[i],nuis[1]) ;
                    l2=one_log_logistic(zj,mean2[i],nuis[1])  ;
                    bl=2*log(biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1]))
                     -(l1+l2);
                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);

  *res+= weights*bl;
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}








/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL Gaussian COPULA *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/



/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL Gaussian COPULA *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Cond_GaussGCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=1; int model=1; int cond=0;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}


    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}




void Comp_Cond_BetaGCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=1; int model=28; int cond=0;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Beta2GCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=1; int model=50; int cond=0;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                       *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_KumaraswamyGCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=1; int model=33; int cond=0;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
void Comp_Cond_Kumaraswamy2GCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=1; int model=42; int cond=0;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL Clayton COPULA *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/



// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Cond_GaussCCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=2; int model=1; int cond=1;
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){

                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                       bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                        *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
void Comp_Cond_BetaCCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=2; int model=28; int cond=1;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Cond_Beta2CCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=2; int model=50; int cond=1;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                      *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}



void Comp_Cond_KumaraswamyCCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=2; int model=33; int cond=1;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                         *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
void Comp_Cond_Kumaraswamy2CCop2mem(int *cormod, double *data1,double *data2,int *NN,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int type_cop=2; int model=42; int cond=1;
    /*############*/
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
           bl=biv_cop(corr,type_cop,cond,data1[i],data2[i],mean1[i],mean2[i],nuis,model,NN[0]);
                    *res+= bl*weights;
                    }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
