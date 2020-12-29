#include "header.h"
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_Gauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double  weights=1.0,sill,nugget,corr,bl;
    sill=nuis[1];nugget=nuis[0];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
  //Rprintf("%f %f %f  %d\n",lags[i],data1[i],data2[i],*npairs);
                      corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      bl=log_biv_Norm((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,0);
                      *res+= bl*weights;
                    }}                          
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Diff_Gauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
  int i=0;
  double vario=0.0,u,v,weights=1.0;
     
 double nugget=nuis[0];
    double sill=nuis[1];

      if(sill<0 || nugget<0||nugget>=1){*res=LOW; return;}

 for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
  vario=Variogram(cormod,lags[i],0,nuis[0],nuis[1],par);
         u=data1[i];v=data2[i];
            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
    *res+= -0.5*(log(2*M_PI)+log(vario)+
                   R_pow(u-v,2)/(2*vario))*weights;}}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_WrapGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double  u=0.0,v=0.0,weights=1.0,corr=0.0;
    double wrap_gauss;
    double alfa=2.0;
     double nugget=nuis[0];
    double sill=nuis[1];

      if(sill<0 || nugget<0||nugget>=1){*res=LOW; return;}
 for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                u=data1[i];
                v=data2[i];
                corr=CorFct(cormod,lags[i],0,par,0,0);
                wrap_gauss=biv_wrapped(alfa,u,v,mean1[i],mean2[i],nuis[0],nuis[1],corr);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    *res+=log(wrap_gauss)*weights ;
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_SinhGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;double corr,zi,zj,bb=0.0,weights=1.0;
           if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}

   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bb=log(biv_sinh((1-nuis[0])*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]));
                    *res+= weights*bb;
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_SkewGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
      double sill=nuis[1];double nugget=nuis[0];
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}

    int i=0;double corr,zi,zj,weights=1.0;
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                             if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                     *res+= weights*log(biv_skew(corr,zi,zj,mean1[i],mean2[i],nuis[1],nuis[2],nuis[0]));
                 }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gamma2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i;double corr,zi,zj,weights=1.0,bl=1.0;
    double nugget=nuis[0];
    if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    bl=biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);

                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);       
                           
  *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Weibull2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double corr,zi,zj,weights=1.0,bl=0.0;
      double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);

                    bl=biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                    //Rprintf("%f\n",bl);
                     *res+= weights*log(bl);
                
                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Kumaraswamy2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i;double corr,zi,zj,weights=1.0,bl;
    double nugget=nuis[0]; 
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);             
                  bl=biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
       
        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_Kumaraswamy22mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i;double corr,zi,zj,weights=1.0,bl;
    double nugget=nuis[0]; 
         double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);             
                  bl=biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
       
        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Beta2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i;double corr,zi,zj,weights=1.0,bl;
    double nugget=nuis[0]; 
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=(data1[i]); zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);             
                  bl=biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
       
        *res+= weights*log(bl);
                  }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_LogGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;double corr,zi,zj,weights=1.0,bl=0.0;
    double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl=d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],(1-nugget)*corr); 
                    *res+= weights*log(bl);
                    }}
    
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_PoisbinnegGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double bl,u,v,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
     double nugget=nuis[0];
     if(nugget<0||nugget>=1){*res=LOW; return;}

        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                   ai=mean1[i];aj=mean2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);

                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u;  vv=(int) v; 
                        bl=biv_poisbinneg(NN[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_PoisbinGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
 double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                   ai=mean1[i];aj=mean2[i];
                   corr=CorFct(cormod,lags[i],0,par,0,0);
                    p11=pbnorm22(ai,aj,(1-nugget)*corr);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data1[i];v=data2[i];
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u; vv=(int) v; 
                        bl=biv_poisbin(NN[0],uu,vv,p1,p2,p11);
                      
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_BinomnegGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
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
                        bl=biv_binomneg (NN[0],uu,vv,p1,p2,p11);   
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_BinomnegGaussZINB2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,  uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
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
                        bl=biv_binomnegZINB(NN[0],corr,uu,vv,ai,aj,nugget1,nugget2,mup);  
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_BinomGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,vv=0;
    double u,v,bl=0.0,weights=1.0,ai=0.0,aj=0.0,corr=0.0;
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
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          uu=(int) u; vv=(int) v; 
                        bl=biv_binom (NN[0],uu,vv,p1,p2,p11);
                    *res+= weights*log(bl);
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_LogLogistic2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i;double bl,corr,zi,zj,weights=1.0,nugget=0.0;
     nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<=2) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl=biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                    *res+= weights*log(bl);
                }}
    
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Logistic2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,nugget=0.0;
        nugget=nuis[0];
    if(nugget>=1||nugget<0.0 ) {*res=LOW;  return;}

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=(data1[i]);zj=(data2[i]);
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl= biv_Logistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[1]);
                    *res+= weights*log(bl);
                }}
    
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Pois2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl;
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
                      bl=biv_Poisson((1-nugget)*corr,uu,ww,mui, muj);
                      *res+= log(bl)*weights;
                    }}        
    
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_PoisZIP2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu,ww;
    double weights=1.0,corr,mui,muj,bl;
   double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];

      
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                    // if(fabs(corr)>1|| !R_FINITE(corr)) {*res=LOW; return;}
                        if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      uu=(int) data1[i];  ww=(int) data2[i];
                      bl=biv_PoissonZIP(corr,uu,ww,mui, muj,mup,nugget1,nugget2);
                      *res+= log(bl)*weights;
                    }}        
    
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_Pois2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,N=2;
    double  weights=1.0,corr,corr1,mui,muj,bl;
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
                      
                      bl=dNnorm(N,M,dat);
                      *res+= log(bl)*weights;
                    }}  
   for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);         
    
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_Gauss_misp_PoisZIP2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double weights=1.0,corr,mui,muj,bl;
    double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];

      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}

  // Rprintf("%d   \n",npairs[0]);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                    mui=exp(mean1[i]);muj=exp(mean2[i]);
                     corr=CorFct(cormod,lags[i],0,par,0,0);
                      if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                      bl=biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2);
                  //    Rprintf("%f %f\n",bl,mup);
                      *res+= log(bl)*weights;
                    }}        
    
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Gauss_misp_SkewT2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double weights=1.0,sill,nugget,skew,corr,corr2,df,bl;
  

    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];

    if(df<2||fabs(skew)>1||sill<0||nugget<0||nugget>=1){*res=LOW; return;}
    //auxuliary variables
    double D1=(df-1)/2;
    double D2=df/2;
    //double delta=skew/sqrt(1-skew*skew);
    double MM=sqrt(df)*gammafn(D1)*skew/(sqrt(M_PI)*gammafn(D2));
    double FF=(df/(df-2)-MM*MM);

     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=CorFct(cormod,lags[i],0,par,0,0)*(1-nugget);
                     corr2=corr_skewt(corr,df,skew);
                     if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                          bl=log_biv_Norm(corr2,data1[i],data2[i],mean1[i]+sqrt(sill)*MM,
                                                                 mean2[i]+sqrt(sill)*MM,sill*FF,0);
                        *res+= bl*weights;
                    }}            
    if(!R_FINITE(*res))  *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_T2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double weights=1.0,corr,df=0.0,bl;

     double sill=nuis[2];
    double nugget=nuis[1];

    if( sill<0||nugget<0||nugget>=1||nuis[0]<0||nuis[0]>0.5){*res=LOW; return;}
    df=1/nuis[0];
     for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){

                     corr=(1-nugget)*CorFct(cormod,lags[i],0,par,0,0);
        //if(df<170) corr=0.5*(df-2)*R_pow(gammafn((df-1)/2),2)/(R_pow(gammafn(df/2),2))* corr *hypergeo(0.5,0.5,df/2,R_pow(corr,2));
        if(fabs(corr)>0)  corr=exp(log(df-2)+2*lgammafn(0.5*(df-1))-(log(2)+2*lgammafn(df/2))+log(hypergeo(0.5,0.5, df/2,corr*corr))+log(corr*(1-nugget)));

           if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                     bl=log_biv_Norm(corr,data1[i],data2[i],mean1[i],mean2[i],sill*df/(df-2),0);
                       *res+= bl*weights;
                    }}           
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_Tukeyhh2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                   bl=biv_tukey_hh((1-nugget)*corr,zi,zj,mean1[i],mean2[i],sill,h1,h2);
                             *res+= weights*log(bl);
                }}
    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}





/*********************************************************/
void Comp_Pair_Tukeyh2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}
        for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                   bl=biv_tukey_h((1-nugget)*corr,zi,zj,mean1[i],mean2[i],tail,sill);
                             *res+= weights*log(bl);
                }}
    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gauss_misp_Tukeygh2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,corr2,zi,zj,weights=1.0,eta,tail,sill,nugget,u,eta2,mu,vv;
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
            //  Rprintf("%f %f-- %f %f \n",mean1[i],mean2[i],zi,zj);
                    *res+= weights*bl;
                }}
    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_T2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0;
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
      if( sill<0||nugget<0||nugget>=1||df<0||df>0.5){*res=LOW; return;}

   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                   bl=biv_T(corr,(zi-mean1[i])/sqrt(sill),
                                            (zj-mean2[i])/sqrt(sill),df,nugget)/sill;
                             *res+= weights*log(bl);
                }}

    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_TWOPIECETukeyh2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,tail,qq,sill,nugget;
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
                      p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]); 
                    bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]);  
                    *res+= weights*log(bl);
                }}
    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_TWOPIECET2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,qq;

    double eta=nuis[3];  //skewness parameter
    double df=nuis[0];
    double sill=nuis[2];
    double nugget=nuis[1];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}

       qq=qnorm((1-eta)/2,0,1,1,0);
      for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                      p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    /********************************************************/
                    bl=biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget);
  //Rprintf("%f- %f- %f %f %f %f %f %f %f %f %f  \n",lags[i],corr,eta,sill,zi,zj,df,eta,p11,mean1[i],mean2[i]);
                    /********************************************************/
                           *res+= weights*log(bl);
                }}
    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_TWOPIECEGauss2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{

    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,nugget;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;} 
   for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                      //p11=pbnorm(cormod,lags[i],0,qq,qq,nugget,1,par,0);

                      p11=pbnorm22(qq,qq,corr);
                       if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    bl=biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]);
                    *res+= weights*log(bl);
                }}
    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_TWOPIECEBIMODAL2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i;double bl,corr,zi,zj,weights=1.0,p11,eta,qq,sill,df,nugget,delta;

    eta=nuis[4];  //skewness parameter
    delta=nuis[3];
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
 if( fabs(eta)>1||df<0||nugget<0||nugget>=1||delta<0||sill<0) {*res=LOW;  return;} 
    
    qq=qnorm((1-eta)/2,0,1,1,0);    

    for(i=0;i<npairs[0];i++){
if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                zi=data1[i];zj=data2[i];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags[i],0,par,0,0);
                        p11=pbnorm22(qq,qq,corr);
                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0]);
                    /********************************************************/
                   bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]);
          // Rprintf("%f %f  --%f  %f %f %f %f  -%f %f \n",bl,lags[i],df,delta,eta,sill,corr,par[0],par[1]);
                    /********************************************************/
                           *res+= weights*log(bl);
                }}}
    
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATI0-TEMPORAL CASE ***********************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

void Comp_Pair_Gauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double corr,weights=1.0,bl;
    double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}


    for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
             corr=CorFct(cormod,lags[i], lagt[i],par,0,0);
             bl=log_biv_Norm((1-nugget)*corr,data1[i],data2[i],mean1[i],mean2[i],sill,0);  
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
             *res+= bl*weights;                                  
}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_WrapGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double  u=0.0, w=0.0,weights=1.0;
    double bl=0.0,corr=0.0;
    double alfa=2.0;  double nugget=nuis[0];  double sill =nuis[1];
   if(sill<0 || nugget<0||nugget>=1){*res=LOW; return;}
   //    if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                u=data1[i];//-2*atan(mean[i])-M_PI;
                                w=data2[i];//-2*atan(mean[i])-M_PI;
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    bl=biv_wrapped(alfa,u,w,mean1[i],mean2[i],nuis[0],nuis[1],(1-nugget)*corr);
                                             if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);  
                             *res+= weights*log(bl);
                                }}
             
    if(!R_FINITE(*res))*res = LOW;
    return;
} 
/******************************************************************************************/
void Comp_Pair_T_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double corr,zi,zj,weights=1.0,bl;
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
   bl=biv_T(corr,(zi-mean1[i])/sqrt(sill),(zj-mean2[i])/sqrt(sill),df,nugget)/sill;
                             *res+= weights*log(bl);
 
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_T_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double  corr,df=0.0,u=0.0, w=0.0,weights=1.0,bl;
     double sill=nuis[2];
    double nugget=nuis[1];


    if( sill<0||nugget<0||nugget>=1||nuis[0]<0||nuis[0]>0.5){*res=LOW; return;}

      df=1/nuis[0];

    for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            corr=(1-nugget)*CorFct(cormod,lags[i],lagt[i],par,0,0);
                 corr=exp(log(df-2)+2*lgammafn(0.5*(df-1))-(log(2)+2*lgammafn(df/2))+log(hypergeo(0.5,0.5, df/2,corr*corr))+log(corr*(1-nugget)));
                               u=data1[i];      
                                w=data2[i];   
                  if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
bl=log_biv_Norm(corr,u,w,mean1[i],mean2[i],sill*df/(df-2),0);
                *res+= bl*weights;   
                                    }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_Pois_stmem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, N=2;
     double weights=1.0,corr,corr2,mui,muj,bl,u=0.0, w=0.0;
   double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
    double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);} 
    double *dat;
    dat=(double *) Calloc(N,double);
     for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                          u=data1[i];      w=data2[i];
                      corr= (1-nugget) * CorFct(cormod,lags[i],lagt[i],par,0,0);
                            mui=exp(mean1[i]);muj=exp(mean2[i]);
                             corr2=corr_pois(corr,mui, muj);
                                   
                            M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr2;M[1][0]= M[0][1];
                           dat[0]=u-mui;dat[1]=w-muj;
                                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                              bl=dNnorm(N,M,dat);
                              *res+= log(bl)*weights;  
                                    }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Pois_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)

{
    int i=0,uu,ww;
     double weights=1.0,corr,mui,muj,bl,u=0.0, w=0.0;
    double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}
    // Computes the log-likelihood:
  for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                          u=data1[i];      w=data2[i];
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                     mui=exp(mean1[i]);
                     muj=exp(mean2[i]);
                          uu=(int) u;  ww=(int) w;
                      bl=biv_Poisson((1-nugget)*corr,uu,ww,mui, muj); 
                //   Rprintf("%d %d--%f %f %f  %f \n",uu,ww,lags[i],lagt[i],corr,bl);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                       *res+= log(bl)*weights;

                                    }}
              
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)

{
    int i=0,uu,ww;
     double weights=1.0,corr,mui,muj,bl,u=0.0, w=0.0;


       double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];

      
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    // Computes the log-likelihood:
  for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                          u=data1[i];      w=data2[i];
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                     mui=exp(mean1[i]);
                     muj=exp(mean2[i]);
                          uu=(int) u;  ww=(int) w;
                      bl=biv_PoissonZIP(corr,uu,ww,mui, muj,mup,nugget1,nugget2);
                //   Rprintf("%d %d--%f %f %f  %f \n",uu,ww,lags[i],lagt[i],corr,bl);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                       *res+= log(bl)*weights;

                                    }}
              
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_Pois_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,N=2;
     double weights=1.0,corr,corr2,mui,muj,bl,u=0.0, w=0.0;
   double nugget=nuis[0];
     if(nugget<0||nugget>=1){*res=LOW; return;}
// ###
    double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);} 
    double *dat;
    dat=(double *) Calloc(N,double);
// ###
    for(i=0;i<npairs[0];i++){
       u=data1[i];w=data2[i];
if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nugget);
                            mui=exp(mean1[i]);muj=exp(mean2[i]);
                            corr2=corr_pois(corr,mui, muj);
                                   
                            M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr2;M[1][0]= M[0][1];
                           dat[0]=u-mui;dat[1]=w-muj;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                              bl=dNnorm(N,M,dat);
                              *res+= log(bl)*weights;  
                }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gauss_misp_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
     double weights=1.0,corr,mui,muj,bl;
     double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];

      
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
    // Computes the log-likelihood:
  for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                     corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                     mui=exp(mean1[i]);
                     muj=exp(mean2[i]);
                         
               bl=biv_Mis_PoissonZIP(corr,data1[i],data2[i],mui, muj,mup,nugget1,nugget2);
                //   Rprintf("%d %d--%f %f %f  %f \n",uu,ww,lags[i],lagt[i],corr,bl);
                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                       *res+= log(bl)*weights;

                                    }}
              
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Tukeyh_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    
    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                    
                                   if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);

 bl=biv_tukey_h((1-nugget)*corr,zi,zj,mean1[i],mean2[i],tail,sill);
                             *res+= weights*log(bl);
 
                         }}
                
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Tukeyhh_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    
    int i=0;
    double corr,zi,zj,weights=1.0,bl;
double sill=nuis[1];
    double nugget=nuis[0];
    double h1=nuis[3];
    double h2=nuis[2];
      if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
bl=biv_tukey_hh((1-nugget)*corr,zi,zj,mean1[i],mean2[i],sill,h1,h2);
                             *res+= weights*log(bl);
                         }}      
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void  Comp_Pair_TWOPIECEGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
       int i=0;
    double corr,zi,zj,weights=1.0,p11,eta,qq,sill,bl,nugget;
   eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;} 
  qq=qnorm((1-eta)/2,0,1,1,0);

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                   if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                        p11=pbnorm22(qq,qq,corr);
  bl=biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean1[i],mean2[i]);
        
                           *res+= weights*log(bl);
                         }}
                 
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void  Comp_Pair_TWOPIECETukeyh_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
       int i=0;
    double corr,zi,zj,weights=1.0,p11,eta,qq,sill,bl,nugget,tail;
     eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

    if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;} 


     qq=qnorm((1-eta)/2,0,1,1,0);
      //   if( fabs(eta)>1  || tail<=0) {*res=LOW;  return;} 

for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                  if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                      // p11=pbnorm(cormod,lags[i],lagt[i],qq,qq,nugget,1,par,0);
                        p11=pbnorm22(qq,qq,corr);
                         bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[i],mean2[i]);  
        
                           *res+= weights*log(bl);

                         }}
                  
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void  Comp_Pair_TWOPIECET_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
       int i=0;
    double corr,zi,zj,weights=1.0,p11,qq,bl;
   double eta=nuis[3];  //skewness parameter
    double df=nuis[0];
    double sill=nuis[2];
    double nugget=nuis[1];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}
       qq=qnorm((1-eta)/2,0,1,1,0);
       //  if( fabs(eta)>1|| sill<0||df >0.5||df<0) {*res=LOW;  return;} 


for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                         if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                        p11=pbnorm22(qq,qq,corr);
                         bl=biv_two_pieceT(corr,zi,zj,sill,df,eta,p11,mean1[i],mean2[i],nugget);  
        
                           *res+= weights*log(bl);

                         }}
            
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_PoisbinGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,ww=0;
    double dens=0.0,weights=1.0,u,w,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
     double nugget=nuis[0];

      if(nugget<0||nugget>=1){*res=LOW; return;}

     for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
              corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                //psj=pbnorm(cormod,lags[i],lagt[i],mean1[i],mean2[i],nuis[0],1,par,0);
                p11=pbnorm22(mean1[i],mean2[i],(1-nugget)*corr);
                p1=pnorm(mean1[i],0,1,1,0);
                p2=pnorm(mean2[i],0,1,1,0);
                u=data1[i];w=data2[i];
                                     uu=(int) u; 
                                     ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    dens=biv_poisbin (NN[0],uu,ww,p1,p2,p11);
                                 *res+=log(dens)*weights;
                               }}
           
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_PoisbinnegGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0,uu=0,ww=0;
    double dens=0.0,weights=1.0,u,w,corr;
    double p1=0.0,p2=0.0;//probability of marginal success
    double p11=0.0;//probability of joint success
 double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                 corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                 p11=pbnorm22(mean1[i],mean2[i],(1-nugget)*corr);
                 p1=pnorm((mean1[i]),0,1,1,0);
                 p2=pnorm((mean2[i]),0,1,1,0);
                                u=data1[i];w=data2[i];
        
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    dens=biv_poisbinneg (NN[0],uu,ww,p1,p2,p11);
                                 *res+=log(dens)*weights;}}
              
    if(!R_FINITE(*res))*res = LOW;
    return;
}


void Comp_Pair_BinomGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success

      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){

   
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                  
                            a=mean1[i];b=mean2[i];
                            psj=pbnorm22(a,b,(1-nugget)*corr);
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_binom (NN[0],uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}


void Comp_Pair_BinomnegGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0,p1=0.0,p2=0.0,psj=0.0;
    
      double nugget=nuis[0];
      if(nugget<0||nugget>=1){*res=LOW; return;}

           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){
                            
                             corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                            psj=pbnorm22(a,b,(1-nugget)*corr);
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
             bl=biv_binomneg (NN[0],uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

void Comp_Pair_BinomnegGaussZINB_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0, uu=0,ww=0;
    double bl=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,corr=0.0;
    
     double nugget1=nuis[0];double nugget2=nuis[1];
    double mup=nuis[2];
      if(nugget1<0||nugget1>=1||nugget2<0||nugget2>=1){*res=LOW; return;}
           for(i=0;i<npairs[0];i++){
                   u=data1[i];w=data2[i];
             if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                            a=mean1[i];b=mean2[i];
                          uu=(int) u;  ww=(int) w;
                          if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                              bl=biv_binomnegZINB(NN[0],corr,uu,ww,a,b,nugget1,nugget2,mup); 

                                 *res+=log(bl)*weights;
                }}
    if(!R_FINITE(*res)) *res = LOW;
    return;
}




/*********************************************************************************/

// Composite marginal (difference) log-likelihood for the spatial-temporal Gaussian model:
void Comp_Diff_Gauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=00;
    double u,w,vario=0.0,weights=1.0;
    double nugget=nuis[0];
    double sill=nuis[1];
    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}

      for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                vario=Variogram(cormod,lags[i],lagt[i],nuis[0],nuis[1],par);
                                      u=data1[i];w=data2[i];
                                              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                        *res+= (-0.5*(log(2*M_PI)+log(vario)+
                                                     R_pow(u-w,2)/(2*vario)))*weights;}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_SkewGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;
    double bl,corr,zi,zj,weights=1.0;
    double nugget=nuis[0];
    double sill=nuis[1];
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_skew(corr,zi,zj,mean1[i],mean2[i],nuis[1],nuis[2],nugget);
                             *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_SinhGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;
    double bl,corr,zi,zj,weights=1.0;
       if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}
   
       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                              
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0)*(1-nuis[0]);
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                    bl=biv_sinh(corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],nuis[1]);
         
                             *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Gamma_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
{
    
    int i=0;
    double corr,zi,zj,weights=1.0,bl=1.0;
        double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
    for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                            if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                  bl=biv_gamma((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                                    *res+= weights*log(bl);
                         }}}}
  
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Kumaraswamy_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    
    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    //double sill=nuis[1];
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
  if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);      
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                   bl= biv_Kumara((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}

/******************************************************************************************/
void Comp_Pair_Kumaraswamy2_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    
    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    //double sill=nuis[1];
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
  if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);      
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                   bl= biv_Kumara2((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Beta_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    int i=0;
    double corr,zi,zj,weights=1.0,bl;
    //double sill=nuis[1];
    double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
  if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
   for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);      
                                    if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                   bl= biv_beta((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2],nuis[3],min,max);
                             *res+= weights*log(bl);
                         }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Weibull_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;
    double corr,zi,zj,weights=1.0,bl=0.0;
  double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                       
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                   // if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                      bl=biv_Weibull((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);

                                      *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_LogGauss_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;
    double corr,zi,zj,weights=1.0,bl=0.0;
     double sill=nuis[1];double nugget=nuis[0];

    if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
     for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                                zi=data1[i];
                                zj=data2[i];
                                    corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                              if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                  bl=d2lognorm(zi,zj,sill,nugget, mean1[i], mean2[i],corr);
                             *res+= weights*log(bl);
                         }}
         
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_LogLogistic_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;
    double bl,corr,zi,zj,weights=1.0,nugget=0.0;
      nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<=2) {*res=LOW;  return;}

       for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];

                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                bl=biv_LogLogistic((1-nugget)*corr,zi,zj,mean1[i],mean2[i],nuis[2]);
                                *res+= weights*log(bl);
                            }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/

void Comp_Pair_Logistic_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;
    double bl,corr=0.0,zi=0.0,zj=0.0,weights=1.0;
    if( nuis[1]<=0)  {*res=LOW;  return;}

    double sill=1-nuis[0];

          for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];
 
                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                bl=biv_Logistic(sill*corr,zi,zj,mean1[i],mean2[i],nuis[1]);
                                *res+= weights*log(bl);
                            }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}



void Comp_Pair_TWOPIECEBIMODAL_st2mem(int *cormod, double *data1,double *data2,int *NN, 
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *GPU,int *local)
{
    
    int i=0;
    double p11,qq,bl,corr,zi,zj,weights=1.0;
   double eta=nuis[4];  //skewness parameter
   double delta=nuis[3];
   double sill=nuis[2];
   double nugget=nuis[1];
   double df=nuis[0];
 if( fabs(eta)>1||df<0||nugget<0||nugget>=1||delta<0||sill<0) {*res=LOW;  return;} 

           qq=qnorm((1-eta)/2,0,1,1,0); 

          for(i=0;i<npairs[0];i++){
             if(!ISNAN(data1[i])&&!ISNAN(data2[i]) ){
                            zi=data1[i];zj=data2[i];
                                corr=CorFct(cormod,lags[i],lagt[i],par,0,0);
                                //p11=pbnorm(cormod,lags[i],lagt[i],qq,qq,nugget,1,par,0);
                                p11=pbnorm22(qq,qq,corr);
                                if(*weigthed) weights=CorFunBohman(lags[i],maxdist[0])*CorFunBohman(lagt[i],maxtime[0]);
                                bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,delta,eta,p11,mean1[i],mean2[i]);
                             *res+= weights*log(bl);
                         }}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* BIVARIATE CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

/* pairwise for bivariate GRF*/
void Comp_Pair_Gauss_biv2mem(int *cormod, double *data1,double *data2,int *NN, 
    double *par, int *weigthed,double *res,double *mean1,double *mean2,
    double *nuis, int *GPU,int *local)
{
    int i=0;
    double u=0.0, w=0.0, dens=0.0,weights=1.0;
    if(  par[0]<0|| par[1]<0|| par[2]<0|| par[3]<0) {*res=LOW;  return;} 

 for(i=0;i<npairs[0];i++){
  if(!ISNAN(u)&&!ISNAN(w) ){
          dens=log_biv2gauss(cormod,lags[i],par, data1[i]-mean1[i], data2[i]-mean2[i],first[i], second[i]);
          *res+= dens*weights;
                                }}

    if(!R_FINITE(*res))*res = LOW;
    return;
}


/* pairwise  skew for bivariate GRF*/
void Comp_Pair_SkewGauss_biv2mem(int *cormod, double *data1,double *data2,int *NN, 
    double *par, int *weigthed,double *res,double *mean1,double *mean2,
    double *nuis, int *GPU,int *local)

{
    int i=0;
    double u=0.0, w=0.0, rhotv=0.0,weights=1.0;
    int N=2;
      double *vari;vari=(double *) Calloc(N,double);vari[0]=par[0];vari[1]=par[1];  /// variances of the skew gaussian
    par[0]=1;par[1]=1;
    if(vari[0]<0||vari[1]<0){*res=LOW; return;}
      weights=1;

        for(i=0;i<npairs[0];i++){
                             u=data1[i];
                             w=data2[i];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    rhotv=CorFct(cormod,lags[i],0,par,first[i],second[i]);
                     *res+= log(biv_skew2(rhotv,u,w,vari[first[i]],vari[second[i]],1,nuis[first[i]],nuis[second[i]]))*weights;
                                }}
            
        Free(vari);
    if(!R_FINITE(*res))*res = LOW;
    return;
}

