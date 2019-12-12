#include "header.h"




// cdf of  a bivariate Gausssian distribution
void cdf_norm_call(double *lim1,double *lim2,double *a11,double *a12, double *res)
{
    *res=   cdf_norm(lim1[0],lim2[0],a11[0],a12[0]);
    //*res=   cdf_norm(0,0,.5,1);
}

void CBessel(double *xxx, double *nuu, int *expscale,double *res, int *tipo)
{
    switch (*tipo) {
        case 1:
            *res = bessel_k(*xxx,*nuu,*expscale);
            break;
        case 2:
            *res = bessel_i(*xxx,*nuu,*expscale);
        default:
            break;
}}




  /******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPACE TIME CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
// Composite marginal (pariwise) log-likelihood for the spatial-temporal Gaussian model:
void Comp_Pair_Gauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns,int *NS,int *GPU,int *local)
{
    int i=0, j=0,  t=0, v=0;
    double corr,lags=0.0, lagt=0.0;
    double  u=0.0, w=0.0,weights=1.0,bl;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr)
    //if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
   if(nuis[1]<0 || nuis[0]<0) {*res=LOW; return;}
    // Set nuisance parameters:
    double sill=nuis[1];
    double nugget=nuis[0];
    // Computes the log-likelihood:
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
           
                        if(lags<=maxdist[0]){
                            corr=CorFct(cormod,lags, 0,par,t,v);
                                u=data[(i+NS[t])];      
                                w=data[(j+NS[v])];   
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
bl=log_biv_Norm(corr,u,w,mean[(i+NS[t])],mean[(j+NS[v])],sill,nugget);
                                     //  if(!R_FINITE(bl)) { bl=1;}
                *res+= bl*weights;   
                                    }}}}
               else {
          lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]&&lagt<=maxtime[0]){
                        corr=CorFct(cormod,lags, lagt,par,t,v);
                                u=data[(i+NS[t])];    
                                w=data[(j+NS[v])];    
                             if(!ISNAN(u)&&!ISNAN(w) ){
                           bl=log_biv_Norm(corr,u,w,mean[(i+NS[t])],mean[(j+NS[v])],sill,nugget);
                *res+= bl*weights;  
                                   
                             }}}}}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
// Composite marginal (pariwise) log-likelihood for the spatial-temporal wrapped Gaussian model:
void Comp_Pair_WrapGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0,  t=0, v=0;
    double  u=0.0, w=0.0,weights=1.0;
    double bl=0.0,lags=0.0, lagt=0.0,corr=0.0;
    double alfa=2.0;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr)
    //if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
    // Set nuisance parameters:
   //    if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}

  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                                u=data[(i+NS[t])];//-2*atan(mean[(i+NS[t])])-M_PI;
                                w=data[(j+NS[v])];//-2*atan(mean[(j+NS[v])])-M_PI;
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    corr=CorFct(cormod,lags,0,par,t,v);
                                    bl=biv_wrapped(alfa,u,w,mean[(i+NS[t])],mean[(j+NS[v])],nuis[0],nuis[1],corr);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);    
                        
                             *res+= weights*log(bl);
                                }}}}
               else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                                u=data[(i+NS[t])];//-2*atan(mean[(i+NS[t])])-M_PI;
                                w=data[(j+NS[v])];//-2*atan(mean[(j+NS[v])])-M_PI;
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     corr=CorFct(cormod,lags,lagt,par,t,v);
                                    bl=biv_wrapped(alfa,u,w,mean[(i+NS[t])],mean[(j+NS[v])],nuis[0],nuis[1],corr);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);    
                               
                             *res+= weights*log(bl);
                                }}
         }}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
} 



void Comp_Pair_Gauss_misp_T_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns,int *NS,int *GPU,int *local)
{
    int i=0, j=0,  t=0, v=0;
    double corr,corr1, lags=0.0, lagt=0.0;
    double  df=0.0,u=0.0, w=0.0,weights=1.0,bl;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr)
    //if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
    if(nuis[1]<0 || nuis[2]<0 || nuis[0]>0.5 || nuis[0]<0){*res=LOW; return;}
    // Set nuisance parameters:
    df=1/nuis[0];
    // Set nuisance parameters:
    double sill=nuis[2];   
    double nugget=nuis[1];
 
    // Computes the log-likelihood:
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
           
                        if(lags<=maxdist[0]){
                            corr=CorFct(cormod,lags, 0,par,t,v);
                               if(df<170)corr1=0.5*(df-2)*R_pow(gammafn((df-1)/2),2)/(R_pow(gammafn(df/2),2))* corr *hypergeo(0.5,0.5,df/2,R_pow(corr,2));
                                else       corr1=corr;//exp(log(df-2)+2*lgammafn(0.5*(df-1))-log(2)-2*lgammafn(df/2)+log(hypergeo(0.5,0.5, df/2,R_pow(corr,2)))+log(corr));
    
                                u=data[(i+NS[t])];      
                                w=data[(j+NS[v])];   
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
bl=log_biv_Norm(corr1,u,w,mean[(i+NS[t])],mean[(j+NS[v])],sill,nugget);
                *res+= bl*weights;   
                                    }}}}
               else {
          lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]&&lagt<=maxtime[0]){
                         corr=CorFct(cormod,lags, lagt,par,t,v);
                          if(df<170) corr1=0.5*(df-2)*R_pow(gammafn((df-1)/2),2)/(R_pow(gammafn(df/2),2))* corr *hypergeo(0.5,0.5,df/2,R_pow(corr,2));
                             else       corr1=corr;//exp(log(df-2)+2*lgammafn(0.5*(df-1))-log(2)-2*lgammafn(df/2)+log(hypergeo(0.5,0.5, df/2,R_pow(corr,2)))+log(corr));
    
                                u=data[(i+NS[t])];    
                                w=data[(j+NS[v])];    
                             if(!ISNAN(u)&&!ISNAN(w) ){
                           bl=log_biv_Norm(corr1,u,w,mean[(i+NS[t])],mean[(j+NS[v])],sill,nugget);
                *res+= bl*weights;  
                                   
                             }}}}}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_Gauss_misp_Pois_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns,int *NS,int *GPU,int *local)
{
    int i=0, j=0,  t=0, v=0,N=2;
     double lags=0.0, lagt=0.0,weights=1.0,nugget,corr,corr2,mui,muj,bl,u=0.0, w=0.0;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
   //if(nuis[1]<0 || nuis[2]<0 || nuis[0]<2 || CheckCor(cormod,par)==-2){*res=LOW; return;}
   if( CheckCor(cormod,par)==-2){*res=LOW; return;} 
    double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);} 
    double *dat;
    dat=(double *) Calloc(N,double);
    nugget=nuis[0];
    // Computes the log-likelihood:
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                          if(!ISNAN(u)&&!ISNAN(w) ){
                            corr=CorFct(cormod,lags, 0,par,t,v);
                            mui=exp(mean[(i+NS[t])]);muj=exp(mean[(j+NS[v])]);
                             corr2=corr_pois(corr,mui, muj);
                                u=data[(i+NS[t])];      
                                w=data[(j+NS[v])];   
                            M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr2;M[1][0]= M[0][1];
                           dat[0]=u-mui;dat[1]=w-muj;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                              bl=dNnorm(N,M,dat);
                              *res+= log(bl)*weights;  
                                    }}}}
               else {
          lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]&&lagt<=maxtime[0]){

                           if(!ISNAN(u)&&!ISNAN(w)){
                        corr=CorFct(cormod,lags, 0,par,t,v);
                            mui=exp(mean[(i+NS[t])]);muj=exp(mean[(j+NS[v])]);
                             corr2=corr_pois(corr,mui, muj);
                                u=data[(i+NS[t])];      
                                w=data[(j+NS[v])];   

                            M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr2;M[1][0]= M[0][1];
                           dat[0]=u-mui;dat[1]=w-muj;
                               
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                              bl=dNnorm(N,M,dat);
                              *res+= log(bl)*weights;  
                                   
                             }}}}}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/*tukey h space time */
void Comp_Pair_Tukeyh_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local)
{
    
    
    int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,bl;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget>=1||nugget<0){*res=LOW; return;}

      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);

 bl=biv_tukey_h((1-nugget)*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],tail,sill);
                             *res+= weights*log(bl);
 
                         }}}}
                    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                               zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
 bl=biv_tukey_h((1-nugget)*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],tail,sill);
                             *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*studen two piece t space time */
void Comp_Pair_T_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local)
{
    int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,bl;
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];

    if( sill<0||df<0||df>0.5||nugget>=1||nugget<0){*res=LOW; return;}

      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
   bl=biv_T((1-nugget)*corr,(zi-mean[(i+NS[t])])/sqrt(sill),(zj-mean[(j+NS[v])])/sqrt(sill),df)/sill;
                             *res+= weights*log(bl);
 
                         }}}}
                    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                               
                               zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
   bl=biv_T((1-nugget)*corr,(zi-mean[(i+NS[t])])/sqrt(sill),(zj-mean[(j+NS[v])])/sqrt(sill),df)/sill;
                             *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/*skewn two piece t space time */
void  Comp_Pair_TWOPIECEGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local)
{
       int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,p11,eta,qq,sill,bl,nugget;
   eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];

     if( fabs(eta)>1|| sill<0||nugget>=1||nugget<0) {*res=LOW;  return;} 
  qq=qnorm((1-eta)/2,0,1,1,0);

      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                       p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
  bl=biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean[(i+NS[t])],mean[(j+NS[v])]);
        
                           *res+= weights*log(bl);

                         }}}}
                    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                               zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                      p11=pbnorm(cormod,lags,lagt,qq,qq,nugget,sill,par,0);
  bl=biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean[(i+NS[t])],mean[(j+NS[v])]);
                           *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/*two piece tukeyh space time */
void  Comp_Pair_TWOPIECETukeyh_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local)
{
       int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,p11,eta,qq,sill,bl,nugget,tail;
     eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];
     qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1 || sill<0 ||nugget>=1 || nugget<0 || tail<=0) {*res=LOW;  return;} 

      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                       p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
                         bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean[(i+NS[t])],mean[(j+NS[v])]);  
        
                           *res+= weights*log(bl);

                         }}}}
                    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                               zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                      p11=pbnorm(cormod,lags,lagt,qq,qq,nugget,sill,par,0);
   bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean[(i+NS[t])],mean[(j+NS[v])]);  
                           *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}



/* two piece t space time */
void  Comp_Pair_TWOPIECET_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local)
{
       int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,p11,eta,qq,sill,bl,nugget,df;
 eta=nuis[3];  //skewness parameter
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||df >0.5||df<0||nugget>=1||nugget<0) {*res=LOW;  return;} 


      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                       p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
                         bl=biv_two_pieceT((1-nugget)*corr,zi,zj,sill,df,eta,p11,mean[(i+NS[t])],mean[(j+NS[v])]);  
        
                           *res+= weights*log(bl);

                         }}}}
                    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                               zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                      p11=pbnorm(cormod,lags,lagt,qq,qq,nugget,sill,par,0);
   bl=biv_two_pieceT((1-nugget)*corr,zi,zj,sill,df,eta,p11,mean[(i+NS[t])],mean[(j+NS[v])]);  
                           *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
// Composite marginal (pariwise) log-likelihood for the poisonspatial-temporal Gaussian model:
void Comp_Pair_PoisbinGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, t=0,v=0,uu=0,ww=0;
    double dens=0.0,lags=0.0,lagt=0.0,weights=1.0,u,w;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    // Checks the validity of the mean[n+kk ],nuis[0],nuis[1]ance and correlation parameters (nugget, sill and corr)
       if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){  
                psj=pbnorm(cormod,lags,0,mean[(i+NS[t])],mean[(j+NS[v])],nuis[0],nuis[1],par,0);
                p1=pnorm((mean[(i+NS[t])]),0,1,1,0);
                p2=pnorm((mean[(j+NS[v])]),0,1,1,0);
                u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; 
                                     ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    dens=biv_poisbin (NN[0],uu,ww,p1,p2,psj);
                                 *res+=log(dens)*weights;}}}}
                else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                           psj=pbnorm(cormod,lags,lagt,mean[(i+NS[t])],mean[(j+NS[v])],nuis[0],nuis[1],par,0);
                             p1=pnorm((mean[(i+NS[t])]),0,1,1,0);
                             p2=pnorm((mean[(j+NS[v])]),0,1,1,0);
                             u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; 
                                     ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                    dens=biv_poisbin (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(dens)*weights;;
                                }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the poisonspatial-temporal Gaussian model:
void Comp_Pair_PoisbinnegGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, t=0,v=0,uu=0,ww=0;
    double dens=0.0,lags=0.0,lagt=0.0,weights=1.0,u,w;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    // Checks the validity of the mean[n+kk ],nuis[0],nuis[1]ance and correlation parameters (nugget, sill and corr)
    //if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
       if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){           
                 psj=pbnorm(cormod,lags,0,mean[(i+NS[t])],mean[(j+NS[v])],nuis[0],nuis[1],par,0);
                 p1=pnorm((mean[(i+NS[t])]),0,1,1,0);
                 p2=pnorm((mean[(j+NS[v])]),0,1,1,0);
                                u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    dens=biv_poisbinneg (NN[0],uu,ww,p1,p2,psj);
                                 *res+=log(dens)*weights;}}}}
                 else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){   
                       psj=pbnorm(cormod,lags,lagt,mean[(i+NS[t])],mean[(j+NS[v])],nuis[0],nuis[1],par,0);
                       p1=pnorm((mean[(i+NS[t])]),0,1,1,0);
                       p2=pnorm((mean[(j+NS[v])]),0,1,1,0);
                                u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                    dens=biv_poisbinneg (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(dens)*weights;;
                                }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/

// Composite conditional log-likelihood for the spatial-temporal Gaussian model:
void Comp_Cond_Gauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN, 
    double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0,  t=0, v=0;
    double s1=0.0, s12=0.0, lags=0.0, lagt=0.0,weights=1.0;
    double det=0.0, u=0.0, u2=0.0, w=0.0, w2=0.0;
      if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
   // if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
    // Set nuisance parameters:
    s1=nuis[0]+nuis[1];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                            s12=nuis[1]*CorFct(cormod,lags,0,par,t,v);
                            det=R_pow(s1,2)-R_pow(s12,2);
                             u=data[(i+NS[t])]-mean[(i+NS[t])];
                                w=data[(j+NS[v])]-mean[(j+NS[v])];
                                 if(!ISNAN(u)&&!ISNAN(w) ){
                                u2=R_pow(u,2);w2=R_pow(w,2);
                                if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                *res+= (-log(2*M_PI)-log(det)+log(s1)+
                                     (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det)*weights;}}}}
                 else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                            s12=nuis[1]*CorFct(cormod,lags,lagt,par,t,v);
                            det=R_pow(s1,2)-R_pow(s12,2);
                                 u=data[(i+NS[t])]-mean[(i+NS[t])];
                                w=data[(j+NS[v])]-mean[(j+NS[v])];
                                 if(!ISNAN(u)&&!ISNAN(w) ){
                                u2=R_pow(u,2);
                                w2=R_pow(w,2);
                                     if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                *res+= (-log(2*M_PI)-log(det)+log(s1)+
                                     (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det)*weights;}}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/*********************************************************************************/

// Composite marginal (difference) log-likelihood for the spatial-temporal Gaussian model:
void Comp_Diff_Gauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0,  t=0, v=0;
    double u,w,vario=0.0,lags=0.0, lagt=0.0,weights=1.0;
       if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    //Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
    //if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lags<=maxdist[0]){
                                vario=Variogram(cormod,lags,0,nuis[0],nuis[1],par);
                                      u=data[(i+NS[t])];w=data[(j+NS[v])];
                                       if(!ISNAN(u)&&!ISNAN(w) ){
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                        *res+= (-0.5*(log(2*M_PI)+log(vario)+
                                                     R_pow(u-w,2)/(2*vario)))*weights;}}}}
                else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lagt<=maxtime[0] && lags<=maxdist[0]){
                                vario=Variogram(cormod,lags,lagt,nuis[0],nuis[1],par);
                                     u=data[(i+NS[t])];w=data[(j+NS[v])];
                                       if(!ISNAN(u)&&!ISNAN(w) ){
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                        *res+= (-0.5*(log(2*M_PI)+log(vario)+
                                                     R_pow(u-w,2)/(2*vario)))*weights;}}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/

void Comp_Pair_SkewGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double bl,corr,zi,zj,lags,lagt,weights=1.0;
      if(nuis[1]<0 || nuis[0]<0|| nuis[0]>1) {*res=LOW;  return;}

        for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    bl=biv_skew((1-nuis[0])*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[1],nuis[2]);
                             *res+= weights*log(bl);
                         }}}}
                    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                               
                               zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
   bl=biv_skew((1-nuis[0])*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[1],nuis[2]);
                             *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/******************************************************************************************/
/******************************************************************************************/
void Comp_Pair_TWOPIECEBIMODAL_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double p11,qq,bl,corr,zi,zj,lags,lagt,weights=1.0;
    double eta=nuis[3];  //skewness parameter
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
    if( fabs(eta)>1|| sill<0||df >0.5||df<0||nugget>=1||nugget<0) {*res=LOW;  return;} 

           qq=qnorm((1-eta)/2,0,1,1,0); 

          for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lags<=maxdist[0]){
                              
                                zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                                p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                     bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,eta,p11,mean[i],mean[j]);
                             *res+= weights*log(bl);
                         }}}}
                     else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                              zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                              p11=pbnorm(cormod,lags,lagt,qq,qq,nugget,sill,par,0);
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                                   bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,eta,p11,mean[i],mean[j]);
                                   
                             *res+= weights*log(bl);
                                        
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
void Comp_Pair_SinhGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double bl,corr,zi,zj,lags,lagt,weights=1.0;
       if(nuis[3]<0||nuis[1]<0) {*res=LOW;  return;}

          for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lags<=maxdist[0]){
                              
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    bl=biv_sinh(corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2],nuis[3],nuis[1]);
         
                             *res+= weights*log(bl);
                         }}}}
                     else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                              zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                                   bl=biv_sinh(corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2],nuis[3],nuis[1]);
                                   
                             *res+= weights*log(bl);
                                        
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/

/******************************************************************************************/
/******************************************************************************************/
void Comp_Pair_Gamma_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,bl=1.0;
 double sill=1-nuis[0];
       if(nuis[2]<1||sill<0||sill>1) {*res=LOW;  return;}

    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                  bl=biv_gamma(sill*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2]);
                                 
                                    *res+= weights*log(bl);
                         }}}}
                     else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                                           bl=biv_gamma(sill*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2]);
                                            
                                            *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/



void Comp_Pair_Kumaraswamy_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local)
{
    
    
    int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,bl;
    //double sill=nuis[1];
    double nugget=nuis[0];
  if(nuis[2]<0||nuis[3]<0) {*res=LOW;  return;}

      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);      
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                   bl= biv_Kumara((1-nugget)*corr,zi,zj,mean[i],mean[j],nuis[2],nuis[3]);
                             *res+= weights*log(bl);
                         }}}}
                    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                               zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                                    bl= biv_Kumara((1-nugget)*corr,zi,zj,mean[i],mean[j],nuis[2],nuis[3]);         
                             *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
void Comp_Pair_Weibull_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,bl=0.0;
  double sill=1-nuis[0];
     if(nuis[2]<=0||sill<0||sill>1) {*res=LOW;  return;}

    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                      bl=biv_Weibull(sill*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2]);

                                      *res+= weights*log(bl);
                         }}}}
                     else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                                         bl=biv_Weibull(sill*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2]);
           
                                      *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
void Comp_Pair_LogGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double corr,zi,zj,lags,lagt,weights=1.0,bl=0.0;
    //if(  CheckCor(cormod,par)==-2)  {*res=LOW;  return;}
  if(nuis[1]<0||nuis[0]<0) {*res=LOW;  return;}
    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
         // lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,0,par,0,0);
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                  bl=d2lognorm(zi,zj,nuis[1],nuis[0], mean[(i+NS[t])], mean[(j+NS[v])],corr);
                             *res+= weights*log(bl);
                         }}}}
                     else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                            if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                                zi=data[(i+NS[t])];
                                zj=data[(j+NS[v])];
                                if(!ISNAN(zi)&&!ISNAN(zj) ){
                                    corr=CorFct(cormod,lags,lagt,par,0,0);
                                           if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
   bl=d2lognorm(zi,zj,nuis[1],nuis[0], mean[(i+NS[t])], mean[(j+NS[v])],corr);
                             *res+= weights*log(bl);
                                }}}}
                }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_BinomGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, t=0,v=0,uu=0,ww=0;
    double bl=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                            a=mean[(i+NS[t])];b=mean[(j+NS[v])];
                            psj=pbnorm(cormod,lags,0,a,b,nuis[0],nuis[1],par,0);
                            p1=pnorm(a,0,1,1,0);
                            p2=pnorm(b,0,1,1,0);
                            u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u;  ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    bl=biv_binom (NN[0],uu,ww,p1,p2,psj);
                                 *res+=log(bl)*weights;}}}}
                 else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                              a=mean[(i+NS[t])];b=mean[(j+NS[v])];
                              psj=pbnorm(cormod,lags,lagt,a,b,nuis[0],nuis[1],par,0);
                              p1=pnorm(a,0,1,1,0);p2=pnorm(b,0,1,1,0);
                              u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                    bl=biv_binom (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;;
                                }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

// Composite marginal pairwise log-likelihood for the binomial two piece  Gaussian model:
void Comp_Pair_BinomTWOPIECEGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, t=0,v=0,uu=0,ww=0;
    double bl=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,ki=0.0,kj=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    double eta=nuis[2];
    if( eta < -1 || eta > 1 || nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                            a=mean[(i+NS[t])];b=mean[(j+NS[v])];
                            ki=pnorm_two_piece(-a,eta); kj=pnorm_two_piece(-b,eta);
                            p1=  1- ki; p2=  1- kj;
                            psj=  1 + pbnorm_two_piece(cormod,lags,0,-a,-b,nuis[0],nuis[1],eta,par) - ki - kj;
                            u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u;  ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    bl=biv_binom (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;;
                               }}}}
                 else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                              a=mean[(i+NS[t])];b=mean[(j+NS[v])];

                               ki=pnorm_two_piece(-a,eta); kj=pnorm_two_piece(-b,eta);
                               p1=  1- ki; p2=  1- kj;
                               psj=  1 + pbnorm_two_piece(cormod,lags,lagt,-a,-b,nuis[0],nuis[1],eta,par) - ki - kj;
                              
                              u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                    bl=biv_binom (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;;
                                }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_Binom2Gauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, t=0,v=0,uu=0,ww=0;
    double bl=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
      int kk=nuis[2];
      if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                          a=mean[(i+NS[t])];b=mean[(j+NS[v])];
               psj=pbnorm(cormod,lags,0,a,b,nuis[0],nuis[1],par,0);
                p1=pnorm((mean[(i+NS[t])]),0,1,1,0);
                p2=pnorm((mean[(j+NS[v])]),0,1,1,0);
                                u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                   bl=biv_binom2 (NN[i],NN[j],kk,uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;;}}}}
                else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                           a=mean[(i+NS[t])];b=mean[(j+NS[v])];
                           psj=pbnorm(cormod,lags,lagt,a,b,nuis[0],nuis[1],par,0);
                            p1=pnorm(a,0,1,1,0);p2=pnorm(b,0,1,1,0);
                            u=data[(i+NS[t])];w=data[(j+NS[v])];
                               if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                   bl=biv_binom2 (NN[i],NN[j],kk,uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;
                                }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_BinomnegGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, t=0,v=0,uu=0,ww=0;
    double bl=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
       if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
              a=mean[(i+NS[t])];b=mean[(j+NS[v])];
              psj=pbnorm(cormod,lags,0,a,b,nuis[0],nuis[1],par,0);
              p1=pnorm(a,0,1,1,0);p2=pnorm(b,0,1,1,0);
          u=data[(i+NS[t])]; w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    bl=biv_binomneg (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;
                                }}}}
                 else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                              a=mean[(i+NS[t])];b=mean[(j+NS[v])];
                              psj=pbnorm(cormod,lags,lagt,a,b,nuis[0],nuis[1],par,0);
                              p1=pnorm(a,0,1,1,0); p2=pnorm(b,0,1,1,0);
                                u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                    bl=biv_binomneg (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;
                                }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_BinomnegTWOPIECEGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, t=0,v=0,uu=0,ww=0;
    double bl=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0,a=0.0,b=0.0,ki=0.0,kj=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
      double eta=nuis[2];
    if( eta < -1 || eta > 1 || nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
              a=mean[(i+NS[t])];b=mean[(j+NS[v])];

               ki=pnorm_two_piece(-a,eta); kj=pnorm_two_piece(-b,eta);
               p1=  1- ki; p2=  1- kj;
               psj=  1 + pbnorm_two_piece(cormod,lags,0,-a,-b,nuis[0],nuis[1],eta,par) - ki - kj;
                       
          u=data[(i+NS[t])]; w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                    bl=biv_binomneg (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;
                                }}}}
                 else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] && lags<=maxdist[0]){
                              a=mean[(i+NS[t])];b=mean[(j+NS[v])];

                               ki=pnorm_two_piece(-a,eta); kj=pnorm_two_piece(-b,eta);
                               p1=  1- ki; p2=  1- kj;
                               psj=  1 + pbnorm_two_piece(cormod,lags,lagt,-a,-b,nuis[0],nuis[1],eta,par) - ki - kj;
                      
                                u=data[(i+NS[t])];w=data[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                     uu=(int) u; ww=(int) w;
                                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
                                    bl=biv_binomneg (NN[0],uu,ww,p1,p2,psj);
                                   *res+=log(bl)*weights;
                                }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

void Comp_Pair_LogLogistic_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                               int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double bl,corr,zi,zj,lags,lagt,weights=1.0,nugget=0.0;
    if( nuis[2]<=2)  {*res=LOW;  return;}
      nugget=nuis[0];

      for(t=0;t<ntime[0];t++){
      for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                            zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                            if(!ISNAN(zi)&&!ISNAN(zj) ){
                                corr=CorFct(cormod,lags,0,par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                bl=biv_LogLogistic((1-nugget)*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2]);
                                *res+= weights*log(bl);
                            }}}}
                 else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                            zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                            if(!ISNAN(zi)&&!ISNAN(zj) ){
                                corr=CorFct(cormod,lags,lagt,par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                                 bl=biv_LogLogistic((1-nugget)*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[2]);
                                *res+= weights*log(bl);
                            }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


void Comp_Pair_Logistic_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                            int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0,t=0,v=0;
    double bl,corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0;
    if( nuis[1]<=0)  {*res=LOW;  return;}

    double sill=1-nuis[0];

    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                            zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                            if(!ISNAN(zi)&&!ISNAN(zj) ){
                                corr=CorFct(cormod,lags,0,par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                                bl=biv_Logistic(sill*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[1]);
                                *res+= weights*log(bl);
                            }}}}
                else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lagt<=maxtime[0] &&lags<=maxdist[0]){
                            zi=data[(i+NS[t])];zj=data[(j+NS[v])];
                            if(!ISNAN(zi)&&!ISNAN(zj) ){
                                corr=CorFct(cormod,lags,lagt,par,0,0);
                                if(*weigthed) weights=CorFunBohman(lags,maxdist[0])*CorFunBohman(lags,maxdist[0]);
                                    bl=biv_Logistic(sill*corr,zi,zj,mean[(i+NS[t])],mean[(j+NS[v])],nuis[1]);
                                *res+= weights*log(bl);
                            }}}}
            }}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
// Composite conditional log-likelihood for the spatial Gaussian model:
void Comp_Cond_Gauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, 
      double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0;
    double s1=0.0, s12=0.0, lags=0.0,weights=1.0;
    double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
    // Checks the validity of the the nuisance and correlation parameters (nugget, sill and corr):
   // if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
       if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    // Set nuisance parameters:
    s1=nuis[0]+nuis[1];//set nugget + sill

    for(i=0; i<(ncoord[0]-1);i++)
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                s12=nuis[1]*CorFct(cormod, lags, 0, par,0,0); //sill * corr
                det=R_pow(s1,2)-R_pow(s12,2);
                    u=data[i]-mean[i]; 
                    v=data[j]-mean[j]; 
                    if(!ISNAN(u)&&!ISNAN(v) ){
                        u2=R_pow(u,2);v2=R_pow(v,2);
                        if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                        *res+= (-log(2*M_PI)-log(det)+log(s1)+
                        (u2+v2)*(0.5/s1-s1/det)+2*s12*u*v/det)*weights;
                    }}}
    // Checks the return values
    if(!R_FINITE(*res))*res = LOW;
    return;
}


void Comp_Diff_Gauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
  int i=0, j=0;
  double vario=0.0,lags=0.0,u,v,weights=1.0;
     if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
  // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
  //if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){ *res=LOW; return;}

  for(i=0; i<(ncoord[0]-1);i++){
    for(j=(i+1);j<ncoord[0];j++){
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
      if(lags<=*maxdist){
	vario=Variogram(cormod,lags,0,nuis[0],nuis[1],par);
	       u=data[i];v=data[j];
        if(!ISNAN(u)&&!ISNAN(v) ){
            if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
	  *res+= -0.5*(log(2*M_PI)+log(vario)+
                   R_pow(u-v,2)/(2*vario))*weights;}}}}
  if(!R_FINITE(*res)|| !*res)
    *res = LOW;
  return;
}

/*  for maximum likelihood  in R^1
void Comp_Diff_Gauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
  double lags=0.0,u,v,r1=0.0,weights=0.0,pp=0.0,qq=0.0;int k;

weights=data[0]*data[0]/nuis[1];
  for(k=1; k<(ncoord[0]-1);k++){
        //lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
        lags=fabs(coordx[k]-coordx[k-1]);
    //vario=Variogram(cormod,lags,0,nuis[0],nuis[1],par);
           u=data[k];v=data[k-1];
           pp=lags/par[0];
           qq=1-exp(-2*pp);
      r1+= R_pow(u-v*exp(-pp),2)/(nuis[1]*(qq) )
        + log(nuis[1]*(qq));
    }
  
   *res=-(r1+weights+log(nuis[1])+ncoord[0]*log(2*M_PI));
  return;
}

*/

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_WrapGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0;
    double  u=0.0,v=0.0,weights=1.0,corr=0.0;
    double wrap_gauss;
    double lags=0.0,alfa=2.0;
       // nuis[0]=0.0;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
   // if(nuis[1]<0 || nuis[0]<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
       if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    // Set nuisance parameters:
 //set nugget + sill
    for(i=0;i<(ncoord[0]-1);i++)
    for(j=(i+1); j<ncoord[0];j++){
        // Pairwise distances
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
        if(lags<=maxdist[0]){
                u=data[i];
                v=data[j];
                if(!ISNAN(u)&&!ISNAN(v) ){
                corr=CorFct(cormod,lags,0,par,0,0);
                wrap_gauss=biv_wrapped(alfa,u,v,mean[i],mean[j],nuis[0],nuis[1],corr);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    *res+=log(wrap_gauss)*weights ;
                }}}
    // Checks the return values
    if(!R_FINITE(*res))*res = LOW;
    return;
}



void Comp_Pair_SinhGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0;double corr,zi,zj,lags,bb=0.0,weights=1.0;
      if(nuis[3]<0||nuis[1]<0) {*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
            for(j=(i+1); j<ncoord[0];j++){
                lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                 if(lags<=maxdist[0]){
                    zi=data[i];zj=data[j];
                      if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    bb=log(biv_sinh(corr,zi,zj,mean[i],mean[j],nuis[2],nuis[3],nuis[1]));
                    *res+= weights*bb;
                 }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}

void Comp_Pair_SkewGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0;double corr,zi,zj,lags,weights=1.0;
         if(nuis[1]<0 || nuis[0]<0|| nuis[0]>1){*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
            for(j=(i+1); j<ncoord[0];j++){
                lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                 if(lags<=maxdist[0]){
                    zi=data[i];zj=data[j];
                      if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                          if(*weigthed) {weights=CorFunBohman(lags,maxdist[0]);}
                     *res+= weights*log(biv_skew((1-nuis[0])*corr,zi,zj,mean[i],mean[j],nuis[1],nuis[2]));
                 }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_2Gamma2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0;double kk,corr,zi,zj,lags,weights=1.0;
    //if(nuis[2]<2 || nuis[3]<=0||  nuis[4]<0
     // || CheckCor(cormod,par)==-2) {*res=LOW;  return;}
  //if(nuis[1]<0 || nuis[0],<=0 || CheckCor(cormod,par)==-2){*res=LOW; return;}
    for(i=0;i<(ncoord[0]-1);i++){
            for(j=(i+1); j<ncoord[0];j++){
                 lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                 if(lags<=maxdist[0]){
                    zi=(data[i]-mean[i]);zj=(data[j]-mean[j]);
                      if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                  kk=biv_gamma2(corr,zi/sqrt(nuis[1]),zj/sqrt(nuis[1]),nuis[2],nuis[3],nuis[4]);
                  *res+= weights*log(kk);
                  }}}}
    if(!R_FINITE(*res)||!res) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Gamma2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i,j;double corr,zi,zj,lags,weights=1.0,bl=1.0;
    double sill=1-nuis[0];
     if(nuis[2]<1||sill<0||sill>1) {*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
            for(j=(i+1); j<ncoord[0];j++){
                lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                 if(lags<=maxdist[0]){
                    zi=(data[i]); zj=(data[j]);
                      if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    bl=biv_gamma(sill*corr,zi,zj,mean[i],mean[j],nuis[2]);
                     if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);        
                           
  *res+= weights*log(bl);
                  }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}

void Comp_Cond_Gamma2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
  double lags=0.0,a1,a2;int k=0;   int nn=ncoord[0];
  double PP1=0.0,PP2=0.0,MM=0.0,sum=0.0,cc=0.0,second=0.0,first=0.0;

  double m=nuis[2];
      if(m<0||par[0]<0){*res=LOW;  return;}

  double x1 = data[0]/exp(mean[0]);
  double xn = data[nn-1]/exp(mean[nn-1]);
  double lags01=fabs(coordx[0]-coordx[1]);
  double lagsn1n=fabs(coordx[nn-2]-coordx[nn-1]);
  /********************************/
  for(k=1; k<(nn-1);k++){
              lags=fabs(coordx[k-1]-coordx[k]); a1=exp(-lags/par[0]);
              lags=fabs(coordx[k]-coordx[k+1]); a2=exp(-lags/par[0]);
      sum= sum+ (data[k]/exp(mean[k]))*(1-R_pow(a1*a2,2))/((1-a1*a1)*(1-a2*a2));
    }
    second=-(m/2)*( x1/(1-R_pow(exp(-lags01 /par[0]),2))  
                   +xn/(1-R_pow(exp(-lagsn1n/par[0]),2))
                         +sum);
/********************************/

  for(k=0; k<nn;k++) MM=MM-mean[k];
/*******+********************************************/
  for(k=0; k<(nn-1);k++){
        lags=fabs(coordx[k]-coordx[k+1]);
           cc=exp(-lags/par[0]);
            PP1=PP1 + log((1-cc*cc)*R_pow(fabs(cc),m/2-1));
            PP2=PP2 + log(bessel_i( m*fabs(cc)*sqrt(
                                 (data[k]/exp(mean[k]))*
                                 (data[k+1]/exp(mean[k+1])))/(1-cc*cc) , m/2-1 ,1));
    }
  first= (m/2-1+nn)*log(m/2) + (m/4-1/2)*log(x1*xn) - ( log(gammafn(m/2)) + PP1);

  *res=MM+first+second+ PP2;
    if(!R_FINITE(*res)|| !*res) *res = LOW;
  return;
}

/*********************************************************/
void Comp_Pair_Kumaraswamy2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i,j;double corr,zi,zj,lags,weights=1.0,bl;
    double sill=1-nuis[0]; 
     if(nuis[2]<0||nuis[3]<0) {*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
            for(j=(i+1); j<ncoord[0];j++){
                lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                 if(lags<=maxdist[0]){
                    zi=(data[i]); zj=(data[j]);
                      if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                     if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);             
                  bl=biv_Kumara(sill*corr,zi,zj,mean[i],mean[j],nuis[2],nuis[3]);
       
        *res+= weights*log(bl);
                  }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}
/*********************************************************/
void Comp_Pair_Weibull2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i,j;double corr,zi,zj,lags,weights=1.0,bl=0.0;
      double sill=1-nuis[0];
     if(nuis[2]<=0||sill<0||sill>1) {*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=(data[i]);zj=(data[j]);
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    bl=biv_Weibull(sill*corr,zi,zj,mean[i],mean[j],nuis[2]);

                     *res+= weights*log(bl);
                
                }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}

void Comp_Cond_Weibull2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
  double lags=0.0,a1,a2;int k=0;   int nn=ncoord[0];
  double PP0=0.0, PP1=0.0,PP2=0.0,MM=0.0,sum=0.0,cc=0.0,second=0.0,first=0.0;

  double m=nuis[2]; 
  if(m<0||par[0]<0){*res=LOW;  return;}
  double nu=1/gammafn(1+1/m);
  double x1 = data[0]/exp(mean[0]);
  double xn = data[nn-1]/exp(mean[nn-1]);
  double lags01=fabs(coordx[0]-coordx[1]);
  double lagsn1n=fabs(coordx[nn-2]-coordx[nn-1]);

  /********************************/
  for(k=1; k<(nn-1);k++){
              lags=fabs(coordx[k-1]-coordx[k]); a1=exp(-lags/par[0]);
              lags=fabs(coordx[k]-coordx[k+1]); a2=exp(-lags/par[0]);
      sum= sum+ R_pow(data[k]/exp(mean[k]),m)*(1-R_pow(a1*a2,2))/((1-a1*a1)*(1-a2*a2));
    }
    second=-R_pow(nu,-m)*( R_pow(x1,m)/(1-R_pow(exp(-lags01 /par[0]),2))  
                          +R_pow(xn,m)/(1-R_pow(exp(-lagsn1n/par[0]),2))
                          +sum);
/********************************/
  for(k=0; k<nn;k++) {MM=MM-mean[k];PP0=PP0+ (m-1)*log(data[k]/exp(mean[k]));}
/*******+********************************************/
  for(k=0; k<(nn-1);k++){
        lags=fabs(coordx[k]-coordx[k+1]);
           cc=exp(-lags/par[0]);
            PP1=PP1 + log(1-cc*cc);
            PP2=PP2 + log(bessel_i( 2*fabs(cc)*R_pow(
                                 (data[k]/exp(mean[k]))*
                                 (data[k+1]/exp(mean[k+1])),m/2)/((1-cc*cc)*R_pow(nu,m)) , 0 ,1));
    }
  first= nn*log(m) + PP0 - ( nn*m*log(nu) + PP1);
  *res=MM+first+second+ PP2;
    if(!R_FINITE(*res)|| !*res) *res = LOW;
  return;
}

void Comp_Pair_LogGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
    int i=0,j=0;double corr,zi,zj,lags,weights=1.0,bl=0.0;
      // if(nuis[1]<0 || nuis[2]<0|| CheckCor(cormod,par)==-2) {*res=LOW;  return;}
  if(nuis[1]<0||nuis[0]<0) {*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
            for(j=(i+1); j<ncoord[0];j++){
                lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                 if(lags<=maxdist[0]){
                    zi=data[i];zj=data[j];
                      if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    bl=d2lognorm(zi,zj,nuis[1],nuis[0], mean[i], mean[j],corr);
                     ///   
                    *res+= weights*log(bl);
                    }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_Gauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0;
    double lags=0.0, weights=1.0,sill,nugget,corr,bl;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
   if(nuis[1]<0 || nuis[0]<0){*res=LOW; return;}
    // Set nuisance parameters:
    sill=nuis[1];nugget=nuis[0];
    //if(sill<0||nugget<0) {*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
 			lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  if(!ISNAN(data[i])&&!ISNAN(data[j]) ){
                     corr=CorFct(cormod,lags,0,par,0,0);
                      if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                      bl=log_biv_Norm(corr,data[i],data[j],mean[i],mean[j],sill,nugget);
                        *res+= bl*weights;
                    }}}}            
    // Checks the return values
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

void Comp_Pair_PoisbinnegGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, uu=0,vv=0;
    double bl,u,v,lags=0.0,weights=1.0,ai=0.0,aj=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
       if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    for(i=0; i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                   ai=mean[i];aj=mean[j];
                psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data[i];v=data[j];
                    if(!ISNAN(u)&&!ISNAN(v) ){
                        if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          uu=(int) u;  vv=(int) v; 
                        bl=biv_poisbinneg(NN[0],uu,vv,p1,p2,psj);
                           
                    *res+= weights*log(bl);
                }}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

// Composite marginal pairwise log-likelihood for the binomial spatial Gaussian model:
void Comp_Pair_PoisbinGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, uu=0,vv=0;
    double u,v,bl=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
       if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    //compute the composite log-likelihood:
    for(i=0; i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                   ai=mean[i];aj=mean[j];
                psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data[i];v=data[j];
                    if(!ISNAN(u)&&!ISNAN(v) ){
                        if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          uu=(int) u; vv=(int) v; 
                        bl=biv_poisbin(NN[0],uu,vv,p1,p2,psj);
                      
                    *res+= weights*log(bl);
                }}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
// Composite marginal pairwise log-likelihood for the binomial spatial Gaussian model:
void Comp_Pair_BinomnegGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, uu=0,vv=0;
    double u,v,bl=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
       if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    //compute the composite log-likelihood:

    for(i=0; i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  ai=mean[i];aj=mean[j];
                psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data[i];v=data[j];
                    if(!ISNAN(u)&&!ISNAN(v) ){
                         if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          uu=(int) u; 
                         vv=(int) v; 
                        bl=biv_binomneg (NN[0],uu,vv,p1,p2,psj);
                           
                    *res+= weights*log(bl);
                }}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

// Composite marginal pairwise log-likelihood for the binomial spatial Gaussian model:
void Comp_Pair_BinomGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, uu=0,vv=0;
    double u,v,bl=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    //compute the composite log-likelihood:
    for(i=0; i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                 ai=mean[i];aj=mean[j];
                psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
                p1=pnorm(ai,0,1,1,0);
                p2=pnorm(aj,0,1,1,0);
                u=data[i];v=data[j];
                    if(!ISNAN(u)&&!ISNAN(v) ){
                        if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          uu=(int) u; vv=(int) v; 
                        bl=biv_binom (NN[0],uu,vv,p1,p2,psj);
                           
                    *res+= weights*log(bl);
                }}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

// Composite marginal pairwise log-likelihood for the binomial two piece spatial Gaussian model:
void Comp_Pair_BinomTWOPIECEGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, uu=0,vv=0;
    double u,v,bl=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0,ki=0.0,kj=0.0;
    double pi=0.0,pj=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    double eta=nuis[2];
    if( eta < -1 || eta > 1 || nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0]; // nuis[0] is the nugget
    //compute the composite log-likelihood:
    for(i=0; i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                 ai=mean[i];aj=mean[j];
  ki=pnorm_two_piece(-ai,eta); kj=pnorm_two_piece(-aj,eta);
  pi=  1- ki; pj=  1- kj;
  psj=  1 + pbnorm_two_piece(cormod,lags,0,-ai,-aj,nuis[0],nuis[1],eta,par) - ki - kj;
  u=data[i];v=data[j];
                if(!ISNAN(u)&&!ISNAN(v) ){
                        if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          uu=(int) u; vv=(int) v; 
                        bl=biv_binom (NN[0],uu,vv,pi,pj,psj);  
                     
                    *res+= weights*log(bl);
                }}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

// Composite marginal pairwise log-likelihood for the binomial spatial Gaussian model:
void Comp_Pair_BinomnegTWOPIECEGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, uu=0,vv=0;
    double u,v,bl=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0,ki=0.0,kj=0.0;
    double pi=0.0,pj=0.0,psj=0.0;//probability of marginal success
    double eta=nuis[2];
    if( eta < -1 || eta > 1 || nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    //compute the composite log-likelihood:
    for(i=0; i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  ai=mean[i];aj=mean[j];
  ki=pnorm_two_piece(-ai,eta); kj=pnorm_two_piece(-aj,eta);
  pi=  1- ki; pj=  1- kj;
  psj=  1 + pbnorm_two_piece(cormod,lags,0,-ai,-aj,nuis[0],nuis[1],eta,par) - ki - kj;
                    u=data[i];v=data[j];
                    if(!ISNAN(u)&&!ISNAN(v) ){
                         if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          uu=(int) u; 
                         vv=(int) v; 
                        bl=biv_binomneg (NN[0],uu,vv,pi,pj,psj);
                               
                    *res+= weights*log(bl);
                }}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}
// Composite marginal pairwise log-likelihood for the binomial spatial Gaussian model:
void Comp_Pair_Binom2Gauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0, uu=0,vv=0;
    double u,v,dens=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    int kk=nuis[2];
        if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    nuis[1]=1-nuis[0];
    for(i=0; i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  ai=mean[i];aj=mean[j];
                psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
                p1=pnorm(ai,0,1,1,0);p2=pnorm(aj,0,1,1,0);
                    u=data[i];v=data[j];
                    if(!ISNAN(u)&&!ISNAN(v) ){
                        if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          uu=(int) u; vv=(int) v; 
                        dens=biv_binom2 (NN[i],NN[j],kk,uu,vv,p1,p2,psj);
                         *res+=log(dens)*weights;
                     
                }}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_LogLogistic2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                            int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    
  
    int i,j;double bl,corr,zi,zj,lags,weights=1.0,nugget=0.0;
    nugget=nuis[0];
    if( nuis[2]<=2||nugget>=1||nugget<0.0 ) {*res=LOW;  return;}
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=data[i];zj=data[j];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    bl=biv_LogLogistic((1-nugget)*corr,zi,zj,mean[i],mean[j],nuis[2]);
                    *res+= weights*log(bl);
                }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_Logistic2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i,j;double bl,corr,zi,zj,lags,weights=1.0,nugget=0.0;
        nugget=nuis[0];
    if(nugget>=1||nugget<0.0 ) {*res=LOW;  return;}

    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=(data[i]);zj=(data[j]);
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    bl= biv_Logistic((1-nugget)*corr,zi,zj,mean[i],mean[j],nuis[1]);
                    *res+= weights*log(bl);
                }}}}
    // Checks the return values
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;
}


// Composite marginal (pariwise) log-likelihood for poisson model
void Comp_Pair_Pois2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0,uu,ww;
    double lags=0.0, weights=1.0,nugget,corr,mui,muj,bl;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
   //if(nuis[1]<0 || nuis[2]<0 || nuis[0]<2 ){*res=LOW; return;}
   //if( CheckCor(cormod,par)==-2){*res=LOW; return;} 
    nugget=nuis[0];
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){

      lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  if(!ISNAN(data[i])&&!ISNAN(data[j]) ){
             //***********/
                    mui=exp(mean[i]);muj=exp(mean[j]);

                     corr=CorFct(cormod,lags,0,par,0,0);

                      if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                      uu=(int) data[i];  ww=(int) data[j];
                      bl=biv_Poisson(corr,uu,ww,mui, muj);
                      Rprintf("%f %f %f %f \n",log(bl),corr,data[i],data[j]);

                     // Rprintf("%f %f %f \n",bl,corr,corr1);
                      *res+= log(bl)*weights;
                    }}}}          
    // Checks the return values
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


// Composite marginal (pariwise) log-likelihood for the spatial  Gaussian misspecification model:
void Comp_Pair_Gauss_misp_Pois2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0,N=2;
    double lags=0.0, weights=1.0,nugget,corr,corr1,mui,muj,bl;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
   //if(nuis[1]<0 || nuis[2]<0 || nuis[0]<2 ){*res=LOW; return;}
   //if( CheckCor(cormod,par)==-2){*res=LOW; return;} 
double **M;
        M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}
    
    double *dat;
    dat=(double *) Calloc(N,double);
    nugget=nuis[0];
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
      lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  if(!ISNAN(data[i])&&!ISNAN(data[j]) ){
             //***********/
                    mui=exp(mean[i]);muj=exp(mean[j]);

                     corr=CorFct(cormod,lags,0,par,0,0);
                    // corr2=corr*corr;
                    // z=2*sqrt(mui*muj)/(1-corr2);
                     
                    // if(z<700)corr1=corr2*(1-exp(-z)*(bessel_i(z,0,1)+bessel_i(z,1,1)));
                    // else corr1=corr2*(1-(sqrt(2/(M_PI*z)))*(1+1/(8*z))); /// approx for large z
                      corr1=corr_pois(corr,mui, muj);
                      if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                        M[0][0]=mui; M[1][1]=muj;M[0][1]=sqrt(mui*muj)*corr1;M[1][0]= M[0][1];
                        dat[0]=data[i]-mui;dat[1]=data[j]-muj;
                      
                      bl=dNnorm(N,M,dat);
                     // Rprintf("%f %f %f \n",bl,corr,corr1);
                      *res+= log(bl)*weights;
                    }}}}   
   for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);         
    // Checks the return values
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

// Composite marginal (pariwise) log-likelihood for the spatial  Gaussian misspecification model:
void Comp_Pair_Gauss_misp_SkewT2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0;
    double lags=0.0, weights=1.0,sill,nugget,skew,corr,corr2,corr3,df,bl;
    if(nuis[1]<0 || nuis[2]<0 || nuis[0]>0.5 || nuis[0]<0||fabs(nuis[3])>1){*res=LOW; return;}

    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];
    //auxuliary variables
    double D1=(df-1)/2;
    double D2=df/2;
    double skew2=R_pow(skew,2);
    double w=1-skew2;
    double KK=2*skew2/M_PI;

    ///double MM=sqrt(df)*gammafn(D1)*skew/(sqrt(M_PI)*gammafn(D2));
   // double FF=sill*(df/(df-2)-MM*MM);

    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
      lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  if(!ISNAN(data[i])&&!ISNAN(data[j]) ){
                     corr=CorFct(cormod,lags,0,par,0,0);
         corr2=(2*skew2/(M_PI*w+skew2*(M_PI-2)))*(sqrt(1-corr*corr)+corr*asin(corr)-1)+w*corr/(w+skew2*(1-2/M_PI));
corr3=(M_PI*(df-2)*R_pow(gammafn(D1),2)/(2*(M_PI*R_pow(gammafn(D2),2)-skew2*(df-2)*R_pow(gammafn(D1),2))))*
(hypergeo(0.5,0.5,D2,R_pow(corr,2))*((1-KK)*corr2+KK)-KK);

                      if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                          bl=log_biv_Norm(corr3,data[i],data[j],mean[i],mean[j],sill,nugget);
                        *res+= bl*weights;
                    }}}}            
    if(!R_FINITE(*res))  *res = LOW;
    return;
}

// Composite marginal (pariwise) log-likelihood for the spatial  Gaussian misspecification model:
void Comp_Pair_Gauss_misp_T2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0;
    double lags=0.0, weights=1.0,corr,df=0.0,bl,sill;
    // Checks the validity of the nuisance and correlation parameters (nugget, sill and corr):
    if(nuis[1]<0 || nuis[2]<0 || nuis[0]>0.5 || nuis[0]<0){*res=LOW; return;}
    // Set nuisance parameters:
    df=1/nuis[0];
    sill=nuis[2];
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
      lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                  if(!ISNAN(data[i])&&!ISNAN(data[j]) ){
                     corr=CorFct(cormod,lags,0,par,0,0);
        if(df<170) corr=0.5*(df-2)*R_pow(gammafn((df-1)/2),2)/(R_pow(gammafn(df/2),2))* corr *hypergeo(0.5,0.5,df/2,R_pow(corr,2));
       // else     corr=exp(log(df-2)+2*lgammafn(0.5*(df-1))-log(2)-2*lgammafn(df/2)+log(hypergeo(0.5,0.5, df/2,R_pow(corr,2)))+log(corr));
         if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                     bl=log_biv_Norm(corr,data[i],data[j],mean[i],mean[j],sill,nuis[1]);
                       *res+= bl*weights;
                    }}}}            
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


/* bivariate spatial tukey */
void Comp_Pair_Tukeyh2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i,j;double bl,corr,zi,zj,lags,weights=1.0;
    double sill=nuis[1];
    double nugget=nuis[0];
    double tail=nuis[2];
      if( sill<0||tail<0||tail>0.5||nugget>=1||nugget<0){*res=LOW; return;}
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=data[i];zj=data[j];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                   bl=biv_tukey_h((1-nugget)*corr,zi,zj,mean[i],mean[j],tail,sill);
                             *res+= weights*log(bl);
                }}}}

    // Checks the return values
    if(!R_FINITE(*res)) *res = LOW;
    return;
}

/*********************************************************/
void Comp_Pair_T2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i,j;double bl,corr,zi,zj,lags,weights=1.0;
    double sill=nuis[2];
    double nugget=nuis[1];
    double df=nuis[0];
      if( sill<0||df<0||df>0.5||nugget>=1||nugget<0){*res=LOW; return;}
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=data[i];zj=data[j];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                   bl=biv_T((1-nugget)*corr,(zi-mean[i])/sqrt(sill),
                                            (zj-mean[j])/sqrt(sill),df)/sill;
                             *res+= weights*log(bl);
                }}}}

    // Checks the return values
    if(!R_FINITE(*res)) *res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_TWOPIECETukeyh2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,
                         double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i,j;double bl,corr,zi,zj,lags,weights=1.0,p11,eta,tail,qq,sill,nugget;
    eta  = nuis[2];  //skewness parameter
    tail = nuis[3];  //tail parameter
    sill =nuis[1];
    nugget=nuis[0];

       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1 || sill<0 ||nugget>=1 || nugget<0 || tail<=0) {*res=LOW;  return;} 
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=data[i];zj=data[j];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                      p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]); 
                    bl=biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean[i],mean[j]);  
                    *res+= weights*log(bl);
                }}}}
    // Checks the return values
    if(!R_FINITE(*res)) *res = LOW;
    return;
}




/*********************************************************/
void Comp_Pair_TWOPIECEBIMODAL2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,
                         double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i,j;double bl,corr,zi,zj,lags,weights=1.0,p11,eta,qq,sill,df,nugget;

    eta=nuis[3];  //skewness parameter
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
    if( fabs(eta)>1|| sill<0||df >0.5||df<0||nugget>=1||nugget<0) {*res=LOW;  return;} 
    
    qq=qnorm((1-eta)/2,0,1,1,0);         
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=data[i];zj=data[j];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                      p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    /********************************************************/
                   bl=biv_two_piece_bimodal((1-nugget)*corr,zi,zj,sill,df,eta,p11,mean[i],mean[j]);
                   ///Rprintf("%f %f %f %f %f %f -- %f %f \n",zi,zj,bl,eta,sill,mean[i],corr,p11);
                    /********************************************************/
                           *res+= weights*log(bl);
                }}}}
    // Checks the return values
    if(!R_FINITE(*res)) *res = LOW;
    return;
}



/*********************************************************/
void Comp_Pair_TWOPIECET2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,
                         double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i,j;double bl,corr,zi,zj,lags,weights=1.0,p11,eta,qq,sill,df,nugget;

    eta=nuis[3];  //skewness parameter
    sill=nuis[2];
    nugget=nuis[1];
    df=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||df >0.5||df<0||nugget>=1||nugget<0) {*res=LOW;  return;} 
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=data[i];zj=data[j];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                      p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    /********************************************************/
                    bl=biv_two_pieceT((1-nugget)*corr,zi,zj,sill,df,eta,p11,mean[i],mean[j]);
                    /********************************************************/
                           *res+= weights*log(bl);
                }}}}
    // Checks the return values
    if(!R_FINITE(*res)) *res = LOW;
    return;
}


/*********************************************************/
void Comp_Pair_TWOPIECEGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,
                         double *nuis,int *ns,int *NS, int *GPU,int *local)
{

    int i,j;double bl,corr,zi,zj,lags,weights=1.0,p11,eta,qq,sill,nugget;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       qq=qnorm((1-eta)/2,0,1,1,0);
         if( fabs(eta)>1|| sill<0||nugget>=1||nugget<0) {*res=LOW;  return;} 
    for(i=0;i<(ncoord[0]-1);i++){
        for(j=(i+1); j<ncoord[0];j++){
            lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
            if(lags<=maxdist[0]){
                zi=data[i];zj=data[j];
                if(!ISNAN(zi)&&!ISNAN(zj) ){
                    corr=CorFct(cormod,lags,0,par,0,0);
                      p11=pbnorm(cormod,lags,0,qq,qq,nugget,sill,par,0);
                    if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
                    bl=biv_two_pieceGaussian((1-nugget)*corr,zi,zj,sill,eta,p11,mean[i],mean[j]);
                    *res+= weights*log(bl);
                }}}}
    // Checks the return values
    if(!R_FINITE(*res)) *res = LOW;
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
void Comp_Pair_Gauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
    int i=0, j=0,  t=0, v=0;
    double det=0.0,u=0.0, w=0.0, dens=0.0, rhott=0.0,rhovv=0.0,rhotv=0.0,lags=0.0,weights=1.0;
    //if(CheckCor(cormod,par)==-2){*res=LOW; return;}
    if(  par[0]<0|| par[1]<0|| par[2]<0|| par[3]<0) {*res=LOW;  return;} 

      //// Rprintf(" %f %f %f %f  \n",dista[1][1],dista[0][1],dista[1][0],dista[0][0]);

    // Computes the log-likelihood:
      weights=1;
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
                    
      if(t==v){

         for(j=i+1;j<ns[v];j++){

          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);   
          //lags=dist(type[0],coordx[(i)],coordx[(j)],coordy[(i)],coordy[(j)],*REARTH); 
     
                        if(lags<=dista[t][v]){
                                  
                            rhott=CorFct(cormod,0,0,par,t,t);
                            rhovv=CorFct(cormod,0,0,par,v,v);
                            rhotv=CorFct(cormod,lags,0,par,t,v);
                            det=rhott*rhovv-R_pow(rhotv,2);
                            u=data[(i+NS[t])]-mean[(i+NS[t])];
                            w=data[(j+NS[v])]-mean[(j+NS[v])];
     
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    if(*weigthed)   weights=CorFunBohman(lags,dista[t][v]);
                                    dens=-0.5*(2*log(2*M_PI)+log(det)+(rhovv*R_pow(u,2)+rhott*R_pow(w,2)-2*(u*w)*rhotv)/det);
                            
                                    *res+= dens*weights;

                                }}}
        }
            else {  
         for(j=0;j<ns[v];j++){
         lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=dista[t][v]){
                            rhott=CorFct(cormod,0,0,par,t,t);
                            rhovv=CorFct(cormod,0,0,par,v,v);
                            rhotv=CorFct(cormod,lags,0,par,t,v);
                               det=rhott*rhovv-R_pow(rhotv,2);
                            u=data[(i+NS[t])]-mean[(i+NS[t])];
                            w=data[(j+NS[v])]-mean[(j+NS[v])];
                         
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    if(*weigthed)   weights=CorFunBohman(lags,dista[t][v]);
                                     dens=-0.5*(2*log(2*M_PI)+log(det)+(rhovv*R_pow(u,2)+rhott*R_pow(w,2)-2*(u*w)*rhotv)/det);
                                    *res+= dens*weights;
                                }}}}}}}
    if(!R_FINITE(*res))*res = LOW;
    return;
}


/* pairwise  skew for bivariate GRF*/
void Comp_Pair_SkewGauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)

{
    int i=0, j=0,  t=0, v=0;
    double u=0.0, w=0.0, rhotv=0.0,lags=0.0,weights=1.0;
    int N=2;
      double *vari;vari=(double *) Calloc(N,double);vari[0]=par[0];vari[1]=par[1];  /// variances of the skew gaussian
    par[0]=1;par[1]=1;
    //if(CheckCor(cormod,par)==-2){*res=LOW; return;}
    if(vari[0]<0||vari[1]<0){*res=LOW; return;}
    // Computes the log-likelihood:
      weights=1;
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
             lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=dista[t][v]){
                            rhotv=CorFct(cormod,lags,0,par,t,v);
                             u=data[(i+NS[t])]-mean[(i+NS[t])];
                             w=data[(j+NS[v])]-mean[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                         if(*weigthed)   weights=CorFunBohman(lags,dista[t][v]);
                     *res+= log(biv_skew2(rhotv,u,w,vari[t],vari[v],1,nuis[v],nuis[v]))*weights;
                                }}}}
            else {  
           for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=dista[t][v]){
                            rhotv=CorFct(cormod,lags,0,par,t,v);
                                u=data[(i+NS[t])]-mean[(i+NS[t])];
                            w=data[(j+NS[v])]-mean[(j+NS[v])];
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    if(*weigthed)   weights=CorFunBohman(lags,dista[t][v]);
                *res+= log(biv_skew2(rhotv,u,w,vari[t],vari[v],1,nuis[t],nuis[v]))*weights;      
                                }}}}}}}
        Free(vari);
    if(!R_FINITE(*res))*res = LOW;
    return;
}





/* pairwise for bivariate wrapped GRF*/
void Comp_Pair_WrapGauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, 
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)

{
    int i=0, j=0,  t=0, v=0;
    double det=0.0,u=0.0, w=0.0, rhott=0.0,rhovv=0.0,rhotv=0.0,lags=0.0,weights=1.0;
    double quadr=0.0,wrap_gauss;
    double alfa=2.0;
    // Checks the validity of the nuisance and correlation parameters (nuggets, sills and corr):
    //if(CheckCor(cormod,par)==-2){*res=LOW; return;}
    // Computes the log-likelihood:
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
        lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=dista[t][v]){
                            rhott=CorFct(cormod,0,0,par,t,t);
                            rhovv=CorFct(cormod,0,0,par,v,v);
                            rhotv=CorFct(cormod,lags,0,par,t,v);
                            det=rhott*rhovv-R_pow(rhotv,2);
                            u=data[(i+NS[t])]-2*atan(mean[(i+NS[t])])+M_PI;
                            w=data[(j+NS[v])]-2*atan(mean[(j+NS[v])])+M_PI;
                                if(!ISNAN(u)&&!ISNAN(w) ){
                                    if(*weigthed)   weights=CorFunBohman(lags,dista[t][v]);
                                    double k1=-alfa,k2=-alfa;wrap_gauss = 0;
                                    while(k1<=alfa){
                                        while(k2<=alfa){
        quadr = -0.5*(1.0/det)*(rhovv*R_pow((w+2*k1*M_PI),2.0)+rhott*R_pow((u+2*k2*M_PI),2.0)-2.0*rhotv*(u+2*k2*M_PI)*(w+2*k1*M_PI));
        wrap_gauss = wrap_gauss +  (1/(2.0*M_PI))*(1/sqrt(det)*exp(quadr) ) ;
                                            k2 = k2+1;}
                                        k1 = k1+1;k2 = -alfa;}
                                    *res+=log(wrap_gauss)*weights ;
                                }}}}
           else {  
         for(j=0;j<ns[v];j++){
            lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=dista[t][v]){
                            rhott=CorFct(cormod,0,0,par,t,t);
                            rhovv=CorFct(cormod,0,0,par,v,v);
                            rhotv=CorFct(cormod,lags,0,par,t,v);
                            det=rhott*rhovv-R_pow(rhotv,2);
                            u=data[(i+NS[t])]-2*atan(mean[(i+NS[t])])+M_PI;
                            w=data[(j+NS[v])]-2*atan(mean[(j+NS[v])])+M_PI;
                            if(!ISNAN(u)&&!ISNAN(w) ){
                                    if(*weigthed)   weights=CorFunBohman(lags,dista[t][v]);
                                    double k1=-alfa,k2=-alfa;wrap_gauss = 0;
                                    while(k1<=alfa){
                                        while(k2<=alfa){
                                            quadr = -0.5*(1.0/det)*(rhovv*R_pow((w+2*k1*M_PI),2.0)+rhott*R_pow((u+2*k2*M_PI),2.0)
                                                                    -2.0*rhotv*(u+2*k2*M_PI)*(w+2*k1*M_PI));
                                            wrap_gauss = wrap_gauss +  (1/(2.0*M_PI))*(1/sqrt(det)*exp(quadr) ) ;
                                            k2 = k2+1;
                                        }
                                        k1 = k1+1;k2 = -alfa;}
                                    *res+=log(wrap_gauss)*weights ;
                                    }}}}}}}
    if(!R_FINITE(*res)|| !*res) *res = LOW;
    return;  
}

/*  Pair of pairwise for bivariate GRF */
void Comp_Cond_Gauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  
    double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
       
    int N=4,i=0, j=0, t=0, v=0;
    double weights=1.0,lags=0.0;
    double **M; 
    M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}
    
    double *dat;
    dat=(double *) Calloc(N,double);
    // Checks the validity of the nuisance and correlation parameters (nuggets, sills and corr):
 //   if(CheckCor(cormod,par)==-2){*res=LOW; return;}
    // Computes the log-likelihood:
    for(t=0;t<*ntime-1;t++){
        for(v=t+1;v<*ntime;v++){
            for(i=0;i<ncoord[0]-1;i++){
                for(j=i+1;j<ncoord[0];j++){
                    lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                    if(lags<=dista[t][v]){
                        M[0][0]=CorFct(cormod,0,0,par,t,t);M[0][1]=CorFct(cormod,lags,0,par,t,t);  M[0][2]=CorFct(cormod,0,0,par,t,v);     M[0][3]=CorFct(cormod,lags,0,par,t,v);
                        M[1][0]=M[0][1];                   M[1][1]=CorFct(cormod,0,0,par,t,t);                        M[1][2]=CorFct(cormod,lags,0,par,v,t);  M[1][3]=CorFct(cormod,0,0,par,v,t);
                        M[2][0]=M[0][2];                   M[2][1]= M[1][2];                       M[2][2]=CorFct(cormod,0,0,par,v,v);     M[2][3]=CorFct(cormod,lags,0,par,v,v);
                        M[3][0]=M[0][3];                   M[3][1]= M[1][3];                       M[3][2]=M[2][3];                        M[3][3]= CorFct(cormod,0,0,par,v,v);
                      
                            dat[0]=data[(i+NS[t])]-mean[(i+NS[t])];//z_1(s_i)
                            dat[1]=data[(j+NS[t])]-mean[(j+NS[t])];//z_1(s_j)
                            dat[2]=data[(i+NS[v])]-mean[(i+NS[v])];//z_2(s_i)
                            dat[3]=data[(j+NS[v])]-mean[(j+NS[v])];//z_2(s_j)
                            
                            if(!ISNAN(dat[0])&&!ISNAN(dat[1])&&!ISNAN(dat[2])&&!ISNAN(dat[3]) ){

                                if(*weigthed) weights=CorFunBohman(lags,dista[t][v]);
                                *res+=log(dNnorm(N,M,dat))*weights;  // pair of pairwise likelihood
                            }}
                }}
        }}
    if(!R_FINITE(*res))*res = LOW;
    Free(dat);
    for(i=0;i<N;i++)  {Free(M[i]);}
    Free(M);
    return;
}
