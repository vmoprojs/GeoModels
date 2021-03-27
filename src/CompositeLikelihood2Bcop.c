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
void Comp_Pair_Gauss2memBcop(int *cormod, double *data1,double *data2,int *NN, 
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
