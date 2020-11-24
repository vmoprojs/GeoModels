#include "header36.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_Gauss_misp_T2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    
    double lags,weights=1.0, sum=0.0;
    double zi, zj, bl,corr;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    
    double df=1/nuis0;
    double sill=nuis2;
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                if(!isnan(zi)&&!isnan(zj) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0)*(1-nuis1);
                    
                    if(df<=170) corr=0.5*(df-2)*pow(tgamma((df-1)/2),2)/(pow(tgamma(df/2),2))* corr *hypergeo(0.5,0.5,df/2,pow(corr,2));
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    bl=log_biv_Norm(corr,data[gid+j],data[j],mean[gid+j],mean[j],sill,0);
                    sum+= weights*bl;
                }
            }
        }
        else
            continue;
    }
    res[gid] = sum;
}
