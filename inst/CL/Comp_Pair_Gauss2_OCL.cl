#include "header36.h"

__kernel void Comp_Pair_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    //entretanto(x);
    int j, gid = get_global_id(0);
    
    double  corr=0.0, lags=0.0,weights=1.0, sum=0.0;
    //double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];//nugget
    double nuis1 = dou_par[5]; //sill
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod  = int_par[0];
    int ncoord  = int_par[1];
    int weigthed    = int_par[2];
    int type    = int_par[3];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[gid+j],coordx[j],coordy[gid+j],coordy[j],REARTH);
            if(lags<=maxdist){
                
                if(!isnan(data[gid+j])&&!isnan(data[j]) )
                {
                    corr=CorFct(cormod, lags, 0, par0,par1,par2,par3,0,0);
                    if(weigthed) weights=CorFunBohman(lags,maxdist);
                    sum+=log_biv_Norm((1-nuis0)*corr,data[gid+j],data[j],mean[gid+j],mean[j],nuis1,0)*weights;
                }}}
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}
