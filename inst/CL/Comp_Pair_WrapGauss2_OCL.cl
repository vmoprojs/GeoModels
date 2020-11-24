#include "header36.h"

__kernel void Comp_Pair_WrapGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double u=0.0,v=0.0,weights=1.0;
    double lags=0.0, corr=0.0,sum=0.0,wrap_gauss,alfa=2.0;
    
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                u=data[gid+j];//-2*atan(mean[gid+j])-M_PI;
                v=data[j];//-2*atan(mean[j])-M_PI;;
                //printf("u,v: %f\t%f\n",u,v);
                if(!isnan(u)&&!isnan(v) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    wrap_gauss=biv_wrapped(alfa,u,v,mean[gid+j],mean[j],nuis0,nuis1,corr);
                    
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    sum+=  log(wrap_gauss)*weights ;
                }
                
            }
        }
        
        else
            continue;
    }
    //barrier(CLK_GLOBAL_MEM_FENCE);
    res[gid] = sum;
    
}
