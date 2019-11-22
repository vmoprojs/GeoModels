#include "header.h"

__kernel void Comp_Diff_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double lags=0.0,weights=1.0, sum=0.0,vario=0.0;
    double u, v;
    
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
    
    
    for (j = 0; j < (ncoord); j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[gid+j],coordx[j],coordy[gid+j],coordy[j],REARTH);
            if(lags<=maxdist){
                
                vario=Variogram(cormod,lags,0,nuis0,nuis1,par0,par1,par2,par3);
                u=data[gid+j];
                v=data[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    if(weigthed) weights=CorFunBohman(lags,maxdist);
                    sum+=  -0.5*(log(2*M_PI)+log(vario)+
                                 pow(u-v,2)/(2*vario))*weights;
                }
                
            }
        }
        
        else
            continue;
    }
    //barrier(CLK_GLOBAL_MEM_FENCE);
    //printf("Ey!\n");
    
    res[gid] = sum;
    
}
