#include "header32.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Cond_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double s1=0.0, s12=0.0, lags=0.0,weights=1.0, sum=0.0;
    double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int ncoord  = int_par[1];
    int cormod  = int_par[0];
    int type    = int_par[3];
    
    
    
    
    s1=nuis0+nuis1;
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            
            if(lags<=maxdist){
            
                s12=nuis1*CorFct(cormod, lags, 0, par0,par1,par2,par3,0,0);
                det=pow(s1,2)-pow(s12,2);
            
                u=data[gid+j]-mean[gid+j];
                v=data[j]-mean[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    u2=pow(u,2);
                    v2=pow(v,2);
                    sum+= (-log(2*M_PI)-log(det)+log(s1)+
                           (u2+v2)*(0.5/s1-s1/det)+2*s12*u*v/det)*weights;
                }
                
            }
            
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}
