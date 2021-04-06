#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Cond_Pois2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0),uu,ww;
    double corr, bl, l1,l2,weights=1.0, sum=0.0,zi,zj,mui,muj;
    double nugget = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    int cormod  = int_par[0];


    zi = data1[gid];
    zj = data2[gid];
   if(!isnan(zi)&&!isnan(zj) )
    {
    
        mui=exp(mean1[gid]);muj=exp(mean2[gid]);

        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
        
        uu=(int) data1[gid];  ww=(int) data2[gid];
        l1=dpois(uu,mui,1);
        l2=dpois(ww,muj,1);
        //printf("l1: %f, l2: %f\n",l1,l2);
        
        bl=2*log(biv_Poisson((1-nugget)*corr,uu,ww,mui, muj))
          - (l1+l2);
        sum+= (bl)*weights;
    
        
        res[gid] = sum;
    }
}
