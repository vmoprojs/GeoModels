#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_Kumaraswamy2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    int j, gid = get_global_id(0);
    double corr, bl, l1,weights=1.0, sum=0.0,zi,zj;
    double l2=0.0;
    //double maxdist = dou_par[6];
    double nugget = dou_par[4];//nuis0
    double nuis1 = dou_par[5];//nuis1
    double nuis2 = dou_par[9];//nuis2
    double nuis3 = dou_par[10];//nuis3
    
    double min=dou_par[11]; //nuis4
    double max=dou_par[12]; //nuis5

    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    //double REARTH = dou_par[8];

    int cormod  = int_par[0];
    
    zi = data1[gid];zj =data2[gid];
   if(!isnan(zi)&&!isnan(zj) )
    {
        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
        
        bl=biv_Kumara((1-nugget)*corr,zi,zj,mean1[gid],mean2[gid],nuis2,nuis3,min,max);
        
        sum+= log(bl)*weights;
        
        res[gid] = sum;
    }
    
    
}
