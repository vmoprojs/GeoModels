#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Cond_Gamma2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    int j, gid = get_global_id(0);
    double corr, bl, l1,weights=1.0, sum=0.0,zi,zj;
    double l2=0.0;
    //double maxdist = dou_par[6];
    double nugget = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    //double REARTH = dou_par[8];
    
    //int ncoord  = int_par[1];
    int cormod  = int_par[0];
    //int type    = int_par[3];
    zi = data1[gid];zj =data2[gid];
   if(!isnan(zi)&&!isnan(zj) )
    {
        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
        
        l1=one_log_gamma(zi,mean1[gid],nuis2);
        l2=one_log_gamma(zj,mean2[gid],nuis2);
        bl=2*log(biv_gamma((1-nugget)*corr,zi,zj,mean1[gid],mean2[gid],nuis2))-(l1+l2);
        sum+= (bl)*weights;
        
        res[gid] = sum;
        
    }
    
    
}
