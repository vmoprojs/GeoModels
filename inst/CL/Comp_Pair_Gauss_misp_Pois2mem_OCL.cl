#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_Gauss_misp_Pois2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    int j, gid = get_global_id(0);
    double corr,corr2, bl, l1,weights=1.0, sum=0.0,zi,zj,mui,muj,corr1;
    double l2=0.0;
    //double maxdist = dou_par[6];
    double nugget = dou_par[4];
    double sill = dou_par[5];
    double eta = dou_par[9];
    double tail = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    //double REARTH = dou_par[8];
    
    //int ncoord  = int_par[1];
    int cormod  = int_par[0];
    //int type    = int_par[3];

    zi = data1[gid];
    zj = data2[gid];
   if(!isnan(zi)&&!isnan(zj) )
    {
        mui=exp(mean1[gid]);muj=exp(mean2[gid]);
        
        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0)*(1-nugget);
         corr1=corr_pois(corr,mui, muj);
        

        //l1=dnorm(zi ,mui,sqrt(mui),1);
        //l2=dnorm(zj ,muj,sqrt(muj),1);;
        bl=dNnorm(zi-mui,zj-muj,mui,muj,sqrt(mui*muj)*corr1);
        sum+= log(bl)*weights;
        
        res[gid] = sum;
        
    }
    
    
}
