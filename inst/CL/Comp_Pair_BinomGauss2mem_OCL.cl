#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_BinomGauss2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    int j, gid = get_global_id(0);
    double corr, bl, l1,weights=1.0, sum=0.0,ai=0,aj=0,p11=0.0;
    double l2=0.0,p1=0.0,p2=0.0,u,v;
    //double maxdist = dou_par[6];
    double nugget = dou_par[4];// nuis0
    double nuis1 = dou_par[5];// nuis1
    double nuis2 = dou_par[9];// nuis2
    double nuis3 = dou_par[10];// nuis3
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    int NN = int_par[6],uu=0,vv=0;
    int cormod  = int_par[0];

    
    u=data1[gid];v=data2[gid];
   if(!isnan(u)&&!isnan(v) )
    {
        
         ai=mean1[gid];aj=mean2[gid];
        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
        p11=pbnorm22(ai,aj,(1-nugget)*corr);
        
        p1=pnorm(ai,0,1,1,0);
        
        p2=pnorm(aj,0,1,1,0);
        //u=data1[gid];v=data2[gid];
        
        uu=(int) u; vv=(int) v;
        //l1=dbinom(uu,NN,p1,1);
       
        //printf("l1:%f uu:%d NN:%d p1:%f \n\n",l1,uu,NN,p1);
          //l2=dbinom(vv,NN,p2,1);
        bl=biv_binom (NN,uu,vv,p1,p2,p11);
        //bl=2*log(0.1)-(l1+l2);
        sum+= log(bl)*weights;
        
        res[gid] = sum;
        
    }
    
    
}
