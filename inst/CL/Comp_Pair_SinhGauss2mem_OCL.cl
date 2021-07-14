#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_SinhGauss2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    int j, gid = get_global_id(0);
    double corr, bb, l1,weights=1.0, sum=0.0,ai=0,aj=0,p11=0.0;
    double l2=0.0,p1=0.0,p2=0.0,zi,zj;
    //double maxdist = dou_par[6];
    double nuis0 = dou_par[4];// nuis0
    double nuis1 = dou_par[5];// nuis1
    double nuis2 = dou_par[9];// nuis2
    double nuis3 = dou_par[10];// nuis3
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    int NN = int_par[6],uu=0,vv=0;
    int cormod  = int_par[0];

    zi=data1[gid];zj=data2[gid];
   if(!isnan(zi)&&!isnan(zj) )
    {
    
         //ai=mean1[gid];aj=mean2[gid];
        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
        
        //l1=one_log_sas(zi,mean1[gid],nuis2,nuis3,nuis1);
        //l2=one_log_sas(zj,mean2[gid],nuis2,nuis3,nuis1);
        
        bb=log(biv_sinh((1-nuis0)*corr,zi,zj,mean1[gid],mean2[gid],nuis2,nuis3,nuis1));
        
        sum+= (bb)*weights;
        res[gid] = sum;
        
    }
    
    
}
