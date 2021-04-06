#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_TWOPIECETukeyh2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    int j, gid = get_global_id(0);
    double corr, bl, l1,weights=1.0, sum=0.0,zi=0.0,zj=0.0;
    double l2=0.0, p11,qq;
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
    zi=data1[gid];zj=data2[gid];
    qq=qnorm55((1-eta)/2,0,1,1,0);
   if(!isnan(zi)&&!isnan(zj) )
    {

        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
    
       //l1=one_log_two_pieceTukey(zi,mean1[gid],sill,tail,eta);
       //l2=one_log_two_pieceTukey(zj,mean2[gid],sill,tail,eta);


        p11=pbnorm22(qq,qq,corr);
        bl=(biv_two_pieceTukeyh((1-nugget)*corr,zi,zj,sill,eta,tail,p11,mean1[gid],mean2[gid]));
        
        sum+= weights*log(bl);
        
        res[gid] = sum;
    }
}
