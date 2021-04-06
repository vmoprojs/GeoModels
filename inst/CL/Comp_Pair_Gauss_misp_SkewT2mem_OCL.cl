#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_Gauss_misp_SkewT2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0),uu,ww;
    double corr,corr2, bl, l1,l2,weights=1.0, sum=0.0,zi,zj,mui,muj;
    double df = 1/dou_par[4];//nuis0
    double nugget = dou_par[5];//nuis1
    double sill = dou_par[9];//nuis2
    double skew = dou_par[10];//nuis3
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    int cormod  = int_par[0];
    
    double D1=(df-1)/2;
    double D2=df/2;
    //double delta=skew/sqrt(1-skew*skew);
    double MM=sqrt(df)*tgamma(D1)*skew/(sqrt(M_PI)*tgamma(D2));
    double FF=(df/(df-2)-MM*MM);


    zi = data1[gid];
    zj = data2[gid];
   if(!isnan(zi)&&!isnan(zj) )
    {

        corr = CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
        
        corr2=corr_skewt(corr,df,skew);
        
        bl=log_biv_Norm(corr2,zi,zj,mean1[gid]+sqrt(sill)*MM,mean2[gid]+sqrt(sill)*MM,sill*FF,0);
        
        //l1=dnorm(zi,mean1[gid]+sqrt(sill)*MM,sqrt(sill*FF),1);
        //l2=dnorm(zj,mean2[gid]+sqrt(sill)*MM,sqrt(sill*FF),1);
        sum+= (bl)*weights;
    
        
        res[gid] = sum;
    }
}
