#include "header45.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Cond_Gauss_misp_T2mem_OCL(__global const double *data1,__global const double *data2,__global const double *mean1,__global const double *mean2, __global const double *lags, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    int j, gid = get_global_id(0);
    double corr,corr1, bl, l1,weights=1.0, sum=0.0,zi,zj,df,var;
    double l2=0.0;
    //double maxdist = dou_par[6];
    double nuis0 = dou_par[4]; //nuis0
    double nugget = dou_par[5];//nuis1
    double sill = dou_par[9];//nuis2
    double nuis3 = dou_par[10];//nuis3
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    //double REARTH = dou_par[8];
    
    
    //int ncoord  = int_par[1];
    int cormod  = int_par[0];
    //int type    = int_par[3];
    
    df = 1/nuis0;
    var = sill*df/(df-2);
    double mu,vv,u,eta2;
    
    
    zi = data1[gid];
    zj = data2[gid];
   if(!isnan(zi)&&!isnan(zj) )
    {
        
        corr = (1-nugget)*CorFct(cormod, lags[gid], 0, par0,par1,par2,par3,0,0);
        corr1=exp(log(df-2)+2*lgamma(0.5*(df-1))-(log(double (2))+2*lgamma(df/2))+log(hypergeo(0.5,0.5, df/2,corr*corr))+log(corr*(1-nugget)));
        
        l1=dnorm(zi,mean1[gid],sqrt(sill*df/(df-2)),1);
        l2=dnorm(zj,mean2[gid],sqrt(sill*df/(df-2)),1);
        bl=2*log_biv_Norm(corr1,data1[gid],data2[gid],mean1[gid],mean2[gid],sill*df/(df-2),0)-(l1+l2);
        
        sum+= bl*weights;
        
        res[gid] = sum;
        
    }
    
    
}
