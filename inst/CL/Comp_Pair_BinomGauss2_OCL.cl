#include "header36.h"

__kernel void Comp_Pair_BinomGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0),uu=0,vv=0;
    double u,v,dens=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0, sum=0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[5];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                ai=mean[gid+j];
                aj=mean[j];
                
                psj=pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);//pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);
                p1=pnorm_OCL(ai,0,1);//pnorm_OCL(ai,0,1)
                p2=pnorm_OCL(aj,0,1);//pnorm_OCL(aj,0,1)
                
                u=data[gid+j];
                v=data[j];
                if(!isnan(u)&&!isnan(v) )
                {
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    uu=(int) u;
                    vv=(int) v;
                    dens=biv_binom(NN,uu,vv,p1,p2,psj);//biv_binom(NN,uu,vv,p1,p2,psj)
                    sum+=  log(dens)*weights; //
                    
                } } }
        else
            continue;
    }
    res[gid] = sum;
}
