#include "header36.h"


/******************************************************************************************/
/********************* SPACE TIME CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_BinomnegGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
{
    
    double maxdist = dou_par[6];
    double maxtime	=	dou_par[11];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[6];
    
    
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    int uu=0,ww=0;
    double dens=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0, a=0.0,b=0.0,sum=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    lags=dist(type,coordx[(l+NS[t])],coordx[(m+NS[v])],coordy[(l+NS[t])],coordy[(m+NS[v])],REARTH);
                    if(lags<=maxdist)
                    {
                        a = mean[(l+NS[t])];
                        b = mean[(m+NS[v])];
                        
                        psj=pbnorm_st(cormod,lags,0,a,b,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(a,0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(b,0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            
                            uu=(int) u;
                            ww=(int) w;
                            //if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            dens=biv_binomneg(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                        //printf("GPU: %d\t%d\t%d\t%d\t%f\n",l,t,v,m,sum);
                    }}}
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[(l+NS[t])],coordx[(m+NS[v])],coordy[(l+NS[t])],coordy[(m+NS[v])],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        a = mean[(l+NS[t])];
                        b = mean[(m+NS[v])];
                        
                        psj=pbnorm_st(cormod,lags,lagt,a,b,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(mean[(l+NS[t])],0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(mean[(m+NS[v])],0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        
                        if(!isnan(u)&&!isnan(w) ){
                            uu=(int) u;
                            ww=(int) w;
                           // if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            dens=biv_binomneg(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                    }}}}
        res[i] = sum;
    }
}
