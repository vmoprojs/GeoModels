#include "header36.h"


/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_Gauss_misp_T_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
{
    double maxdist = dou_par[6];
    double maxtime    =    dou_par[11];
    double nuis0 = dou_par[4];//nugget
    double nuis1 = dou_par[5];//sill
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
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    int m=0,v =0;
    double  lags=0.0, lagt=0.0,weights=1.0,sum=0.0, corr=0.0, bl;
    double u=0.0,  w=0.0;
    
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
    
    double df=1/nuis0;
    // Set nuisance parameters:
    double sill=nuis2;
    double nugget=nuis1;
   
    
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
                         corr=CorFct_st(cormod,lags,0,par0,par1,par2,par3,par4,par5,par6,t,v)*(1-nugget);
                        
                        if(df<=170) corr=0.5*(df-2)*pow(tgamma((df-1)/2),2)/(pow(tgamma(df/2),2))* corr *hypergeo(0.5,0.5,df/2,pow(corr,2));
                        
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            bl= log_biv_Norm(corr,u,w,mean[(l+NS[t])],mean[(m+NS[v])],sill,0);
                            sum+= bl*weights;
                            
                        }
                        
                    }
                }
            }
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[(l+NS[t])],coordx[(m+NS[v])],coordy[(l+NS[t])],coordy[(m+NS[v])],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        
                        corr=CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,t,v)*(1-nugget);
                        
                        if(df<=170) corr=0.5*(df-2)*pow(tgamma((df-1)/2),2)/(pow(tgamma(df/2),2))* corr *hypergeo(0.5,0.5,df/2,pow(corr,2));
                        
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            bl= log_biv_Norm(corr,u,w,mean[(l+NS[t])],mean[(m+NS[v])],sill,0);
                            sum+= bl*weights;
                            
                        }
                    }
                }
            }
        }
        res[i] = sum;
        //sum = 0.0;
    }
}
