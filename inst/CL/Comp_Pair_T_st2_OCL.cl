#include "header.h"

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Pair_T_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0,bl=0.0;
    
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
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[(l+NS[t])],coordx[(m+NS[v])],coordy[(l+NS[t])],coordy[(m+NS[v])],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            //if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                                if(bl<0||bl>9999999999999999||!isfinite(bl)||isnan(bl)||isinf(bl)) bl=1;
                    sum+= weights*log(bl);
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[(l+NS[t])],coordx[(m+NS[v])],coordy[(l+NS[t])],coordy[(m+NS[v])],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            //if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                        bl=biv_T(corr*(1-nuis1),(zi-mean[(l+NS[t])])/sqrt(nuis2), (zj-mean[(m+NS[v])])/sqrt(nuis2),nuis0)/nuis2;
                                        if(bl<0||bl>9999999999999999||!isfinite(bl)||isnan(bl)||isinf(bl)) bl=1;
                    sum+= weights*log(bl);
                        }
                    }}}
        }
        res[i] = sum;
    }
}