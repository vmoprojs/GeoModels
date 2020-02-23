#include "header32.h"


/******************************************************************************************/
/********************* SPACE TIME CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Cond_Gauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    
    //int ls = get_global_size(0);
    //int ms = get_global_size(1);
    
    int m=0,v =0;
    double s1=0.0, s12=0.0, lags=0.0, lagt=0.0,weights=1.0,sum=0.0;
    double det=0.0, u=0.0, u2=0.0, w=0.0, w2=0.0;
    
    //int gid = (npts*t+l);
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    //int j = (ntime*gidx+gidy);
    
    //for (j = 0; j < ncoord; j++) {
    //    if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    //if(l >= ncoord) isValid = false;
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        s1=nuis0+nuis1;
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                
                for(m=l+1;m<ns[t];m++)
                {
                    
                    lags=dist(type,coordx[(l+NS[t])],coordx[(m+NS[v])],coordy[(l+NS[t])],coordy[(m+NS[v])],REARTH);
                    if(lags<=maxdist)
                    {
                        s12=nuis1*CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,t,v);
                        det=pow(s1,2)-pow(s12,2);
                        u=data[(l+NS[t])]-mean[(l+NS[t])];
                        w=data[(m+NS[v])]-mean[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            u2=pow(u,2);
                            w2=pow(w,2);
                           // if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= (-log(2*M_PI)-log(det)+log(s1)+
                                   (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det)*weights;}
                    }}
            }
            else
            {
                
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    
                    lags=dist(type,coordx[(l+NS[t])],coordx[(m+NS[v])],coordy[(l+NS[t])],coordy[(m+NS[v])],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        
                        s12=nuis1*CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,t,v);
                        det=pow(s1,2)-pow(s12,2);
                        
                        u=data[(l+NS[t])]-mean[(l+NS[t])];
                        w=data[(m+NS[v])]-mean[(m+NS[v])];
    
                        if(!isnan(u)&&!isnan(w) ){
                            u2=pow(u,2);
                            w2=pow(w,2);
                            
                            //if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= (-log(2*M_PI)-log(det)+log(s1)+
                                   (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det)*weights;}
                    }
                }
            }
        }
        res[i] = sum;
    }
}
