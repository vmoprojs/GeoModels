#include "header.h"

void Comp_Pair_Gauss_misp_SkewT2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_Gauss_misp_SkewT2mem_OCL";
    int *int_par;
    double *dou_par;
    
    double sill,nugget,skew,df;
    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];

     if(df<2||fabs(skew)>1||sill<0||nugget<0||nugget>=1){*res=LOW; return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}
void Comp_Cond_Gauss_misp_SkewT2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_Gauss_misp_SkewT2mem_OCL";
    int *int_par;
    double *dou_par;
    
    double sill,nugget,skew,df;
    df=1/nuis[0];
    nugget=nuis[1];
    sill=nuis[2];
    skew=nuis[3];

     if(df<2||fabs(skew)>1||sill<0||nugget<0||nugget>=1){*res=LOW; return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Pois2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_Pois2mem_OCL";
    int *int_par;
    double *dou_par;
    
   double nugget=nuis[0];

    if(nugget<0||nugget>=1){*res=LOW; return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Cond_Pois2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_Pois2mem_OCL";
    int *int_par;
    double *dou_par;
    
   double nugget=nuis[0];

    if(nugget<0||nugget>=1){*res=LOW; return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Kumaraswamy22mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_Kumaraswamy22mem_OCL";
    int *int_par;
    double *dou_par;
    
    double min=nuis[4];
    double max=nuis[5];
    if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Cond_Kumaraswamy22mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_Kumaraswamy22mem_OCL";
    int *int_par;
    double *dou_par;
    
    double min=nuis[4];
    double max=nuis[5];
    if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Kumaraswamy2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_Kumaraswamy2mem_OCL";
    int *int_par;
    double *dou_par;
    
    double min=nuis[4];
    double max=nuis[5];
    if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_Kumaraswamy2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_Kumaraswamy2mem_OCL";
    int *int_par;
    double *dou_par;
    
    double min=nuis[4];
    double max=nuis[5];
    if(nuis[2]<0||nuis[3]<0||min>max) {*res=LOW;  return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Beta2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_Beta2mem_OCL";
    int *int_par;
    double *dou_par;
    
    //double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Cond_Beta2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_Beta2mem_OCL";
    int *int_par;
    double *dou_par;
    
    //double nugget=nuis[0];
    double min=nuis[4];
     double max=nuis[5];
     if(nuis[2]<0||nuis[3]<0||min>max)  {*res=LOW;  return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Gauss_misp_T2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_Gauss_misp_T2mem_OCL";
    int *int_par;
    double *dou_par;
    double sill=nuis[2];
      double nugget=nuis[1];
      if( sill<0||nugget<0||nugget>=1||nuis[0]<0||nuis[0]>0.5){*res=LOW; return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_Gauss_misp_T2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_Gauss_misp_T2mem_OCL";
    int *int_par;
    double *dou_par;
    double sill=nuis[2];
      double nugget=nuis[1];
      if( sill<0||nugget<0||nugget>=1||nuis[0]<0||nuis[0]>0.5){*res=LOW; return;}
    
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Gauss_misp_Pois2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_Gauss_misp_Pois2mem_OCL";
    int *int_par;
    double *dou_par;
   double nugget=nuis[0];
        
    if(nugget<0||nugget>=1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Cond_Gauss_misp_Pois2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_Gauss_misp_Pois2mem_OCL";
    int *int_par;
    double *dou_par;
   double nugget=nuis[0];
        
    if(nugget<0||nugget>=1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_SinhGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_SinhGauss2mem_OCL";
    int *int_par;
    double *dou_par;
   //double nugget=nuis[0];
        
    if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Cond_SinhGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_SinhGauss2mem_OCL";
    int *int_par;
    double *dou_par;
   //double nugget=nuis[0];
        
    if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_BinomnegGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_BinomnegGauss2mem_OCL";
    int *int_par;
    double *dou_par;
   double nugget=nuis[0];
        
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_BinomnegGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_BinomnegGauss2mem_OCL";
    int *int_par;
    double *dou_par;
   double nugget=nuis[0];
        
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_BinomGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_BinomGauss2mem_OCL";
    int *int_par;
    double *dou_par;
   double nugget=nuis[0];
        
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_BinomGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{
    char *f_name = "Comp_Cond_BinomGauss2mem_OCL";
    int *int_par;
    double *dou_par;
   double nugget=nuis[0];
        
    if( nugget>=1 || nugget<0){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_T2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_T2mem_OCL";
    int *int_par;
    double *dou_par;
       double sill=nuis[2];
        double nugget=nuis[1];
        double df=nuis[0];
        
          if( sill<0||nugget<0||nugget>=1||df<0||df>0.5){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_T2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_T2mem_OCL";
    int *int_par;
    double *dou_par;
       double sill=nuis[2];
        double nugget=nuis[1];
        double df=nuis[0];
        
          if( sill<0||nugget<0||nugget>=1||df<0||df>0.5){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_LogGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_LogGauss2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1];double nugget=nuis[0];
        if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_LogGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_LogGauss2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1];double nugget=nuis[0];
        if(sill<0 || nugget<0||nugget>1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_Weibull2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_Weibull2mem_OCL";
    int *int_par;
    double *dou_par;
        double nugget=nuis[0];
        if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Weibull2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_Weibull2mem_OCL";
    int *int_par;
    double *dou_par;
        double nugget=nuis[0];
        if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_Gamma2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_Gamma2mem_OCL";
    int *int_par;
    double *dou_par;
        double nugget=nuis[0];
        if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Gamma2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_Gamma2mem_OCL";
    int *int_par;
    double *dou_par;
        double nugget=nuis[0];
        if(nugget<0||nugget>=1||nuis[2]<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}
void Comp_Cond_Gauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_Gauss2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1]; double nugget=nuis[0];
    if(nugget>=1||nugget<0||sill<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Cond_TWOPIECEGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                                     int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_TWOPIECEGauss2mem_OCL";
    int *int_par;
    double *dou_par, eta, sill, nugget;
        eta=nuis[2];  //skewness parameter
        sill=nuis[1];
        nugget=nuis[0];
             if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_TWOPIECEGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                                     int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_TWOPIECEGauss2mem_OCL";
    int *int_par;
    double *dou_par, eta, sill, nugget;
        eta=nuis[2];  //skewness parameter
        sill=nuis[1];
        nugget=nuis[0];
             if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_TWOPIECETukeyh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                                     int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_TWOPIECETukeyh2mem_OCL";
    int *int_par;
    double *dou_par, eta, sill, nugget,tail;
        eta=nuis[2];  //skewness parameter
        tail = nuis[3];  //tail parameter
        sill=nuis[1];
        nugget=nuis[0];
             if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_TWOPIECETukeyh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                                     int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_TWOPIECETukeyh2mem_OCL";
    int *int_par;
    double *dou_par, eta, sill, nugget,tail;
        eta=nuis[2];  //skewness parameter
        tail = nuis[3];  //tail parameter
        sill=nuis[1];
        nugget=nuis[0];
             if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}

    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_TWOPIECET2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                                     int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_TWOPIECET2mem_OCL";
    int *int_par;
    double *dou_par, eta, sill, nugget,df;
        sill=nuis[2];  //skewness parameter
        eta = nuis[3];  //tail parameter
        nugget=nuis[1];
        df=nuis[0];
             if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}

    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}
void Comp_Pair_TWOPIECET2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                                     int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_TWOPIECET2mem_OCL";
    int *int_par;
    double *dou_par, eta, sill, nugget,df;
        sill=nuis[2];  //skewness parameter
        eta = nuis[3];  //tail parameter
        nugget=nuis[1];
        df=nuis[0];
             if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}

    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_Tukeyh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_Tukeyh2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1];
        double nugget=nuis[0];
        double tail=nuis[2];
          if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_Tukeyhh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_Tukeyhh2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1];
        double nugget=nuis[0];
        double h1=nuis[3];
        double h2=nuis[2];
          if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Tukeyhh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_Tukeyhh2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1];
        double nugget=nuis[0];
        double h1=nuis[3];
        double h2=nuis[2];
          if( sill<0||h1<0||h1>0.5||h2<0||h2>0.5||nugget<0||nugget>=1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_Gauss_misp_Tukeygh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_Gauss_misp_Tukeygh2mem_OCL";
    int *int_par;
    double *dou_par;
        double eta  = nuis[2];  //skewness parameter
        double tail = nuis[3];  //tail parameter
        double sill =nuis[1];
        double nugget=nuis[0];
    double mu,vv,u,eta2;

        eta2=eta*eta;
        u=1-tail;
        mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
        vv=((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*
                               sqrt(1-2*tail))-mu*mu);
        if(fabs(eta)<1e-5)
               {
               mu=0.0;
               vv=R_pow(1-2*tail,-3/2);
               }
             if(sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}
void Comp_Pair_Gauss_misp_Tukeygh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    
    char *f_name = "Comp_Pair_Gauss_misp_Tukeygh2mem_OCL";
    int *int_par;
    double *dou_par;
        double eta  = nuis[2];  //skewness parameter
        double tail = nuis[3];  //tail parameter
        double sill =nuis[1];
        double nugget=nuis[0];
    double mu,vv,u,eta2;

         eta2=eta*eta;
           u=1-tail;
           mu=(exp(eta2/(2*u))-1)/(eta*sqrt(u));
           vv=((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*
                                  sqrt(1-2*tail))-mu*mu);
               if(fabs(eta)<1e-5)
                  {
                  mu=0.0;
                  vv=R_pow(1-2*tail,-3/2);
                  }
                if(sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
    
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Tukeyh2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_Tukeyh2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1];
        double nugget=nuis[0];
        double tail=nuis[2];
          if( sill<0||tail<0||tail>0.5||nugget<0||nugget>=1){*res=LOW; return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Gauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_Gauss2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1]; double nugget=nuis[0];
   if(nugget>=1||nugget<0||sill<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_SkewGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Pair_SkewGauss2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1]; double nugget=nuis[0];
    if(nugget>=1||nugget<0||sill<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Cond_SkewGauss2mem_OCL(int *cormod, double *data1, double *data2, int *NN,
                    double *par,  int *weigthed,double *res,double *mean1,double *mean2,double *nuis,
                    int *local_wi, int *dev)
{

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_SkewGauss2mem_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1]; double nugget=nuis[0];
    if(nugget>=1||nugget<0||sill<0) {*res=LOW;  return;}
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    //Rprintf("++++++++ npairs: %d \n",npairs[0]);
    param_OCL_mem(cormod,NN,npairs,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_mem(data1,data2,mean1,mean2, lags, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

// Composite conditional log-likelihood for the spatial Gaussian model:
void Comp_Cond_Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                    double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                    int *local_wi, int *dev)
{ 

    //printf("Modelo de correlacion: %d\n",*cormod);
    char *f_name = "Comp_Cond_Gauss2_OCL";
    int *int_par;
    double *dou_par;
        double sill=nuis[1]; double nugget=nuis[0];
    if(nugget<0||nugget>=1||sill<0) {*res=LOW;  return;}
      int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy,mean, data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                          double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                          int *local_wi, int *dev)
{
  
    double sill,nugget;
    sill=nuis[1];nugget=nuis[0];
    if(nugget<0||nugget>=1||sill<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Gauss2_OCL";
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

// Composite marginal (difference) log-likelihood for the spatial Gaussian model:
void Comp_Diff_Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
                      double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Diff_Gauss2_OCL";
    int *int_par;
    double *dou_par;
     double   sill=nuis[1];double nugget=nuis[0];
    if(nugget<0||nugget>=1||sill<0) {*res=LOW;  return;}
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}



void Comp_Pair_Logistic2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{

    double sill=nuis[1]; double nugget=nuis[0];
    if(nugget<0||nugget>=1||sill<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Logistic2_OCL";
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res)) *res = LOW;
}




void Comp_Pair_Gamma2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                          double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
     double nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Gamma2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res)) *res = LOW;    
}



void Comp_Pair_LogGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_LogGauss2_OCL";
        double sill,nugget;
    sill=nuis[1];nugget=nuis[0];
    if(nugget<0||nugget>=1||sill<0) {*res=LOW;  return;}
    
    int *int_par;
    double *dou_par;
         int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Weibull2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                                int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    double  nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Weibull2_OCL";
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res)) *res = LOW;
}

void Comp_Pair_LogLogistic2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    char *f_name = "Comp_Pair_LogLogistic2_OCL";
    int *int_par;
    double *dou_par;
      double nugget=nuis[0];
     if(nugget<0||nugget>=1||nuis[2]<=2) {*res=LOW;  return;}
         int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res)) *res = LOW;
}


// Composite marginal (pariwise) log-likelihood for the spatial Gaussian model:
void Comp_Pair_T2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                          double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                          int *local_wi, int *dev)
{
    double df=nuis[0]; double nugget=nuis[1];
    if( df<0||df>0.5||nugget<0||nugget>=1){*res=LOW; return;}
    char *f_name = "Comp_Pair_T2_OCL";
    int *int_par;
    double *dou_par;
        int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}






void Comp_Pair_WrapGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                          double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_WrapGauss2_OCL";
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_SinhGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
      if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_SinhGauss2_OCL";
    
    int *int_par;
    double *dou_par;
         int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}



void Comp_Pair_SkewGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    double sill=nuis[1];double nugget=nuis[0];
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}
    char *f_name = "Comp_Pair_SkewGauss2_OCL";
    int *int_par;
    double *dou_par;
         int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}




void Comp_Pair_PoisbinnegGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                          double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
   double nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}
    
    char *f_name = "Comp_Pair_PoisbinnegGauss2_OCL";
    int *int_par;
    double *dou_par;
         int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_PoisbinGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                    double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
   double nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}
    
    char *f_name = "Comp_Pair_PoisbinGauss2_OCL";
    
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_BinomGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                    double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
   double nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_BinomGauss2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_BinomnegGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                               double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
  double nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}
    
    char *f_name = "Comp_Pair_BinomnegGauss2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Binom2Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
  double nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Binom2Gauss2_OCL";
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
      Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}











void Comp_Pair_TWOPIECET2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                      double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                      int *local_wi, int *dev)
{
        double eta=nuis[3];  //skewness parameter
    double df=nuis[0];
    double sill=nuis[2];
    double nugget=nuis[1];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_TWOPIECET2_OCL";
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_TWOPIECEGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                              double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                              int *local_wi, int *dev)
{
   double eta,sill,nugget;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
       //qq=qnorm((1-eta)/2,0,1,1,0);*/
      //   if( fabs(eta)>1) {*res=LOW;  return;}
    if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;} 
    char *f_name = "Comp_Pair_TWOPIECEGauss2_OCL";
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Kumaraswamy2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                                  double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                                  int *local_wi, int *dev)
{
    //if(nuis[2]<0||nuis[3]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Kumaraswamy2_OCL";
    int *int_par;
    double *dou_par;
      int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_Tukeyh2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                           double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                           int *local_wi, int *dev)
{

    double tail=nuis[2];double nugget= nuis[0];double sill =nuis[1];
    if( tail>0.5||tail<0||nugget<0||nugget>=1||sill<0 ){*res=LOW; return;}
    char *f_name = "Comp_Pair_Tukeyh2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res)||!*res)*res = LOW;
}


void Comp_Pair_Gauss_misp_T2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                           double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                           int *local_wi, int *dev)
{
    
    if( nuis[0]<0 || nuis[0]>0.5|| nuis[1]<0 || nuis[1]>=1 ){*res=LOW; return;}
    char *f_name = "Comp_Pair_Gauss_misp_T2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res)||!*res)*res = LOW;
}

void Comp_Pair_TWOPIECETukeyh2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                      double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                      int *local_wi, int *dev)
{
     double eta,tail,sill,nugget;
       eta  = nuis[2];  //skewness parameter
       tail = nuis[3];  //tail parameter
       sill =nuis[1];
       nugget=nuis[0];
           // if( fabs(eta)>1 || tail<=0) {*res=LOW;  return;} 

    if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1||tail<0||tail>0.5) {*res=LOW;  return;} 

    char *f_name = "Comp_Pair_TWOPIECETukeyh2_OCL";
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


// Composite marginal (pariwise) log-likelihood for poisson model
void Comp_Pair_Pois2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
int *local_wi, int *dev)
{
    double nugget=nuis[0];
    if(nugget<0||nugget>=1) {*res=LOW;  return;}

    char *f_name = "Comp_Pair_Pois2_OCL";
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
     Free(int_par);
    Free(dou_par);
    
    // Checks the return values
    if(!R_FINITE(*res))  *res = LOW;
    return;
}


void Comp_Pair_Gauss_misp_Pois2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                           double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                           int *local_wi, int *dev)
{
    

    double nugget=nuis[0];
    if(nugget<0||nugget>=1){*res=LOW; return;}
    char *f_name = "Comp_Pair_Gauss_misp_Pois2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel(coordx,coordy, mean,data, int_par, dou_par, local_wi,dev,res,f_name);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res)||!*res)*res = LOW;
}

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/********************* SPACE TIME CASE *****************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

// Composite marginal (pariwise) log-likelihood for the spatial-temporal Gaussian model:

void Comp_Pair_Gauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
   // if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Gauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}
void Comp_Pair_WrapGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_WrapGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
      int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);       
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_PoisbinGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
   // if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    char *f_name = "Comp_Pair_PoisbinGauss_st2_OCL";
    int *int_par;
    double *dou_par;
   int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}

void Comp_Pair_PoisbinnegGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                    double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    char *f_name = "Comp_Pair_PoisbinnegGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;

}


void Comp_Cond_Gauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Cond_Gauss_st2_OCL";
    int *int_par;
    double *dou_par;
      int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Diff_Gauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if(nuis[1]<0 || nuis[0]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Diff_Gauss_st2_OCL";
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_SkewGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
     double sill=nuis[1];double nugget=nuis[0];
     if(nugget<0|| nugget>=1||sill<0){*res=LOW;  return;}
    char *f_name = "Comp_Pair_SkewGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}



void Comp_Pair_SinhGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    if(nuis[3]<0||nuis[1]<0||nuis[0]<0||nuis[0]>=1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_SinhGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Gamma_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if(nuis[2]<1||sill<0||sill>1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Gamma_st2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_LogGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    // if(nuis[1]<0||nuis[0]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_LogGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
    
}


void Comp_Pair_BinomGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    char *f_name = "Comp_Pair_BinomGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
    
}

void Comp_Pair_BinomnegGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
  //  if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    char *f_name = "Comp_Pair_BinomnegGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
     int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
    
}

void Comp_Pair_LogLogistic_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                     double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if( nuis[0]>1 || nuis[0]<0){*res=LOW; return;}
    char *f_name = "Comp_Pair_LogLogistic_st2_OCL";
    
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
    
}


void Comp_Pair_Logistic_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                   double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if( nuis[2]<=0)  {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Logistic_st2_OCL";
    
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
    
}


void Comp_Pair_Weibull_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
   // if(nuis[2]<=0||sill<0||sill>1) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Weibull_st2_OCL";
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_T_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                               double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    double df=nuis[0];
    if( df<0||df>0.5){*res=LOW; return;}
    char *f_name = "Comp_Pair_T_st2_OCL";
    int *int_par;
    double *dou_par;
   int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}





void Comp_Pair_TWOPIECET_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                         double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    double eta=nuis[3];  //skewness parameter
    double df=nuis[0];
    double sill=nuis[2];
    double nugget=nuis[1];
    if(sill<0||nugget<0||nugget>=1 ||fabs(eta)>1|| df >0.5||df<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_TWOPIECET_st2_OCL";
    
    int *int_par;
    double *dou_par;
   int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}



void Comp_Pair_TWOPIECEGauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
       double eta,sill,nugget;
    eta=nuis[2];  //skewness parameter
    sill=nuis[1];
    nugget=nuis[0];
    if( fabs(eta)>1|| sill<0||nugget<0||nugget>=1) {*res=LOW;  return;} 
    char *f_name = "Comp_Pair_TWOPIECEGauss_st2_OCL";
    
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}




void Comp_Pair_Kumaraswamy_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                         double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if(nuis[2]<0||nuis[3]<0) {*res=LOW;  return;}
    char *f_name = "Comp_Pair_Kumaraswamy_st2_OCL";
    
    int *int_par;
    double *dou_par;
       int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
     Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}



void Comp_Pair_Tukeyh_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                   double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //double sill=nuis[1];
    //double nugget=nuis[0];
    double tail=nuis[2];
   if( tail>0.5||tail<0 ){*res=LOW; return;}

    char *f_name = "Comp_Pair_Tukeyh_st2_OCL";
    
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}
void Comp_Pair_Gauss_misp_T_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    if( nuis[0]>0.5 || nuis[0]<0||nuis[1]>=1 || nuis[1]<0){*res=LOW; return;}
    char *f_name = "Comp_Pair_Gauss_misp_T_st2_OCL";
    
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Gauss_misp_Pois_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
      double nugget=nuis[0];
    if(nugget<0||nugget>=1){*res=LOW; return;}

    char *f_name = "Comp_Pair_Gauss_misp_Pois_st2_OCL";
    
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}


void Comp_Pair_Pois_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev)
{
    //if( CheckCor(cormod,par)==-2){*res=LOW; return;} 

     double nugget=nuis[0];
    if(nugget<0||nugget>=1){*res=LOW; return;}
    char *f_name = "Comp_Pair_Pois_st2_OCL";
    
    int *int_par;
    double *dou_par;
    int_par = (int*)Calloc((50), int *);
    dou_par = (double*)Calloc((50), double *);
    param_st_OCL(cormod,NN,par,weigthed,nuis,int_par,dou_par);
    exec_kernel_st_dyn(coordx,coordy, coordt,mean,data, int_par, dou_par, local_wi,dev,res,f_name,ns,NS);
    Free(int_par);
    Free(dou_par);
    if(!R_FINITE(*res))*res = LOW;
}
