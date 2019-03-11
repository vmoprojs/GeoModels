

void Comp_Cond_Weibull2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
 double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local)
{
  double lags=0.0,a1,a2;int k=0;   int nn=ncoord[0];
  double PP1=0.0,PP2=0.0,MM=0.0,sum=0.0,cc=0.0,second=0.0,first=0.0;

  double m=nuis[2]; if(m<0||par[0]<0){*res=LOW;  return;}
  double nu=1/gammafn(1+1/m);

  double x1 = data[0]/exp(mean[0]);
  double xn = data[nn-1]/exp(mean[nn-1]);
  double lags01=fabs(coordx[0]-coordx[1]);double lagsn1n=fabs(coordx[nn-1]-coordx[nn]]);
  /********************************/
  for(k=1; k<(nn-1);k++){
              lags=fabs(coordx[k-1]-coordx[k]); a1=exp(-lags/par[0]);
              lags=fabs(coordx[k]-coordx[k+1]); a2=exp(-lags/par[0]);
      sum= sum+ R_pow(data[k]/exp(mean[k]),m)*(1-R_pow(a1*a2,2))/((1-a1*a1)*(1-a2*a2));
    }
    second=-R_pow(nu,-m)*( R_pow(x1,m)/(1-R_pow(exp(-lags01 /par[0]),2))  
                   +R_pow(xn,m)/(1-R_pow(exp(-lagsn1n/par[0]),2))
                         +sum);
/********************************/
  for(k=0; k<nn;k++) {MM=MM-mean[k];PP0=PP0+ (m-1)*log(data[k]/exp(mean[k]));}
/*******+********************************************/
  for(k=0; k<(nn-1);k++){
        lags=fabs(coordx[k]-coordx[k+1]);
           cc=exp(-lags/par[0]);
            PP1=PP1 + log(1-cc*cc);
            PP2=PP2 + log(bessel_i( 2*fabs(cc)*R_pow(
                                 (data[k]/exp(mean[k]))*
                                 (data[k+1]/exp(mean[k+1])),m)/((1-cc*cc)*R_pow(nu,m)) , 0 ,1));
    }
  first= nn*log(m) + PP0 - ( nn*m*log(nu) + PP1);

  //Rprintf(" %f %f %f %f \n", MM,first,second,PP2);

  *res=MM+first+second+ PP2;
    if(!R_FINITE(*res)) *res = LOW;
  return;
}
