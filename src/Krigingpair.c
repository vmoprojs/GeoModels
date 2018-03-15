
#include "header.h"
/*************************************************************************************************/
/********** functions for pairwise binary  predictor ****************************/
/*************************************************************************************************/
double cond_exp_bin(int *cormod,double data_i,double data_j,double lags_i,double lags_j,double lags,double *nuis,double *par,double psm)
{
double   a,b,p010,p111,p110,p011,p001,p101,p100,p000;
p111=ptnorm(1,cormod,lags_i,lags_j,lags,0,0,0,nuis,par,0);
                   p110=ptnorm(2,cormod,lags_i,lags_j,lags,0,0,0,nuis,par,0);
                   p101=ptnorm(3,cormod,lags_i,lags_j,lags,0,0,0,nuis,par,0);
                   p011=ptnorm(4,cormod,lags_i,lags_j,lags,0,0,0,nuis,par,0); 
                   p010=psm-p111-p110-p011;
                   p001=psm-p111-p101-p011;
                   p100=psm-p111-p101-p110;
                   p000=2*p111-3*psm+p101+p011+p110+1;
                   a=pow(p100,(1-data_i)*(1-data_j))*
                       pow(p101,(1-data_i)*(  data_j))*
                       pow(p110,(  data_i)*(1-data_j))*
                       pow(p111,(  data_i)*(  data_j));
                   b=pow(p000,(1-data_i)*(1-data_j))*
                       pow(p001,(1-data_i)*(  data_j))*
                       pow(p010,(  data_i)*(1-data_j))*
                       pow(p011,(  data_i)*(  data_j));
                   return(a/(a+b));    

}


/*************************************************************************************************/
/********** functions for pairwise skew gaussian  predictor ****************************/
/*************************************************************************************************/
/* integrand  in  skew pairwise predictoor*/
double int_gen_skew(double x,double data_i, double data_j,double c_0i,double c_0j,double rho,double *nuis)
{
    double res=0.0,triv,biv;
    triv=triv_skew(x,c_0i,c_0j,rho,data_i,data_j,nuis);  // x iz z_0
    biv=biv_skew(rho,data_i,data_j,0,0,nuis[2],nuis[3]); 
    res=x*triv/biv;
    return(res);///(pow(2,alpha-1)*gamma(alpha)*pow(supp,2*alpha)));
}
// function skew predictor  to integrate
void integr_gen_skew(double *x, int n, void *ex){
    int i;double data_i,data_j,c_0i,c_0j,rho, *nuis;
    nuis=(double *) Calloc(4,double);
    data_i =    ((double*)ex)[0];  data_j = ((double*)ex)[1];  
    c_0i =     ((double*)ex)[2]; c_0j =     ((double*)ex)[3];  
    rho =     ((double*)ex)[4];  
    nuis[0]=((double*)ex)[5];   nuis[1]=((double*)ex)[6]; nuis[2]=((double*)ex)[7]; nuis[3]=((double*)ex)[8];                 
    for (i=0;i<n;i++) {x[i]=int_gen_skew(x[i],data_i,data_j,c_0i,c_0j,rho,nuis);}
    Free(nuis);
    return;
}
// function computing pairwise skew gaussian predictor
double cond_exp_skew(double c_0i,double c_0j,double rho, double data_i,double data_j,double *nuis) { 
    double ex[9], bound, epsabs, epsrel, result, abserr, *work;
    int inf,neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;         /* as instructed in WRE */
    iwork =   (int *) Calloc(subdiv, int); 
    work = (double *) Calloc(lenw, double); 
    ex[0] = data_i; ex[1] = data_j; 
    ex[2] = c_0i;ex[3]=c_0j;
    ex[4]=rho; 
    ex[5]=nuis[0];ex[6]=nuis[1];ex[7]=nuis[2];ex[8]=nuis[3];   
    bound=0;
    inf= 2;
    Rdqagi(integr_gen_skew, (void *) &ex,
               &bound, &inf, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
    Free(iwork);Free(work);
    return(result);
}

/*************************************************************************************************/
/*************************************************************************************************/
/*************************************************************************************************/





/*kriging based on pairwise */
double uni_kri_space(double *coordx, double *coordy, int *cormod, double *data,double locx,  double locy,int *n,int *ncoord, 
  double *maxdist,double *nuis,int *model,double *par, double *radius,int *type_krig,int *type_dist,int *lin_opt,int *weighted)
{
    int i=0;
    double mean=0.0,var=nuis[2],nugget=nuis[1],var1=0.0,psm=0.0;
    double c_0i=0.0,sum1=0.0,sum2=0.0,lags_i=0.0,weights=1.0,
       weights_i=1,res=0.0; 
   double t1,t2=0.0,sk=0.0,sk2=0.0,tm,consta=0.0,tl=0.0;   
      
     

 /*############### setting means ###################################################*/   
      switch(*model) 
              {
                 case 1: // gaussian
                     mean=nuis[0];
                     var1=var;
                     break;
                case 2: // binary o binomial
                case 11:
                     psm=pnorm((nuis[0]-0)/sqrt(var+nugget),0,1,1,0);
                     mean=n[0]*psm;
                     var1=n[0]*psm*(1-psm);
                     break;
                case 14: //geometric
                     psm=pnorm((nuis[0]-0)/sqrt(var+nugget),0,1,1,0); //geometric case
                     mean=(1-psm)/psm; 
                     var1=(1-psm)/pow(psm,2);
                     break; 

                case 9:  //tukey
                 tl=nuis[4]; sk=nuis[3];
                 t1=1-tl;t2=pow(t1,2)-pow(tl,2);sk2=pow(sk,2);
                 if(sk==0)  {tm=0; consta= ( 2/(1-tl*(2))   -1)/sqrt(pow(t1,2) - pow(tl,2));}                        
                 else            {tm=(exp(sk2/(2*t1))-1)/(sk*sqrt(t1));
                     consta= ((exp(sk2 * 2/(1-2*tl)) -
                        2* exp( sk2 *0.5)+1))/(sk2*sqrt(t2)) - pow(tm,2);}
                 mean=nuis[0] + sqrt(var) * tm;
                 var1=var;

                break;
                case 10:  //skew
                     sk=nuis[3];
                     mean=  nuis[0]+sqrt(2/M_PI)*sk;
                     var1=var+(1-2/M_PI)*pow(sk,2);
                     break;
              }
/*##################################################################*/


    for(i=0;i<(ncoord[0]-1);i++){
        lags_i=dist(*type,coordx[i],locx,coordy[i],locy,radius[0]);
          if(lags_i<=*maxdist){
          if(*weighted) {
             /*weights_i=CorFunStable(lags_i, 2, *maxdist/3);*/
             weights_i=CorFunBohman(lags_i,*maxdist);
          }

 


      if(*lin_opt){ //// pairwise LINEAR
         /*************************************************************/ 
              switch(*model) 
              {
                 case 1: // gaussian
                  c_0i=CorFct(cormod,lags_i,0,par,0,0);
                 
                  break;
                case 2: // binary o binomial
                case 11:
                   break;
                  case 14: // geometric
                    break;

             case 9:  //tukey gaussian
                        break;

             case 10:  //skew gaussian
                 break;
               }
              /*************************************************************/ 

              if(*type_krig==0)  {  //simple kriging
                  sum1=sum1 + (c_0i * (data[i]-mean))*weights;
                  sum2=sum2 + weights;
                }
            //if(*type_krig==1)  ordinary kriging
            }
     else{   // pairwise COND Expectation
               switch(*model) //
              {
           // case 1: // gaussian
             case 2: // binary 
                    
                  break; 
             case 9:  //tukey 
            // c_0i=CorFct(cormod,lags_i,0,par,0,0);c_0j=CorFct(cormod,lags_j,0,par,0,0);;
            // rho=CorFct(cormod,lags,0,par,0,0);
            // det=pow(var1+nugget,2)-pow(rho*var1,2);   

            //  aij=(var1*c_0i*(var1+nugget)-c_0j*rho*pow(var1,2))/det;
            //  aji=(var1*c_0j*(var1+nugget)-c_0i*rho*pow(var1,2))/det;

            //  b=aij*inv_tuk((data[i]- nuis[0])/sqrt(var)) +aji*inv_tuk((data[j]- nuis[0])/sqrt(var));
             // vart=1-(aij*c_0i+aji*c_0j*c_0j);
             // a=1-tl*vart;
             // sum1=sum1 + (mean+ sqrt(var)/(sk*sqrt(a))*exp(0.5*tl*pow(b,2)/a)*
             //           (exp(0.5*(sk2*vart+2*sk*b)/a)-1))*weights;
             // sum2=sum2 + weights;
              break;
            case 10:// skew gaussian
              
               break; 
            case 11:
             break;
             case 14: //geometric case
               break;      

    
           }}}}

    if(*lin_opt)   res=mean+(sum1/sum2);//+sum3/sum4)/2;
    if(!*lin_opt)  res=      sum1/sum2;
//Rprintf("qua\n");
    return(res);
}

/*#################################################################################*/
/*#################################################################################*/
/*#################################################################################*/


/*kriging based on pairwise */
double pair_kri_space(double *coordx, double *coordy, int *cormod, double *data,double locx,  double locy,int *n,int *ncoord, 
  double *maxdist,double *nuis,int *model,double *par, double *radius,int *type_krig,int *type_dist,int *lin_opt,int *weighted)
{
    int i=0,j=0;
    double mean=0.0,var=nuis[2],nugget=nuis[1],var1=0.0,psm=0.0,psi=0.0,psj=0.0,psij=0.0,k=0.0,a=0.0,b=0.0;
    double aij=0.0,aji=0.0,c_0i=0.0,c_0j=0.0,det=0.0,
       sum1=0.0,sum2=0.0,lags=0.0,lags_i=0.0,lags_j=0.0,rho=0.0,weights=1.0,
       weights_i=1,weights_j=1,res=0.0; 
   double t1,t2=0.0,sk=0.0,sk2=0.0,tm,consta=0.0,tl=0.0;   
      
     

 /*############### setting means ###################################################*/   
      switch(*model) 
              {
                 case 1: // gaussian
                     mean=nuis[0];
                     var1=var;
                     break;
                case 2: // binary o binomial
                case 11:
                     psm=pnorm((nuis[0]-0)/sqrt(var+nugget),0,1,1,0);
                     mean=n[0]*psm;
                     var1=n[0]*psm*(1-psm);
                     break;
                case 14: //geometric
                     psm=pnorm((nuis[0]-0)/sqrt(var+nugget),0,1,1,0); //geometric case
                     mean=(1-psm)/psm; 
                     var1=(1-psm)/pow(psm,2);
                     break; 

                case 9:  //tukey
                 tl=nuis[4]; sk=nuis[3];
                 t1=1-tl;t2=pow(t1,2)-pow(tl,2);sk2=pow(sk,2);
                 if(sk==0)  {tm=0; consta= ( 2/(1-tl*(2))   -1)/sqrt(pow(t1,2) - pow(tl,2));}                        
                 else            {tm=(exp(sk2/(2*t1))-1)/(sk*sqrt(t1));
                     consta= ((exp(sk2 * 2/(1-2*tl)) -
                        2* exp( sk2 *0.5)+1))/(sk2*sqrt(t2)) - pow(tm,2);}
                 mean=nuis[0] + sqrt(var) * tm;
                 var1=var;

                break;
                case 10:  //skew
                     sk=nuis[3];
                     mean=  nuis[0]+sqrt(2/M_PI)*sk;
                     var1=var+(1-2/M_PI)*pow(sk,2);
                     break;
              }
/*##################################################################*/


    for(i=0;i<(ncoord[0]-1);i++){
        lags_i=dist(*type,coordx[i],locx,coordy[i],locy,radius[0]);
          if(lags_i<=*maxdist){
          if(*weighted) {
            /* weights_i=CorFunStable(lags_i, 2, *maxdist/3);*/
            weights_i=CorFunBohman(lags_i,*maxdist);
           }

     for(j=i+1; j<ncoord[0];j++){
          lags_j=dist(*type,coordx[j],locx,coordy[j],locy,radius[0]);
          if(lags_j<=*maxdist){

          if(*weighted) {
            /* weights_j=CorFunStable(lags_j, 2, *maxdist/3);*/
            weights_j=CorFunBohman(lags_j,*maxdist);

                        weights=weights_i*weights_j;}
          lags=dist(*type,coordx[i],coordx[j],coordy[i],coordy[j],radius[0]);
 


      if(*lin_opt){ //// pairwise LINEAR
         /*************************************************************/ 
              switch(*model) 
              {
                 case 1: // gaussian
                  c_0i=CorFct(cormod,lags_i,0,par,0,0);c_0j=CorFct(cormod,lags_j,0,par,0,0);;
                  rho=CorFct(cormod,lags,0,par,0,0);
                  break;
                case 2: // binary o binomial
                case 11:
                  k=psm*(1-psm);
                  c_0i=(pbnorm(cormod,lags_i,0,mean,mean,nugget,var,par,0)-pow(psm,2))/k;
                  c_0j=(pbnorm(cormod,lags_j,0,mean,mean,nugget,var,par,0)-pow(psm,2))/k;
                  rho=(pbnorm(cormod,lags,0,mean,mean,nugget,var,par,0)-pow(psm,2))/k;
                break;
                  case 14: // geometric
                 k=2*psm;psi=pbnorm(cormod,lags_i,0,mean,mean,nugget,var,par,0);
                 c_0i=(psi*(1-2*psm+psi)-pow(psm-psi,2))*((1-psm)/(k-psi));
                 psj=pbnorm(cormod,lags_j,0,mean,mean,nugget,var,par,0);
                 c_0j=(psj*(1-2*psm+psj)-pow(psm-psj,2))*((1-psm)/(k-psj));
                 psij=pbnorm(cormod,lags,0,mean,mean,nugget,var,par,0);
                 rho=(psij*(1-2*psm+psij)-pow(psm-psij,2))*((1-psm)/(k-psij));
               break;

             case 9:  //tukey gaussian
               psi=CorFct(cormod,lags_i,0,par,0,0);
               a=(1+psi)/(1-tl*(1+psi));b=0.5* (1-tl*(1-psi*psi));
               if(sk==0)   c_0i=(( a -2* b)/sqrt(t2) )/consta;
               else             c_0i=((exp(sk2 * a) -2* exp( sk2 *b)+1)/(sk2*sqrt(t2)) - pow(tm,2))/consta;
               psj=CorFct(cormod,lags_j,0,par,0,0);
               a=(1+psj)/(1-tl*(1+psj));b=0.5* (1-tl*(1-psj*psj));
               if(sk==0)   c_0j=(( a -2* b)/sqrt(t2) )/consta;
               else             c_0j=((exp(sk2 * a) -2* exp( sk2 *b)+1)/(sk2*sqrt(t2)) - pow(tm,2))/consta;
               k=CorFct(cormod,lags,0,par,0,0);
               a=(1+k)/(1-tl*(1+k));b=0.5* (1-tl*(1-k*k));
               if(sk==0)   rho=(( a -2* b)/sqrt(t2) )/consta;
               else             rho=((exp(sk2 * a) -2* exp( sk2 *b)+1)/(sk2*sqrt(t2)) - pow(tm,2))/consta;
               break;

             case 10:  //skew gaussian
               psi=CorFct(cormod,lags_i,0,par,0,0);psj=CorFct(cormod,lags_j,0,par,0,0);
               k=CorFct(cormod,lags,0,par,0,0);
               a=2*pow(sk,2)/M_PI; b=pow(sk,2)*(1-2/M_PI);
               c_0i=(a*(sqrt(1-psi*psi) + psi*asin(psi)-1) + psi*var)/(var+b);
               c_0j=(a*(sqrt(1-psj*psj) + psj*asin(psj)-1) + psj*var)/(var+b);
               rho =(a*(sqrt(1-k*k) + k*asin(k)-1) + k*var)/(var+b);
               break;
               }
              /*************************************************************/ 

              det=pow(var1+nugget,2)-pow(rho*var1,2);   
              aij=(var1*c_0i*(var1+nugget)-c_0j*rho*pow(var1,2))/det;
              aji=(var1*c_0j*(var1+nugget)-c_0i*rho*pow(var1,2))/det;
              if(*type_krig==0)  {  //simple kriging
                  sum1=sum1 + (aij * (data[i]-mean) + aji * (data[j]-mean))*weights;
                  sum2=sum2 + weights;
                }
            //if(*type_krig==1)  ordinary kriging
            }
     else{   // pairwise COND Expectation
               switch(*model) //
              {
           // case 1: // gaussian
             case 2: // binary 
                   sum1=sum1 + cond_exp_bin(cormod,data[i],data[j],lags_i,lags_j,lags,nuis,par,psm)*weights;
                   sum2=sum2 + weights;
                  
                  break; 
             case 9:  //tukey 
            // c_0i=CorFct(cormod,lags_i,0,par,0,0);c_0j=CorFct(cormod,lags_j,0,par,0,0);;
            // rho=CorFct(cormod,lags,0,par,0,0);
            // det=pow(var1+nugget,2)-pow(rho*var1,2);   

            //  aij=(var1*c_0i*(var1+nugget)-c_0j*rho*pow(var1,2))/det;
            //  aji=(var1*c_0j*(var1+nugget)-c_0i*rho*pow(var1,2))/det;

            //  b=aij*inv_tuk((data[i]- nuis[0])/sqrt(var)) +aji*inv_tuk((data[j]- nuis[0])/sqrt(var));
             // vart=1-(aij*c_0i+aji*c_0j*c_0j);
             // a=1-tl*vart;
             // sum1=sum1 + (mean+ sqrt(var)/(sk*sqrt(a))*exp(0.5*tl*pow(b,2)/a)*
             //           (exp(0.5*(sk2*vart+2*sk*b)/a)-1))*weights;
             // sum2=sum2 + weights;
              break;
            case 10:// skew gaussian
              
              
               c_0i=CorFct(cormod,lags_i,0,par,0,0);c_0j=CorFct(cormod,lags_j,0,par,0,0);;
               rho=CorFct(cormod,lags,0,par,0,0);
               sum1=sum1 + cond_exp_skew(c_0i,c_0j,rho,data[i],data[j],nuis)*weights;
               sum2=sum2 + weights; 
               break; 
            case 11:
             break;
             case 14: //geometric case
               break;      

    
           }}}}}}

    if(*lin_opt)   res=mean+(sum1/sum2);//+sum3/sum4)/2;
    if(!*lin_opt)  res=      sum1/sum2;
//Rprintf("qua\n");
    return(res);
}



void pair_k(double *coordx, double *coordy, double *coordt,int *cormod, double *data,double *locx,  
  double *locy, double *loct,int *n,int *ncoord,int *nloc,int *tloc,int *ntime, double *maxdist,
  double *maxtime,double *nuis,int *model,double *par,int *pair,double *radius, int *type_krig,int *type_dist,int *spt,
  int *lin_opt,int *weighted,double *res)
{
int i,h=0;


if(pair[0])
{
if(!*spt){  //spatial
for(i=0;i<*nloc;i++){
    res[h]=pair_kri_space(coordx, coordy, cormod, data,locx[i],locy[i],n,ncoord,maxdist,
      nuis,model,par,radius,type_krig,type_dist,lin_opt,weighted);
    h++;}}

/*
if(*spt)  //spacetime
{
for(j=0;j<*tloc;j++){
for(i=0;i<*nloc;i++){
    res[h]=pair_kri_spacetime(coordx, coordy,coordt, cormod, data,locx[i],locy[i],loct[j],n,ncoord,ntime,
     maxdist,maxtime,nuis,model,par,radius,type_krig,type_dist,weighted);
    h++;
}}}*/
  }
  if(!pair[0])
    {

if(!*spt){  //spatial
for(i=0;i<*nloc;i++){
    res[h]=uni_kri_space(coordx, coordy, cormod, data,locx[i],locy[i],n,ncoord,maxdist,
      nuis,model,par,radius,type_krig,type_dist,lin_opt,weighted);
    h++;}}

/*
if(*spt)  //spacetime
{
for(j=0;j<*tloc;j++){
for(i=0;i<*nloc;i++){
    res[h]=uni_kri_spacetime(coordx, coordy,coordt, cormod, data,locx[i],locy[i],loct[j],n,ncoord,ntime,
     maxdist,maxtime,nuis,model,par,radius,type_krig,type_dist,weighted);
    h++;
}}}*/


    }
}

/*  theoretical MSE pairwise kriging*/
double MSE_pair_space(double *coordx, double *coordy, int *cormod,double locx,  double locy,int *n,int *ncoord,
              double *maxdist,double *nuis,int *model,double *par, double *radius,int *type_krig,int *type_dist,
              int *lin_opt, int *weighted)
{
    int i=0,j=0,l=0,k=0;// check the paramaters range:
    double lagsij=0.0,lagslk=0.0,lagsjk=0.0,lagsik=0.0,lagsil=0.0,lagsjl=0.0,lags_i=0.0,lags_j=0.0,lags_l=0.0,lags_k=0.0;
    double c_0i=0.0,c_0j=0.0,c_0k=0.0,c_0l=0.0,cij=0.0,clk=0.0,cil=0.0,cjk=0.0,cjl=0.0,cik=0.0;
    double sum2ij=0.0,sum3ij=0.0,sum2ijlk=0.0;
    double wi=1,wj=1,wij=1,wl=1,wk=1,wlk=1;
    double detij,detlk,aijlk=0.0;
    double aij=0.0,aji=0.0,alk,akl;
    double mean,var,nugget,varpred,cov,res=0.0,psi=0.0,psj=0.0,a=0.0,b=0.0,kk=0.0,psm=0.0,p11=0.0;
    mean=nuis[0];var=nuis[2];nugget=nuis[1];
    if(*model==10){ a=2*pow(nuis[3],2)/M_PI; b=pow(nuis[3],2)*(1-2/M_PI); } //skew gaussian
    if(*model==2||*model==11||*model==14)  psm=pnorm((nuis[0]-0)/sqrt(nuis[1]+nuis[2]),0,1,1,0);
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    for(i=0;i<(ncoord[0]-1);i++){
        lags_i=dist(*type,coordx[i],locx,coordy[i],locy,radius[0]);
        if(lags_i<=*maxdist){

            for(j=(i+1);j<ncoord[0];j++){
                lags_j=dist(*type,coordx[j],locx,coordy[j],locy,radius[0]);
                if(lags_j<=*maxdist){

                    lagsij=dist(*type,coordx[i],coordx[j],coordy[i],coordy[j],radius[0]);
                    if(weighted[0]) { /*wi=CorFunStable(lags_i, 2, *maxdist/3);
                                      wj=CorFunStable(lags_j, 2, *maxdist/3);*/
                                     wi=CorFunBohman(lags_i,*maxdist);
                                     wj=CorFunBohman(lags_j,*maxdist);wij=wi*wj;}



if(*lin_opt){ 
    
               if(*model==1)     {cij=CorFct(cormod,lagsij,0,par,0,0);c_0i=CorFct(cormod,lags_i,0,par,0,0);c_0j=CorFct(cormod,lags_j,0,par,0,0);}
               if(*model==10)   {psi=CorFct(cormod,lags_i,0,par,0,0);psj=CorFct(cormod,lags_j,0,par,0,0); kk=CorFct(cormod,lagsij,0,par,0,0);
                                 c_0i=(a*(sqrt(1-psi*psi) + psi*asin(psi)-1) + psi*var)/(var+b);
                                 c_0j=(a*(sqrt(1-psj*psj) + psj*asin(psj)-1) + psj*var)/(var+b);
                                 cij =(a*(sqrt(1-kk*kk) + kk*asin(kk)-1) + kk*var)/(var+b);}   
               if(*model==2||*model==11)    { c_0i=  pbnorm(cormod,lags_i,0,mean,mean,nugget,var,par,0)-psm*psm;

                                              c_0j=  pbnorm(cormod,lags_j,0,mean,mean,nugget,var,par,0)-psm*psm;           
                                              cij=  pbnorm(cormod,lagsij,0,mean,mean,nugget,var,par,0)-psm*psm;}
               if(*model==14)  { p11=pbnorm(cormod,lagsij,0,mean,mean,nugget,var,par,0);
                                cij=p11*(1-2*psm+p11 - pow(psm-p11,2)) ;}                          

                    detij=pow(1+nugget,2)-pow(cij,2);
                    aij=(c_0i-c_0j*cij)/detij; aji=(c_0j-c_0i*cij)/detij;
                   
                    sum2ij=sum2ij + wij;
                    sum3ij=sum3ij + (aij*c_0i+aji*c_0j)*wij;
                    /************************************************************************************************************/
                    for(l=0;l<(ncoord[0]-1);l++){
                        lags_l=dist(*type,coordx[l],locx,coordy[l],locy,radius[0]);
                        if(lags_l<=*maxdist){

                            for(k=(l+1);k<ncoord[0];k++){
                                    lags_k=dist(*type,coordx[k],locx,coordy[k],locy,radius[0]);
                                if(lags_k<=*maxdist) {

                                    lagslk=dist(*type,coordx[l],coordx[k],coordy[l],coordy[k],radius[0]);
                                    if(weighted[0]) {
                                      /*wl=CorFunStable(lags_l, 2, *maxdist/3);
                                      wk=CorFunStable(lags_k, 2, *maxdist/3);*/
                                      wl=CorFunBohman(lags_l,*maxdist);
                                      wk=CorFunBohman(lags_k,*maxdist);wlk=wl*wk;
                                      }
                                    
          if(*model==1)  {c_0l=CorFct(cormod,lags_l,0,par,0,0);c_0k=CorFct(cormod,lags_k,0,par,0,0);clk=CorFct(cormod,lagslk,0,par,0,0);}
          if(*model==10) {psi=CorFct(cormod,lags_l,0,par,0,0);psj=CorFct(cormod,lags_k,0,par,0,0); kk=CorFct(cormod,lagslk,0,par,0,0);
                                 c_0l=(a*(sqrt(1-psi*psi) + psi*asin(psi)-1) + psi*var)/(var+b);
                                 c_0k=(a*(sqrt(1-psj*psj) + psj*asin(psj)-1) + psj*var)/(var+b);
                                 clk =(a*(sqrt(1-kk*kk) + kk*asin(kk)-1) + kk*var)/(var+b);}  
          if(*model==2||*model==11) {c_0l=pbnorm(cormod,lags_l,0,mean,mean,nugget,var,par,0)-psm*psm; 
                                     c_0k=pbnorm(cormod,lags_k,0,mean,mean,nugget,var,par,0)-psm*psm;
                                     clk=pbnorm(cormod,lagslk,0,mean,mean,nugget,var,par,0)-psm*psm;}    
             if(*model==14)  { p11=pbnorm(cormod,lags_l,0,mean,mean,nugget,var,par,0);c_0l=p11*(1-2*psm+p11 - pow(psm-p11,2) );
                              p11=pbnorm(cormod,lags_k,0,mean,mean,nugget,var,par,0);c_0k=p11*(1-2*psm+p11 - pow(psm-p11,2)) ;
                              p11=pbnorm(cormod,lagslk,0,mean,mean,nugget,var,par,0);clk=p11*(1-2*psm+p11 - pow(psm-p11,2)) ;}                                                   

                                    detlk=pow(1+nugget,2)-pow(clk,2);     
                                    alk=(c_0l-c_0k*clk)/detlk;akl=(c_0k-c_0l*clk)/detlk;
    
                                    lagsjl=dist(*type,coordx[j],coordx[l],coordy[j],coordy[l],radius[0]);
                                    lagsjk=dist(*type,coordx[j],coordx[k],coordy[j],coordy[k],radius[0]);
                                    lagsik=dist(*type,coordx[i],coordx[k],coordy[i],coordy[k],radius[0]);
                                    lagsil=dist(*type,coordx[i],coordx[l],coordy[i],coordy[l],radius[0]);


              if(*model==1)       {cjl= CorFct(cormod,lagsjl,0,par,0,0);
                                   cjk= CorFct(cormod,lagsjk,0,par,0,0);
                                   cik= CorFct(cormod,lagsik,0,par,0,0);
                                   cil= CorFct(cormod,lagsil,0,par,0,0);}
              if(*model==10)      {psi=CorFct(cormod,lagsjl,0,par,0,0);cjl=(a*(sqrt(1-psi*psi) + psi*asin(psi)-1) + psi*var)/(var+b);
                                   psi=CorFct(cormod,lagsjk,0,par,0,0);cjk=(a*(sqrt(1-psi*psi) + psi*asin(psi)-1) + psi*var)/(var+b); 
                                   psi=CorFct(cormod,lagsik,0,par,0,0);cik=(a*(sqrt(1-psi*psi) + psi*asin(psi)-1) + psi*var)/(var+b);
                                   psi=CorFct(cormod,lagsil,0,par,0,0);cil=(a*(sqrt(1-psi*psi) + psi*asin(psi)-1) + psi*var)/(var+b);

                                    }  
              if(*model==2||*model==11) {cjl=pbnorm(cormod,lagsjl,0,mean,mean,nugget,var,par,0)-psm*psm;
                                         cjk=pbnorm(cormod,lagsjk,0,mean,mean,nugget,var,par,0)-psm*psm;
                                         cik=pbnorm(cormod,lagsik,0,mean,mean,nugget,var,par,0)-psm*psm;
                                         cil=pbnorm(cormod,lagsil,0,mean,mean,nugget,var,par,0)-psm*psm;}
                   
              if(*model==14)  { p11=pbnorm(cormod,lagsjl,0,mean,mean,nugget,var,par,0);cjl=p11*(1-2*psm+p11 - pow(psm-p11,2) );
                                p11=pbnorm(cormod,lagsjk,0,mean,mean,nugget,var,par,0);cjk=p11*(1-2*psm+p11 - pow(psm-p11,2) );
                                p11=pbnorm(cormod,lagsik,0,mean,mean,nugget,var,par,0);cik=p11*(1-2*psm+p11 - pow(psm-p11,2) );
                                p11=pbnorm(cormod,lagsil,0,mean,mean,nugget,var,par,0);cil=p11*(1-2*psm+p11 - pow(psm-p11,2) );}

                                    aijlk= aij*alk*cil   +  aij*akl*cik  + aji*alk*cjl  +  aji*akl*cjk;
                                    sum2ijlk = sum2ijlk + aijlk*(wij*wlk);
                                    
              }}}}}}}}
             }    
   if(!*lin_opt){ 



          if(*model==2||*model==11)    { p11=1;
            sum2ij=sum2ij + wij;
            sum3ij=sum3ij + (aij*c_0i+aji*c_0j)*wij;}

                  for(l=0;l<(ncoord[0]-1);l++){
                        lags_l=dist(*type,coordx[l],locx,coordy[l],locy,radius[0]);
                        if(lags_l<=*maxdist){

                            for(k=(l+1);k<ncoord[0];k++){
                                    lags_k=dist(*type,coordx[k],locx,coordy[k],locy,radius[0]);
                                if(lags_k<=*maxdist) {

                                    lagslk=dist(*type,coordx[l],coordx[k],coordy[l],coordy[k],radius[0]);
                                    if(weighted[0]) {wl=CorFunBohman(lags_l,*maxdist);
                                                     wk=CorFunBohman(lags_k,*maxdist);
                                                     wlk=wl*wk;}
         
             if(*model==2||*model==11)    { p11=1;}
        }}}} 
      }

    cov = sum3ij/sum2ij;
    varpred = sum2ijlk/pow(sum2ij,2);
      if(*model==1)  res=var*(1 - 2 * cov + varpred);
      if(*model==10) res=(var+b)*(1 - 2 * cov + varpred);
    return(res);
    
} 

/* main function for  Theoretical MSE  pairwise kriging*/
void mse_pair_k(double *coordx, double *coordy, double *coordt,int *cormod,double *locx,
            double *locy, double *loct,int *n,int *ncoord,int *nloc,int *tloc,int *ntime, double *maxdist,
            double *maxtime,double *nuis,int *model,double *par,double *radius, int *type_krig,int *type_dist,int *spt,
           int *lin_opt,int *weighted,double *res)
{
    int i=0,j=0,h=0;
  
    if(!*spt){
        for(i=0;i<*nloc;i++){
            res[h]=MSE_pair_space(coordx, coordy, cormod,locx[i],locy[i],n,ncoord,maxdist,
                                  nuis,model,par,radius,type_krig,type_dist,lin_opt,weighted);
            h++;
        }}
    if(*spt)
    {
        for(j=0;j<*tloc;j++){
            for(i=0;i<*nloc;i++){
                //res[h]=MSE_pair_spacetime(coordx, coordy,coordt, cormod, data,locx[i],locy[i],loct[j],n,ncoord,ntime,
                  //                        maxdist,maxtime,nuis,model,par,radius,type_krig,type_dist,weighted);
                h++;
            }}}
}
