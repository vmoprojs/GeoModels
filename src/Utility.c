#include "header.h"

 
/***********************Legendre polinomial****************************************/

void lgnd (int *lmax,double *x, double *p){

    int l;
    p[0] = 1;
    p[1] = x[0];
    for (l = 1; l < *lmax; ++l)
    {
        p[l+1] = ((2*l+1)*x[0]*p[l]-l*p[l-1])/(l+1);
    }
  }
/*************************************************************************************/
/**************** matrix functions ***************************************************/
/*************************************************************************************/


void ludcmp(double **a, int n, int *indx, double *dd)
{
    int i,imax=0,j,k;
    double   big,dum,sum,temp,*vv;
    vv=(double *) Calloc(n,double);
    *dd = 1.0;
    for (i=0;i<n;i++)
    {
        big = 0.0;
        for (j=0;j<n;j++)
        {if ((temp=fabs(a[i][j])) > big) big = temp;}
        vv[i] = 1.0/big;}
    for (j=0;j<n;j++){
        for (i=0;i<j;i++){
            sum = a[i][j];
            for (k=0;k<i;k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;}
        big = 0.0;
        for (i=j;i<n;i++){
            sum = a[i][j];
            for (k=0;k<j;k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum=vv[i]*fabs(sum)) >= big){
                big = dum;
                imax = i;}}
        if (j != imax){
            for (k=0;k<n;k++){
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;}
            *dd = -(*dd);
            vv[imax] = vv[j];}
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = MAXERR;
        if (j != n-1){
            dum = 1.0 / a[j][j];
            for (i=j+1;i<n;i++) a[i][j] *= dum;}
    }
    Free(vv);
}


void lubksb(double **a, int n, int *indx, double *b)
{
    int i,ip,j,ii=-1;
    double   sum;
    for (i=0;i<n;i++){
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii>=0)
            for (j=ii;j<i;j++) sum -= a[i][j] * b[j];
        else if (sum) ii = i;
        b[i] = sum;}
    for (i=n-1;i>=0;i--){
        sum = b[i];
        for (j=i+1;j<n;j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];}
}




void Matrix_prod_vec(double **a,double *b,double *c,int n)
{
    int i,j,k;double sum=0;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            for(k=0; k<n; k++) {sum=sum+ a[i][k]*b[k];}
            c[j]=sum;sum=0;
    }}
}



void Matrix_prod(double **a,double **b,double **c,int n)
{
    int i,j,k;double sum=0;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            for(k=0; k<n; k++) {sum=sum+ a[i][k]*b[k][j];}
            c[i][j]=sum;sum=0;
    }}
}
/*compute the trace of  square matrix*/
double Trace(double **A,int n)
{
    int i;double x = 0;
    for(i=0;i<n;i++) x=x+A[i][i];
    return(x);
}
/*compute the quadratic form given a simmetric square matrix*/
double QFORM(double **A,double *x,double *y,int n)
{
    int i,j;double qf = 0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            qf=qf + x[i]*A[i][j]*y[j];}}
    return(qf);
}

/*compute the quadratic form given a simmetric square matrix*/
double QFORM2(double **A,double *x,double *y,int n, int m)
{
    int i,j;double qf = 0;
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            qf=qf + x[i]*A[i][j]*y[j];}}
    return(qf);
}

/*Transpose matrix of a square matrix*/
void Transpose(double **a,int n,double k)
{
    int i,j;double tmp;
    for (i=1;i<n;i++) {
        for (j=0;j<i;j++) {
            tmp = a[i][j]/k;
            a[i][j] = a[j][i]/k;
            a[j][i] = tmp;}}
}
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/




// Computes the Geodesic distance between to coordinates:
double Dist_geodesic(double loni, double lati, double lonj, double latj,double radius)
{
  double ai, bi, aj, bj, val=0.0,val2=0.0;
 if (loni == lonj && lati == latj) return val;
  ai = (lati)*M_PI/180;
  bi = (loni)*M_PI/180;
  aj = (latj)*M_PI/180;
  bj = (lonj)*M_PI/180;
  val = sin(ai) * sin(aj) + cos(ai) * cos(aj) * cos(bi - bj);
  if(val<= -1)  {val2=M_PI*radius;return(val2);}
  if(val>=1)    {val2=0;;return(val2);}
  val2 = acos(val)*radius; 
  return(val2);
}

double Dist_chordal(double loni, double lati, double lonj, double latj,double radius)
 {
   double val=0.0;/*,val1=0.0,ai, bi, aj, bj ;
   if (loni == lonj && lati == latj) return val;
   ai = (lati)*M_PI/180;
   bi = (loni)*M_PI/180;
   aj = (latj)*M_PI/180;
   bj = (lonj)*M_PI/180;
 val1=radius  *sqrt(R_pow(cos(ai) * cos(bi)-cos(aj)  *cos(bj) ,2) +
                    R_pow(cos(ai) * sin(bi)-cos(aj) * sin(bj) ,2)+
                         R_pow(sin(ai)-sin(aj) ,2));
 */
        val=radius*2*sin(Dist_geodesic(loni,lati,lonj,latj,1)/2);
 return(val);
 }



double dist(int type_dist,double coordx,double locx,double coordy,double locy,double radius)
{
double lags=0.0;

if(type_dist==0) lags=hypot(coordx-locx,coordy-locy);                        /*euclidean*/
if(type_dist==2) lags=Dist_geodesic(coordx,coordy,locx,locy,radius);           /*great circle*/
if(type_dist==1) lags=Dist_chordal(coordx,coordy,locx,locy,radius);      /*chordal*/

return(lags);
}


void Maxima_Minima_dist(double *res,double *coordx,double *coordy,int *nsize,int *type_dist,double *radius)
{
  double res1=0.0,res2=-LOW,lags=0.0;
  int i=0,j=0;
    for(i=0; i<(*nsize-1);i++){
        for(j=(i+1); j<*nsize;j++){
          lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*radius);
          res1 = fmax(res1, lags);
          res2 = fmin(res2, lags);
}}
res[0]=res2;res[1]=res1;
return;
}


void Maxima_Minima_time(double *res,double *coordt,int *nsize)
{
  double res1=0.0,res2=-LOW,lagt=0.0;
  int i=0,j=0;
    for(i=0; i<(*nsize-1);i++){
        for(j=(i+1); j<*nsize;j++){
               lagt=fabs(coordt[i]-coordt[j]);
          res1 = fmax(res1, lagt);
          res2 = fmin(res2, lagt);
}}
res[0]=res2;res[1]=res1;
return;
}



/*void RangeDist(double *max, double *min)
{
  *max=*maximdista;
  *min=*minimdista;
  return;
}*/


double aux_euv_binomneg (int N, double p1,double p2,double p11)
{

 int i=0;
 double a=0.0,b=0.0,kk1=0.0,kk2=0.0,kk3=0.0,kk4=0.0,kk5=0.0,kk6=0.0,kk7=0.0,euv1=0.0,euv2=0.0,euv3=0.0;

 double p00=1+p11-(p1+p2);
 double p10=p1-p11;double p01=p2-p11;
 double P=1-p00;

for(i=1;i<=fmin(N-1,2*N-3);i++){
kk1=exp(lgammafn(2*N-i-3+1)-(lgammafn(N-i)+lgammafn(N-i)+lgammafn(i)));
a=(R_pow(N,2)*p00+2*R_pow(N,2)-2*i*N-i-1)*p00 + (N-2*i-2)*N + R_pow(i,2) + 2*i+1;
b=R_pow(p01/P,N-i)*pow(p10/P,N-i)*R_pow(p11/P,i);

kk2=a*p1-(N*p00-i-1)*p00 + N-i-1;
kk3=P*p1*p10;
euv1+=kk1*(kk2/kk3)*b;
                     

kk4=a*p2-(N*p00-i-1)*p00 + N-i-1;
kk5=P*p2*p01;
euv2+=kk1*(kk4/kk5)*b;
                      
kk6=p11*a;
kk7=P*p01*p10;
euv3+=kk1*(kk6/kk7)*b;
} 


 return(euv1+euv2+euv3);
}


double corr_binomneg (int N, double p1,double p2,double p11)
{
double corr=0.0;
 corr=(p1*p2*aux_euv_binomneg(N,p1,p2,p11)-R_pow(N,2)*(1-p1)*(1-p2))/(N*R_pow((1-p1),1/2)*R_pow((1-p2),1/2));
 return(corr);
}
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
// Computes the spatial distances 
void Space_Dist(double *coordx,double *coordy,int *ia,int *idx,
		int *ismal,int *ja,int *colidx,int *rowidx ,double thres) 

{
  int i=0,h=0,j=0;
  double dij=0.0;


  /******************************************************************************/
  if(*istap){   // tapering case
      ia[0]=1;
	for(i=0;i<ncoord[0];i++)
	  for(j=0;j<ncoord[0];j++){
	    dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
	    if(dij<= thres){
	      tlags[h]=dij;
	      ja[h]=j+1;
	      idx[h]=i*(ncoord[0])+j+1;
	      ia[i+1]=ia[i+1]+1;
	      h++;}}
	for(i=0;i<ncoord[0];i++)
	  ia[i+1]=ia[i+1]+ia[i];
}  //end tapering
  /******************************************************************************/
else{  //no tapering
     h=0;


	  for(i=0;i<(ncoord[0]-1);i++){
	    for(j=(i+1);j<ncoord[0];j++){

      dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
  
     if(dij<= thres){

        tlags[h]=dij;
        colidx[h]=i;  rowidx[h]=j; 
        h++;}
	      }}

   }
 // saving  spatial distances   
    *npairs=h;
    lags= (double *) Calloc(*npairs,double);
    for(i=0;i<*npairs;i++)  lags[i]=tlags[i];  
     Free(tlags);
	return;
}



// Computes the spatial-temporal distances on regular and irregular grid:
void SpaceTime_Dist(double *coordx,double *coordy,double *coordt,int *ia,int *idx,int *ismal,int *ja,
                    int *tapmodel,int *ns, int  *NS,int *colidx,int *rowidx ,double *thres,double *thret)
{
  int i=0,cc=0,j=0,h,k=0,t=0,v=0;
  double dij=0.0,dtv=0.0;


  if (*istap) {// start tapering case

  double *thre,*c_supp;
  c_supp=(double *) R_alloc(2, sizeof(double));   // vector of compact support in space time tapering
  thre=(double *)   R_alloc(2, sizeof(double));
  thre[0]=thres[1];thre[1]=thret[1];
  
  h=0;

 if(isst[0]){  // space time case
        ia[0] = 1;
        for(t=0;t<*ntime;t++){
        for(i=0;i<ncoord[0];i++){
          cc=0;
           for(v=0;v<*ntime;v++){
               dtv=fabs(coordt[t]-coordt[v]);
          for(j=0;j<ncoord[0];j++){
          dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                Comp_supp(c_supp,tapmodel, dij, dtv,thre);
                  if((dij<c_supp[0]||is_equal(dij,c_supp[0]))&&(dtv<c_supp[1] ||  is_equal(dtv,c_supp[1]))){     
                               tlags[h]=dij;
                               tlagt[h]=dtv;
                               idx[h] =(t * (ncoord[0]) * (ncoord[0]) * ntime[0]) +  (i*  ntime[0] *  *ncoord) +  (1+j+ *ntime * v);
                               ja[h]=1+j+(*ncoord) * v;
                               cc=cc+1;
                               h = h +1 ;
                }}}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}
   }   // end space time case
}    // end tapering case

  else {   // no tapering
  h=0;
        
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
           dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                      if(dij<=thres[1]){
                          tlags[h]=dij;tlagt[h]=0;
                 
                           colidx[h]=i+NS[t];  rowidx[h]=j+NS[v];
                          h++;
                                    }}}
               else {
          dtv=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
           dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                          if(dij<=thres[1] && dtv<=thret[1]){
                            tlags[h]=dij;tlagt[h]=dtv;
                           colidx[h]=i+NS[t];  rowidx[h]=j+NS[v];  
                           h++; 
                             }}}}}}    
}  // end no tapering case

    *npairs=h;
    lags= (double *) Calloc(*npairs,double);
    lagt= (double *) Calloc(*npairs,double);
    for(i=0;i<*npairs;i++)  
        {lags[i]=tlags[i];lagt[i]=tlagt[i];  }
     Free(tlags); Free(tlagt);   

  return;
}


/*********************************************************************************************/
/*********************************************************************************************/

// Computes the spatial  distances bivariate case
void SpaceBiv_Dist(double *coordx,double *coordy,double *coordt,int *ia,int *idx,int *ismal,int *ja,
                    int *tapmodel,int *ns, int  *NS,int *colidx,int *rowidx ,double *thres)
{
  int i=0,cc=0,j=0,k=0,t=0,v=0;
  double dij=0.0;
  
      int count=0;

  if (*istap) {// start tapering case


       ia[0] = 1;
   for(t=0;t<*ntime;t++){
   for(i=0;i<ncoord[0];i++){  //for(i=0;i<ns[t];i++){
          cc=0;
           for(v=0;v<*ntime;v++){  
             for(j=0;j<ncoord[0];j++){   //  for(j=0;j<ns[v];j++){
                    dij=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                    if(dij<=dista[t][v]){
                  if(j<i)  {tfirst[count]=v;tsecond[count]=t;}
                  else     {tfirst[count]=t;tsecond[count]=v;}
                               tlags[count]=dij;
                               idx[count] =(t * ncoord[0] * ncoord[0] * ntime[0]) +  (i*  ntime[0] *  ncoord[0]) +  (1+j+ ncoord[0] * v);
                               ja[count]=1+j+ncoord[0] * v;
                               cc=cc+1;
                               count = count +1 ;
                }}}
                ia[k+1]=ia[k]+cc;
                k=k+1;}}

}    // end tapering case

  else {   // no tapering
        count=0;

  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){    
      if(t==v){
         for(j=i+1;j<ns[v];j++){

          dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);        
                        if(dij<=dista[t][v]){

                           tlags[count]=dij;
                           colidx[count]=i+NS[t];  rowidx[count]=j+NS[v];
                           tfirst[count]=t;tsecond[count]=v;

                           count++;
                                }
        }}
             else {  
         for(j=0;j<ns[v];j++){
         dij=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(dij<=dista[t][v]){

                            tlags[count]=dij;
                           colidx[count]=i+NS[t];  rowidx[count]=j+NS[v];
                           tfirst[count]=t;tsecond[count]=v;
                           count++;
                                 }}}}}}
}  // end no tapering case


  *npairs=count;
  first =(int *)  Calloc(count,int);
  second=(int *)  Calloc(count,int);
  lags=(double *) Calloc(count,double);
  for(i=0;i<count;i++) { lags[i]=tlags[i];
                         first[i]=tfirst[i];
                         second[i]=tsecond[i];}
   Free(tlags);Free(tfirst);Free(tsecond);
  return;
}

/*********************************************************************************************/
/*********************************************************************************************/
/*********************************************************************************************/

// comp support function:
void Comp_supp(double *c_supp,int *cormod, double h,double u, double *par)
{
    switch(*cormod) // computing comp supports
    {
        
        case 30:// Wendland1
        case 32:// Wendland2
        case 34:// Wendland3
        c_supp[0]=par[0];
        c_supp[1]=-LOW;
        break;
        case 200: // separable spacetime wendland
        case 202:
        case 204:
        case 206:
        case 208:
        case 210:
        case 212:
        case 214:
        case 216:
        c_supp[0]=par[0];
        c_supp[1]=par[1];
        break;
        case  218:   // quasi taper in time
        case  220:
        case  222:
        case 64:
        case 66:
        case 68:
        c_supp[0]=-LOW;
        c_supp[1]=par[1]*pow(1+pow(h/par[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));// arg=pow(1+pow(h/scale_s,power_s/2),-sep/(power_s/2));
        break;
        case  224: // quasi taper in space
        case  226:
        case  228:
        case 63:
        case 65:
        case 67:

        c_supp[0]=tapsep[2]*pow(1+pow(u/tapsep[1],tapsep[1]/2),-tapsep[4]/(tapsep[1]/2));
        c_supp[1]=-LOW;
        break;
        //case  226: // quasi taper in space x quasi taper in time;
        //c_supp[0]=par[0]*pow(1+pow(u/par[1],1),-tapsep[0]);
        //c_supp[1]=par[1]*pow(1+pow(h/par[0],1),-tapsep[1]);
        break;
        case 230:
        case 69:
        case 70:
        case 71:
        case 72:
        case 73:
        case 74:
        case 75:
        case 76:
        case 77:
            c_supp[0]=par[0];
            c_supp[1]=par[1];
            break;
            
    }
}




// check if 'val1' is equal to 'val2'  when both are double
int is_equal(double val1, double val2)
{
  return fabs(val1-val2)<MAXERR;
}


double Maxima_i(int *x, int size)
{
  int res=0;
  int i=0;

  res = x[0];

  for(i = 1; i < size; i++)
    res = imax2(res, x[i]);

  return res;
}


double Maxima(double *x, int *size)
{
  double res=0.0;
  int i=0;

  res = x[0];

  for(i = 1; i < *size; i++)
    res = fmax2(res, x[i]);

  return res;
}

double Minima(double *x, int *size)
{
  double res=0.0;
  int i=0;

  res = x[0];

  for(i = 1; i < *size; i++)
    res = fmin2(res, x[i]);

  return res;
}


//return x*(x-1)*...j  j<x
int fmax_int(int u,int v)
{
  if(u>=v) return u;
  else return(v);
}
/**********************************************/
int fmin_int(int u,int v)
{
  if(u<=v) return u;
  else return(v);
}

// sterling approx
double stirling(double x)
{
    return (x*log(x)-x);
}


/*******  stirling appro***************************************/
double logfac(int n)
{
    int i;
    double factorial = 1.0;
    if(n<=170)
    {
    for(i=1; i<=n; i++) factorial *= i; 
    return(log(factorial));
    }
    else 
      {
  
        return(stirling(n));}
}



/**********************************************/
double fac(int n,int j)
{
    int i;
    double factorial = 1.0;
    for(i=j; i<=n; ++i) factorial *= i;              
    return(factorial);
}
/**********************************************/
void Range(double *x, double *ran, int *size)
{
  int i=0;

  ran[0] = x[0];
  ran[1] = x[0];

  for(i = 1; i < *size; i++)
    {
      ran[0] = fmin(ran[0], x[i]);
      ran[1] = fmax(ran[1], x[i]);
    }

  return;
}
// define a sequence of points from x[0] to x[1] of length 'len'
void Seq(double *x, int len, double *res)
{
  double delta=0.0;
  int i=0;

  res[0] = x[0];
  delta = (x[1] - x[0]) / (len - 1);

  for(i = 1; i < len; i++)
    res[i] = res[i - 1] + delta;

  return;
}
// define a sequence of 'len' points of of 'delta' steps from the starting point x[0]
void SeqStep(double *x, int len, double step, double *res)
{
  int i=0;
  res[0]=x[0];
  for(i=0;i<len;i++) res[i+1]=res[i]+step;
  return;
}


// Determine (for the sub-sampling procedure) the sub-coordinates and
// the sub-data given spatial data, coordinates and an interval:
void SetSampling(double *coordx, double *coordy, double *data, int n, int *npts, int nbetas,
		 double *scoordx, double *scoordy, double *sdata, double xmax,
		 double xmin, double ymax, double ymin, double **sX,double **X)
{
  int i=0, j=0,k=0;
  for(i=0;i<ncoord[0];i++)
    if((xmin<coordx[i]||is_equal(xmin,coordx[i]))&&
       (xmax>coordx[i]||is_equal(xmax,coordx[i]))&&
       (ymin<coordy[i]||is_equal(ymin,coordy[i]))&&
       (ymax>coordy[i]||is_equal(ymax,coordy[i]))){
	scoordx[j]=coordx[i];
	scoordy[j]=coordy[i];
  
  for(k=0;k<nbetas;k++)  sX[j][k]=X[i][k];
	
  sdata[j]=data[i];
	
  j++;}
  *npts = j;
  return;
}



// subsampling in space
void SetSampling_s(double *coordx, double *coordy, double *data, int *npts, int nbetas,
                   double *scoordx, double *scoordy, double *sdata, double xmax,
                   double xmin, double ymax, double ymin,double **sX,double **X,int *ns,int *NS, int *NS_sub, double *res_sub, double *coordt,int *ns_sub)
{
    int i=0, j=0, f=0,k=0;
    int NTOT=(NS[ntime[0]-1]+ns[ntime[0]-1]);
        for(i=0;i<NTOT;i++)
        {
            if((xmin<coordx[i]||is_equal(xmin,coordx[i]))&&
               (xmax>coordx[i]||is_equal(xmax,coordx[i]))&&
               (ymin<coordy[i]||is_equal(ymin,coordy[i]))&&
               (ymax>coordy[i]||is_equal(ymax,coordy[i])))
            {
                scoordx[j]=coordx[i];
                scoordy[j]=coordy[i];
                sdata[j] = data[i];
                
                for(f=0;f<ntime[0];f++)
                {
                    if(res_sub[i]==coordt[f]){ns_sub[f]++;}
                }
          for(k=0;k<nbetas;k++) {sX[j][k]=X[i][k];}
                j++;   
            }
        }
        cumvec(ns_sub,NS_sub,ntime[0]);
        *npts = j;
    return;
}


// Determine (for the sub-sampling procedure) the sub-coordinates and
// the sub-data given spatial data, coordinates and an interval:
void SetSampling_st(double *data,double *sdata,int *ncoordss,int *ntime, int nbetas,
        int wint,int k, double **sX,double **X)
{
  int i=0,j=0,p=0,f=0;
  for(i=0;i<(*ncoordss);i++){
    for(j=(k+(ntime[0]*i));j<(k+wint+(ntime[0]*i));j++) {
        sdata[p]=data[j];
       for(f=0;f<nbetas;f++) sX[p][f]=X[j][f];

        p++;


    }}
  return;
}  

void SetSampling_t(double *data,double *sdata, int nbetas,int npts,
                   int nt,int wint,int k,double **sX,double **X,int *ns_sub,int *NS_sub,int nsub_t, int *ntimeS, double *s2cx, double *s2cy, double *scoordx, double *scoordy)
{
    int i=0,j=0,p=0,f=0;
    for(j=(k*wint);j<(k*wint+wint);j++)
    {
        for(i=NS_sub[j];i<NS_sub[j]+ns_sub[j];i++)
        {
            sdata[p]=data[i];
            s2cx[p] = scoordx[i];
            s2cy[p] = scoordy[i];
            for(f=0;f<nbetas;f++) {sX[p][f]=X[i][f];}
            p++;
        }
    }
    *ntimeS = p;
    return;
}

// Determine (for the sub-sampling procedure) the sub-coordinates and
// the sub-data given spatial  bivariate data, coordinates and an interval:
void SetSampling_biv(double *coordx, double *coordy, double *data, int n, int *npts,
                 double *scoordx, double *scoordy, double *sdata, double xmax,
                 double xmin, double ymax, double ymin)
{
    int i=0, j=0,f=0;
    for(i=0;i<ncoord[0];i++){
        if((xmin<coordx[i]||is_equal(xmin,coordx[i]))&&
           (xmax>coordx[i]||is_equal(xmax,coordx[i]))&&
           (ymin<coordy[i]||is_equal(ymin,coordy[i]))&&
           (ymax>coordy[i]||is_equal(ymax,coordy[i]))){
            scoordx[j]=coordx[i];scoordy[j]=coordy[i];
            sdata[f]=data[n     + *ntime*i];
            f++;
            sdata[f]=data[n+ 1 + *ntime*i];
            f++;
            j++; }}
    *npts = j;
    return;
}



// Set the global variables for the spatial and spatial-temporal fitting:
void SetGlobalVar(int *biv,double *coordx,double *coordy,double *coordt,int *grid,int *ia,
		  int *idx,int *ismal,int *ja,int *mem, int *nsite,int *nsitex,int *nsitey,
		  int *npair,double *radius,double *srange, double *sep,int *st, int *times,double *trange,
		  int *tap,int *tapmodel,int *tp,int *weighted, int *colidx,int *rowidx, 
      int *ns, int *NS, int *dyn)
{
  //Spatial settings: //Spatial settings:
  maxdist=(double *) Calloc(3,double);//spatial threshould
  if(maxdist==NULL) {*ismal=0; return;}

  maxtime=(double *) Calloc(1,double);//temporal threshold
    if(maxtime==NULL) {*ismal=0; return;}

   ntime=(int *) Calloc(1,int);//number of times
    if(ntime==NULL) {*ismal=0; return;}
    *ntime=*times;

/**********************************/
  ncoord=(int *) Calloc(1,int);//number of total spatial coordinates

   // ncoord=(int *) R_alloc(1, sizeof(int));

  if(ncoord==NULL) {*ismal=0; return;}
  ncoord[0]=*nsite;

  ncoordx=(int *) Calloc(1,int);//number of the first spatial coordinates
  if(ncoordx==NULL) {*ismal=0; return;}
  *ncoordx=*nsitex;
  ncoordy=(int *) Calloc(1,int);//number of the second spatial coordinates
  if(ncoordy==NULL) {*ismal=0; return;}
  *ncoordy=*nsitey;
/************************/
  npairs=(int *) Calloc(1,int);//effective number of pairs
  //npairs=(int *) R_alloc(1, sizeof(int));
  if(npairs==NULL) {*ismal=0; return;}

  isbiv=(int *) Calloc(1,int);//is a bivariate random field?
  if(isbiv==NULL) {*isbiv=0; return;}
  isbiv[0]=biv[0];  

  ismem=(int *) Calloc(1,int);//is distances computed using memory allocation
  if(ismem==NULL) {*ismal=0; return;}
  *ismem=mem[0];

  isst=(int *) Calloc(1,int);//is a spatio-temporal random field?
  if(isst==NULL) {*ismal=0; return;}
  isst[0]=st[0]; 
    
  cdyn=(int *) Calloc(1,int);//dynamic coords
  if(dyn==NULL) {*cdyn=0; return;}
  cdyn[0]=dyn[0]; 

  istap=(int *) Calloc(1,int);//is tapering?
  if(istap==NULL) {*ismal=0; return;}
  istap[0]=tap[0];

  type=(int *) Calloc(1,int);//type of distance
  if(type==NULL) {*ismal=0; return;}

  REARTH=(double *) Calloc(1,double);//radius of hearth
    if(REARTH==NULL) {*ismal=0; return;}

  tapsep=(double *) Calloc(5,double);
           if(tapsep==NULL){*ismal=0; return;}

  *REARTH=*radius;
  *type=*tp;
  
 // case of not saving distances
  if(!ismem[0]) {
                   if(srange[1]) maxdist[0]=srange[1];
                   else maxdist[0]=-LOW;

                   if(isst[0]){ if(trange[1]) maxtime[0]=trange[1];
                                else maxtime[0]=-LOW;
                       }
                   if(isbiv[0]){  
                                   int i=0;
                                   dista= (double **) Calloc(ntime[0],double *);
                                   if(dista==NULL) {*ismal=0; return;}
                                   for(i=0;i<ntime[0];i++){
                                   dista[i]=(double *) Calloc(ntime[0],double);
                                   if(dista[i]==NULL) {*ismal=0; return;}}
                                  
                                  if(srange[1]) maxdist[0]=srange[1];
                                  else maxdist[0]=-LOW;
                                  if(srange[2]) maxdist[1]=srange[2];
                                  else maxdist[1]=-LOW;
                                  if(srange[3]) maxdist[2]=srange[3];
                                  else maxdist[2]=-LOW;
                                  dista[0][0]=maxdist[0];dista[0][1]=maxdist[1];dista[1][0]=dista[0][1];dista[1][1]=maxdist[2];  
                    }
                  return;    // return here for the not  saving distance  
                } // end not saving distances
  else 
  {     // start  saving distances
  /***********************************************************/  
if(!isst[0]&&!isbiv[0]) {// spatial case
           // settting compact support
        if(srange[1]) maxdist[0]=srange[1];
        else maxdist[0]=-LOW;    
      if(istap[0])  // tapering case
           {
              npairs[0]=(int)(ncoord[0]*ncoord[0]);
              tlags=(double *) Calloc(*npairs,double);
              if(tlags==NULL){*ismal=0; return;}
           }  // end tapering case
      else { // distances composite likelihood  


           *npairs=*ncoord *(int) ((*ncoord-1)/2);

   //Rprintf("sdds22 %ld  %d %d %d %d\n",*npairs,npairs[0],*ncoord * (int)((*ncoord-1)/2),ncoord[0], *ncoord);

           tlags= (double *) Calloc(*npairs,double *);
           //  tlags=(double *) R_alloc(*npairs, sizeof(double));
           if(tlags==NULL) {*ismal=0; return;}
               //   Rprintf("hhk%d",ncoord[0]);
           } // end  no tapering case

 // computing spatial distances and indexes      
 Space_Dist(coordx,coordy,ia,idx,ismal,ja,colidx,rowidx,srange[1]);
      if(!ismal[0]) return;
  /***********************************************************/  
}  // end spatial case
else { //spatio temporal case or bivariate case
       int qq=ncoord[0]*ntime[0];
       
    // setting compact supports for space-time and bivariate case
       if(isst[0]){ 
            if(srange[1]) maxdist[0]=srange[1];
                 else maxdist[0]=-LOW;
            if(trange[1]) maxtime[0]=trange[1];
                else maxtime[0]=-LOW; 
        }
       if(isbiv[0])              {
                                   int i=0;
                                   dista= (double **) Calloc(ntime[0],double *);
                                   if(dista==NULL) {*ismal=0; return;}
                                   for(i=0;i<ntime[0];i++){
                                   dista[i]=(double *) Calloc(ntime[0],double);
                                   if(dista[i]==NULL) {*ismal=0; return;}}
                                  
                                  if(srange[1]) maxdist[0]=srange[1];
                                  else maxdist[0]=-LOW;
                                  if(srange[2]) maxdist[1]=srange[2];
                                  else maxdist[1]=-LOW;
                                  if(srange[3]) maxdist[2]=srange[3];
                                  else maxdist[2]=-LOW;
                                  dista[0][0]=maxdist[0];dista[0][1]=maxdist[1];dista[1][0]=dista[0][1];dista[1][1]=maxdist[2]; 
                                  }
     if(istap[0])  // tapering case
        {

          // allocating vectors
           npairs[0]=(int)(qq*qq);
           tlags=(double *) Calloc(*npairs,double);
           if(tlags==NULL){*ismal=0; return;}
         
         if(isst[0])    { //spatio temporal case
                           tlagt=(double *) Calloc(*npairs,double);
                           if(tlagt==NULL){*ismal=0; return;}
                           }
         if(isbiv[0])         { //bivariate case
                           tfirst=(int *) Calloc(*npairs,int);
                           if(tfirst==NULL){*ismal=0; return;}
                           tsecond=(int *) Calloc(*npairs,int);
                           if(tsecond==NULL){*ismal=0; return;}
                              }
           tapsep[0]=sep[0];tapsep[1]=sep[1];tapsep[2]=sep[2];tapsep[3]=sep[3];tapsep[4]=sep[4];
           if(tapsep[0]==1) tapsep[0]=0.99999999;
       }  // end tapering
else {  // distance for composite likelihood
              
               if(isst[0])  npairs[0]=qq*(int)((qq-1)*0.5);
               if(isbiv[0]) npairs[0]=qq*(int)((qq-1)*0.5);

             tlags= (double *) Calloc(*npairs,double *);
            if(tlags==NULL) {*ismal=0; return;}
          // allocates the matrix of temporal distances:
          if(isst[0]) {
             //memory allocation of matrix temporal distances
             tlagt= (double *) Calloc(*npairs,double *);
             if(tlagt==NULL) {*ismal=0; return;}
            
                     }
          if(isbiv[0]) {
            tfirst=(int *) Calloc(*npairs,int);
            if(tfirst==NULL){*ismal=0; return;}
            tsecond=(int *) Calloc(*npairs,int);
            if(tsecond==NULL){*ismal=0; return;}
                      }
       }
if(isst[0])  SpaceTime_Dist(coordx,coordy,coordt,ia,idx,ismal,ja,tapmodel,
                                        ns,NS,colidx,rowidx,srange,trange);
if(isbiv[0]) SpaceBiv_Dist(coordx,coordy,coordt,ia,idx,ismal,ja,tapmodel,
                                     ns,NS,colidx,rowidx,srange);
  if(!ismal[0]) return;
 } //end spatio temporal case or bivariate case


  npair[0]=npairs[0];
  return;
     }
}

void DeleteGlobalVar()
{

  int i=0;
  // Delete all the global variables:
  Free(maxdist);Free(maxtime);
  Free(ncoordx);Free(ncoordy); 
  Free(ncoord);
  Free(npairs);
  Free(type);Free(REARTH);
  Free(tapsep);
  if(isbiv[0])for(i=0;i<ntime[0];i++)  Free(dista[i]);
  Free(dista);
  Free(ntime);
  if(ismem[0]) {
           Free(lags);
           if(isst[0])    {Free(lagt);}
           if(isbiv[0])   {Free(first);Free(second);}
  }

  Free(isbiv); Free(istap);
  Free(isst);Free(ismem);
  Free(cdyn);
  return;
}






void Rep(double *coordt,int *ns, double *res)
{
    int i, j,ppb=0;
    for(i =0;i<ntime[0];i++)
    {for(j =0;j<ns[i];j++){
            res[ppb] = coordt[i];
            ppb++;}}
}


void cumvec(int *ns,int *res,int len)
{
    int i=0;
    res[0] =0;
    int sum=0;
    for(i =1;i<len;i++)
    {
        sum += ns[i-1];
        res[i] = sum;
    }
}




// ============= TEST QQNORM



#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
if (log_p) {					\
if(p > 0)					\
return NAN;				\
if(p == 0) /* upper bound*/			\
return lower_tail ? _RIGHT_ : _LEFT_;	\
if(p == -INFINITY)				\
return lower_tail ? _LEFT_ : _RIGHT_;	\
}							\
else { /* !log_p */					\
if(p < 0 || p > 1)				\
return NAN;				\
if(p == 0)					\
return lower_tail ? _LEFT_ : _RIGHT_;	\
if(p == 1)					\
return lower_tail ? _RIGHT_ : _LEFT_;	\
}

#define R_D_Cval(p)	(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */

#define R_DT_CIv(p)	(log_p ? (lower_tail ? -expm1(p) : exp(p)) \
: R_D_Cval(p))

#define R_D_Lval(p)	(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */


#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - expm1(p)) \
: R_D_Lval(p))


double qnorm55(double p, double mu, double sigma, int lower_tail, int log_p)
{
    double p_, q, r, val;
    
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
        return p + mu + sigma;
#endif
    R_Q_P01_boundaries(p, -INFINITY, INFINITY);
    
    if(sigma  < 0)	return NAN;
    if(sigma == 0)	return mu;
    
    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;
    
#ifdef DEBUG_qnorm
    REprintf("qnorm(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
             p,mu,sigma, lower_tail, log_p, q);
#endif
    
    
    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
     
     Produces the normal deviate Z corresponding to a given lower
     tail area of P; Z is accurate to about 1 part in 10**16.
     
     (original fortran code used PARAMETER(..) for the coefficients
     and provided hash codes for checking them...)
     */
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
        q * (((((((r * 2509.0809287301226727 +
                   33430.575583588128105) * r + 67265.770927008700853) * r +
                 45921.953931549871457) * r + 13731.693765509461125) * r +
               1971.5909503065514427) * r + 133.14166789178437745) * r +
             3.387132872796366608)
        / (((((((r * 5226.495278852854561 +
                 28729.085735721942674) * r + 39307.89580009271061) * r +
               21213.794301586595867) * r + 5394.1960214247511077) * r +
             687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */
        
        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = R_DT_CIv(p);/* 1-p */
        else
            r = p_;/* = R_DT_Iv(p) ^=  p */
        
        r = sqrt(- ((log_p &&
                     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
                    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
#ifdef DEBUG_qnorm
        REprintf("\t close to 0 or 1: r = %7g\n", r);
#endif
        
        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                         .0227238449892691845833) * r + .24178072517745061177) *
                       r + 1.27045825245236838258) * r +
                      3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                   1.42343711074968357734)
            / (((((((r *
                     1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                   .14810397642748007459) * r + .68976733498510000455) *
                 r + 1.6763848301838038494) * r +
                2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                         2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                      r + .29656057182850489123) * r +
                     1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772)
            / (((((((r *
                     2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                    r + 1.8463183175100546818e-5) * r +
                   7.868691311456132591e-4) * r + .0148753612908506148525)
                 * r + .13692988092273580531) * r +
                .59983220655588793769) * r + 1.);
        }
        
        if(q < 0.0)
            val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}






void qnorm55_call(double *p, double *mu, double *sigma, int *lower_tail, int *log_p, double *res)
{
    *res=   qnorm55( *p,  *mu,  *sigma,  *lower_tail,  *log_p);
    
}
