#include "header.h"

// get all pair between certain distance 
void pairs(int *ncoords,double *data,double *coordx, double *coordy, double *numbins, double *bins, double *v0,double *v1, double *v2,double *maxdist)
{
  int ncrd,numbin,h=0,k=0,i,j;
  double max_dist;
  
  ncrd     = *ncoords; //printf("num coords =  %d \n",ncrd);
  numbin   = *numbins; //printf("num bins =  %d \n",numbin);
  max_dist = *maxdist; //printf("max distance =  %f \n",max_dist);

  double distance=0.0;
  for(h=0;h<=numbin;h++){
      for(i=0; i<(ncrd-1);i++){
        for(j=(i+1);j<ncrd;j++){
          distance = dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
          if(distance <= max_dist){
            if((bins[h] < distance) && (distance <= bins[(h+1)])){
              v0[k] = bins[(h)];
              v1[k] = data[i];
              v2[k] = data[j];
              k = k+1;
         } }  
      }
  }}

return;
}

// binned spatial variogram:
void Binned_Variogram2(double *bins, double *coordx, double *coordy, double *coordt,double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0;
  double x,y,lags=0.0,step=0.0,*mm;
  //Set the binnes step:
  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, ncoord,type,REARTH);
  if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  step=(mm[1]-mm[0])/(*nbins-1);
  bins[0]= mm[0];
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned moments:
  for(i=0;i<(ncoord[0]-1);i++)
    for(j=(i+1);j<ncoord[0];j++){
                         lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
      if(lags<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=lags) && (lags<bins[h+1])){
            x=data[n+i ];
            y=data[n+j ];
            if(!(ISNAN(x)||ISNAN(y)))
              moms[h]+=0.5*pow(x-y,2);lbins[h]+=1;
            }}}
  return;
}
/***********************************************************************************************************************************/
// binned spatial variogram:
void Binned_Variogram2new(double *bins, int *np,double *data1, double *data2, 
    double *vdist, int *lbins, double *moms, int *nbins,double *mm)
{
  int h=0,  k=0;
  double x,y,step=0.0;
  //Set the binnes step:
  //if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  step=(mm[1]-mm[0])/(*nbins-1);
  bins[0]= mm[0];
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned moments:
  for(k=0;k<*np;k++){
  
      if(vdist[k]<=*maxdist){
  for(h=0;h<(*nbins-1);h++)
    if((bins[h]<=vdist[k]) && (vdist[k]<bins[h+1])){
            x=data1[k];   y=data2[k];
            if(!(ISNAN(x)||ISNAN(y)))
              moms[h]+=0.5*pow(x-y,2);lbins[h]+=1;
            }}
        k=k+1;
          }
  return;
}
/***********************************************************************************************************************************/


void Binned_Variogram_22(double *bins, double *coordx, double *coordy, double *coordt,double *data, int *lbins, double *moms, int *nbins)
{
  int h=0, i=0, j=0, n=0,p=0;
  double x,y,step=0.0,*mm;
  //Set the binnes step:
  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, ncoord,type,REARTH);
  if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  step=(mm[1]-mm[0])/(*nbins-1);
  bins[0]= mm[0];
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned moments:
  for(i=0;i<(ncoord[0]-1);i++){
    for(j=(i+1);j<ncoord[0];j++){
      if(lags[p]<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=lags[p]) && (lags[p]<bins[h+1])){
            x=data[n+i ]; y=data[n+j ];
            if(!(ISNAN(x)||ISNAN(y))){
	      moms[h]+=0.5*pow(x-y,2);
	      lbins[h]+=1;}}
	      p++;}}}
  return;
}



// binned spatial-temporal variogram:

void Binned_Variogram_st2(double *bins, double *bint, double *coordx, double *coordy, double *coordt,double *data, int *lbins, int *lbinst,
       int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint, int *ns,int *NS)
{
int h=0, i=0, j=0;
  int q=0, t=0, u=0, v=0;
  double x,y,lags=0.0,lagt=0.0,step=0.0,*mm,*tt;
  //defines the spatial bins:
  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, ncoord,type,REARTH); // computing max and min distances
  if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  //Set the binnes step:
  step=(mm[1])/(*nbins-1);
  bins[0]= mm[0];
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
   tt=(double *) R_alloc(2, sizeof(double));
 Maxima_Minima_time(tt,coordt,ntime);
  bint[0]=0;
  for(u=1;u<*nbint;u++)
    bint[u]=bint[u-1]+tt[0];


   for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
  if(t==v){// computes the marginal spatial variogram:
             for(j=i+1;j<ns[v];j++){
                                lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                         if(lags<=*maxdist){
                            for(h=0;h<(*nbins-1);h++){
                             if((bins[h]<=lags) && (lags<bins[h+1])){
                             x=data[(i+NS[t])];y=data[(j+NS[v])];
                             if(!(ISNAN(x)||ISNAN(y))){
                             moms[h]+=0.5*pow(x-y, 2);
                             lbins[h]+=1;}
                           }}}}
          } 
     else {
         lagt=fabs(coordt[t]-coordt[v]);
          for(j=0;j<ns[v];j++){
                if(i==j){// computes the marginal temp variogram:
                    if(lagt<=*maxtime)
                    {
                    for(u=0;u<*nbint;u++){
                       if(is_equal (bint[u],lagt)){

                     x=data[(i+NS[t])];y=data[(i+NS[v])];
                    if(!(ISNAN(x)||ISNAN(y))){
                    momt[u]+=0.5*pow(x-y, 2);
                    lbint[u]+=1;}
                  }}}}
          else{// computes the spatial-temporal variogram:
                 lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                  if(lags<=*maxdist && lagt<=*maxtime){
                    q=0;
                     for(h=0;h<(*nbins-1);h++){
                      for(u=0;u<*nbint;u++){
        if((bins[h]<=lags) && (lags<bins[h+1]) && (is_equal(bint[u],lagt))){
                 x=data[(i+NS[t])];y=data[(j+NS[v])];
               if(!(ISNAN(x)||ISNAN(y)))   {
                                            momst[q]+=0.5*pow(x-y, 2);
                                            lbinst[q]+=1;
                                          }
                                        }
        q++;}}

      }

      }}}}}}
 }    


// binned spatial-temporal variogram:
void Binned_Variogram_st2_dyn(double *bins, double *bint, double *coordx, double *coordy, double *coordt,double *data, int *lbins, int *lbinst,
       int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint, int *ns,int *NS)
{
int h=0, i=0, j=0;
  int q=0, t=0, u=0, v=0;
  double x,y,lags=0.0,lagt=0.0,step=0.0,*mm,*tt;
  //defines the spatial bins:
  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, ncoord,type,REARTH); // computing max and min distances
  if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  //Set the binnes step:
  step=(mm[1])/(*nbins-1);
  bins[0]= mm[0];
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
   tt=(double *) R_alloc(2, sizeof(double));
 Maxima_Minima_time(tt,coordt,ntime);
  bint[0]=0;
  for(u=1;u<*nbint;u++)
    bint[u]=bint[u-1]+tt[0];


   for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
  if(t==v){// computes the marginal spatial variogram:
             for(j=i+1;j<ns[v];j++){
                                lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                         if(lags<=*maxdist){
                            for(h=0;h<(*nbins-1);h++){
                             if((bins[h]<=lags) && (lags<bins[h+1])){
                             x=data[(i+NS[t])];y=data[(j+NS[v])];
                             if(!(ISNAN(x)||ISNAN(y))){
                             moms[h]+=0.5*pow(x-y, 2);
                             lbins[h]+=1;}
                           }}}}
          } 
     else {
         lagt=fabs(coordt[t]-coordt[v]);
          for(j=0;j<ns[v];j++){

            lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                

    // a "marginal" temporal semivariogram
             if((bins[0]/2<=lags) && (lags<bins[1]/2)  && lagt<=*maxtime){
                    {
                    for(u=0;u<*nbint;u++){
                       if(is_equal (bint[u],lagt)){
                     x=data[(i+NS[t])];y=data[(j+NS[v])];
                    if(!(ISNAN(x)||ISNAN(y))){ momt[u]+=0.5*pow(x-y, 2);lbint[u]+=1;}
                  }}}}
          

          if(lags<=*maxdist && lagt<=*maxtime){

                     q=0;
                     for(h=0;h<(*nbins-1);h++){
                      for(u=0;u<*nbint;u++){
        if((bins[h]<=lags) && (lags<bins[h+1]) && (is_equal(bint[u],lagt))){
                 x=data[(i+NS[t])];y=data[(j+NS[v])];
               if(!(ISNAN(x)||ISNAN(y)))   {
                                           
                               //  if(h>0)         {  
                                  momst[q]+=0.5*pow(x-y, 2);lbinst[q]+=1; //}
                                // if(h==0)        {     momt[u]+=0.5*pow(x-y, 2);lbint[u]+=1;}
                                          }
                                        }q++;
                                      }}}}}}}}
 }    




/*

// binned spatial-temporal variogram:
void Binned_Variogram_st2(double *bins, double *bint, double *coordx, double *coordy, double *coordt,double *data, int *lbins, int *lbinst,
       int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint, int *ns,int *NS)
{
int h=0, i=0, j=0;
  int q=0, t=0, u=0, v=0;
  double x,y,lags=0.0,lagt=0.0,step=0.0,*mm,*tt;
  //defines the spatial bins:
  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, ncoord,type,REARTH); // computing max and min distances
  if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  //Set the binnes step:
  step=(mm[1])/(*nbins-1);
  bins[0]= mm[0];
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //defines the temporal bins:
   tt=(double *) R_alloc(2, sizeof(double));
 Maxima_Minima_time(tt,coordt,ntime);
  bint[0]=0;
  for(u=1;u<*nbint;u++)
    bint[u]=bint[u-1]+tt[0];

    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){       //-1
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]){
                            for(h=0;h<(*nbins-1);h++){
                             if((bins[h]<=lags) && (lags<bins[h+1])){
                             x=data[(i+NS[t])];y=data[(j+NS[v])];
                             if(!(ISNAN(x)||ISNAN(y))){
                             moms[h]+=0.5*pow(x-y, 2);
                             lbins[h]+=1;}
                           }}}}
          } 
     else {
            
         lagt=fabs(coordt[t]-coordt[v]);
        for(j=0;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                //if(i==j){// computes the marginal temp variogram:
                 if(lags==0||lags<=bins[1]){
                    if(lagt<=*maxtime)
                    {
                    for(u=0;u<*nbint;u++){
                        if(is_equal (bint[u],lagt)){
                    x=data[(i+NS[t])];y=data[(j+NS[v])];
                    if(!(ISNAN(x)||ISNAN(y))){
                    momt[u]+=0.5*pow(x-y, 2);
                    lbint[u]+=1;}
                  }}}}
          else{// computes the spatial-temporal variogram:
                 // lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                   if(lags<=*maxdist && lagt<=*maxtime){
                    q=0;
                     for(h=0;h<(*nbins);h++){
                      for(u=0;u<*nbint;u++){
        if((bins[h]<=lags) && (lags<bins[h+1]) && (is_equal(bint[u],lagt))){
                 x=data[(i+NS[t])];y=data[(j+NS[v])];
               if(!(ISNAN(x)||ISNAN(y)))   {
                                            momst[q]+=0.5*pow(x-y, 2);
                                            lbinst[q]+=1;
                                          }}
        q++;}}}}}}}}}
 }    
*/

void Binned_Variogram_biv2(double *bins, double *coordx, double *coordy, double *coordt,double *data, int *cross_lbins, double *cross_moms, int *nbins,
                          int *marg_lbins, double *marg_moms,int *ns, int *NS)
{
int h=0, i=0, j=0;
  int t=0, v=0;
  double x,y,a,b,lags=0.0,step=0.0,*mm,md;
    //Set the binnes step:
  //Set the binnes step:

  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, ncoord,type,REARTH);
  md=fmax(dista[0][1],fmax(dista[1][1],dista[0][0])); // we consider the max of the ditances and we build bins on [0,mm]
  if(md<mm[1]) mm[1]=md;
  step=mm[1]/(*nbins-1);
  bins[0]=0;
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //computes the empirical variogram:
  // Computes the log-likelihood:
    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
           
              if(lags<=dista[t][v]) {
                     for(h=0;h<(*nbins-1);h++)     {
      if((bins[h]<=lags) && (lags<bins[h+1])){
               x=data[(i+NS[t])];
               y=data[(j+NS[v])];
              if(!(ISNAN(x)||ISNAN(y))){
                  marg_moms[h+t*(*nbins-1)]+=0.5*pow(x-y,2);
                  marg_lbins[h+t*(*nbins-1)]+=1;
                  
              }}}}}}
   else {
         for(j=0;j<ns[v];j++){
             lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
         
            if(lags<=dista[t][v]) {
             for(h=0;h<(*nbins-1);h++){
        if((bins[h]<=lags) && (lags<bins[h+1])){
                x=data[(i+NS[t])];y=data[(j+NS[t])];
                a=data[(i+NS[v])];b=data[(j+NS[v])];
                 if(!(ISNAN(x)||ISNAN(y)||ISNAN(a)||ISNAN(b))){
                     cross_moms[h+(v-t-1)*(*nbins-1)]+=0.5*(x-y)*(a-b);
                     cross_lbins[h+(v-t-1)*(*nbins-1)]+=1;
                 }}}}}
        }
               }}}
  return;
}

/***********************************************************************************************************************************/
// variogram cloud:
void Cloud_Variogram2(double *bins, double *coordx, double *coordy, double *coordt,double *data, int *lbins, double *moms, int *nbins)
{
  int  h=0,i=0, j=0, n=0;double lags=0.0,x,y;
 //Computes the cloud moments:
  for(i=0;i<(ncoord[0]-1);i++)
    for(j=(i+1);j<ncoord[0];j++){
          dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
      bins[h]=lags;
        x=data[n+i ];  y=data[n+j ];
        if(!(ISNAN(x)||ISNAN(y))){
	        moms[h]+=0.5*pow(x-y,2);
        lbins[h]=1;
        h++;}}
  return;
}
/***********************************************************************************************************************************/
/***********************************************************************************************************************************/




// Least square method for Gaussian spatial-temporal random field:
void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0, i=0, u=0;
  double vario=0.0, varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis[1],nuis[2],par); 
      *res=*res-pow(varhat-vario,2);// Computes the least squares
      i++;}
  return;
}
// Weighted least square method for Gaussian spatial-temporal random field:
void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		       int *nbins, int *nbint, double *nuis, double *par, double *res)
{
  int h=0,i=0,u=0;
  double vario=0.0,varhat=0.0;
  //Checks the nuisance parameters (nugget, sill and corr):
  if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2) {
    *res=LOW; return;}
  // Computes the least squares:
  for(u=0;u<*nbint;u++)
    for(h=0;h<(*nbins-1);h++){
      vario=moms[i]/lbins[i];// Computes the empirical variogram
      varhat=Variogram(cormod,0.5*(bins[h]+bins[h+1]),bint[u],nuis[1],nuis[2],par);
      if(vario) *res=*res-pow(varhat-vario,2)*(lbins[i]/pow(vario,2));// Computes the weighted least squares
      i++;}
  return;
}


