#include "header.h"


/***********************************************************************************************************************************/
/***********************************************************************************************************************************/
// binned spatial lorelogram:

void Binned_Lorelogram2(double *bins, double *coordx, double *coordy, double *coordt,double *data, int *lbins, double *moms,int *nbins)
{
  int h=0, i=0, j=0, n=0;
  double lags=0.0,step=0.0,*n11,*n10,*n01,*n00,*mm;

    
    n11=(double *) Calloc(*nbins-1,double);n10=(double *) Calloc(*nbins-1,double);
    n01=(double *) Calloc(*nbins-1,double);n00=(double *) Calloc(*nbins-1,double);

  mm=(double *) R_alloc(2, sizeof(double));
  Maxima_Minima_dist(mm, coordx, coordy, ncoord,type,REARTH); // computing max and min distances
  if(maxdist[0]<mm[1]) mm[1]=maxdist[0];
  //Set the binnes step:
  step=(mm[1]-mm[0])/(*nbins-1);
  bins[0]= mm[0];
  //define bins:
  for(h=1;h<*nbins;h++)
    bins[h]=bins[h-1]+step;
  //Computes the binned statistics:
  for(i=0;i<(ncoord[0]-1);i++)
    for(j=(i+1);j<ncoord[0];j++){
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
      if(lags<=*maxdist){
	for(h=0;h<(*nbins-1);h++)
	  if((bins[h]<=lags) && (lags<bins[h+1])){
	      if(data[n+i ] && data[n+j ]) n11[h]++;if(data[n+i ] && !data[n+j ]) n10[h]++;
	      if(!data[n+i ] && data[n+j ]) n01[h]++;if(!data[n+i ] && !data[n+j ]) n00[h]++;
	      }}}
// computing log odds ration in each bin
 for(h=0;h<(*nbins-1);h++){
   if(n11[h]&&n10[h]&&n01[h]&&n00[h]){
     moms[h]=log((n11[h]*n00[h])/(n01[h]*n10[h]));
     lbins[h]=1;}
   else{
     moms[h]=1;
     lbins[h]=0;}}
        Free(n11);Free(n10);Free(n01);Free(n00);
  return;
}
/***********************************************************************************************************************************/
/***********************************************************************************************************************************/
// binned spatial-temporal variogram:
void Binned_Lorelogram_st2(double *bins, double *bint, double *coordx, double *coordy, double *coordt,double *data, int *lbins, int *lbinst,
			 int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint)
{
  int h=0, i=0, j=0;
  int q=0, t=0, u=0, v=0;
  double x,y,lags=0.0,lagt=0.0,step=0.0,*n11s,*n10s,*n01s,*n00s,*n11t;
  double *n10t,*n01t,*n00t,*n11,*n10,*n01,*n00,*mm,*tt;
  n11s=(double *) Calloc((*nbins-1), double);
  n10s=(double *) Calloc((*nbins-1), double);
  n01s=(double *) Calloc((*nbins-1), double);
  n00s=(double *) Calloc((*nbins-1), double);
  n11t=(double *) Calloc((*nbint), double);
  n10t=(double *) Calloc((*nbint), double);
  n01t=(double *) Calloc((*nbint), double);
  n00t=(double *) Calloc((*nbint), double);
  n11=(double *) Calloc((*nbins-1)*(*nbint), double);
  n10=(double *) Calloc((*nbins-1)*(*nbint), double);
  n01=(double *) Calloc((*nbins-1)*(*nbint), double);
  n00=(double *) Calloc((*nbins-1)*(*nbint), double);
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
  //computes the empirical variogram:


  for(i=0;i<ncoord[0];i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<ncoord[0];j++){
	if(i==j){//computes the marignal temporal lorelogram:
	  for(v=(t+1);v<*ntime;v++){
	       lagt=fabs(coordt[t]-coordt[v]);
	    if(lagt<=*maxtime){
	      for(u=0;u<*nbint;u++){
		if(is_equal (bint[u],lagt)){
              x=data[(t+i * *ntime)];
              y=data[(v+i * *ntime)];
              if(!(ISNAN(x)||ISNAN(y))){
		    if( x&&y ) n11t[u]++;if( x&&!y) n10t[u]++;if(!x&&y ) n01t[u]++;if(!x&&!y) n00t[u]++;
		    }}}}}}
	else{
                    lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
	  for(v=0;v<*ntime;v++){
	    if(t==v){// computes the marginal spatial lorelogram:
	      if(lags<=*maxdist){
          
		for(h=0;h<(*nbins-1);h++){
		  if((bins[h]<=lags) && (lags<bins[h+1])){
                x=data[(t+i * *ntime)];
                y=data[(t+j * *ntime)];
                if(!(ISNAN(x)||ISNAN(y))){
		      if(x&&y) n11s[h]++;if(x&&!y) n10s[h]++;if(!x&&y) n01s[h]++;if(!x&&!y) n00s[h]++;}}}}}
	    else{// computes the spatial-temporal lorelogram:
	         lagt=fabs(coordt[t]-coordt[v]);
	      if(lags<=*maxdist && lagt<=*maxtime){
		q=0;
		for(h=0;h<(*nbins-1);h++){
		  for(u=0;u<*nbint;u++){
		    if((bins[h]<=lags) && (lags<bins[h+1]) && (is_equal (bint[u],lagt))){
                  x=data[(t+i * *ntime)];
                  y=data[(v+j * *ntime)];
                  if(!(ISNAN(x)||ISNAN(y))){
			if(x&&y)  n11[q]++;
			if(x&&!y)  n10[q]++;
			if(!x&&y)  n01[q]++;
			if(!x&&!y)  n00[q]++;
                  }}
			q++;}}}}}}}
    }}
  // computing  space time log odds ratio in each spacetime bin
  q=0;
  for(h=0;h<(*nbins-1);h++){
   if(!n10s[h]||!n01s[h]||!n11s[h]||!n00s[h]||((n11s[h]*n00s[h])==(n10s[h]*n01s[h])))
     moms[h]=0;
   else moms[h]=log((n11s[h]*n00s[h])/(n01s[h]*n10s[h]));}
   for(u=0;u<(*nbint);u++){
     if(!n11t[u]||!n01t[u]||!n10t[u]||!n00t[u]||((n11t[u]*n00t[u])==(n10t[u]*n01t[u])))
       momt[u]=0;
     else momt[u]=log((n11t[u]*n00t[u])/(n01t[u]*n10t[u]));}
   for(q=0;q<((*nbins-1)*(*nbint));q++){
     if(!n11[q]||!n01[q]||!n10[q]||!n00[q]||((n11[q]*n00[q])==(n10[q]*n01[q])))
       momst[q]=0;
     else momst[q]=log((n11[q]*n00[q])/(n01[q]*n10[q]));}
  ;
  Free(n10s);Free(n01s);Free(n00s);Free(n11t);Free(n10t);Free(n01t);
  Free(n00t);Free(n11); Free(n10);Free(n01);Free(n00);
  return;
}
/***********************************************************************************************************************************/

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
/*
// binned spatial-temporal variogram:
void Binned_Variogram_st2(double *bins, double *bint, double *coordx, double *coordy, double *coordt,double *data, int *lbins, int *lbinst,
       int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint)
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


//Rprintf("%d %d %f %f \n", ntime[0],ncoord[0],*maxtime,*maxdist);
for(t=0;t<ntime[0];t++){
    for(i=0;i<ncoord[0];i++){
      for(v=t;v<ntime[0];v++){

      if(t==v){
         for(j=i+1;j<ncoord[0];j++){
          lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
      if(lags<=*maxdist){
  // marginal spatial
        for(h=0;h<(*nbins-1);h++){
      if((bins[h]<=lags) && (lags<bins[h+1])){
                   x=data[i+ncoord[0]*t];  // x=data[t+ntime[0]*i];          
                   y=data[j+ncoord[0]*v];    //y=data[v+ntime[0]*j];  
        if(!(ISNAN(x)||ISNAN(y))){
          moms[h]+=0.5*pow(x-y, 2);
          lbins[h]+=1;}}}
    }}}
  else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ncoord[0];j++){

if(i==j){            // marginal temporal
  if(lagt<=*maxtime){
            for(u=0;u<*nbint;u++){
     if(is_equal (bint[u],lagt)){
             x=data[i+ncoord[0]*t];  // x=data[t+ntime[0]*i];          
             y=data[j+ncoord[0]*v];    //y=data[v+ntime[0]*j];      
              if(!(ISNAN(x)||ISNAN(y))){
        momt[u]+=0.5*pow(x-y, 2);
        lbint[u]+=1;
        }}}}}

else{
         lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);

       if(lags<=*maxdist && lagt<=*maxtime){
    q=0;
    for(h=0;h<(*nbins-1);h++)
      for(u=0;u<*nbint;u++){
        if((bins[h]<=lags) && (lags<bins[h+1]) && (is_equal(bint[u],lagt))){
             x=data[i+ncoord[0]*t];  // x=data[t+ntime[0]*i];          
             y=data[j+ncoord[0]*v];    //y=data[v+ntime[0]*j];  
               if(!(ISNAN(x)||ISNAN(y))){
      momst[q]+=0.5*pow(x-y, 2);
      lbinst[q]+=1;}
        q++;}}}}
    }}
  }}}
  return;
}

*/



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
/*******************************************/
/*******************************************/
   for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
  if(t==v){// computes the marginal spatial variogram:
             for(j=i+1;j<ns[v];j++){
                         lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
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
                    x=data[(i+NS[t])];y=data[(j+NS[v])];
                    if(!(ISNAN(x)||ISNAN(y))){
                    momt[u]+=0.5*pow(x-y, 2);
                    lbint[u]+=1;}
                  }}}}
          else{// computes the spatial-temporal variogram:
                  lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                  if(lags<=*maxdist && lagt<=*maxtime){
                    q=0;
                     for(h=0;h<(*nbins-1);h++){
                      for(u=0;u<*nbint;u++){
        if((bins[h]<=lags) && (lags<bins[h+1]) && (is_equal(bint[u],lagt))){
                 x=data[(i+NS[t])];y=data[(j+NS[v])];
               if(!(ISNAN(x)||ISNAN(y)))   {
                                            momst[q]+=0.5*pow(x-y, 2);
                                            lbinst[q]+=1;
                                          } }
        q++;}}

                     }
                  }
            }
       }
     }}}
 }    


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

  for(t=0;t<ntime[0];t++){
 for(v=t;v<ntime[0];v++){
   if(t==v){  // computes the marginal spatial variograms:
            for(i=0;i<ns[t];i++){
          for(j=i+1;j<ns[v];j++){
                  lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
        if(lags<=dista[t][v]) {


    for(h=0;h<(*nbins-1);h++)     {
      if((bins[h]<=lags) && (lags<bins[h+1])){
               x=data[(i+NS[t])];
               y=data[(j+NS[v])];
              if(!(ISNAN(x)||ISNAN(y))){
                  marg_moms[h+t*(*nbins-1)]+=0.5*pow(x-y,2);
                  marg_lbins[h+t*(*nbins-1)]+=1;
                  
              }}}


            }}}}
   else{// computes the   cross  variogram:
             for(i=0;i<ns[t];i++){
          for(j=i+1;j<ns[v];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
        if(lags<=dista[t][v]) {
    for(h=0;h<(*nbins-1);h++){
        if((bins[h]<=lags) && (lags<bins[h+1])){
                x=data[(i+NS[t])];y=data[(j+NS[t])];
                a=data[(i+NS[v])];b=data[(j+NS[v])];
                 if(!(ISNAN(x)||ISNAN(y)||ISNAN(a)||ISNAN(b))){
                     cross_moms[h+(v-t-1)*(*nbins-1)]+=0.5*(x-y)*(a-b);
                     cross_lbins[h+(v-t-1)*(*nbins-1)]+=1;
                 }}}}}}}
  }}
  return;
}

/*

void Binned_Variogram_biv2(double *bins, double *coordx, double *coordy, double *coordt,double *data, int *cross_lbins, double *cross_moms, int *nbins,
                          int *marg_lbins, double *marg_moms)
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

   for(t=0;t<*ntime;t++){
   for(v=t;v<*ntime;v++){
   if(t==v){  // computes the marginal spatial variograms:
      for(i=0;i<(ncoord[0]-1);i++){
      for(j=i+1;j<ncoord[0];j++){
                  lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
	      if(lags<=dista[t][v]) {


		for(h=0;h<(*nbins-1);h++)     {
		  if((bins[h]<=lags) && (lags<bins[h+1])){
              x=data[(t+*ntime*i)] ; y=data[(t+*ntime*j)];
              if(!(ISNAN(x)||ISNAN(y))){
                  marg_moms[h+t*(*nbins-1)]+=0.5*pow(x-y,2);
                  marg_lbins[h+t*(*nbins-1)]+=1;
                  
              }}}


            }}}}
   else{// computes the   cross  variogram:
	      for(i=0;i<ncoord[0];i++){
           for(j=i;j<ncoord[0];j++){
             lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
	      if(lags<=dista[t][v]) {
		for(h=0;h<(*nbins-1);h++){
		    if((bins[h]<=lags) && (lags<bins[h+1])){
                x=data[(t+*ntime*i)];y=data[(t+*ntime*j)];
                a=data[(v+*ntime*i)];b=data[(v+*ntime*j)];
                 if(!(ISNAN(x)||ISNAN(y)||ISNAN(a)||ISNAN(b))){
                     cross_moms[h+(v-t-1)*(*nbins-1)]+=0.5*(x-y)*(a-b);
                     cross_lbins[h+(v-t-1)*(*nbins-1)]+=1;
                 }}}}}}}
  }}
  return;
}
*/
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


