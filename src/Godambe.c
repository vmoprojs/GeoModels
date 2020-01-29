#include "header.h"

// Computing the gradient and  the components of J matrix for Godmabe matrix computation

void GodambeMat(double *betas,int *biv,double *coordx, double *coordy, double *coordt, int *cormod, double *data, int *dst,
		      double *eps,int *flagcor, int *flagnuis, int *grid, int *like, double *mean,int *model,double *NN, int *nbetas,
		      int *npar, int *nparc,int *nparcT, double *parcor, double *nuis, double *score,
		      double *sensmat, int *spt, int *type_lik, double *varimat,
		      int *vartype, double *winc, double *winstp,double *winct,double *winstp_t,int *weigthed, double *X,int *ns,int *NS)
{
   // Rprintf("winc: %f winstp:%f\n",*winc,*winstp);
  //---------- COMPUTATION OF THE GODAMBE MATRIX ---------//
  int *np;np=(int *) R_alloc(1, sizeof(int));  //numbers of effective pairs

  switch(*vartype)
    {
    case 1://------------ START EMPIRICAL ESTIMATION ------------//
   //   GodambeMat_emp(coordx,coordy,coordt,cormod,data,eps,flagcor,flagnuis,like,
   //		     model,npar,nparc,parcor,nuis,score,
   //		     sensmat,spt,varimat,type_lik,weigthed);
      break;//------------ END EMPIRICAL ESTIMATION ------------//
    case 2://------------ START SUB-SAMPLE ESTIMATION ------------//
 
            //Rprintf("START: Sensitivity!!!!!\n");

     Sensitivity(betas,biv,coordx,coordy,coordt,cormod,data,eps,flagcor,flagnuis,like,mean,model,NN,nbetas,
	  npar,nparc,nparcT,parcor,nuis,np,score,sensmat,spt,type_lik,weigthed,X,ns,NS);
            //Rprintf("END: Sensitivity!!!!!\n");
   
            if(!*spt&&!*biv) Vari_SubSamp(betas,coordx,coordy,coordt,cormod,data,dst,eps,flagcor,flagnuis,
                                          grid,like,mean,model,NN,nbetas,npar,nparc,nparcT,nuis,np,parcor, type_lik, 
                                          varimat,winc,winstp,weigthed,X);
            else {
                if(*spt)   Vari_SubSamp_st2(betas,coordx,coordy,coordt,cormod,data,dst,eps,flagcor,flagnuis,
                                         like,mean,model,NN,nbetas,npar,nparc,nparcT,nuis,np,parcor,type_lik,
                                         varimat,winc,winstp,winct,winstp_t,weigthed,X,ns,NS);
                      //Vari_SubSamp_st(betas,coordx,coordy,coordt,cormod,data,dst,eps,flagcor,flagnuis,
                      //                   like,mean,model,NN,nbetas,npar,nparc,nuis,np,parcor,type_lik,
                        //                 varimat,winc,winstp,weigthed,X);


                if(*biv) Vari_SubSamp_biv(betas,coordx,coordy,coordt,cormod,data,dst,eps,flagcor,flagnuis,
                                          grid,like,model,NN,npar,nparc,nparcT,nuis,np,parcor, type_lik, varimat,winc,winstp,weigthed,X,ns,NS);

            }
      break;//------------ END SUB-SAMPLE ESTIMATION ------------//
    
    }
  return;
}



// Compute the Sensitivity matrix of a random field:
void Sensitivity(double *betas,int *biv,double *coordx,double *coordy,double *coordt,int *cormod,  double *data, double *eps, int *flagcor, int *flagnuis, int *like,
		 double *mean,int *model, double *NN,int *nbetas, int *npar, int *nparc,int *nparcT, double *parcor, double *nuis, int *np,double *score,
		 double *sensmat, int *spt,  int *type_lik,int *weigthed,double *Z,int *ns, int *NS)
{
  double *grad,*gradcor;
  grad=(double *) R_alloc(*npar,sizeof(double));// gradient of the ijth composite log-likelihood
  gradcor=(double *) R_alloc(*nparc, sizeof(double));// gradient of the correlation
  switch(*model)// Compute the Sensitivity matrix
    {
    /***********************************************/    
    case 1: //gaussian
    case 2: //binary
    case 11: //binomial
    case 12: // studentT
    case 13: // wrapped
    case 14: //geom
    case 16: //negtivebinomial
    case 10: // skew
    case 20: //sinh
    case 21: // gamma
    case 22: //loggaussian
    case 24: //loglogistic
    case 25: // logistic
    case 26: //Weibull
    case 27: //twopieceT
    case 30: //poisson
    case 34: //tukeyh
    // more........

        if( (!*spt) && (!*biv)) 
          Sens_Pair(betas,coordx,coordy,coordt,cormod,data,eps,flagcor,flagnuis,NN,nuis,np,nbetas,
				   npar,nparc,nparcT,mean,model,parcor,score,sensmat,weigthed,Z);
         else {
             if(*spt) 
          Sens_Pair_st(betas,coordx,coordy,coordt,cormod,data,eps,flagcor,flagnuis,NN,nuis,np,nbetas,
            npar,nparc,nparcT,mean,model,parcor,score,sensmat,weigthed,Z,ns,NS);
             if(*biv)
           Sens_Pair_biv(betas,coordx,coordy,coordt,cormod,data,eps,flagcor,flagnuis,NN,nuis,np,
                                         npar,nparc,nparcT,mean,model,parcor,score,sensmat,weigthed,Z,ns,NS);
         }
    break;
    /***********************************************/  
       /***********************************************/ 
    }

  return;
}

 

// Compute the Sensitivity matrix for the pairwise composite  likelihood:
void Sens_Pair(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, 
    int *flagcor, int *flagnuis, double *NN,double *nuis, int *np, int *nbetas,int *npar, int *nparc,int *nparcT, 
             double *mean, int *model,double *parcor,double *score, double *sensmat,int *weigthed,double *Z)
{
  // Initialization variables:
  int b=0,i=0,j=0,d=0,qq=0,o=0,nsens=0;
  double lags=0.0,rho=0.0, *grad, *gradcor, *sens,weights;
      double *Xi,*Xj,**X;
     Xi=(double *) Calloc(*nbetas,double);
     Xj=(double *) Calloc(*nbetas,double);
     X= (double **) Calloc(ncoord[0],double *);
     for(i=0;i<ncoord[0];i++) X[i]=(double *) Calloc(nbetas[0],double);

  // matrix of covariates 
      qq=0; 
     for(i=0;i<ncoord[0];i++){
     for(j=0;j<nbetas[0];j++){       
         X[i][j]=Z[qq];qq++;
    
     }}
  nsens=*npar * (*npar+1)/2;
  gradcor=(double *) Calloc(*nparc,double);// Correlation gradient  (double *) Calloc(*nparc,double);
  grad=(double *) Calloc(*npar,double);// ijth component of the gradient
  sens=(double *) Calloc(nsens,double);// One sensitive contribute
//Rprintf("%d %d \n ",ncoord[0],*model);
      for(i=0; i<(ncoord[0]-1);i++){
    for(j=(i+1); j<ncoord[0];j++){
    /*****************************************************/    
    /*****************************************************/  
         if(!ISNAN(data[i])&&!ISNAN(data[j]) ){
    //Compute the correlation function for the elements i,j
    lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
if(lags<maxdist[0]){
    rho=CorFct(cormod,lags,0,parcor,0,0);
    //Compute the gradient for the given correlation
   // if(*model==1) GradCorrFct(rho,cormod,eps[0],flagcor,gradcor,lags,0,0,0,parcor);
        /***/
        for(o=0;o<*nbetas;o++) {Xi[o]=X[i][o];
                                Xj[o]=X[j][o];}
/******************************************************/  
      switch(*model){
      case 1: // gaussian
     //  Grad_Pair_Gauss(rho,flagnuis,flagcor,gradcor,grad,npar,nparc,nbetas[0],nuis,
		   //   data[i]-mean[i],data[j]-mean[j],Xi,Xj);
        Grad_Pair_Gauss2(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas);  
      break;

      case 2:
      case 11: // binomial
       Grad_Pair_Binom(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas);   
       break;

     case 10://
       Grad_Pair_Skewgauss(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
     break;

        case 12:
      Grad_Pair_StudenT(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;

      case 13://wrapped
       Grad_Pair_Wrapped(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
      break;

      case 14:
      case 16: // negative binomial
        Grad_Pair_Binomneg(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
      break;

       case 20: //
       Grad_Pair_Sinh(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;
            case 34: //
       Grad_Pair_Tukeyh(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;


       case 21: // gamma
           Grad_Pair_Gamma(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;

         case 30: // poisson
           Grad_Pair_Poisson(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;

       case 22: // loggaus
           Grad_Pair_LogGauss(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;

       case 24: // loglogistic
           Grad_Pair_LogLogistic(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;

       case 25: // logistic
           Grad_Pair_Logistic(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;

       case 26: // Weibull
           Grad_Pair_Weibull(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i],data[j],mean[i],mean[j],Xi,Xj,X,i,j,betas); 
       break;
       //case... more
       //break;
   }
/*************************************************************/
      //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
     // Sens_Pair_Gauss_ij(rho,flagnuis,gradcor,npar,nbetas[0],nparc,nuis,sens,Xi,Xj);
        if(*weigthed) weights=CorFunBohman(lags,maxdist[0]);
        else          weights=1;
        //Rprintf("%f %f %f\n",score[0],score[1],score[2]);
     // for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l]*weights;
      for(d=0;d<*npar;d++) {score[d]=score[d]+grad[d]*weights;}
      b++;
     }}}}
     // Rprintf("ddfdfdfdf. %f %f %f %f \n",score[0],score[1],score[2],score[3]);
    *np=b;
    Free(gradcor);Free(grad);Free(sens);
    Free(Xi);Free(Xj);
    for(i=0;i<ncoord[0];i++) { Free(X[i]);}Free(X);
  return;
}





// Comgaufte the Sensitivity matrix for the space time pairwise composite  likelihood:
void Sens_Pair_st(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, 
         double *eps, int *flagcor, int *flagnuis, double *NN,double *nuis, int *np,
			int *nbetas, int *npar,int *nparc,int *nparcT, double *mean, int *model,
            double *parcor, double *score, double *sensmat,  int *weigthed,double *Z, int *ns, int *NS)
{
  // Initialization variables:
  int  b=0,d=0,i=0,nsens=0,j=0,t=0,v=0,o=0,qq=0;
  double lags=0.0,lagt=0.0,rho=0.0, *grad, *gradcor,weigths=1.0;
   // Rprintf("%d %d  %f %f \n",ncoord[0],*ntime,maxdist[0],maxtime[0]);
    nsens=*npar * (*npar+1)/2;
    gradcor=(double *) Calloc(*nparc,double);// Correlation gradient  (double *) Calloc(*nparc,double);
    grad=(double *) Calloc(*npar,double);// ijth component of the gradient
    //sens=(double *) Calloc(nsens,double);// One sensitive contribute
    double *Xi,*Xj,**X;
     Xi=(double *) Calloc(*nbetas,double);
     Xj=(double *) Calloc(*nbetas,double);
     X=(double **) Calloc(ncoord[0]*ntime[0],double *);
    if(cdyn[0]==0)
    {
        for(i=0;i<(ncoord[0]*ntime[0]);i++) {X[i]=(double *) Calloc(nbetas[0],double);  }
        /***/
        for(i=0;i<(ncoord[0]*ntime[0]);i++){
            for(j=0;j<nbetas[0];j++){
                X[i][j]=Z[qq];
                qq++;
            }}
    }
    else
    {
        for(i=0;i<(ncoord[0]);i++) {X[i]=(double *) Calloc(nbetas[0],double);  }
        /***/
        for(i=0;i<(ncoord[0]);i++){
            for(j=0;j<nbetas[0];j++){
                X[i][j]=Z[qq];
                qq++;
            }}
    }
    
    

      /********************************************************/
      /********************************************************/
for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){

         if(!ISNAN(data[i+NS[t]]-mean[i+NS[t]])&&!ISNAN(data[j+NS[v]]-mean[j+NS[v]]) ){
                   lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
	             if(lags<=maxdist[0]){

                b++;
                                        for(o=0;o<*nbetas;o++) {
                                            Xi[o]=X[i+NS[t]][o];
                                            Xj[o]=X[j+NS[v]][o];
                                        }
	             // Compute the temporal correlation function for for a given pair:
	             rho=CorFct(cormod,lags,0,parcor,0,0);
	             // Compute the gradient for the temporal correlation function:
	     //  if(*model==1)      GradCorrFct(rho,cormod,eps[0],flagcor,gradcor,lags,0,0,0,parcor);
       /******************************************************/        
	switch(*model){
      case 1:
		  //       Grad_Pair_Gauss(rho,flagnuis,flagcor,gradcor,grad,npar,nparc,nbetas[0],nuis,
		//		 data[i+NS[t]]-mean[i+NS[t]],
				// data[j+NS[v]]-mean[j+NS[v]],
       //          Xi,Xj);
             Grad_Pair_Gauss2(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
      case 2:
      case 11:
       Grad_Pair_Binom(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;

       case 10:
       Grad_Pair_Skewgauss(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
    case 12:
       Grad_Pair_StudenT(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],Xi,Xj,X , i+NS[t],j+NS[v],betas);
     break;
      case 13:
       Grad_Pair_Wrapped(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
        case 14:
        case 16:
       Grad_Pair_Binomneg(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
         case 20:
       Grad_Pair_Sinh(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
            case 34:
       Grad_Pair_Tukeyh(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
        case 21:
       Grad_Pair_Gamma(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;

        case 30:
       Grad_Pair_Poisson(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;

        case 22:
       Grad_Pair_LogGauss(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
          case 24:
       Grad_Pair_LogLogistic(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
           case 25:
       Grad_Pair_Logistic(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
             case 26:
       Grad_Pair_Weibull(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
  }
        /********************************************************/

		         //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
		         //Sens_Pair_Gauss_ij(rho,flagnuis,gradcor,npar,nbetas[0],nparc,nuis,sens,Xi,Xj);
                 if(*weigthed) weigths=CorFunBohman(lags,maxdist[0]);
                 else          weigths=1;
		         //for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l]*weigths;
                 for(d=0;d<*npar;d++) {
                     if(R_FINITE(grad[d]))
                  score[d]=score[d]+grad[d]*weigths;}
	             }}}}

else {
             lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          if(!ISNAN(data[i+NS[t]]-mean[i+NS[t]])&&!ISNAN(data[j+NS[v]]-mean[j+NS[v]])){
            lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=maxdist[0]&&lagt<=maxtime[0]){
                           b++;
                                   for(o=0;o<*nbetas;o++) {
                                            Xi[o]=X[i+NS[t]][o];
                                            Xj[o]=X[j+NS[v]][o];
                                        }
	      // Compute the spatial-temporal correlation function for a given pair:
	      rho=CorFct(cormod,lags,lagt,parcor,0,0);
	      // Compute the gradient of the patial-temporal correlation function:
	    //  if(*model==1)GradCorrFct(rho,cormod,eps[0],flagcor,gradcor,lags,lagt,0,0,parcor);
		     /******************************************************/        
    switch(*model){
      case 1:
                // Grad_Pair_Gauss(rho,flagnuis,flagcor,gradcor,grad,npar,nparc,nbetas[0],nuis,
                 //data[i+NS[t]]-mean[i+NS[t]],
                 //data[j+NS[v]]-mean[j+NS[v]],
                 //Xi,Xj);
                  Grad_Pair_Gauss2(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
      case 2:
      case 11:
       Grad_Pair_Binom(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
       case 10:
       Grad_Pair_Skewgauss(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
     case 12:
        Grad_Pair_StudenT(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
         case 13:
       Grad_Pair_Wrapped(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
        case 14:
        case 16:
       Grad_Pair_Binomneg(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
        case 20:
       Grad_Pair_Sinh(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
           case 34:
       Grad_Pair_Tukeyh(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
        case 21:
       Grad_Pair_Gamma(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
         case 30:
       Grad_Pair_Poisson(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
        case 22:
       Grad_Pair_LogGauss(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
             case 24:
       Grad_Pair_LogLogistic(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
           case 25:
       Grad_Pair_Logistic(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
               case 26:
       Grad_Pair_Weibull(rho,cormod,flagnuis,flagcor,gradcor,grad,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                 data[i+NS[t]],data[j+NS[v]],mean[i+NS[t]],mean[j+NS[v]],
                 Xi,Xj,X , i+NS[t],j+NS[v],betas); 
      break;
  }
        /********************************************************/
		  //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
		 // Sens_Pair_Gauss_ij(rho,flagnuis,gradcor,npar,nbetas[0],nparc,nuis,sens,Xi,Xj);
              if(*weigthed) weigths=CorFunBohman(lags,maxdist[0])*CorFunBohman(lagt,maxtime[0]);
              else          weigths=1;
		 // for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l]*weigths;
              //for(d=0;d<*npar;d++) score[d]=score[d]+grad[d]*weigths;

                    for(d=0;d<*npar;d++) {
                     if(R_FINITE(grad[d]))
                  score[d]=score[d]+grad[d]*weigths;}
		}}}}}}}
          *np=b;
    Free(gradcor);// Correlation gradient  (double *) Calloc(*nparc,double);
    Free(grad);// ijth component of the gradient
    //Free(sens);// One sensitive contribute
    Free(Xi);Free(Xj);
    
    if(cdyn[0] ==0)
    {
        for(i=0;i<(ncoord[0]*ntime[0]);i++) { Free(X[i]);}
        
    }else
    {
        for(i=0;i<(ncoord[0]);i++) { Free(X[i]);}
    }
    
    Free(X);
  return;
}




// Compute the Sensitivity matrix for the space time pairwise composite Gaussian likelihood for bivariate GRF:
void Sens_Pair_biv(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, 
    int *flagcor, int *flagnuis, double *NN,double *nuis, int *np,int *npar, int *nparc,int *nparcT,double *mean, int *model,
     double *parcor, double *score, double *sensmat,int *weigthed,double *Z, int *ns, int *NS)
{
    // Initialization variables:
    int  b=0,d=0,i=0,l=0,nsens=0,j=0,t=0,v=0;
    double lags=0.0,rhott=0.0,rhovv=0.0,rhotv=0.0,rhovt=0.0, *grad, *gradcortt,*gradcortv,*gradcorvt,*gradcorvv, *sens,weigths;
    nsens=*npar * (*npar+1)/2;
    gradcortt=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcortv=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvt=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvv=(double *) Calloc(*nparc,double);// Correlation gradient
    grad=(double *) Calloc(*npar,double);// ijth component of the gradient
    sens=(double *) Calloc(nsens,double);// One sensitive contribute
      for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=dista[t][v]){
                            // Compute the  correlation function for for a given pair:
                            // Compute the l correlation function for a given pair:
                            rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                            rhotv=CorFct(cormod,lags,0,parcor,t,v);rhovt=rhotv;

                            // Compute the gradient of the patial-temporal correlation function:
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lags,0,t,v,parcor);
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);
/**********************************************************/
                                 switch(*model){
      case 1: 
                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,
                                                    gradcortt,gradcortv,gradcortv,gradcorvv,grad,npar,nuis,
                                                    data[(i+NS[t])],
                                                    data[(j+NS[v])]);
                                //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
                                Sens_Pair_Gauss_biv_ij(rhott,rhotv,rhovt,rhovv,gradcortt,gradcortv,gradcortv,gradcorvv,
                                 flagnuis,npar,nparc,nuis,sens);

      break;
      }    
/**********************************************************/                      
                                weigths=1;
                                if(*weigthed&&t>v) weigths=CorFunBohman(lags,dista[t][v]);
                                for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l]*weigths;
                                for(d=0;d<*npar;d++) score[d]=score[d]+grad[d]*weigths;
                            
                            b++;}}}
            else {
                for(j=0;j<ns[v];j++){
         lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        if(lags<=dista[t][v]){
    
                        // Compute the l correlation function for a given pair:
                            rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                            rhotv=CorFct(cormod,lags,0,parcor,t,v);rhovt=rhotv;
                            // Compute the gradient of the patial-temporal correlation function:
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lags,0,t,v,parcor);
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);
                  switch(*model){
      case 1:
                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,
                                                    gradcortt,gradcortv,gradcortv,gradcorvv,grad,npar,nuis,
                                                   data[(i+NS[t])],
                                                    data[(j+NS[v])]);
                                //ADD TO THE SENSITIVITY MATRIX THE CONTRIBUTE OF THE GIVEN PAIR
                               Sens_Pair_Gauss_biv_ij(rhott,rhotv,rhovt,rhovv,gradcortt,gradcortv,gradcortv,gradcorvv,
                                                         flagnuis,npar,nparc,nuis,sens);
      break;
  }
                                if(*weigthed) weigths=CorFunBohman(lags,dista[t][v]);
                                else          weigths=1;
                                for(l=0;l<nsens;l++) sensmat[l]=sensmat[l]-sens[l]*weigths;
                                for(d=0;d<*npar;d++) score[d]=score[d]+grad[d]*weigths;
                            b++;}}}
                        }}}
    *np=b;
    Free(gradcortt);// Correlation gradient
    Free(gradcortv);// Correlation gradient
    Free(gradcorvt);// Correlation gradient
    Free(gradcorvv);// Correlation gradient
    Free(grad);// ijth component of the gradient
    Free(sens);// One sensitive contribute
    return;
}





// Compute the Sensitivity matrix for the space time Bipairwise composite Gaussian likelihood for bivariate GRF:
void Sens_Cond_Gauss_biv(double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
                         int *npar, int *nparc, double *parcor, double *score, double *sensmat,int *weigthed)
{
    // Initialization variables:
    int  N=4,b=0,d=0,i=0,n=0,m=0,nsens=0,j=0,t=0,v=0,r=0;
    double dd,rhottii=0.0,rhovvii=0.0, rhotvij=0.0, rhovtij=0.0, rhottij=0.0, rhovvij=0.0, rhotvii=0.0,rhovtii=0.0,weigths;
    double lags=0.0, *grad, *sens;
    nsens=*npar * (*npar+1)/2;
    double *gradcorttii,*gradcorvvii ,*gradcorvtii  ,*gradcortvii ,*gradcorttij ,*gradcorvvij ,*gradcorvtij,*gradcortvij;
    gradcorttii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvvii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvtii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcortvii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorttij=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvvij=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvtij=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcortvij=(double *) Calloc(*nparc,double);// Correlation gradient
    int *indx;
    indx = (int *) R_alloc(N, sizeof(int));
    double *col;
    col =(double *) R_alloc(N, sizeof(double));
    double **M;
    M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}
    
    double **inverse;
    inverse= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){inverse[i]=(double *) Calloc(N,double);}
    double *dat; //data
    dat=(double *) Calloc(N,double);
     grad=(double *) Calloc(*npar,double);
     sens=(double *) Calloc(nsens,double);
    
    for(t=0;t<*ntime-1;t++){
        for(v=t+1;v<*ntime;v++){
                for(i=0;i<ncoord[0]-1;i++){
                    for(j=i+1;j<ncoord[0];j++){
                          lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                        if(lags<=dista[t][v]){
                        rhottii=CorFct(cormod,0,0,parcor,t,t);rhovvii=CorFct(cormod,0,0,parcor,v,v);
                        rhottij=CorFct(cormod,lags,0,parcor,t,t);rhovvij=CorFct(cormod,lags,0,parcor,v,v);
                        rhotvii=CorFct(cormod,0,0,parcor,t,v);rhovtii=CorFct(cormod,0,0,parcor,v,t);
                        rhotvij=CorFct(cormod,lags,0,parcor,t,v);rhovtij=CorFct(cormod,lags,0,parcor,v,t);
                        M[0][0]=rhottii;           M[0][1]=rhottij;       M[0][2]=rhotvii;        M[0][3]=rhotvij;
                        M[1][0]=M[0][1];           M[1][1]=rhottii;       M[1][2]=rhovtij;        M[1][3]=rhovtii;
                        M[2][0]=M[0][2];           M[2][1]= M[1][2];      M[2][2]=rhovvii;        M[2][3]=rhovvij;
                        M[3][0]=M[0][3];           M[3][1]= M[1][3];      M[3][2]=M[2][3];        M[3][3]=rhovvii;
                            ludcmp(M,N,indx,&dd);    //Lu decomposition
                            for(n=0;n<N;n++) {
                                for(m=0;m<N;m++) col[m]=0.0;
                                col[n]=1.0;
                                lubksb(M,N,indx,col);
                                for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                            }
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcorttii,0,0,t,t,parcor);    GradCorrFct(0,cormod,eps[0],flagcor,gradcorvvii,0,0,v,v,parcor);
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcorttij,lags,0,t,t,parcor); GradCorrFct(0,cormod,eps[0],flagcor,gradcorvvij,lags,0,v,v,parcor);
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcortvii,0,0,t,v,parcor);  GradCorrFct(0,cormod,eps[0],flagcor,gradcorvtii,0,0,v,t,parcor);
                            GradCorrFct(0,cormod,eps[0],flagcor,gradcortvij,lags,0,t,v,parcor);  GradCorrFct(0,cormod,eps[0],flagcor,gradcorvtij,lags,0,v,t,parcor);
                        for(r=0;r<*nrep;r++){
                            dat[0]=data[(t+*ntime*i)+r* *nrep];dat[1]=data[(t+*ntime*j)+r* *nrep];
                            dat[2]=data[(v+*ntime*i)+r* *nrep];dat[3]=data[(v+*ntime*j)+r* *nrep];
                                //Compute the gradient of the log pairwise likelihood
                            Grad_Cond_Gauss_biv(gradcorttii,gradcorvvii,gradcorvtii,gradcortvii ,gradcorttij ,gradcorvvij, gradcorvtij,gradcortvij,
                                                  inverse,flagnuis,grad,npar,nuis,dat,N);
                            Sens_Cond_Gauss_biv_ij(gradcorttii,gradcorvvii,gradcorvtii,gradcortvii,gradcorttij ,gradcorvvij,gradcorvtij,gradcortvij,
                                                   inverse,flagnuis,npar,nparc,nuis,N,sens);
                            if(*weigthed) weigths=CorFunBohman(lags,dista[t][v]);
                            else          weigths=1;
                            for(d=0;d<nsens;d++) {sensmat[d]=sensmat[d]-sens[d]*weigths;  }
                            for(d=0;d<*npar;d++) score[d]=score[d]+grad[d]*weigths;
                            }b++;}}}}}
    *np=b;
    Free(gradcorttii);// Correlation gradient
    Free(gradcorvvii);// Correlation gradient
    Free(gradcorvtii);// Correlation gradient
    Free(gradcortvii);// Correlation gradient
    Free(gradcorttij);// Correlation gradient
    Free(gradcorvvij);// Correlation gradient
    Free(gradcorvtij);// Correlation gradient
    Free(gradcortvij);// Correlation gradient
    Free(dat);
    Free(grad);
    Free(sens);
    for(i=0;i<N;i++)  {Free (inverse[i]);Free(M[i]);}
    Free(M);Free(inverse);
    return;
}



// Compute the Sensitivity matrix for the Gaussian pairwise likelihood for bivariate GRF:
void Sens_Pair_Gauss_biv_ij(double rhott,double rhotv,double rhovt,double rhovv,double *gradcortt,double *gradcortv,double *gradcorvt,double *gradcorvv,int *flag, int *npar,int *nparc, double *par, double *sensmat)
{
    // Initialization variables:
    double det,a,b;
    int  i=0, l=0, j=0, p=0;
    // Defines useful quantities:
    det=rhott*rhovv-rhotv*rhovt;

    //---- START COMPUTATION THE SENSITIVITY MATRIX-----//
    // Derivatives of the PLLik respect with the  first mean
    if(flag[0]==1)
    {
        sensmat[p]=rhovv/det; // Second derivative
        i++; p++; l++;
        // Mixed derivatives:
        if(flag[1]==1){ sensmat[p]=-0.5*(rhovt+rhotv)/det; p++; l++; }
        for(j=l;j<*npar;j++){ sensmat[p]=0; p++;}
    }
    //  // Derivatives of the PLLik respect with the  second mean
    if(flag[1]==1)
    {
        l = i; i++;
        sensmat[p]=rhott/det;
        p++; l++;
        for(j=l;j<*npar;j++){sensmat[p]=0;p++;}
    }

    //Derivative of the PLLik respect with the correlation parameters
     for(j=0;j<*nparc;j++)
     for(l=j;l<*nparc;l++)
     {
         a=2*det*(gradcortt[j]*gradcorvv[l]-gradcortv[j]*gradcorvt[l])   ;
         b=(rhott*gradcorvv[j]+rhovv*gradcortt[j]-2*rhotv*gradcortv[j])*(rhott*gradcorvv[l]+rhovv*gradcortt[l]-2*rhotv*gradcortv[l]);
     sensmat[p]=pow(2*det*det,-1)*(a-b);
     p++;
     }
    //---- END COMPUTATION OF THE SENSITIVITY MATRIX-----//
    return;
}



// Compute the Sensitivity matrix for the Gaussian pair of pairwise likelihood for bivariate GRF:
void Sens_Cond_Gauss_biv_ij(double *gradcorttii,double *gradcorvvii,double *gradcorvtii,double *gradcortvii ,
                            double *gradcorttij ,double *gradcorvvij, double *gradcorvtij,double *gradcortvij,
                            double **inverse,int *flag,int *npar, int * nparc, double *par,int N,double *sens)
{
    // Initialization variables:
    int i=0,l=0, j=0, p=0;double *uni1,*uni2;
    uni1=(double *) R_alloc(N, sizeof(double));
    uni2=(double *) R_alloc(N, sizeof(double));
    uni1[0]=1;uni1[1]=0;uni1[2]=1;uni1[3]=0;
    uni2[0]=0;uni2[1]=1;uni2[2]=0;uni2[3]=1;
    
    double **D1;
    D1= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){D1[i]=(double *) Calloc(N,double);}
    
    
    double **D2;
    D2= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){D2[i]=(double *) Calloc(N,double);}

    double **P1;
    P1= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){P1[i]=(double *) Calloc(N,double);}
    
    
    double **P2;
    P2= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){P2[i]=(double *) Calloc(N,double);}
    
    
    double **P3;
    P3= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){P3[i]=(double *) Calloc(N,double);}
    
    
    //for(m=0;m<N;m++){for(n=0;n<N;n++) { P3[m][n]=0; P2[m][n]=0; P1[m][n]=0;}}
    i=0;
    //---- START COMPUTATION THE SENSITIVITY MATRIX-----//
    // Derivatives of the PLLik respect with the  first mean
    if(flag[0]==1)
    {
        sens[p]=-QFORM(inverse,uni1,uni1,N); // Second derivative
        i++; p++; l++;
        // Mixed derivatives:
        if(flag[1]==1){ sens[p]=-QFORM(inverse,uni1,uni2,N);p++;l++ ;}
        for(j=l;j<*npar;j++){ sens[p]=0; p++;}
    }
    //  // Derivatives of the PLLik respect with the  second mean
    if(flag[1]==1)
    {
        l = i; i++;
        sens[p]=-QFORM(inverse,uni2,uni2,N);
        p++; l++;
     for(j=l;j<*npar;j++) {sens[p]=0; p++;}
    }
    for(j=0;j<*nparc;j++){
        for(l=j;l<*nparc;l++){
        D1[0][0]=gradcorttii[j];     D1[0][1]=gradcorttij[j];        D1[0][2]=gradcortvii[j];        D1[0][3]=gradcortvij[j];
        D1[1][0]=D1[0][1];           D1[1][1]=gradcorttii[j];        D1[1][2]=gradcorvtij[j];        D1[1][3]=gradcorvtii[j];
        D1[2][0]=D1[0][2];           D1[2][1]= D1[1][2];             D1[2][2]=gradcorvvii[j];        D1[2][3]=gradcorvvij[j];
        D1[3][0]=D1[0][3];           D1[3][1]= D1[1][3];             D1[3][2]=D1[2][3];              D1[3][3]=gradcorvvii[j];

        D2[0][0]=gradcorttii[l];     D2[0][1]=gradcorttij[l];        D2[0][2]=gradcortvii[l];        D2[0][3]=gradcortvij[l];
        D2[1][0]=D2[0][1];           D2[1][1]=gradcorttii[l];        D2[1][2]=gradcorvtij[l];        D2[1][3]=gradcorvtii[l];
        D2[2][0]=D2[0][2];           D2[2][1]= D2[1][2];             D2[2][2]=gradcorvvii[l];        D2[2][3]=gradcorvvij[l];
        D2[3][0]=D2[0][3];           D2[3][1]= D2[1][3];             D2[3][2]=D2[2][3];              D2[3][3]=gradcorvvii[l];
        Matrix_prod(inverse,D1,P1,N);
        Matrix_prod(inverse,D2,P2,N);Matrix_prod(P1,P2,P3,N);
        sens[p]=-0.5*Trace(P3,N);
        p++;}}
    
    for(i=0;i<N;i++)  {Free (D1[i]);Free(D2[i]);Free(P1[i]);Free(P2[i]);Free(P3[i]);}
    Free(D1);Free(D2);Free(P1);Free(P2);Free(P3);
    return;
}




/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
void Vari_SubSamp(double *betas,double *coordx, double *coordy, double *coordt,int *cormod, double *data,
                  int *dst,double *eps, int *flagcor, int *flagnuis, int *grid,
                  int *like, double *mean,int *model, double *NN, int *nbetas,int *npar, int *nparc,int *nparcT, double *nuis, int *np,
                  double *parcor, int *type_lik, double *varimat,
                  double *winc, double *winstp,int *weigthed,double *Z)
{
    double meanl,meanm,rho=0.0, lag=0.0, *gradcor,*gradient, *rangex, *rangey;
    double *ecoordx,*ecoordy,*scoordx,*scoordy,*sdata,*sumgrad,*subvari,*xgrid,*ygrid;
    double delta=0.0, deltax=0.0, deltay=0.0, dimwinx=0.0, dimwiny=0.0,weigths;
    double winstx=winstp[0], winsty=winstp[0];
    int *npts, numintx=0, numinty=0;
    int h=0,i=0,l=0,m=0,n=0,nsub=0,nvari=0,nwpair=0,p=0,q=0,qq=0,j=0,o=0;


    ////creating matrix
     double *Xl,*Xm,**X,**sX;
     Xl=(double *) Calloc(*nbetas,double);
     Xm=(double *) Calloc(*nbetas,double);
    
     //
     X= (double **) Calloc(ncoord[0],double *);
     for(i=0;i<ncoord[0];i++){
     X[i]=(double *) Calloc(nbetas[0],double);
             }
   //          
      sX= (double **) Calloc(ncoord[0],double *);
     for(i=0;i<ncoord[0];i++){
     sX[i]=(double *) Calloc(nbetas[0],double);
             }        
     /*********/  

    // matrix of covariates        
       qq=0; 
     for(i=0;i<ncoord[0];i++){
     for(j=0;j<nbetas[0];j++){       
         X[i][j]=Z[qq];
         qq++;
     }}
    /***********/   
    nvari=*npar * (*npar+1)/2;
    
    gradcor=(double *) Calloc(*nparc,double);
    gradient=(double *) Calloc(*npar,double);
    sumgrad=(double *) Calloc(*npar,double);
    subvari=(double *) Calloc(nvari,double);
    
    npts=(int *) R_alloc(1, sizeof(int));
    rangex=(double *) R_alloc(2, sizeof(double));
    rangey=(double *) R_alloc(2, sizeof(double));
    scoordx=(double *) Calloc(ncoord[0],double);
    scoordy=(double *) Calloc(ncoord[0],double);
    sdata=(double *) Calloc(ncoord[0],double);
    
    if(*grid){
        ecoordx=(double *) Calloc(ncoord[0],double);
        ecoordy=(double *) Calloc(ncoord[0],double);
        for(i=0;i<*ncoordx;i++)
            for(j=0;j<*ncoordy;j++){
                ecoordx[h]=coordx[i];
                ecoordy[h]=coordy[j];
                h++;}
        coordx=ecoordx;
        coordy=ecoordy;}
    
    Range(coordx,rangex,ncoord);// range of the x-coordinate
    Range(coordy,rangey,ncoord);// range of the y-coordinate
    // set the sub-sampling window based on prototype unit window (R_0)
    // and scaling factor (lambda_n)
    deltax=rangex[1]-rangex[0];// R_n = lambda_n * R_0
    deltay=rangey[1]-rangey[0];
   if(!winc[0]){ 
        delta=fmin(deltax,deltay);  // I set this rule if no constant
        winc[0]=sqrt(delta)/4;
        winc[1]=sqrt(delta)/4; }
    else{ 
      if(!winc[1]) winc[1]=winc[0];
    } 
    if(!winstp[0]) winstp[0]=0.5;
    ///Rprintf("%f %f %f %f \n",deltax,deltay,winc[0],winc[1]);
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[1] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows
    numinty=floor((deltay-dimwiny)/winsty+1);   //number of overlapping sub-windows
    
    xgrid=(double *) R_alloc(numintx, sizeof(double));
    ygrid=(double *) R_alloc(numinty, sizeof(double));
    
    SeqStep(rangex,numintx,winstx,xgrid);
    SeqStep(rangey,numinty,winsty,ygrid);
   
        for(h=0;h<nvari;h++) subvari[h]=0;//initialize the  variance on the subwindows
        nsub=0;
        
        for(i=0;i<=numintx;i++){
            for(j=0;j<=numinty;j++){
                
                
                *npts=0;
                for(h=0;h<*npar;h++){ sumgrad[h]=0;gradient[h]=0;} //initialize the gradient of the sub-window
                for(h=0;h<*nparc;h++)gradcor[h]=0;
                SetSampling(coordx,coordy,data,n,npts,nbetas[0],scoordx,scoordy,
                            sdata,xgrid[i]+dimwinx,xgrid[i],ygrid[j]+dimwiny,ygrid[j],sX,X);//  create data and coordinates of the sub-windows
                if(*npts>4){
                    nsub++;
                    nwpair=0;//initialize the number of pairs in the window
                     for(l=0;l<(*npts-1);l++){
                                for(m=(l+1);m<*npts;m++){
                            if(!ISNAN(sdata[l])&&!ISNAN(sdata[m]) ){
                                     lag=dist(type[0],scoordx[l],scoordx[m],scoordy[l],scoordy[m],*REARTH);
                                    if(lag<maxdist[0]){

                              meanl=0.0;meanm=0.0;
                                        for(o=0;o<*nbetas;o++) {
                                            Xl[o]=sX[l][o];
                                            Xm[o]=sX[m][o];
                                            meanl=meanl+sX[l][o]*betas[o];
                                            meanm=meanm+sX[m][o]*betas[o];}
                                            nwpair++;
                         rho=CorFct(cormod,lag,0,parcor,0,0); 
                       //  if(*model==1) GradCorrFct(rho,cormod,eps[0],flagcor,gradcor,lag,0,0,0,parcor);
               
                    switch(*model){
                        case 1: // Gaussian random field:
                            Grad_Pair_Gauss2(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas);
                            break;
                        case 2:  // binary binomial case
                        case 11: 

                            Grad_Pair_Binom(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas);                                           
                        break;
                        case 10://  skewgaussian random field
                            Grad_Pair_Skewgauss(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas);
                        break;
                        case 14:
                        case 16: // negative binomial case
                            Grad_Pair_Binomneg(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas);
                        break;
                      case 12:
                            Grad_Pair_StudenT(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                       break;
                        case 13: //  wrapped gaussian random field
                                  Grad_Pair_Wrapped(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;
                         case 20: // sinh
                                      Grad_Pair_Sinh(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;
                           case 34: // sinh
                                      Grad_Pair_Tukeyh(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;
                        case 21: // gamma
                                      Grad_Pair_Gamma(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;

                         case 30: // gamma
                                      Grad_Pair_Poisson(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;

                        case 22: // loggauss
                                      Grad_Pair_LogGauss(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas);
                        break;
                        case 24: // 
                                      Grad_Pair_LogLogistic(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;
                        case 25: // 
                                      Grad_Pair_Logistic(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;
                            case 26: // 
                                      Grad_Pair_Weibull(rho,cormod,flagnuis,flagcor,gradcor,gradient,lag,0,NN[0],npar,nparc,nparcT,nbetas[0],
                                       nuis,parcor, sdata[l],sdata[m],meanl,meanm,Xl,Xm,sX,l,m,betas); 
                        break;
                    }
                    if(*weigthed) weigths=CorFunBohman(lag,maxdist[0]);
                                        else          weigths=1;
                                        for(h=0;h<*npar;h++) {
                                          if(R_FINITE(gradient[h]))
                                        sumgrad[h]=sumgrad[h]+gradient[h]*weigths;}// sum the gradient in the subwindow
                }}}}
                    h=0;
                    for(p=0;p<*npar;p++)//update the sub-variance in the subwindow
                        for(q=p;q<*npar;q++){
                            subvari[h]=subvari[h]+sumgrad[p]*sumgrad[q]/nwpair;h++;}
                }
        }}
    
    for(h=0;h<nvari;h++)   varimat[h]= np[0] * subvari[h]/(nsub);//update variability matrix
    Free(gradcor);Free(gradient);
    Free(sumgrad);Free(subvari);
    Free(scoordx); Free(scoordy);
    Free(sdata);
    Free(Xl);Free(Xm);
    for(i=0;i<ncoord[0];i++) { Free(X[i]); Free(sX[i]);}
    Free(X);Free(sX);
    if(*grid){Free(ecoordx);Free(ecoordy);}
    return;
}
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/





/********************************************************************************/

void Vari_SubSamp_st2(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, int *dst, double *eps, int *flagcor,
                      int *flagnuis,int *like,double *mean,int *model,double *NN,int *nbetas,int *npar, int *nparc, int *nparcT, double *nuis,int *np, double *parcor,
                      int *type_lik,double *varimat, double *winc, double *winstp,double *winc_t,double *winstp_t,int *weigthed,double *Z,
                      int *ns, int *NS)

{
    double meanl,meanm,beta,rho=0,lags=0.0,lagt=0.0,weigths=1.0, *rangex, *rangey,*gradient;
    double *sdata,*s2data,*xgrid,*ygrid,*scoordx,*scoordy,*sublagt,*subvari, *s2cx,*s2cy;
    double *sumgrad,*gradcor;
    double delta=0.0, deltax=0.0, deltay=0.0, dimwinx=0.0, dimwiny=0.0;
    double winstx=winstp[0], winsty=winstp[0];
    double step=coordt[1]-coordt[0];  // subsampling works in the regular temporal sampling setting
    
    int *npts, numintx=0, numinty=0,nstime=0,*ntimeS,nnc=0;
    int n_win=0,nwpair=0,qq=0,t=0,v=0,h=0,i=0,l=0,m=0,o=0,nsub,nsub_t=0,j=0,f=0,p=0,q=0;
    int nvari=*npar * (*npar+1)/2;

    int NTOT=(NS[ntime[0]-1]+ns[ntime[0]-1]);
    //==========================================/
    
    //==========================================/
    npts=(int *) R_alloc(1, sizeof(int));          //number of location sites in the subwindow
    ntimeS=(int *) R_alloc(1, sizeof(int));          // number of spatial points for each time in the block
    rangex=(double *) R_alloc(2, sizeof(double));
    rangey=(double *) R_alloc(2, sizeof(double));
    
    Range(coordx,rangex,ncoord);// range of the x-coordinate (min and max)
    Range(coordy,rangey,ncoord);// range of the y-coordinate  (min and max)
    // set the sub-sampling window based on prototype unit window (R_0)
    // and scaling factor (lambda_n)
    deltax=rangex[1]-rangex[0];//  x window length
    deltay=rangey[1]-rangey[0];//  y window length
    
    if(!winc[0]){  delta=fmin(deltax,deltay);  // I set this rule if no constant
        winc[0]=sqrt(delta)/4;
        winc[1]=sqrt(delta)/4; }
    else{ 
      if(!winc[1]) winc[1]=winc[0];
    } 
    
    if(!winstp[0]) winstp[0]=1.0;   //proportion of the overlapping  0< winstp <=1  OJO!!!
    
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[1] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
 
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    
  
    

    // what is this ??
        
    // Start conditions for valid space subwindows
    if(cdyn[0] ==0) // no dym coords
    {
        nnc = ncoord[0];
    }
    if(cdyn[0] ==1)
    {
        nnc = ncoord[0]/ntime[0];
    }
    double Cspace1 = (deltax-dimwinx);
    double Cspace2 = (deltay-dimwiny);
    if(Cspace1>0&&Cspace2>0)
    {
          
        numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
        numinty=floor((deltay-dimwiny)/winsty+1);
    }
     if(Cspace1>0&&Cspace2<0)
    {
          
        numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows is  numintx+1 * numinty+1
        numinty=nnc;
    }
       if(Cspace1<0&&Cspace2>0)
    {
          
        numintx=nnc;   //number of overlapping sub-windows is  numintx+1 * numinty+1
        numinty=floor((deltay-dimwiny)/winsty+1);
    }
    if(Cspace1<0&&Cspace2<0)
    {
        numintx = nnc;
        numinty = nnc;
    }
    // End conditions for valid space subwindows
   
    
    ygrid=(double *) R_alloc(numinty, sizeof(double));
    xgrid=(double *) R_alloc(numintx, sizeof(double));
    
    //Rprintf("winc: %f winstp: %f deltax:%f deltay:%f dimwinx:%f dimwiny:%f winstx:%f winsty:%f numintx:%d numinty:%d nnc:%d\n",*winc,*winstp,deltax,deltay,dimwinx,dimwiny,winstx,winsty,numintx,numinty,nnc);
    
    SeqStep(rangex,numintx,winstx,xgrid);
    SeqStep(rangey,numinty,winsty,ygrid);
    
    //Rprintf("winc: %f winstp: %f deltax:%f deltay:%f dimwinx:%f dimwiny:%f winstx:%f winsty:%f",*winc,*winstp,deltax,deltay,dimwinx,dimwiny,winstx,winsty);
    //==========================================/
    //==========================================/
    //default sub window temporal length
    if(!(*winc_t)){
        
        beta=CorFct(cormod,0,1,parcor,0,0);
        // rule given in Genton paper
        *winc_t=R_pow(2*beta/(1-R_pow(beta,2)),(double) 2/3)*R_pow(*ntime*1.5,(double)1/3);
        
        if(*winc_t<4*step) *winc_t=2*step;// if the length is too small
        if(*winc_t>=*ntime) *winc_t=*ntime-step;
        
    } // if the length is too big
     // Rprintf("%f %f\n",winc[0],*winc_t);
    //set the spatial-temporal windows:
    //double wint = (double) *winc_t; // OJO!!!
    int wint = (int) *winc_t; // OJO!!!
    if(!winstp_t) *winstp_t=wint; //defualt for the forward step:the minimum distance
    //else *winstp=*winstp*wint;   //otherwise a proportion of the temporal window
      // Rprintf(" %f %f %f %f %f %f %d--\n",deltax,deltay, dimwinx, dimwiny,winc[0],winc[1],wint);
    sublagt=(double *) R_alloc(wint, sizeof(double));
    sublagt[0]=step;
    //==========================================/
    double **X,*Xl,*Xm,**sX,**s2X;
    Xl=(double *) Calloc(*nbetas,double);
    Xm=(double *) Calloc(*nbetas,double);
    
    
        X= (double **) Calloc(NTOT,double *);
        sX= (double **) Calloc(NTOT,double *);
        for(i=0;i<(NTOT);i++){
            X[i]=(double *) Calloc(nbetas[0],double);
            sX[i]=(double *) Calloc(nbetas[0],double);
         }

        qq=0;
        for(i=0;i<(NTOT);i++){
            for(j=0;j<nbetas[0];j++){
                X[i][j]=Z[qq];
                qq++;
            }}

    //set the temporal distances for the sub-sample:
    for(i=0;i<wint;i++){sublagt[i+1]=sublagt[i]+step;nstime++;}
    
    // Start conditions for valid time subwindows
    double Ctime = (ntime[0]-wint); //condition for time subsampling
     Rprintf("%f %f %f %d\n",winc[0],winc[1],Cspace1,wint);
    if(Ctime>0)
    {
        nsub_t=floor(((ntime[0]-wint)/(winstp[0]*wint)+1));
    }else
    {
        nsub_t = ntime[0];
    }
    // End conditions for valid time subwindows
  
    if( (Cspace1<=0 || Cspace2<=0) && Ctime<=0 ) {return;} // Conditions for valid SPACE && TIME subwindows
    
    n_win=(numintx+1)*(numinty+1);   //number of spatial windows
    nsub=0;
    
    int *ns_sub,*NS_sub;
    int nsub1 =0;
    
    double *res_sub;
    res_sub=(double *) Calloc(NTOT,double);

    Rep(coordt,ns, res_sub);
    gradcor=(double *) Calloc(*nparc,double);
    gradient=(double *) Calloc(*npar,double);
    subvari=(double *) Calloc(nvari,double);
    
    for(i=0;i<=numintx;i++){
        for(j=0;j<=numinty;j++){  // cycle for each block∫∫
            *npts=0;   // number of points in the block
            ns_sub=(int *) Calloc(ntime[0],int);
            NS_sub=(int *) Calloc(ntime[0],int);
          
                scoordx=(double *) Calloc(NTOT,double);
                scoordy=(double *) Calloc(NTOT,double);
                sdata=(double *) Calloc(NTOT ,double);
            
            SetSampling_s(coordx,coordy,data,npts,nbetas[0],scoordx,scoordy,
                          sdata,xgrid[i]+dimwinx,xgrid[i],ygrid[j]+dimwiny,ygrid[j],
                          sX,X,ns,NS,  NS_sub,res_sub,coordt,ns_sub);//  create data and coordinates of the sub-windows
            //  scoordx    scoordy are the spatial coords of the subwindow // sdata the associated data // npts number of loc sites in the window
            if(       //Anche le finestre che hanno "pochi" punti(griglia irregolare)
               (npts[0]>5))   {//OJO
                nsub1++;
                          nwpair=0;//initialize the number of pairs in the window
                for(f=0;f<nsub_t;f++){//loop for the number of tmporal sub-sampling:
                    
                    sumgrad=(double *) Calloc(npar[0],double);
                    s2data=(double *) Calloc(npts[0] ,double);
                    s2cx=(double *) Calloc(npts[0] ,double);
                    s2cy=(double *) Calloc(npts[0] ,double);
                    s2X= (double **) Calloc(npts[0] ,double *);
                    int iii = 0;
                    for(iii=0;iii<npts[0];iii++) 
                       {s2X[iii]=(double *) Calloc(nbetas[0],double);}
                    // set the sub-sample of the data:
                    
                    *ntimeS=0;// number of spatial points for each time in the block
                    SetSampling_t(sdata,s2data,nbetas[0],npts[0],nstime,wint,f,s2X,sX,ns_sub,NS_sub,nsub_t,ntimeS,s2cx,s2cy,scoordx,scoordy);
                    
                    //======================================/
                    //computing gradient in the window/
                    //======================================/


                    if(ntimeS[0]>5)
                    {
                        for(t=0;t<nstime;t++){
                            for(l=0;l<ns_sub[t+f*nstime];l++){       // what is ns_sub[t]
                                for(v=t;v<nstime;v++){
                                    if(t==v){
                                        for(m=l+1;m<ns_sub[t+f*nstime];m++){      // what is ns_sub[t]
              lags=dist(type[0],s2cx[(l+NS_sub[t])],s2cx[(m+NS_sub[v])],s2cy[(l+NS_sub[t])],s2cy[(m+NS_sub[v])],*REARTH);
            
              // lags=dist(type[0],s2cx[l],s2cx[m],s2cy[l],s2cy[m],*REARTH);

                                            if(lags<=maxdist[0]){
                                                meanl=0.0;meanm=0.0;
                                                for(o=0;o<nbetas[0];o++) {
                                                    Xl[o]=s2X[l+NS_sub[t]][o];   // what is ns_sub[t].
                                                    Xm[o]=s2X[m+NS_sub[v]][o];   // what is ns_sub[t]
                                                   // Rprintf("%f, %f\n",Xl[o],Xm[o]);
                                                    meanl=meanl+Xl[o]*betas[o];
                                                    meanm=meanm+Xm[o]*betas[o];
                                                }
                           // Rprintf("%d --%f %f %f \n",f,lags,s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])]);
                                                if(!ISNAN(s2data[(l+NS_sub[t])])&&!ISNAN(s2data[(m+NS_sub[v])]) ){
                                                    nwpair++;
                                                    rho=CorFct(cormod,lags,0,parcor,0,0); //
                                                    
                                                    //===================   starting cases ================/
                                                    
                                                   // if(*model==1) GradCorrFct(rho,cormod,eps[0],flagcor,gradcor,lags,0,0,0,parcor);
                                                    
                                                    switch(*model){
                                                        case 1:   // gaussian case
                                                      //  Grad_Pair_Gauss(rho,flagnuis,flagcor,gradcor,gradient,npar,nparc,nbetas[0],nuis, //gradient is the output
                                                       //                 s2data[l]-meanl,
                                                         //               s2data[m]-meanm,Xl,Xm);
                                                         Grad_Pair_Gauss2(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                        s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                        betas);
                                                        break;
                                                        case 2:     //binomial  gaussian case
                                                        case 11:
                                                        Grad_Pair_Binom(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 12:
                                                        Grad_Pair_StudenT(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 13:
                                                       Grad_Pair_Wrapped(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 14:
                                                        case 16:
                                                   Grad_Pair_Binomneg(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 20:
                                                        Grad_Pair_Sinh(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                          case 34:
                                                        Grad_Pair_Tukeyh(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 21:
                                                
                                                         Grad_Pair_Gamma(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;

                                                              case 30:
                                                
                                                         Grad_Pair_Poisson(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;

                                                        case 22:
                                                       Grad_Pair_LogGauss(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 24:
                                                       Grad_Pair_LogLogistic(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 25:
                                                        Grad_Pair_Logistic(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                        break;
                                                        case 26:
                                                        Grad_Pair_Weibull(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,0,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                          s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],
                                                                          betas);
                                                         //Grad_Pair_Weibull(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                           //             s2data[l],s2data[m],meanl,meanm,Xl,Xm,s2X,l,m,betas);
                                                        break;
                                                    }
                                                    //==============   end cases ============/
                                                    if(*weigthed) weigths=CorFunBohman(lagt,maxtime[0]);
                                                    else          weigths=1;
                                                    for(h=0;h<*npar;h++) {
                                                        if(R_FINITE(gradient[h]))
                                                        sumgrad[h]=sumgrad[h]+gradient[h]*weigths;}}}}
                                    }  //================================================/
                                    else{
                                        lagt=fabs(sublagt[t]-sublagt[v]);
                                        for(m=0;m<ns_sub[t+f*nstime];m++){
                                            
                                            //lags=dist(type[0],s2cx[l],s2cx[m],s2cy[l],s2cy[m],*REARTH);
                                             lags=dist(type[0],s2cx[(l+NS_sub[t])],s2cx[(m+NS_sub[v])],s2cy[(l+NS_sub[t])],s2cy[(m+NS_sub[v])],*REARTH);
                                           //    Rprintf("%d--%f %f  %f %f \n",f,lags,lagt,s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])]);
                                            if(lagt<=maxtime[0] && lags<=maxdist[0]){
                                                
                                                meanl=0.0;meanm=0.0;
                                                for(o=0;o<nbetas[0];o++) {
                                                    Xl[o]=s2X[l+NS_sub[t]][o];   
                                                    Xm[o]=s2X[m+NS_sub[v]][o];   
                                                    meanl=meanl+Xl[o]*betas[o];
                                                    meanm=meanm+Xm[o]*betas[o];
                                                }
                                                
                                                if(!ISNAN(s2data[(l+NS_sub[t])])&&!ISNAN(s2data[(m+NS_sub[v])]) ){
                                                    nwpair++;
                                                    rho=CorFct(cormod,lags,lagt,parcor,0,0);
                                                    
                                                    //===================   starting cases ================/
                                                    //if(*model==1)   GradCorrFct(rho,cormod,eps[0],flagcor,gradcor,lags,lagt,0,0,parcor);
                                                    switch(*model){
                                                        case 1:   // gaussian case
                                                         Grad_Pair_Gauss2(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                        s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 2:     //binomial  gaussian case
                                                        case 11:
                                                        Grad_Pair_Binom(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                             case 12:
                                                        Grad_Pair_StudenT(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                          break;
                                                        case 13:
                                                        Grad_Pair_Wrapped(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 14:
                                                        case 16:
                                                        Grad_Pair_Binomneg(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 20:
                                                        Grad_Pair_Sinh(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                         case 34:
                                                        Grad_Pair_Tukeyh(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 21:
                                                             Grad_Pair_Gamma(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 30:
                                                             Grad_Pair_Poisson(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 22:
                                                        Grad_Pair_LogGauss(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                  
                                                           break;
                                                        case 24:
                                                        Grad_Pair_LogLogistic(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 25:
                                                        Grad_Pair_Logistic(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                        case 26:
                                                        Grad_Pair_Weibull(rho,cormod,flagnuis,flagcor,gradcor,gradient,lags,lagt,NN[0],npar,nparc,nparcT,nbetas[0],nuis,parcor,
                                                                         s2data[(l+NS_sub[t])],s2data[(m+NS_sub[v])],meanl,meanm,Xl,Xm,s2X,l+NS_sub[t],m+NS_sub[v],betas);
                                                        break;
                                                    }
                                                    //==============   end cases ============/
                                                    if(*weigthed) weigths=CorFunBohman(lagt,maxtime[0]);
                                                    else          weigths=1;
                                                    for(h=0;h<*npar;h++) {
                                                           if(R_FINITE(gradient[h]))
                                                        sumgrad[h]=sumgrad[h]+gradient[h]*weigths;
                                                        
                                                    }}}}}
                                }}}//END t loop
                        h=0;//update the sub-variance in the subwindow
                        for(p=0;p<*npar;p++){
                               for(q=p;q<*npar;q++){subvari[h]=subvari[h]+sumgrad[p]*sumgrad[q]/nwpair;h++;}}
                        
                        Free(sumgrad); Free(s2data);Free(s2cx);Free(s2cy);
                        int sss =0;
                        for(sss=0;sss<(npts[0] );sss++) {Free(s2X[sss]);}
                        Free(s2X);
                        nsub++;}//END f loop
                }
            }  //END ID xgrid[i]
            Free(scoordx); Free(scoordy);Free(sdata);
        }//END j loop
    }//END i loop
    
    Free(ns_sub);Free(NS_sub);Free(res_sub);
    for(h=0;h<nvari;h++)  varimat[h]= np[0]*subvari[h]/(nsub); //update variability matrix
    
    
    Free(gradcor);
    Free(gradient);
    Free(subvari);
    Free(Xl);
    Free(Xm);
    for(i=0;i<NTOT;i++) { Free(X[i]); Free(sX[i]);}
    Free(X);Free(sX);
  
    return;
}



// Computes the variability matrix based on the sub-sampling method (bivariate):
void Vari_SubSamp_biv(double *betas,double *coordx, double *coordy, double *coordt,int *cormod, double *data,
                  int *dst,double *eps, int *flagcor, int *flagnuis, int *grid,
                  int *like, int *model, double *NN,int *npar, int *nparc,int *nparcT, double *nuis, int *np,
                  double *parcor, int *type_lik, double *varimat,
                  double *winc, double *winstp,int *weigthed,double *Z,
                      int *ns, int *NS)

{
    double dd, lag=0.0,*gradient, *rangex, *rangey,rhott=0.0,rhovv=0.0,rhotv=0.0,rhovt=0.0;
    double rhottii=0.0,rhovvii=0.0, rhotvij=0.0, rhovtij=0.0, rhottij=0.0, rhovvij=0.0, rhotvii=0.0,rhovtii=0.0;
    double  *gradcortt,*gradcortv,*gradcorvt,*gradcorvv,*ecoordx,*ecoordy,*scoordx,*scoordy,*sdata,*score,*subvari,*xgrid,*ygrid;
    double delta=0.0, deltax=0.0, deltay=0.0, dimwinx=0.0, dimwiny=0.0,weigths;
    double winstx=winstp[0], winsty=winstp[0];
    int *npts, numintx=0, numinty=0;
    int N=4,h=0,i=0,l=0,m=0,mm=0,n=0,nsub=0,nvari=0,nwpair=0,p=0,q=0,j=0,t=0,v=0,k=0,kk=0;

    nvari=*npar * (*npar+1)/2;

    gradcortt=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcortv=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvt=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvv=(double *) Calloc(*nparc,double);// Correlation gradient
    gradient=(double *) Calloc(*npar,double);
    score=(double *) Calloc(*npar,double);
    subvari=(double *) Calloc(nvari,double);

    double *gradcorttii,*gradcorvvii ,*gradcorvtii  ,*gradcortvii ,*gradcorttij ,*gradcorvvij ,*gradcorvtij,*gradcortvij;
    gradcorttii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvvii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvtii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcortvii=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorttij=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvvij=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcorvtij=(double *) Calloc(*nparc,double);// Correlation gradient
    gradcortvij=(double *) Calloc(*nparc,double);// Correlation gradient
    
    int *indx;
    indx = (int *) R_alloc(N, sizeof(int));
    double *col;
    col =(double *) R_alloc(N, sizeof(double));
    double **M;
    M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}
    
    double **inverse;
    inverse= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){inverse[i]=(double *) Calloc(N,double);}
    double *dat; //dataÁÁÁ
    dat=(double *) R_alloc(N, sizeof(double));

    npts=(int *) R_alloc(1, sizeof(int));

    rangex=(double *) R_alloc(2, sizeof(double));
    rangey=(double *) R_alloc(2, sizeof(double));

    
    scoordx=(double *) Calloc(ncoord[0],double);
    scoordy=(double *) Calloc(ncoord[0],double);
    sdata=(double *) Calloc(ncoord[0],double);
    
    if(*grid){
        ecoordx=(double *) Calloc(ncoord[0],double);
        ecoordy=(double *) Calloc(ncoord[0],double);
        for(i=0;i<*ncoordx;i++)
            for(j=0;j<*ncoordy;j++){
                ecoordx[h]=coordx[i];
                ecoordy[h]=coordy[j];
                h++;}
        coordx=ecoordx;
        coordy=ecoordy;}
    
    
    
    
    
    Range(coordx,rangex,ncoord);// range of the x-coordinate
    Range(coordy,rangey,ncoord);// range of the y-coordinate
    // set the sub-sampling window based on prototype unit window (R_0)
    // and scaling factor (lambda_n)
    deltax=rangex[1]-rangex[0];// R_n = lambda_n * R_0
    deltay=rangey[1]-rangey[0];
    if(fabs(winc[0])<LOW || !winc[0]){
        delta=fmin(deltax,deltay);
        winc[0]=(delta/sqrt(delta))/2;}
    if(fabs(winstp[0])<LOW || !winstp[0]) winstp[0]=0.5;
    dimwinx=winc[0] * sqrt(deltax);// sub-window x length depends on a constant: deafault??
    dimwiny=winc[0] * sqrt(deltay);// sub-window y length depends on a constant: deafault??
    winstx=*winstp * dimwinx;     // x step is a  proportion of sub-window x length (deafult is 0.5)
    winsty=*winstp * dimwiny;     // y step is a  proportion of sub-window y length (deafult is 0.5)
    numintx=floor((deltax-dimwinx)/winstx+1);   //number of overlapping sub-windows
    numinty=floor((deltay-dimwiny)/winsty+1);   //number of overlapping sub-windows

    


    xgrid=(double *) R_alloc(numintx, sizeof(double));
    ygrid=(double *) R_alloc(numinty, sizeof(double));

    SeqStep(rangex,numintx,winstx,xgrid);
    SeqStep(rangey,numinty,winsty,ygrid);

    for(kk=0;kk<*nrep;kk++){
        for(h=0;h<nvari;h++) subvari[h]=0;//initialize the  variance on the subwindows
        nsub=0;

        for(mm=0;mm<=numintx;mm++){
            for(l=0;l<=numinty;l++){
                *npts=0;
                for(h=0;h<*npar;h++){ score[h]=0;gradient[h]=0;} //initialize the gradient of the sub-window
                for(h=0;h<*nparc;h++) {gradcortt[h]=0;gradcorvv[h]=0;gradcorvt[h]=0;gradcortv[h]=0;}
                SetSampling_biv(coordx,coordy,data,kk,npts,scoordx,scoordy,
                            sdata,xgrid[mm]+dimwinx,xgrid[mm],ygrid[l]+dimwiny,ygrid[l]);//  setting data and coordinates of the sub-windows
                
                if(*npts>5){
                    nsub++;nwpair=0;//initialize the number of pairs in the window
                    switch(*model){
                        case 1: // Gaussian random field:
                                 switch(*like){//select the type of composite likelihood
                                 case 1:// Conditional likelihood (bipairwise)
                                     switch(*type_lik){
                                         case 2:   //Marginal pair of pairwise likelihood:
                                         for(t=0;t<*ntime-1;t++){
                                             for(v=t+1;v<*ntime;v++){
                                                 for(i=0;i<*npts-1;i++){
                                                     for(j=i+1;j<*npts;j++){
                                                                 lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                                              if(lag<=dista[t][v]){
                                                             nwpair++;
                                                             rhottii=CorFct(cormod,0,0,parcor,t,t);rhovvii=CorFct(cormod,0,0,parcor,v,v);
                                                             rhottij=CorFct(cormod,lag,0,parcor,t,t);rhovvij=CorFct(cormod,lag,0,parcor,v,v);
                                                             rhotvii=CorFct(cormod,0,0,parcor,t,v);rhovtii=CorFct(cormod,0,0,parcor,v,t);
                                                             rhotvij=CorFct(cormod,lag,0,parcor,t,v);rhovtij=CorFct(cormod,lag,0,parcor,v,t);
                                                             M[0][0]=rhottii;           M[0][1]=rhottij;       M[0][2]=rhotvii;        M[0][3]=rhotvij;
                                                             M[1][0]=M[0][1];           M[1][1]=rhottii;       M[1][2]=rhovtij;        M[1][3]=rhovtii;
                                                             M[2][0]=M[0][2];           M[2][1]= M[1][2];      M[2][2]=rhovvii;        M[2][3]=rhovvij;
                                                             M[3][0]=M[0][3];           M[3][1]= M[1][3];      M[3][2]=M[2][3];        M[3][3]=rhovvii;
                                                             ludcmp(M,N,indx,&dd);    //Lu decomposition
                                                             for(n=0;n<N;n++) {
                                                                 for(m=0;m<N;m++) col[m]=0.0;
                                                                 col[n]=1.0;
                                                                 lubksb(M,N,indx,col);
                                                                 for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                                                             }
                                                             GradCorrFct(0,cormod,eps[0],flagcor,gradcorttii,0,0,t,t,parcor);    GradCorrFct(0,cormod,eps[0],flagcor,gradcorvvii,0,0,v,v,parcor);
                                                             GradCorrFct(0,cormod,eps[0],flagcor,gradcorttij,lag,0,t,t,parcor); GradCorrFct(0,cormod,eps[0],flagcor,gradcorvvij,lag,0,v,v,parcor);
                                                             GradCorrFct(0,cormod,eps[0],flagcor,gradcortvii,0,0,t,v,parcor);  GradCorrFct(0,cormod,eps[0],flagcor,gradcorvtii,0,0,v,t,parcor);
                                                             GradCorrFct(0,cormod,eps[0],flagcor,gradcortvij,lag,0,t,v,parcor);  GradCorrFct(0,cormod,eps[0],flagcor,gradcorvtij,lag,0,v,t,parcor);
                                                             dat[0]=sdata[(t+*ntime*i)+kk* *nrep];dat[1]=sdata[(t+*ntime*j)+kk* *nrep];
                                                             dat[2]=sdata[(t+*ntime*i+1)+kk* *nrep];dat[3]=sdata[(t+*ntime*j+1)+kk* *nrep];
                                                             //Compute the gradient of the log pairwise likelihood
                                                             Grad_Cond_Gauss_biv(gradcorttii,gradcorvvii,gradcorvtii,gradcortvii ,gradcorttij ,gradcorvvij, gradcorvtij,gradcortvij,inverse,flagnuis,gradient,npar,nuis,dat,N);
                                                             if(*weigthed) weigths=CorFunBohman(lag,dista[t][v]);
                                                             else          weigths=1;
                                                             for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                                         }}}}}
                                         break;
                                     }
                                     break;
                                  case 3: //Marginal pairwise likelihood:
                                     switch(*type_lik){
                                     case 2:
                                         for(t=0;t<*ntime;t++){
                                             for(v=0;v<*ntime;v++){
                                                 if(t==v||t>v){
                                                     for(i=0;i<*npts-1;i++){
                                                         for(j=i+1;j<*npts;j++){
                                                                     lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                                               if(lag<=dista[t][v]){
                                                                 nwpair++;
                                                                 // Compute the  correlation function for for a given pair:
                                                                 rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                                 rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                                 // Compute the gradient for the temporal correlation function:
                                                                 GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                                 GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                                 //Compute the gradient of the log pairwise likelihood
                                                                 Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,
                                                                                     gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                                     sdata[(t+*ntime*i)+kk* *nrep],sdata[(v+*ntime*j)+kk* *nrep]);
                                                                  weigths=1;
                                                                 if(*weigthed&&t>v) weigths=CorFunBohman(lag,dista[t][v]);
                                                                 for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                                             }}}}
                                                 else {
                                                     for(i=0;i<*npts;i++){
                                                         for(j=i;j<*npts;j++){
                                                                  lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                                                  if(lag<=dista[t][v]){
                                                                 nwpair++;
                                                                 // Compute the l correlation function for a given pair:
                                                                 rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                                 rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                                 // Compute the gradient of the patial-temporal correlation function:
                                                                 GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                                 GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                                 //Compute the gradient of the log pairwise likelihood
                                                                 Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,
                                                                                     gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                                     sdata[(t+*ntime*i)+kk* *nrep],
                                                                                     sdata[(v+*ntime*j)+kk* *nrep]);
                                                                 if(*weigthed) weigths=CorFunBohman(lag,dista[t][v]);
                                                                 else          weigths=1;
                                                                 for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                                             }}}}}}
                                         break;}
                                         break;}
                        break;
                        case 2: // Binary gaussian
                        for(t=0;t<*ntime;t++){
                            for(v=0;v<*ntime;v++){
                                if(t==v||t>v){
                                    for(i=0;i<*npts-1;i++){
                                        for(j=i+1;j<*npts;j++){
                                                lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the  correlation function for for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient for the temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],sdata[(v+*ntime*j)+kk* *nrep]);
                                                weigths=1;
                                                if(*weigthed&&t>v) weigths=CorFunBohman(lag,dista[t][v]);
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}
                                else {
                                    for(i=0;i<*npts;i++){
                                        for(j=i;j<*npts;j++){
                                                  lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the l correlation function for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient of the patial-temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],
                                                                    sdata[(v+*ntime*j)+kk* *nrep]);
                                                if(*weigthed) weigths=CorFunBohman(lag,dista[t][v]);
                                                else          weigths=1;
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}}}
                        break;
                        case 10:  // skew gaussian
                        for(t=0;t<*ntime;t++){
                            for(v=0;v<*ntime;v++){
                                if(t==v||t>v){
                                    for(i=0;i<*npts-1;i++){
                                        for(j=i+1;j<*npts;j++){
                                                 lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the  correlation function for for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient for the temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],sdata[(v+*ntime*j)+kk* *nrep]);
                                                weigths=1;
                                                if(*weigthed&&t>v) weigths=CorFunBohman(lag,dista[t][v]);
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}
                                else {
                                    for(i=0;i<*npts;i++){
                                        for(j=i;j<*npts;j++){
                                                 lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the l correlation function for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient of the patial-temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],
                                                                    sdata[(v+*ntime*j)+kk* *nrep]);
                                                if(*weigthed) weigths=CorFunBohman(lag,dista[t][v]);
                                                else          weigths=1;
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}}}
                        break;
                        case 11:  // binomial gaussian
                        for(t=0;t<*ntime;t++){
                            for(v=0;v<*ntime;v++){
                                if(t==v||t>v){
                                    for(i=0;i<*npts-1;i++){
                                        for(j=i+1;j<*npts;j++){
                                                 lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the  correlation function for for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient for the temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],sdata[(v+*ntime*j)+kk* *nrep]);
                                                weigths=1;
                                                if(*weigthed&&t>v) weigths=CorFunBohman(lag,dista[t][v]);
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}
                                else {
                                    for(i=0;i<*npts;i++){
                                        for(j=i;j<*npts;j++){
                                             lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the l correlation function for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient of the patial-temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],
                                                                    sdata[(v+*ntime*j)+kk* *nrep]);
                                                if(*weigthed) weigths=CorFunBohman(lag,dista[t][v]);
                                                else          weigths=1;
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}}}
                        break;
                        case 12:  // chis gaussian
                        for(t=0;t<*ntime;t++){
                            for(v=0;v<*ntime;v++){
                                if(t==v||t>v){
                                    for(i=0;i<*npts-1;i++){
                                        for(j=i+1;j<*npts;j++){
                                                  lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the  correlation function for for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient for the temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],sdata[(v+*ntime*j)+kk* *nrep]);
                                                weigths=1;
                                                if(*weigthed&&t>v) weigths=CorFunBohman(lag,dista[t][v]);
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}
                                else {
                                    for(i=0;i<*npts;i++){
                                        for(j=i;j<*npts;j++){
                                                  lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the l correlation function for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient of the patial-temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],
                                                                    sdata[(v+*ntime*j)+kk* *nrep]);
                                                if(*weigthed) weigths=CorFunBohman(lag,dista[t][v]);
                                                else          weigths=1;
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}}}
                        break;
                        case 13:  // wrap gaussian
                        for(t=0;t<*ntime;t++){
                            for(v=0;v<*ntime;v++){
                                if(t==v||t>v){
                                    for(i=0;i<*npts-1;i++){
                                        for(j=i+1;j<*npts;j++){
                                                  lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the  correlation function for for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient for the temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],sdata[(v+*ntime*j)+kk* *nrep]);
                                                weigths=1;
                                                if(*weigthed&&t>v) weigths=CorFunBohman(lag,dista[t][v]);
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}
                                else {
                                    for(i=0;i<*npts;i++){
                                        for(j=i;j<*npts;j++){
                                                lag=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                                            if(lag<=dista[t][v]){
                                                nwpair++;
                                                // Compute the l correlation function for a given pair:
                                                rhott=CorFct(cormod,0,0,parcor,t,t);rhovv=CorFct(cormod,0,0,parcor,v,v);
                                                rhotv=CorFct(cormod,lag,0,parcor,t,v);rhovt=rhotv;
                                                // Compute the gradient of the patial-temporal correlation function:
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcortt,0,0,t,t,parcor);GradCorrFct(0,cormod,eps[0],flagcor,gradcortv,lag,0,t,v,parcor);
                                                GradCorrFct(0,cormod,eps[0],flagcor,gradcorvv,0,0,v,v,parcor);//GradCorrFct(0,cormod,eps[0],flagcor,gradcorvt,lag,0,t,v,parcor);
                                                //Compute the gradient of the log pairwise likelihood
                                                Grad_Pair_Gauss_biv(rhott,rhotv,rhovt,rhovv,flagnuis,           ////// TO CHANGE HERE
                                                                    gradcortt,gradcortv,gradcortv,gradcorvv,gradient,npar,nuis,
                                                                    sdata[(t+*ntime*i)+kk* *nrep],
                                                                    sdata[(v+*ntime*j)+kk* *nrep]);
                                                if(*weigthed) weigths=CorFunBohman(lag,dista[t][v]);
                                                else          weigths=1;
                                                for(k=0;k<*npar;k++) score[k]=score[k]+gradient[k]*weigths;
                                            }}}}}}
                        break;
                    }
                    h=0;
                    for(p=0;p<*npar;p++)//update the sub-variance in the subwindow
                        for(q=p;q<*npar;q++){subvari[h]=subvari[h]+score[p]*score[q]/nwpair;h++;}
                }
            
            
            
            }}}
    //Rprintf("%f %f  %f %f %d  %d %d\n",subvari[0],subvari[1],subvari[2],subvari[3],nsub ,*nrep,nwpair);
    for(h=0;h<nvari;h++) varimat[h]= np[0] * subvari[h]/(nsub* *nrep);//update variability matrix
    
    Free(gradcortt);// Correlation gradient
    Free(gradcortv);// Correlation gradient
    Free(gradcorvt);// Correlation gradient
    Free(gradcorvv);// Correlation gradient
    Free(gradcorttii);// Correlation gradient
    Free(gradcorvvii);// Correlation gradient
    Free(gradcorvtii);// Correlation gradient
    Free(gradcortvii);// Correlation gradient
    Free(gradcorttij);// Correlation gradient
    Free(gradcorvvij);// Correlation gradient
    Free(gradcorvtij);// Correlation gradient
    Free(gradcortvij);// Correlation gradient
    Free(gradient);
    Free(score);
    Free(subvari);
    for(i=0;i<N;i++)  {Free (inverse[i]);Free(M[i]);}
    Free(M);Free(inverse);
    
    Free(scoordx);
    Free(scoordy);
    Free(sdata);
    if(*grid){
        Free(ecoordx);
        Free(ecoordy);
    }
    
    return;
}
