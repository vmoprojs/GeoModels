#include "header.h"

      
// check the validity of the parameters' range:
double CheckCor(int *cormod, double *par)
{
  double col=0.0, R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0.0, R_power_t=0.0, var11=0.0, var22=0.0;
  double rho=0.0, sep=0, scale=0.0, smooth=0.0,smooth_s=0.0,smooth_t=0.0, scale_s=0.0, scale_t=0;
  double scale11=0.0, scale22=0, scale12=0.0, smoo11=0.0, smoo12=0.0, smoo22=0,
    R_power11=0.0, R_power12=0.0, R_power22=0,nug11=0.0,nug22=0.0;
 

  switch(*cormod) // Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      R_power2=par[0];
      scale=par[1];
      if(scale<=0 || R_power2<=0) rho=-2;
     break;
    case 2:
    case 3:
    case 4:// Exponential correlation function
    //case 6:// Gaussian correlation function
    case 10:// skarofsky
    case 16://wave correlation function
      scale=par[0];
      if(scale<=0) rho=-2;
      break;
    case 8: // Generalised Cuachy correlation function
    case 5: // Dagum
      R_power1=par[0];
      R_power2=par[1];
      scale=par[2];
      if(scale<=0 || R_power1<=0 || R_power1>2 || R_power2<=0) rho=-2;
      break;
    case 12:// Stable correlation function
    case 17: //multiquadric sphere
      R_power=par[0];
      scale=par[1];
      if(scale<=0 || R_power<0 || R_power>2) rho=-2;
       break;
    case 11:// wen0 correlation function
        R_power=par[0];
        scale=par[1];
        if(scale<=0 || R_power<1.5) rho=-2;
       break;
    case 13:// wen1 correlation function
        R_power=par[0];
        scale=par[1];
        if(scale<=0 || R_power<2.5) rho=-2;
      break;
    case 15:// wen1 correlation function
        R_power=par[0];
        scale=par[1];
        if(scale<=0 || R_power<3.5) rho=-2;
        break; 
    case 19: // Generalised wendland
    case 6:
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
       // if(scale<=0 ||  R_power1<(1.5+smooth) ||smooth<0) rho=-2;
            if(scale<=0 ||smooth<0) rho=-2;
      break;
    case 18://sinR_power valid on sphere
            R_power=par[0];
            if( R_power<0 ) rho=-2;
            break;        
    case 14://  Whittle-Matern correlation function
    case 20:
      scale=par[0];
      smooth=par[1];
      if(scale<=0 || smooth<=0) rho=-2;
      break;
      // START non-separable correlation functions:
    case 42: //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
    case 46://Porcu model as in (4) of Porcu Bevilacqua et al (2010), with beta_1=beta_2=1
    case 50:
    case 60:
    case 52:
    case 54:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_t>2  || sep<0 || sep>1) rho=-2;  //||
      break;
             case 61:
        R_power_s=par[0];
        scale_s=par[2];
        scale_t=par[3];
        R_power=par[1];
        sep=par[4];
        smooth=par[5];
        if(scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 ||  sep<0 || sep>1 || smooth<=0||R_power<1) rho=-2;  //||
        break;
        case 62:
         R_power_t=par[0];
        scale_s=par[2];
        scale_t=par[3];
        R_power=par[1];
        sep=par[4];
        smooth=par[5];
        if(scale_t<=0 || scale_s<=0 || R_power_t<0 || R_power_t>2 ||  sep<0 || sep>1 || smooth<=0||R_power<1) rho=-2;  //||
        break;
  case 44:// Iaco-Cesare model as in (14) of Gneitint (2006): note that R_power parameters are in [0,1]
      R_power2=par[0];
      R_power_s=par[1];
      R_power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      if(R_power2<=0||scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_t>2) rho=-2;
      break;
    case 48:// Stein model as in (16) of JASA (2005) with epsilon=0*/
      R_power_t=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(scale_s<=0 || scale_t<=0 || R_power_t<0 || R_power_t>2 || smooth<=0) rho=-2;
      break;
    case 58:  //only sphere
    case 56:    //only sphere1
        R_power_s=par[0];
        R_power_t=par[1];
        scale_s=par[2];
        scale_t=par[3];
        if(scale_s<=0  || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_s>2) rho=-2;
        break;
    case 63:  // non separable temporal wendloand
    case 65:
    case 67:
        R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
  if(*cormod==63&&(scale_s<=0||scale_t<=0 || R_power<2.5  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1|| R_power_s<3.5)) rho=-2;
  if(*cormod==65&&(scale_s<=0||scale_t<=0 || R_power<4.5  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1|| R_power_s<4.5)) rho=-2;
  if(*cormod==67&&(scale_s<=0||scale_t<=0 || R_power<6.5  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1||R_power_s<5.5)) rho=-2;
    break;
        case 87:
        R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
  if(scale_s<=0||scale_t<=0 || R_power<(2.5+2*smooth)  || R_power_t>2|| R_power_t<0||sep<0 ||sep >1|| R_power_s<(3.5+smooth)||smooth<0) rho=-2;
      break;
   case 64:  // non separable temporal wendloand
    case 66:
    case 68:
       R_power_s=par[0];
        R_power_t=par[2];
        R_power=par[1];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
   if(*cormod==64&&(scale_s<=0||scale_t<=0 || R_power<2.5  || R_power_s>2|| R_power_s<0||sep<0 ||sep >1)) rho=-2;
  if(*cormod==66&&(scale_s<=0||scale_t<=0 || R_power<3.5  || R_power_s>2|| R_power_s<0||sep<0 ||sep >1)) rho=-2;
  if(*cormod==68&&(scale_s<=0||scale_t<=0 || R_power<4.5  || R_power_s>2|| R_power_s<0||sep<0 ||sep >1)) rho=-2;
    break;
        case 88:
       R_power_s=par[0];
        R_power_t=par[2];
        R_power=par[1];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
   if(scale_s<=0||scale_t<=0 || R_power<(2.5+2*smooth) || R_power_s>2|| R_power_s<0||sep<0 ||sep >1|| R_power_s<(3.5+smooth)||smooth<0)  rho=-2;
    break;
    case 69:  // separable wendlands
    case 70:
    case 71:
    case 72:
    case 73:
    case 74:
    case 75:
    case 76:
    case 77: 
          R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
        if(*cormod==69 &&(scale_s<=0  || scale_t<=0  || R_power_s<2.5||R_power_t<2.5)) rho=-2; 
        if(*cormod==70 &&(scale_s<=0  || scale_t<=0  || R_power_s<2.5||R_power_t<3.5)) rho=-2;
        if(*cormod==71 &&(scale_s<=0  || scale_t<=0  || R_power_s<2.5||R_power_t<4.5)) rho=-2;
        if(*cormod==72 &&(scale_s<=0  || scale_t<=0  || R_power_s<3.5||R_power_t<2.5)) rho=-2;
        if(*cormod==72 &&(scale_s<=0  || scale_t<=0  || R_power_s<3.5||R_power_t<3.5)) rho=-2;
        if(*cormod==74 &&(scale_s<=0  || scale_t<=0  || R_power_s<3.5||R_power_t<4.5)) rho=-2;
        if(*cormod==75 &&(scale_s<=0  || scale_t<=0  || R_power_s<4.5||R_power_t<2.5)) rho=-2;
        if(*cormod==76 &&(scale_s<=0  || scale_t<=0  || R_power_s<4.5||R_power_t<3.5)) rho=-2;
        if(*cormod==77 &&(scale_s<=0  || scale_t<=0  || R_power_s<4.5||R_power_t<4.5)) rho=-2;
        break;          
      // END non-separable correlation functions
      // START separable correlation functions:
    case 82:// Exp-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      if(scale_s<=0 || scale_t<=0 || R_power2<=0) rho=-2;
      break;
    case 84:// Double exp:
    /*case 88:// Exp-cos:
      scale_s=par[0];
      scale_t=par[1];
      if(scale_s<=0 || scale_t<=0) rho=-2;
      break;*/
    case 86:
      scale_s=par[0];
      scale_t=par[1];
      smooth_s=par[2];
      smooth_t=par[3];
      if(scale_s<=0 || scale_t<=0|| smooth_t<=0|| smooth_t<=0) rho=-2;
    break;
    case 90:// Matern-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(scale_s<=0 || scale_t<=0 || R_power2<=0 || smooth<=0) rho=-2;
      break;
    case 92:// Matern-exp:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(scale_s<=0 || scale_t<=0 || smooth<=0) rho=-2;
      break;
    case 94:// Stable-stable:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      if(scale_s<=0 || scale_t<=0 || R_power_s<0 || R_power_s>2 || R_power_t<0 || R_power_s>2) rho=-2;
      break;
      // END separable correlation functions:
    case 111:  //wend sep  k=0
    case 113:   //wend sep  k=1
    case 115:   //wend sep  k=2
       var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power=par[5];
        scale=par[6];
     if(*cormod==111&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 || 
        R_power<2.5||fabs(col)>=1)) { rho=-2;break; }
      if(*cormod==113&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 || 
        R_power<3.5||fabs(col)>=1)) { rho=-2;break; }
     if(*cormod==115&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 || 
        R_power<4.5||fabs(col)>=1)) { rho=-2;break; }
     break;
    case 129:  //wend  contr sep  k=0
    case 131:   //wend contr sep  k=1
    case 120:   //wend contr sep  k=2
       var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     R_power11=par[5];
     R_power22=par[6];
     scale11=par[7];
     scale22=par[8];
      if(*cormod==129&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
         scale22<=0 ||  R_power11<2.5||  R_power22<2.5||fabs(col)>= R_pow( R_pow(0.5*(scale11+scale22),2)/(scale11*scale22) ,2.5) )) { rho=-2;break; }
        if(*cormod==131&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale22<=0 || R_power11<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
        if(*cormod==120&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale22<=0 ||   R_power11<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
    break; 

    case 112:  //wend full  k=0
    case 114:   //wend full  k=1
    case 116:   //wend full  k=2

    var11=par[0];
    var22=par[1];
    nug11=par[2];
    nug22=par[3];
    col=par[4];
    R_power11=par[5];
    R_power12=par[6];
    R_power22=par[7];
    scale11=par[8];
    scale12=par[9];   
    scale22=par[10];
        if(*cormod==112&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale12<=0 || scale22<=0 ||  
        R_power11<2.5|| R_power12<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
        if(*cormod==114&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale12<=0 || scale22<=0 ||  
        R_power11<2.5|| R_power12<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
        if(*cormod==116&&(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 ||
        scale12<=0 || scale22<=0 ||  
        R_power11<2.5|| R_power12<2.5|| R_power22<2.5||fabs(col)>=1)) { rho=-2;break; }
    break;      
   case 122:  //bivariate sep matern
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale=par[5];
     smooth=par[6];
     if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale<=0 || smooth<0 || fabs(col)>=1) rho=-2;
     break;
     case 124:   /*parsimonious LMC*/
nug11=par[3];
nug22=par[4];
scale11=par[5];
scale22=par[6];
if(nug11<0 || nug22<0 ||  scale11<=0 || scale22<=0   ) rho=-2;
     break;
     case 126:   /*not parsimonious LMC*/
nug11=par[4];
nug22=par[5];
scale11=par[6];
scale22=par[7];
if(nug11<0 || nug22<0 ||  scale11<=0 || scale22<=0   ) rho=-2;
     break;
     case 118:  //bivariate matern with contrainsts
     case 121:
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale11=par[5];
     scale22=par[6];
     smoo11=par[7];
     smoo22=par[8];
     scale12=0.5*(scale11+scale22);
     smoo12=0.5*(smoo11+smoo22);
     if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 || scale12<=0 ||
        scale22<0|| smoo11<0 || smoo12<0 || smoo22<0 || fabs(col)>1) { rho=-2;break; }
        
       /*
     if(smoo12<(0.5*(smoo11+smoo22)))  {if(col!=0) { rho=-2;break; }} //ok
     
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12<=fmin(scale11,scale22)))
  
     { if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)) { rho=-2;break;}} //ok
     
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12> fmin(scale11,scale22)) && (scale12<fmax(scale11,scale22))) 
     {
      if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0) &&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0) &&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,
         sqrt(
         ((2*smoo22+2)*R_pow(scale12*scale11,2)+(2*smoo11+2)*R_pow(scale12*scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale22*scale11,2))/
         ((2*smoo11+2)*R_pow(scale11,2)+(2*smoo22+2)*R_pow(scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale12,2))
         ),0))
         {rho=-2;break;}
     }
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12>=fmax(scale11,scale22)))                                     
     {
         if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,1))
         {rho=-2;break;}
     }
     if(smoo12>(0.5*(smoo11+smoo22)))  {
        if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)&&
        R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0))
        {rho=-2;break;}
     }*/
     break;

   case 128:  //bivariate matern
   case 117: //bivariate smoke
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];  
     col=par[4];
     scale11=par[5];
     scale12=par[6];
     scale22=par[7];
     smoo11=par[8];
     smoo12=par[9];
     smoo22=par[10];
     if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 || scale12<=0 ||
        scale22<=0|| smoo11<=0 || smoo12<=0 || smoo22<=0|| fabs(col)>1) { rho=-2;break; }
        
        
     /*
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12<=fmin(scale11,scale22)))                                     
     { if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)) { rho=-2;break;}} //ok
     
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12> fmin(scale11,scale22)) && (scale12<fmax(scale11,scale22))) 
     {
      if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)&&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0)&&
         R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,
         sqrt(
         ((2*smoo22+2)*R_pow(scale12*scale11,2)+(2*smoo11+2)*R_pow(scale12*scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale22*scale11,2))/
         ((2*smoo11+2)*R_pow(scale11,2)+(2*smoo22+2)*R_pow(scale22,2)-2*(smoo11+smoo22+2)*R_pow(scale12,2))
         ),0))
         {rho=-2;break;}
     }
     if((smoo12==(0.5*(smoo11+smoo22))) && (scale12>=fmax(scale11,scale22)))                                     
     {
         if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,1))
         {rho=-2;break;}
     }
    
    if(smoo12<(0.5*(smoo11+smoo22)))  {if(col!=0) { rho=-2;break; }} //ok
            
     if(smoo12>(0.5*(smoo11+smoo22)))  {
         if(R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,0,0)&&
           R_pow(col,2)>bi_matern_bounds(scale11,scale22,scale12,smoo11,smoo22,smoo12,-LOW,0)) 
         {rho=-2;break;}
     }*/
     break;
   case 134:  //bivariate wendhole 1
   //case 69: //bivariate wendhole 2
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale11=par[5];
     scale12=par[6];
     scale22=par[7];
     smoo11=par[8];
     smoo12=par[9];
     smoo22=par[10];
    if(nug11<0 || nug22<0 || var11<=0 || var22<=0 || scale11<=0 || scale12<=0 || scale22<=0 || fabs(col)>=1) rho=-2;
     break;
    }
  return rho; 
}
// list of spatial and spatial-temporal correlation functions:

double CorFct(int *cormod, double h, double u, double *par, int c11, int c22)
{
  double arg=0.0, col=0.0,R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0.0, R_power_t=0.0, var11=0.0, var22=0.0;
  double rho=0.0, sep=0, scale=0.0, smooth=0.0,smooth_s=0.0,smooth_t=0.0, scale_s=0.0, scale_t=0, x=0, nug11=0.0, nug22=0.0;
  double scale11=0.0, scale22=0.0, scale12=0.0, smoo11=0.0, smoo22=0.0, smoo12=0.0,R_power11=0.0, R_power22=0.0, R_power12=0.0;
  switch(*cormod) // Correlation functions are in alphabetical order
    {
    case 1:// Cauchy correlation function
      R_power1=2;
      R_power2=par[0];
      scale=par[1];
      rho=CorFunCauchy(h, R_power2, scale);
      break;
    case 2:// Matern1
      scale=par[0];
      rho=exp(-h/scale)*(1+h/scale);
      break;
    case 3:// Matern2
      scale=par[0];
       rho=exp(-h/scale)*(1+h/scale+pow(h/scale,2)/3);
      break;
    case 4:// Exponential correlation function
      R_power=1;
      scale=par[0];
      rho=CorFunStable(h, R_power, scale);
      break;
    case 5: // Dagum
    R_power1=par[0];
    R_power2=par[1];
    scale=par[2];
    rho=CorFunDagum(h, R_power1, R_power2, scale);
    break;
    case 8: // Generalised Cuachy correlation function
      R_power1=par[0];
      R_power2=par[1];
      scale=par[2];
      rho=CorFunGenCauchy(h, R_power1, R_power2, scale);
      break;
    case 10:// Skarofski correlation function
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      rho=Shkarofski(h*h, scale_s,scale_t,smooth);
      break;
    case 11://wen0
        R_power=par[0];
        scale=par[1];
    rho=CorFunW0(h,scale,R_power);
        break;
    case 13://wen1
        R_power=par[0];
        scale=par[1];
    rho=CorFunW1(h,scale,R_power);
        break;
    case 12:// Stable correlation function
      R_power=par[0];
      scale=par[1];
      rho=CorFunStable(h, R_power, scale);
      break;
    case 14://  Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      rho=CorFunWitMat(h, scale, smooth);
      break;
    case 15://wen2
        R_power=par[0];
        scale=par[1];
        rho=CorFunW2(h,scale,R_power);
        break;
    case 16: //wave
      scale=par[0];
      rho=CorFunWave(h,scale);
      break;
    case 17://  multiquadric correlation function valid on sphere
        R_power=par[0];
        scale=par[1];
        rho=R_pow(1-R_power/2,2*scale)/R_pow(1+R_pow(R_power/2,2)-R_power*cos(h),scale);
        

    break;
    case 18://  sinsphere correlation function valid on sphere
        R_power=par[0];
        rho=1-R_pow(sin(h/2),R_power);    
    break;
    case 19: // Generalised wend correlation function
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
  rho=CorFunW_gen(h, R_power1, smooth, scale);
        break;
    case 6: // Generalised wend correlation function second paramtrizazion
        R_power1=1/par[0];        
        scale=par[1];
        smooth=par[2];
        sep=exp(  (lgammafn(2*smooth+R_power1+1)-lgammafn(R_power1))/ (1+2*smooth) );
        rho=CorFunW_gen(h, R_power1, smooth,  scale * sep);
        //Rprintf("%f %f \n",rho,R_power1);
        break;
     case 20://  Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      rho=CorFunSmoke(h, scale, smooth);
      break;    

      /*case 20: // Generalised wend correlation function
        R_power1=par[0];
        scale=par[1];
        smooth=par[2];
        rho=(pow(1,2*smooth+1)*CorFunW_gen(h, smooth+3+1.5, smooth, 1)-pow(0.75,2*smooth+1)*CorFunW_gen(h, smooth+3+1.5, smooth, 0.75) )/
           (pow(1,2*smooth+1)-pow(0.75,2*smooth+1));
        break;*/
 /***************** spatial tapers****************************/
   case 28:// Bohman taper
      rho=CorFunBohman(h,maxdist[0]);
      break;
   case 29:// Bohman model
      rho=CorFunBohman(h,par[0]);
        break;
  case 30:// Wendland1 for tap
      rho= CorFunW0(h,maxdist[0],2);
      break;
   case 31:// Wendland1 for model
      rho=CorFunW0(h,par[0],2);
            break;
    case 32:// Wendland1 for tap
      rho=CorFunW1(h,maxdist[0],3);
      break;
    case 33:// Wendland1 for tap
      rho=CorFunW1(h,par[0],3);
      break;
    case 34:// Wendland1 for tap
      rho=CorFunW2(h,maxdist[0],4);
      break;
    case 35:// Wendland1 for tap
     rho=CorFunW2(h,par[0],4);
    break;
        case 36:// unit taper
        case 37:// unit taper
      rho=1;
      break; 
 /***************** end spatial tapers****************************/
      // START non-separable correlation functions:
    case 42:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];//1/(1+exp(-par[4]));
      arg=1+R_pow(u/scale_t, R_power_t);
      rho=exp(-(R_pow(h/scale_s, R_power_s))*R_pow(arg, -0.5*sep*R_power_s))/arg;
      //Rprintf("w\n");
      // R_power=exp(sep)/(1+exp(sep));   // with reparametrization of beta parameter
      // rho=exp(-(R_pow(h/scale_s, R_power_s))*R_pow(arg, -0.5*R_power*R_power_s))/R_pow(arg,1);
      break;
    case 44:// Iaco-Cesare model as in (14) of Gneitint (2006): note that R_power parameters are in [0,1]
      R_power2=par[0];
      R_power_s=par[1];
      R_power_t=par[2];
      scale_s=par[3];
      scale_t=par[4];
      rho=R_pow(1+R_pow(h/scale_s, R_power_s)+R_pow(u/scale_t, R_power_t),-R_power2);
      break;
    case 46://Porcu model as in (4) of Porcu Bevilacqua et al (2010), with beta_1=beta_2=1
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      if(sep) rho=R_pow(0.5*R_pow(1+R_pow(h/scale_s, R_power_s),sep)+0.5*R_pow(1+R_pow(u/scale_t, R_power_t),sep),-1/sep);
      else rho=R_pow((1+R_pow(h/scale_s, R_power_s))*(1+R_pow(u/scale_t,R_power_t)),-1);
      break;
    case 48:// Stein model as in (16) of JASA (2005) with epsilon=0*/
      R_power_t=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      arg=smooth+R_pow(u, 0.5*R_power_t)/scale_t;
      if(h==0) rho=1/(R_pow(2, arg)*gammafn(arg+1));
      else rho=R_pow(h/scale_s, arg)*bessel_k(h/scale_s, arg, 1)/(R_pow(2, arg)*gammafn(arg+1));
      break;
    case 50:   //Gneiting correlation with prac ranges "giusto"
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+R_pow(u, R_power_t)/scale_t;
      rho=exp(-(R_pow(h, R_power_s)/scale_s)*R_pow(arg, 0.5*sep*R_power_s))/R_pow(arg,1.5);
      break;
     case 52:// Gneiting correlation model valid on the sphere (equation 8)
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
      arg=1+R_pow(h/scale_s, 2*R_power_s);
      rho=exp(-R_pow(u/scale_t, 2*R_power_t)/(R_pow(arg, sep*R_power_t)))/arg;
      break;
    case 54:// Gneiting correlation model valid on the sphere  (equazione 9)
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4];
     // arg=1+R_pow(h/scale_s, 2*R_power_s);
     // rho=R_pow(1+R_pow(u/scale_t, 2*R_power_t)/R_pow(arg,R_power_t*sep),-1)/R_pow(arg,1);
            arg=1+R_pow(u/scale_t, 2*R_power_t);
            rho=exp(-R_pow(h/scale_s, R_power_s)*R_pow(arg,R_power_s*sep))/R_pow(arg,1);
            
            
      break;
     case 56:    //st sinR_power
        R_power_s=1;
        R_power_t=par[0];
        scale_s=par[1];
        scale_t=par[2];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        sep=cos(h)*arg;
        rho=(exp(R_power_s*sep/scale_s)*(1+R_power_s*sep/scale_s))/((1+R_power_s/scale_s)*exp(R_power_s/scale_s));


      break;
     case 58:  //st multiquaderic
        R_power_s=par[0];
        R_power_t=par[1];
        scale_s=par[2];
        scale_t=par[3];
      arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
     rho= R_pow(R_pow(1-R_power_s/2,2)/(1+R_pow(R_power_s/2,2)-R_power_s*arg*cos(h)),scale_s);   // model B2 in the paper  (eq 9 right part)
     break;
    case 60:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      sep=par[4]; 
      arg=R_pow((1+(R_pow(u/scale_t,R_power_t))),0.5);
      if(h<(scale_s/R_pow(arg,sep))) rho=R_pow(1-h*R_pow(arg,sep)/scale_s,1.5+R_power_s)/arg;
      else                   rho=0;
      break;
    case 61:  //no sep gneiting  with temporal matern margin
        R_power_s=par[0];
        R_power=par[1];
        scale_s=par[2];
        scale_t=par[3];
        sep=par[4];
        smooth=par[5];
        arg=1+R_pow(h/scale_s, R_power_s);  
        if(u) rho=R_pow(arg,-R_power)*CorFunWitMat(u,scale_t*R_pow(arg,sep/2),smooth);
        else  rho=R_pow(arg,-R_power);

        break;

     case 62:  //no sep gneiting  with spatial matern margin

        R_power_t=par[0];
        R_power=par[1];
        scale_s=par[2];
        scale_t=par[3];   
        sep=par[4];
        smooth=par[5];
        arg=1+R_pow(u/scale_t, R_power_t);
        if(h)  rho=R_pow(arg,-R_power)*CorFunWitMat(h,scale_s*R_pow(arg,sep/2),smooth);
        else  rho=R_pow(arg,-R_power);
        break;
    case 63:  // 
         R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];  
       // Rprintf("%f %f  %f\n",sep,R_power_t,R_power);
        //arg=R_pow(1+R_pow(u/scale_t,R_power_t/2),-1/(R_power_t/2));
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        rho=R_pow(arg,R_power)*CorFunW0(h,scale_s*R_pow(arg,sep),R_power_s);   
        break;
          case 64:
        R_power_s=par[0];
        R_power=par[1];
        R_power_t=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
       // Rprintf("%f %f  %f\n",R_power_s,R_power,R_power_t);
        //arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        //rho=R_pow(arg,R_power)*CorFunW0(u,scale_t*R_pow(arg,sep),R_power_t); 
        
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
       /* Rprintf("f\n,arg ");
           arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
                Rprintf("f\n,arg ");  */ 
        rho=R_pow(arg,R_power)*CorFunW0(u,scale_t*R_pow(arg,sep),R_power_t);  //2.5+2*0
         //2.5+2*0
        break;
    case 65:  
          R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        rho=R_pow(arg,R_power)*CorFunW1(h,scale_s*R_pow(arg,sep),R_power_s);   
        break;

     case 66:
        R_power_s=par[0];
        R_power=par[1];
        R_power_t=par[2];  
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        rho=R_pow(arg,R_power)*CorFunW1(u,scale_t*R_pow(arg,sep),R_power_t); //2.5+2*1
        break;
      case 67:  //
         R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);   
        rho=R_pow(arg,R_power)*CorFunW2(h,scale_s*R_pow(arg,sep),R_power_s);   
        break;

     case 68:
        R_power_s=par[0];
        R_power_t=par[2];
        R_power=par[1];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        rho=R_pow(arg,R_power)*CorFunW2(u,scale_t*R_pow(arg,sep),R_power_t); ////2.5+2*2
        break;
    case 87:  
          R_power_t=par[0];
        R_power_s=par[1];
        R_power=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t),-1);
        rho=R_pow(arg,R_power)*CorFunW_gen(h,R_power_s,smooth,scale_s*R_pow(arg,sep)); 
        break;  
      case 88:
        R_power_s=par[0];
        R_power=par[1];
        R_power_t=par[2];
        scale_s=par[3];
        scale_t=par[4];
        sep=par[5];
        smooth=par[6];
        arg=R_pow(1+R_pow(h/scale_s,R_power_s),-1);
        rho=R_pow(arg,R_power)*CorFunW_gen(u,R_power_t,smooth,scale_t*R_pow(arg,sep));  
        break;
        case 69:
          R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW0(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
        break;
        case 70:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW0(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
        break;
        case 71:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW0(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
        break;
        case 72:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW1(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
        break;
        case 73:
          R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW1(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
        break;
          case 74:
                R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW1(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
        break;
          case 75:
                R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW2(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
        break;
          case 76:
                R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW2(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
        break;
        case 77:
              R_power_s=par[0];
          R_power_t=par[1];
          scale_s=par[2];
          scale_t=par[3];
          rho=CorFunW2(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
        break;
      // END non-separable correlation functions
      // START separable correlation functions:
    case 82:// Exp-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      rho=CorFunStable(h,1,scale_s)*R_pow((1+R_pow(u/scale_t, 2)), -R_power2);
      break;
    case 84:// Double exp:
      scale_s=par[0];
      scale_t=par[1];
      rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
      break;
     case 96:// prove
        R_power_s=par[0];
          R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      //Rprintf("%f %f %f %f\n",R_power_s,R_power_t,scale_s,scale_t);
      rho=exp(-    R_pow(h/scale_s,R_power_s)* exp(-u/scale_t)       );
      break;
    case 86:// Matern Matern
      scale_s=par[0];
      scale_t=par[1];
      smooth_s=par[2];
       smooth_t=par[3];
      rho=CorFunWitMat(h, scale_s,smooth_s)*CorFunWitMat(u, scale_t, smooth_t);
      break;
  /*  case 88://  exp_cos:
        scale_s=par[0];
        scale_t=par[1];
        rho=CorFunStable(h,1,scale_s)*CorFunWave(u,scale_t);
        break;*/
    case 90:// Matern-Cauchy:
      R_power2=par[0];
      scale_s=par[1];
      scale_t=par[2];
      smooth=par[3];
      if(h==0) arg=1;
      else arg=R_pow(2,1-smooth)/gammafn(smooth)*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
      rho=arg*R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
      break;
    case 92:// Matern-exp:
      scale_s=par[0];
      scale_t=par[1];
      smooth=par[2];
      if(h==0) arg=1;
      else arg=R_pow(2,1-smooth)/gammafn(smooth)*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
      rho=arg*exp(-u/scale_t);
      break;
    case 94:// Stable-stab:
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      rho=CorFunStable(h,R_power_s,scale_s)*CorFunStable(u,R_power_t,scale_t);
      break;
      /****************************** bivariate ***************************************/

   case 122:       //bivariate sep matern
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale=par[5];
     smooth=par[6];
     if((c11==0)&&(c22==0))                   {if(h==0)  rho=var11+nug11;
                                              else       rho=var11*CorFunWitMat(h, scale, smooth);
                                              break;}
     if((c11==0&&c22==1)||(c11==1&&c22==0))   {
                                          if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                          else     rho=col*sqrt(var11)*sqrt(var22)*CorFunWitMat(h, scale, smooth); 
                                          break;}
     if((c11==1)&&(c22==1))                   {if(h==0)  rho=var22+nug22;
                                              else       rho=var22*CorFunWitMat(h, scale, smooth);
                                              break;}
        break;

         case 119:       //bivariate smoke sep
     var11=par[0];
     var22=par[1];
     nug11=par[2];
     nug22=par[3];
     col=par[4];
     scale=par[5];
     smooth=par[6];
     if((c11==0)&&(c22==0))                   {if(h==0)  rho=var11+nug11;
                                              else       rho=var11*CorFunSmoke(h, scale, smooth);
                                              break;}
     if((c11==0&&c22==1)||(c11==1&&c22==0))   {
                                          if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                          else     rho=col*sqrt(var11)*sqrt(var22)*CorFunSmoke(h, scale, smooth); 
                                          break;}
     if((c11==1)&&(c22==1))                   {if(h==0)  rho=var22+nug22;
                                              else       rho=var22*CorFunSmoke(h, scale, smooth);
                                              break;}
        break;


        case 124:   /*parsimonious LMC */
        var11=par[0];
        col=par[1];
        var22=par[2];
        nug11=par[3];
        nug22=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=  CorFunStable(h, 1, scale11); //CorFunStable(h, 2, scale11);    //; for environmetrics
        smoo22=CorFunStable(h, 1, scale22); //CorFunWave(h,scale22);        //;      for environmetrics
        if((c11==0)&&(c22==0))                  { if(h==0) rho=rho+nug11;break;
                                                  rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;
                                                 }
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11+var22*col*smoo22;break;}
        if((c11==1)&&(c22==1))                  {if(h==0) rho=rho+nug22;break;
                                                 rho=R_pow(col,2)*smoo11+R_pow(var22,2)*smoo22;
                                                 }
        break;

        case 126:   /*not parsimonious LMC */
        var11=par[0];
        col=par[1];
        var22=par[2];
        smooth=par[3];
        nug11=par[4];
        nug22=par[5];
        scale11=par[6];
        scale22=par[7];
        smoo11=CorFunStable(h, 1, scale11);
        smoo22=CorFunStable(h, 1, scale22);
        if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;
                                                if(h==0) rho=rho+nug11;break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*smooth*smoo11+var22*col*smoo22;break;}
        if((c11==1)&&(c22==1))                  {rho=R_pow(smooth,2)*smoo11+R_pow(var22,2)*smoo22;
                                                if(h==0) rho=rho+nug22;break;}
        break;
        case 111:       // multi wend(k=1) separable
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power=par[5];
        scale=par[6];
        if((c11==0)&&(c22==0))            {if(h==0)  rho=var11+nug11;
                                           else       rho=var11*CorFunW0(h,scale,R_power);
                                           break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {
                                         if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                        else rho=col*sqrt(var11)*sqrt(var22)*CorFunW0(h,scale,R_power);break;}
        if((c11==1)&&(c22==1))            {if(h==0)  rho=var22+nug22;
                                           else      rho=var22*CorFunW0(h,scale,R_power);
                                           break;}
        break;
        case 113:       // multi wend(k=1) separable
            var11=par[0];
            var22=par[1];
            nug11=par[2];
            nug22=par[3];
            col=par[4];
            R_power=par[5];
            scale=par[6];
            if((c11==0)&&(c22==0))            {if(h==0)  rho=var11+nug11;
            else       rho=var11*CorFunW1(h,scale,R_power);
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))  {
                                                     if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunW1(h,scale,R_power);break;}
            if((c11==1)&&(c22==1))            {if(h==0)  rho=var22+nug22;
            else      rho=var22*CorFunW1(h,scale,R_power);
                break;}
            break;

        case 115:       // multi wend(k=1) separable
            var11=par[0];
            var22=par[1];
            nug11=par[2];
            nug22=par[3];
            col=par[4];
            R_power=par[5];
            scale=par[6];
            if((c11==0)&&(c22==0))            {if(h==0)  rho=var11+nug11;
            else       rho=var11*CorFunW2(h,scale,R_power);
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))  {
              if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW2(h,scale,R_power);break;}
            if((c11==1)&&(c22==1))            {if(h==0)  rho=var22+nug22;
            else      rho=var22*CorFunW2(h,scale,R_power);
                break;}
            break;

            
        case 112:       // multi wend(k=0)
                var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=par[7];
        scale11=par[8];
        scale22=par[9];
        scale12=par[10];
       // Rprintf("%f %f %f\n",scale11,scale22,scale12);
        if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
                                  else       rho=var11*CorFunW0(h,scale11,R_power11);
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunW0(h,scale12,R_power12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunW0(h,scale22,R_power22);;
                                  break;}
        
        break;
        case 114:       // multi wend(k=1)
            var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=par[7];
        scale11=par[8];
        scale22=par[9];
        scale12=par[10];
       // Rprintf("%f %f %f\n",scale11,scale22,scale12);
        if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
                                  else       rho=var11*CorFunW1(h,scale11,R_power11);
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunW1(h,scale12,R_power12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunW1(h,scale22,R_power22);;
                                  break;}
        
        break;
        case 116:       // multi wend(k=2)
       var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=par[7];
        scale11=par[8];
        scale22=par[9];
        scale12=par[10];
       // Rprintf("%f %f %f\n",scale11,scale22,scale12);
        if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
                                  else       rho=var11*CorFunW2(h,scale11,R_power11);
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunW2(h,scale12,R_power12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunW2(h,scale22,R_power22);;
                                  break;}
        break;
         case 129:       // multi wend(k=0) contr
          var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);
       // Rprintf("%f %f %f\n",scale11,scale22,scale12);
        if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
                                  else       rho=var11*CorFunW0(h,scale11,R_power11);
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunW0(h,scale12,R_power12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunW0(h,scale22,R_power22);;
                                  break;}
        break;
        case 131:       // // multi wend(k=1) contr
            var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power22=par[6];
        R_power12=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);
            if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
            else       rho=var11*CorFunW1(h,scale11,R_power11);
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW1(h,scale12,R_power12);;break;}
            if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
            else      rho=var22*CorFunW1(h,scale22,R_power22);
                break;}
        break;
        case 120:          // // multi wend(k=2) contr
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power22=par[6];
        R_power12=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);
  
            if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
            else       rho=var11*CorFunW2(h,scale11,R_power11);;
                break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0)) { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW2(h,scale12,R_power12);;break;}
            if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
            else      rho=var22*CorFunW2(h,scale22,R_power22);;
                break;}
        break;
        case 118:       // full bivariate matern with contraists
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        scale12=0.5*(scale11+scale22);
        smoo12=0.5*(smoo11+smoo22);
      
        if((c11==0)&&(c22==0))    {if(h==0)  rho=var11+nug11;
                                   else      rho=var11*CorFunWitMat(h, scale11,  smoo11);
                                   break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                else rho=col*sqrt(var11)*sqrt(var22)*CorFunWitMat(h, scale12,  smoo12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunWitMat(h, scale22,  smoo22);
                                  break;}
        break;
        case 121:       // full bivariate smoke with contraists
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        scale12=0.5*(scale11+scale22);
        smoo12=0.5*(smoo11+smoo22);
      // Rprintf("%f\n",col);
        if((c11==0)&&(c22==0))    {if(h==0)  rho=var11+nug11;
                                   else      rho=var11*CorFunSmoke(h, scale11,  smoo11);
                                   break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                else rho=col*sqrt(var11)*sqrt(var22)*CorFunSmoke(h, scale12,  smoo12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunSmoke(h, scale22,  smoo22);
                                  break;}
        break;
        case 128:       // full bivariate matern
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
     // Rprintf("%f %f %f %f %f %f %f %f %f %f %f  \n",var11,var22,nug11,nug22,col,scale11,scale12,scale22,smoo11,smoo12,smoo22);
        if((c11==0)&&(c22==0))  {if(h==0)  rho=var11+nug11;
                                else      rho=var11*CorFunWitMat(h, scale11,  smoo11);
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
                                               if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11*
                                                         var22)*CorFunWitMat(h, scale12,  smoo12);break;}
        if((c11==1)&&(c22==1))  {if(h==0)  rho=var22+nug22;
                                 else      rho=var22*CorFunWitMat(h, scale22,  smoo22);
                                 break;}
        break;
             case 117:       // full bivariate matern
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];     
        smoo12=par[9];
        smoo22=par[10];
     // Rprintf("%f %f %f %f %f %f %f %f %f %f %f  \n",var11,var22,nug11,nug22,col,scale11,scale12,scale22,smoo11,smoo12,smoo22);
        if((c11==0)&&(c22==0))  {if(h==0)  rho=var11+nug11;
                                else      rho=var11*CorFunSmoke(h, scale11,  smoo11);
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
                                               if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11*
                                                         var22)*CorFunSmoke(h, scale12,  smoo12);break;}
        if((c11==1)&&(c22==1))  {if(h==0)  rho=var22+nug22;
                                 else      rho=var22*CorFunSmoke(h, scale22,  smoo22);
                                 break;}
        break;
        /************************************************/
        case 130:       //biv gen wend sep
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power=par[5];
        scale=par[6];
        smooth=par[7];
        if((c11==0)&&(c22==0))            {if(h==0)  rho=var11+nug11;
                                           else       rho=var11*CorFunW_gen(h, R_power, smooth, scale);
                                           break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunW_gen(h, R_power, smooth, scale);break;}
        if((c11==1)&&(c22==1))            {if(h==0)  rho=var22+nug22;
                                           else      rho=var22*CorFunW_gen(h, R_power, smooth, scale);
                                           break;}
        break;

/************************************************/
              case 132:       //biv gen wend 
            var11=par[0];
            var22=par[1];
            nug11=par[2];
            nug22=par[3];
            col=par[4];
            R_power11=par[5];
            R_power12=par[6];
            R_power22=par[7];
            scale11=par[8];
            scale12=par[9];
            scale22=par[10];
            smoo11=par[11];
            smoo12=par[12];
            smoo22=par[13];
        if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
                                  else       rho=var11*CorFunW_gen(h, R_power11, smoo11, scale11);
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunW_gen(h, R_power12, smoo12, scale12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunW_gen(h, R_power22, smoo22, scale22);
                                  break;}
        break;

/************************************************/

                case 134:       // biv gen wend contr
          var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        R_power11=par[5];
        R_power12=par[6];
        R_power22=0.5*(R_power11+R_power22);
        scale11=par[7];
        scale22=par[8];
        scale12=0.5*(scale11+scale22);
        smoo11=par[9];
        smoo22=par[10];
        smoo12=0.5*(smoo11+smoo22);
        if((c11==0)&&(c22==0))   {if(h==0)  rho=var11+nug11;
                                  else       rho=var11*CorFunW_gen(h, R_power11, smoo11, scale11);
                                  break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  { if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*
                                                         sqrt(var22)*CorFunW_gen(h, R_power12, smoo12, scale12);break;}
        if((c11==1)&&(c22==1))   {if(h==0)  rho=var22+nug22;
                                  else      rho=var22*CorFunW_gen(h, R_power22, smoo22, scale22);
                                  break;}       
        break;

    /************************************************/
        case 136:       // matern cauchy
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        R_power22=par[11];
      //  Rprintf("%f %f %f %f%f %f %f \n",scale11,scale12,scale22,smoo11,smoo12,smoo22,R_power22);
       
        if((c11==0)&&(c22==0))  {if(h==0)  rho=var11+nug11;
                                else      rho=var11*CorFunWitMat1(h*h, scale11,  smoo11);
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                else rho=col*sqrt(var11)*
                                                         sqrt(var22)*Shkarofski(h*h,scale12,scale12,smoo12);
                                                       //CorFunWitMatCau(h,scale12,smoo12);
                                                       break;}
        if((c11==1)&&(c22==1))  {if(h==0)  rho=var22+nug22;
                                 else      rho=var22*CorFunGenCauchy2(h*h,smoo22,R_power22,scale22);
                                 break;}
        break;
           case 137:       // gen matern cauchy
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        R_power11=par[11];
        R_power12=par[12];
        R_power22=par[13];

        if((c11==0)&&(c22==0))  {if(h==0)  rho=var11+nug11;
                                else      rho=var11*CorFunGenWitMatCau(h, scale11,  smoo11,R_power11);
                             //   Rprintf("%f %f %f : %d %d\n",h,rho,c11,c22);
                                break;}
if((c11==0&&c22==1)||(c11==1&&c22==0)){ if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunGenWitMatCau(h, scale12,  smoo12,R_power12);
                                    //  Rprintf("%f %f : %d %d     %f %f %f\n",h,rho,c11,c22,scale12,smoo12,R_power12);
                                break;}
        
        if((c11==1)&&(c22==1))  {if(h==0)  rho=var22+nug22;
                                 else      rho=var22*CorFunGenWitMatCau(h, scale22,  smoo22,R_power22);
                               //  Rprintf("%f %f %f : %d %d\n",h,rho,c11,c22);
                                 break;}
        break;
           case 138:       //  bivariate cauchy
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        R_power11=par[8];   
        R_power12=par[9];
        R_power22=par[10];  
        smoo11=par[11];
        smoo12=par[12];
        smoo22=par[13];
      
        if((c11==0)&&(c22==0))  {if(h==0)  rho=var11+nug11;
                                else      rho=var11*CorFunGenCauchy(h, R_power11, smoo11, scale11);
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
          if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
           else rho=col*sqrt(var11)*sqrt(var22)*CorFunGenCauchy(h, R_power12, smoo12, scale12);break;}
        if((c11==1)&&(c22==1))  {if(h==0)  rho=var22+nug22;
                                 else      rho=var22*CorFunGenCauchy(h, R_power22, smoo22, scale22);
                                 break;}
    
        break;
        case 139:       //  bivariate stable
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        R_power11=par[8];   
        R_power12=par[9];
        R_power22=par[10];
        if((c11==0)&&(c22==0))  {if(h==0)  rho=var11+nug11;
                                else      rho=var11*CorFunStable(h, R_power11, scale11);
                                break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {
          if(h==0) rho=col*(sqrt(var11+nug11)*sqrt(var22+nug22));
                                                     else rho=col*sqrt(var11)*sqrt(var22)*CorFunStable(h, R_power12, scale12);break;}
        if((c11==1)&&(c22==1))  {if(h==0)  rho=var22+nug22;
                                 else      rho=var22*CorFunStable(h, R_power22, scale22);
                                 break;}
    
        break;
        /****************************** bivariate  taper ***************************************/
        case 140:       // wendland-gneiting k=0  
        if((c11==0)&&(c22==0))                   {rho=CorFunW0(h,dista[0][0],2);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=tapsep[0]*CorFunW0(h,dista[0][1],2);break;}
        if((c11==1)&&(c22==1))                   {rho=CorFunW0(h,dista[1][1],2);break;};
        break;
        case 141:       // wendland-gneiting k=0
        if((c11==0)&&(c22==0))                   {rho=CorFunW0(h,par[0],2);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=par[3]*CorFunW0(h,par[1],2);break;}
        if((c11==1)&&(c22==1))                   {rho=CorFunW0(h,par[2],2);break;};
        break;
        case 142:       // wendland-gneiting k=1
        if((c11==0)&&(c22==0))                    {rho=CorFunW1(h,dista[0][0],3);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {rho=tapsep[0]*CorFunW1(h,dista[0][1],3);break;}
        if((c11==1)&&(c22==1))                    {rho=CorFunW1(h,dista[1][1],3);break;}
        break;
        case 143:       // wendland-gneiting k=1
        if((c11==0)&&(c22==0))                    {rho=CorFunW1(h,par[0],3);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))    {rho=par[3]*CorFunW1(h,par[1],3);break;}
        if((c11==1)&&(c22==1))                    {rho=CorFunW1(h,par[2],3);break;}
        break;
        case 144:        // wendland-gneiting k=2
        if((c11==0)&&(c22==0))                   {rho=CorFunW2(h,dista[0][0],4);break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=tapsep[0]*CorFunW2(h,dista[0][1],4);break;}
        if((c11==1)&&(c22==1))                   {rho=CorFunW2(h,dista[1][1],4);break;}
         break;
        case 145:        // wendland-gneiting k=2
            if((c11==0)&&(c22==0))                   {rho=CorFunW2(h,par[0],4);break;}
            if((c11==0&&c22==1)||(c11==1&&c22==0))   {rho=par[3]*CorFunW2(h,par[1],4);break;}
            if((c11==1)&&(c22==1))                   {rho=CorFunW2(h,par[2],4);break;}
            break;
        case 146:        // Asymmetric taper with wendland-gneiting k=1
        if((c11==0)&&(c22==0)) {  rho=CorFunWend1_tap(h,dista[0][0],0);break;}
        if((c11==0&&c22==1))  {
            rho=1*(tapsep[0]*CorFunWend1_tap(h,dista[0][0],0)+(1-tapsep[0])*CorFunWend1_tap(h,dista[1][1],0));break;}
        if((c11==1&&c22==0))   {
            rho=1*(tapsep[0]*CorFunWend1_tap(h,dista[1][1],0)+(1-tapsep[0])*CorFunWend1_tap(h,dista[0][0],0));break;}
        if((c11==1)&&(c22==1)) {  rho=CorFunWend1_tap(h,dista[1][1],0);break;}
         break;
       case 147:        // Unit
        
        if((c11==0)&&(c22==0)) {  rho=1;break;}
        if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=1;break;}
        if((c11==1)&&(c22==1)) {  rho=1;break;}
         break;
       
        /******************************space time taper***************************************/
        case 200:
        rho=CorFunW0(h,maxdist[0],2)*CorFunW0(u,maxtime[0],2);
        break;
        case 201:
        rho=CorFunW0(h,par[0],2)*CorFunW0(u,par[1],2);
        break;
        case 202:
        rho=CorFunW0(h,maxdist[0],2)*CorFunW1(u,maxtime[0],3);
        break;
        case 203:
        rho=CorFunW0(h,par[0],2)*CorFunW1(u,par[1],3);
        break;
        case 204:
        rho=CorFunW0(h,maxdist[0],2)*CorFunW2(u,maxtime[0],4);
        break;
        case 205:
        rho=CorFunW0(h,par[0],2)*CorFunW2(u,par[1],4);
        break;
        case 206:
        rho=CorFunW1(h,maxdist[0],3)*CorFunW0(u,maxtime[0],2);
        break;
        case 207:
        rho=CorFunW1(h,par[0],3)*CorFunW0(u,par[1],2);
        break;
        case 208:
        rho=CorFunW1(h,maxdist[0],3)*CorFunW1(u,maxtime[0],3);
        break;
        case 209:
        rho=CorFunW1(h,par[0],3)*CorFunW1(u,par[1],3);
        break;
        case 210:
        rho=CorFunW1(h,maxdist[0],3)*CorFunW2(u,maxtime[0],4);
        break;
        case 211:
        rho=CorFunW1(h,par[0],3)*CorFunW2(u,par[1],4);
        break;
        case 212:
        rho=CorFunW2(h,maxdist[0],4)*CorFunW0(u,maxtime[0],2);
        break;
        case 213:
        rho=CorFunW2(h,par[0],4)*CorFunW0(u,par[1],2);
        break;
        case 214:
        rho=CorFunW2(h,maxdist[0],4)*CorFunW1(u,maxtime[0],3);
        break;
        case 215:
        rho=CorFunW2(h,par[0],4)*CorFunW1(u,par[1],3);
        break;
        case 216:
        rho=CorFunW2(h,maxdist[0],4)*CorFunW2(u,maxtime[0],4);
        break;
        case 217:
        rho=CorFunW2(h,par[0],4)*CorFunW2(u,par[1],4);
        break;
        
    /*non separable taper*/
        case 218:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/maxdist[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*0)*CorFunW0(u,maxtime[0]*arg,tapsep[2]);
        break;
        case 219:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/par[0], 1),-par[2]);
        x=u/(arg*par[1]);
        rho=R_pow(arg,3.5)*CorFunW0(x,1,3.5);
        break;
        case 220:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/maxdist[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*1)*CorFunW1(u,maxtime[0]*arg,tapsep[2]);
        break;  
        case 221:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(h/par[0], 1),-par[2]);
        x=u/(arg*par[1]);
        rho=R_pow(arg,5.5)*CorFunW1(x,1,4.5);
        break;   
         case 222:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(h/maxdist[0],tapsep[1]/2),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*2)*CorFunW2(u,maxtime[0]*arg,tapsep[2]);
        break;  
           case 223:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(h/par[0], 1),-par[2]);
        x=u/(arg*par[1]);
        rho=R_pow(arg,7.5)*CorFunW2(x,1,5.5);
        break; 


        /*non separable taper*/
        case 224:  /* non separable temporal adaptive  taper */

           R_power_t=par[0];
        R_power=par[1];
        scale_s=par[2];
        scale_t=par[3];
        sep=par[4];
        arg=R_pow(1+R_pow(u/scale_t,R_power_t/2),-sep/(R_power_t/2));
        rho=R_pow(arg,2.5+2*0)*CorFunW0(h,scale_s*arg,R_power);
        break;
         case 225:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(u/par[1], 1),-par[2]);
        x=h/(arg*par[0]);
        rho=R_pow(arg,3.5)*CorFunW0(x,1,3.5);
        break;
        case 226:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(u/scale_t,tapsep[1]),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*1)*CorFunW0(h/maxdist[0],arg,tapsep[2]);
        break;  
        case 227:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(u/par[1], 1),-par[2]);
        x=h/(arg*par[0]);
        rho=R_pow(arg,5.5)*CorFunW1(x,1,4.5);
        break;   
         case 228:  /* non separable temporal adaptive  taper */
        arg=R_pow(1+R_pow(u/scale_t,tapsep[1]),-tapsep[0]/(tapsep[1]/2));
        rho=R_pow(arg,2.5+2*2)*CorFunW0(h/maxdist[0],arg,tapsep[2]);
        break;  
           case 229:  /* non separable temporal adaptive  taper */
              arg=R_pow(1+R_pow(u/par[1], 1),-par[2]);
        x=h/(arg*par[1]);
        rho=R_pow(arg,7.5)*CorFunW2(x,1,5.5);
        break; 
         case 230:  //unit  space time taper
            rho=1;
            break;     
        /******************************end space time taper***************************************/
        // END separable correlation functions:

    }
  return rho;
}
// Cauhcy class of correlation models:
double CorFunCauchy(double lag, double R_power2, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=R_pow((1+R_pow(lag/scale,2)),-R_power2/2);
  return rho;
}
// Dagum:
double CorFunDagum(double lag, double R_power1, double R_power2, double scale)
{
    double rho=0.0;
    rho=1-R_pow(R_pow(lag/scale,R_power1)/(1+R_pow(lag/scale,R_power1)), R_power2/R_power1);
    return rho;
}
// Generalised Cauhcy class of correlation models:
double CorFunGenCauchy(double lag, double R_power1, double R_power2, double scale)
{
  double rho=0.0;
  rho=R_pow((1+R_pow(lag/scale,R_power1)), -R_power2/R_power1);
  return rho;
}

// Generalised Cauhcy class of correlation models:
double CorFunGenCauchy2(double lag, double R_power1, double R_power2, double scale)
{
  double rho=0.0;
  rho=R_pow((1+R_pow(lag,R_power2)/scale), -R_power1);
  return rho;
}


double CorFunWitMatCau(double h, double scale12,double smo12)
{
double C,B,dd=h*h;
B=1/bessel_k(1,smo12,1);
C=pow(1+dd/(scale12),-smo12/2)*bessel_k(sqrt(1+dd/(scale12)),smo12,1);
return(C*B);
}



double CorFunGenWitMatCau(double h, double scale,double smoo,double beta)
{
double C,B,dd=h*h;
B=1/bessel_k(sqrt(beta/scale),smoo,1);
C=pow(1+dd/beta,-smoo/2)*bessel_k(sqrt((dd+beta)/scale),smoo,1);
return(C*B);
}



// Stable class of correlation models:
double CorFunStable(double lag, double R_power, double scale)
{
  double rho=0.0;
  // Computes the correlation:
  rho=exp(-R_pow(lag/scale,R_power));
  return rho;
}
// Double Stable class of correlation models:
/*double CorFunDobStable(double lag, double R_power_s, double R_power_t, double scale_s, double scale_t, double tsep)
{
  double rho=0.0;
  // Computes the correlation:
  rho=exp(-R_pow(lag/scale_s,R_power_s)-R_pow(tsep/scale_t,R_power_t));
  return rho;
  }*/
// Sferical class of correlation models:
double CorFunSferical(double lag, double scale)
{
  double rho=0.0;
  if(lag<=scale) rho=1-1.5*lag/scale+0.5*R_pow(lag/scale, 3);
  else rho=0;
  return rho;
}

// Wave  correlation model:
double CorFunWave(double lag, double scale)
{
  double rho=0.0;
  if(lag==0) { rho=1;}
  else       { rho=(scale/lag)*sin(lag/(scale));}
  return rho;
}

// smoke class of correlation models:
double CorFunSmoke(double lag, double scale, double smooth)
{

  double rho=0.0,a=0.0,kk1=0.0,iscale=0.0;
  iscale=1/scale;
    a=0.5+smooth;
  // Computes the correlation:
  if(lag==0) {rho=1;}
  else  {

  kk1=(lgammafn(iscale+a)+lgammafn(iscale+smooth))-(lgammafn(2/scale+a)+lgammafn(smooth));
  rho=exp(kk1)*  pow(1-cos(lag),smooth)*hypergeo(1/scale+a,1/scale+smooth,2/scale+a,cos(lag));
    }
    
 /*   double *param;
    param=(double *) Calloc(3,double);
    kk1=(lgammafn(iscale+a)+lgammafn(iscale+smooth))-(lgammafn(2/scale+a)+lgammafn(smooth));
    param[0]= iscale;  // a
    param[1]= iscale+0.5; //b
    param[2]= 2/scale+a;  //c
    kk2=lgammafn(param[2])-(lgammafn(param[1])+lgammafn(param[2]-param[1]));
    res=HyperG_integral(cos(lag), param);
    Free(param);
    rho=res*exp(kk1+kk2);}*/
  return(rho);
}

// Whittle=matern class of correlation models:
double CorFunWitMat(double lag, double scale, double smooth)
{
  double rho=0.0;
  // Computes the correlation:
  if(lag==0) {rho=1;return rho;}
  if(smooth==0.5) {rho=exp(-lag/scale);return rho;}
  if(smooth==1.5) {rho=exp(-lag/scale)*(1+lag/scale);return rho;}
  if(smooth==2.5) {rho=exp(-lag/scale)*(1+lag/scale+ pow(lag/scale,2)/3);return rho;}
  rho=(R_pow(lag/scale,smooth)*bessel_k(lag/scale,smooth,1))/(R_pow(2,smooth-1)*gammafn(smooth));
  return rho;
}

double Shkarofski(double lag, double a,double b, double k)
{
double corr=0.0;
if(a==0 && k>0) return( R_pow(1+sqrt(lag/b),-2*k));
if(b==0 && k<0) return( R_pow(2,1+k) * R_pow(gammafn(-k),-1)  * 
                           R_pow(sqrt(lag/a),-k) * bessel_k(sqrt(lag/a),k,1));

corr=R_pow(1+lag/b,-k/2)*bessel_k(sqrt((b+lag)/a),k,1)/bessel_k(sqrt(b/a),k,1);
return(corr); 
}

double CorFunWitMat1(double lag, double scale, double smooth)
{
  double rho=0.0;
  double bb=sqrt(lag/scale);
  // Computes the correlation:
  if(lag==0) rho=1;
  else  rho=(R_pow(2,smooth+1)*R_pow(bb,-smooth)*bessel_k(bb,smooth,1))/(gammafn(-smooth));
  return rho;
}
double CorFunBohman(double lag,double scale)
{
  double rho=0.0,x=0;
  x=lag/scale;
  if(x<=1) {
       if (x>0) rho=(1-x)*(sin(2*M_PI*x)/(2*M_PI*x))+(1-cos(2*M_PI*x))/(2*M_PI*M_PI*x);
       else   rho=1;}
  else rho=0;
  return rho;
}

/* wendland function alpha=0*/
double CorFunW0(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=  R_pow(1-x,smoo);
    else rho=0;
    return rho;
}
/* wendland function alpha=1*/
double CorFunW1(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)  rho=R_pow(1-x,smoo+1)*(1+(smoo+1)*x);
    else rho=0;
    return rho;
}
/* wendland function alpha=2*/
double CorFunW2(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)  rho=R_pow(1-x,smoo+2)*(3+x*(3*smoo+6)+x*x*(R_pow(smoo,2)+4*smoo+3))/3;
    else rho=0;
    return rho;
}

/* generalized wendland function*/
double CorFunW_gen(double lag,double R_power1,double smooth,double scale)  // mu alpha beta
{
    double rho=0.0,x=0.0;

   
    if(lag==0) {rho=1; return(rho);}
   x=lag/scale;

    if(smooth==0) {
         if(x<=1)   rho=R_pow(1-x,R_power1);
         else rho=0;
         return(rho);
    }
    if(smooth==1) {
         if(x<=1) rho=R_pow(1-x,R_power1+1)*(1+x*(R_power1+1));
         else rho=0;
         return(rho);
    }
    if(smooth==2) {
        
         if(x<=1) rho=R_pow(1-x,R_power1+2)*(1+x*(R_power1+2)+x*x*(R_power1*R_power1 +4*R_power1 +3 )/3  );
         else rho=0;
         return(rho);
    }     
     
      /*first version  */     
    if(x<=1)
         {
        rho=exp((lgammafn(smooth)+lgammafn(2*smooth+R_power1+1))-(lgammafn(2*smooth)+lgammafn(smooth+R_power1+1)))
         *R_pow(2,-R_power1-1)*R_pow(1-x*x,smooth+R_power1)*hypergeo(R_power1/2,(R_power1+1)/2,smooth+R_power1+1, 1-x*x);
      }
  else {rho=0;}
   /*/second version
  
        x=lag;
        double *param;
        param=(double *) Calloc(3,double);
        param[0]=R_power1;param[1]=smooth;param[2]=scale;  //mu,alpha //beta
        rho=wendintegral(x,param);
        Free(param);*/
    return(rho);
}

double CorFunWend0_tap(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=  R_pow(1-x,smoo+3);
    else rho=0;
    return rho;
}

double CorFunWend1_tap(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,smoo+5)*(1+(smoo+5)*x);
    else rho=0;
    return rho;
}



double CorFunWend2_tap(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,smoo+7)*(1+(smoo+7)*x+R_pow(smoo+7,2)*R_pow(x,2)/3);
    else rho=0;
    return rho;
}


double CorFunWend1(double lag,double scale)
{
  double rho=0.0,x=0;
    x=lag/scale;
  if(x<=1) rho=R_pow(1-x,2)*(1+0.5*x);
  else rho=0;
  return rho;
}
double CorFunWend2(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,4)*(1+4*x);
  else rho=0;
  return rho;
}
double CorFunWend3(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,6)*(1+6*x+(35/3)*x*x);
  else rho=0;
  return rho;
}

double CorFunWend5(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,7)*(1+7*x+(48/3)*x*x);
    else rho=0;
    return rho;
}


double CorFunWend4(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=R_pow(1-x,5)*(1+5*x);
    else rho=0;
    return rho;
}
double CorFunWendhole1(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,4)*(1-3*x);
  else rho=0;
  return rho;
}

double CorFunWendhole2(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,5)*(1+4*x-18*R_pow(x,2));
  else rho=0;
  return rho;
}


double CorFunWendhole3(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,6)*(1+5*x-3.666667*R_pow(x,2)-70.36111*R_pow(x,3));
  else rho=0;
  return rho;
}

double CorFunWendhole(double lag,double scale)
{
  double rho=0.0,x=0;
     x=lag/scale;
  if(x<=1) rho=R_pow(1-x,5)*(1+5*x-27*R_pow(x,2));
  else rho=0;
  return rho;
}
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIAL CORRELATION MATRIX (upper trinagular) *******************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

// Computation of the upper (lower) triangular spatial correlation matrix: spatial case
void CorrelationMat2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod, 
 double *nuis, double *par,double *radius,int *ns, int *NS)
{
  int i=0,j=0,h=0;// check the paramaters range:
  double lags=0.0;
  //if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){rho[0]=-2;return;}// compute the correlations:
  //if(nuis[1]<0 || nuis[2]<=0){rho[0]=-2;return;}// compute the correlations:
   
     for(i=0;i<(ncoord[0]-1);i++){
	    for(j=(i+1);j<ncoord[0];j++){
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
    rho[h]=CorFct(cormod,lags,0,par,0,0);
       h++;
    }}
  return;
}

void CorrelationMat_poi2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod, double *mean,
 double *nuis, double *par,double *radius,int *ns, int *NS)
{
  int i=0,j=0,h=0;// check the paramaters range:
  double lags=0.0,corr=0.0,mui,muj;
  //if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){rho[0]=-2;return;}// compute the correlations:
  //if(nuis[1]<0 || nuis[2]<=0){rho[0]=-2;return;}// compute the correlations:
   
     for(i=0;i<(ncoord[0]-1);i++){
      for(j=(i+1);j<ncoord[0];j++){
      //   for(i=0;i<(ncoord[0]);i++){
      //for(j=i;j<ncoord[0];j++){
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
    corr=CorFct(cormod,lags,0,par,0,0);
     mui=exp(mean[i]);muj=exp(mean[j]);
    rho[h]=sqrt(mui*muj)*corr_pois((1-nuis[0])*corr,mui, muj);
       h++;
    }}
  return;
}
// Computation of the upper (lower) triangular spatial binomial type 1 covmatrix
void CorrelationMat_bin2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod, double *mean, 
        int *n,double *nuis, double *par,double *radius, int *ns, int *NS)
{
    int i=0,j=0,h=0;// check the paramaters range:
    double psj=0.0,lags=0.0,ai=0.0,aj=0.0,p1=0.0,p2=0.0;
    //if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){rho[0]=-2;return;}// compute the correlations:
    //if(nuis[1]<0 || nuis[2]<=0){rho[0]=-2;return;}// compute the correlations:
      for(i=0;i<(ncoord[0]-1);i++){
      for(j=(i+1);j<ncoord[0];j++){
    //     for(i=0;i<(ncoord[0]);i++){
    //  for(j=i;j<ncoord[0];j++){
        ai=mean[i];aj=mean[j];
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
        psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
              p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
             rho[h]=n[0]*(psj-p1*p2);///sqrt((1-p1)*p2*(1-p2));
            h++;
        }}
    return;
}



// Computation of the upper (lower) triangular spatial binomial type 2 covamtrix
void CorrelationMat_binneg2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod,  
  double *mean,int *n,double *nuis, double *par,double *radius, int *ns, int *NS)
{
    int i=0,j=0,h=0;// check the paramaters range:
    double  lags=0.0,psj=0.0,ai=0.0,aj=0.0,p1=0.0,p2=0.0;
      for(i=0;i<(ncoord[0]-1);i++){
      for(j=(i+1);j<ncoord[0];j++){
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
             ai=mean[i];aj=mean[j];
             psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
             p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
            rho[h]= corr_binomneg(n[0],p1,p2,psj);    
            h++;
        }}
    return;
}



// Computation of the upper (lower) triangular spatial geom cov matrix
void CorrelationMat_geom2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod,  double *mean,
  double *nuis, double *par,double *radius,int *ns, int *NS)
{
    int i=0,j=0,h=0;// check the paramaters range:
    double  lags=0.0,psj=0.0,ai=0.0,aj=0.0,p1=0.0,p2=0.0;
         for(i=0;i<(ncoord[0]-1);i++){
      for(j=(i+1);j<ncoord[0];j++){
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
             ai=mean[i];aj=mean[j];
             psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
             p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
            rho[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);
            h++;
        }}
    return;
}

// Computation of the upper (lower) triangular spatial tukey cov matrix
void CorrelationMat_tukeygh2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod, 
 double *nuis, double *par,double *radius,int *ns, int *NS)
{
  int i=0,j=0,h=0;// check the paramaters range:
  double cc,lags=0.0,a=0.0,b=0.0,c=0.0,tm=0.0,sk2=0.0,consta=0.0,t1;
    t1=1-nuis[4];
    c=R_pow(t1,2) - R_pow(nuis[4],2);
    sk2=R_pow(nuis[3],2);
    if(nuis[3]==0) {consta= ( 2/(1-nuis[4]*(2))   -1)/sqrt(c) ;}
    else {tm=(exp(sk2/(2*t1))-1)/(nuis[3]*sqrt(t1));
          consta= ((exp(sk2 * 2/(1-2*nuis[4])) -
               2* exp( sk2 *0.5)+1))/(sk2*sqrt(c)) - R_pow(tm,2);}

     for(i=0;i<(ncoord[0]-1);i++){
      for(j=(i+1);j<ncoord[0];j++){
        lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
        cc=CorFct(cormod,lags,0,par,0,0);

        a=(1+cc)/(1-nuis[4]*(1+cc));
        b=0.5* (1-nuis[4]*(1-cc*cc));

  if(nuis[3]==0)        rho[h]=(( a -2* b)/sqrt(c) )/consta;
  else                  rho[h]=((exp(sk2 * a) -2* exp( sk2 *b)+1)/(sk2*sqrt(c)) - R_pow(tm,2))/consta;

    h++;
    }}
  return;
}



// Computation of the correlations for spatial tapering:
void CorrelationMat_tap(double *rho,double *coordx, double *coordy, double *coordt, int *cormod,  double *nuis, double *par,double *radius,
  int *ns, int *NS)
{
  int i=0;// check the paramaters range:
  //if(nuis[1]<0 || nuis[2]<=0 ){rho[0]=-2;return;} //CheckCor(cormod,par)==-2
  //if(nuis[1]<0 || nuis[2]<=0){rho[0]=-2;return;}
  for(i=0;i<*npairs;i++) {rho[i]=CorFct(cormod,lags[i],0,par,0,0);}
  return;
}


/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIO-TEMPORAL CORRELATION MATRIX (upper trinagular) ***********************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/*************************************************************************************************/



// Computation of the correlations for spatio-temporal tapering:
void CorrelationMat_st_tap(double *rho,double *coordx, double *coordy, double *coordt, int *cormod, double *nuis, 
  double *par,double *radius, int *ns, int *NS)
{
  int i=0;
//if(nuis[1]<0 || nuis[2]<=0 ){rho[0]=-2;return;}
//if(nuis[1]<0 || nuis[2]<=0 || CheckCor(cormod,par)==-2){rho[0]=-2;return;}

  for(i=0;i<*npairs;i++) {
      rho[i]=CorFct(cormod,lags[i],lagt[i],par,0,0);
  }
return;
}


/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** Dynamic SPATIO-TEMPORAL CORRELATION MATRIX (upper trinagular) ***************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

// Computation of the upper (lower) triangular  correlation matrix: spatial-temporal case
void CorrelationMat_st_dyn2(double *rho, double *coordx, double *coordy, double *coordt,int *cormod,  
  double *nuis, double *par,double *radius, int *ns,int *NS)
{

  int i=0,j=0,t=0,v=0,h=0; double lags=0.0,lagt=0.0;
  for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
           rho[h]=CorFct(cormod,lags,0,par,t,v);
           h++;}}

    else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
        
rho[h]=CorFct(cormod,lags,lagt,par,t,v);
              h++;
              }}
    }}}
  return;
}

// Computation of the upper (lower) triangular  correlation matrix:  geom spatial-temporal case
void CorrelationMat_st_dyn_geom2(double *rho,double *coordx, double *coordy, double *coordt, 
                    int *cormod,  double *mean,double *nuis, double *par,double *radius, int *ns, int *NS)

{
    int i=0,j=0,t=0,v=0,h=0; double lags=0.0,lagt=0.0;
       double psj,ai,aj,p1,p2;
for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
    lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        ai=mean[i+ns[t]*t];aj=mean[j+ns[t]*v];
                        psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
                        p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
                   rho[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);
                        h++;}}
          else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                        ai=mean[i+ns[v]*t];aj=mean[j+ns[v]*v];
                        psj=pbnorm(cormod,lags,lagt,ai,aj,nuis[0],nuis[1],par,0);
                        p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
                  rho[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);
                        h++;}}
            }}}
    return;
}

// Computation of the upper (lower) triangular  correlation matrix: spatial-temporal case
void CorrelationMat_st_dyn_bin2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod,  double *mean,int *n,
  double *nuis, double *par,double *radius, int *ns, int *NS)

{
    int i=0,j=0,t=0,v=0,h=0; double lags=0.0,lagt=0.0;
       double psj,ai,aj,p1,p2;
for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
 lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                         ai=mean[i+ns[t]*t];aj=mean[j+ns[t]*v];
                        psj=pbnorm(cormod,lags,0,ai,aj,nuis[0],nuis[1],par,0);
                        p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
                        rho[h]=n[0]*(psj-p1*p2);
                        h++;}}
            else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
        
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                       ai=mean[i+ns[v]*t];aj=mean[j+ns[v]*v];
                        psj=pbnorm(cormod,lags,lagt,ai,aj,nuis[0],nuis[1],par,0);
                        p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
                        rho[h]=n[0]*(psj-p1*p2);
                        h++;}}
            }}}
    return;
}
// Computation of the upper (lower) triangular  correlation matrix: spatial-temporal case
void CorrelationMat_st_dyn_poi2(double *rho,double *coordx, double *coordy, double *coordt,  int *cormod,  double *mean,int *n,
  double *nuis, double *par,double *radius, int *ns, int *NS)

{
    int i=0,j=0,t=0,v=0,h=0; double lags=0.0,lagt=0.0;
       double mui,muj,corr;
for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[v];j++){
 lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                         mui=exp(mean[i+ns[t]*t]);muj=exp(mean[j+ns[t]*v]);
                        corr=CorFct(cormod,lags,0,par,t,v);
                      rho[h]= sqrt(mui* muj)*corr_pois((1-nuis[0])*corr,mui, muj);
                        h++;}}
               else {  
         lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
        
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                       mui=exp(mean[i+ns[v]*t]);muj=exp(mean[j+ns[v]*v]);
                               corr=CorFct(cormod,lags,lagt,par,t,v);
                     rho[h]=sqrt(mui* muj)*corr_pois((1-nuis[0])*corr,mui, muj);
                        h++;}}
            }}}
    return;
}

//###########################
void CorrelationMat_st_dyn_tukeygh2(double *rho, double *coordx, double *coordy, double *coordt,int *cormod,  double *nuis, double *par,double *radius,
  int *ns, int *NS)
{
  int i=0,j=0,t=0,v=0,h=0; 
   double lagt=0.0,cc,lags=0.0,a=0.0,b=0.0,c=0.0,tm=0.0,sk2=0.0,consta=0.0,t1;
  t1=1-nuis[4];
    c=R_pow(t1,2) - R_pow(nuis[4],2);
    sk2=R_pow(nuis[3],2);
    if(nuis[3]==0) {consta= ( 2/(1-nuis[4]*(2))   -1)/sqrt(c) ;}
    else {tm=(exp(sk2/(2*t1))-1)/(nuis[3]*sqrt(t1));
          consta= ((exp(sk2 * 2/(1-2*nuis[4])) -
               2* exp( sk2 *0.5)+1))/(sk2*sqrt(c)) - R_pow(tm,2);}

for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i+1;j<ns[t];j++){
           lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
         cc=CorFct(cormod,lags,0,par,t,v);
       a=(1+cc)/(1-nuis[4]*(1+cc));
        b=0.5* (1-nuis[4]*(1-cc*cc));

  if(nuis[3]==0)        rho[h]=(( a -2* b)/sqrt(c) )/consta;
  else                  rho[h]=((exp(sk2 * a) -2* exp( sk2 *b)+1)/(sk2*sqrt(c)) - R_pow(tm,2))/consta;
           h++;}}
    else {  
      lagt=fabs(coordt[t]-coordt[v]);
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
     cc=CorFct(cormod,0,lagt,par,t,v);
      a=(1+cc)/(1-nuis[4]*(1+cc));
        b=0.5* (1-nuis[4]*(1-cc*cc));

  if(nuis[3]==0)        rho[h]=(( a -2* b)/sqrt(c) )/consta;
  else                  rho[h]=((exp(sk2 * a) -2* exp( sk2 *b)+1)/(sk2*sqrt(c)) - R_pow(tm,2))/consta;
              h++;}}
    }}}
  return;
}

/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIO-DYNAMICAL BIVARIATE CORRELATION MATRIX (upper trinagular) **********/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/


// Computation of the upper (lower) triangular covariance matrix: bivariate case
void CorrelationMat_biv_dyn2(double *rho,double *coordx, double *coordy, double *coordt, int *cormod,  double *nuis, 
  double *par,double *radius, int *ns,int *NS)
{
  int i=0,j=0,t=0,v=0,h=0;double lags=0.0;
    for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i;j<ns[t];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                rho[h]=CorFct(cormod,lags,0,par,t,v);
                h++;}}
    else {  
         for(j=0;j<ns[v];j++){
          lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                rho[h]=CorFct(cormod,lags,0,par,t,v);
                h++;}}
    }}}
  return;
}


// Computation of the upper (lower) triangular covariance matrix: bivariate case
void CorrelationMat_biv_skew_dyn2(double *rho,double *coordx, double *coordy, double *coordt, int *cormod, 
 double *nuis, double *par,double *radius, int *ns,int *NS)
{
 int i=0,j=0,t=0,v=0,h=0;
  double cc=0.0,lags=0.0;
  int N=2;

 // if(CheckCor(cormod,par)==-2){rho[0]=-2;return;}
  //  if(par[0]<=0 || par[1]<=0 || par[2]<0 || par[3]<0){rho[0]=-2;return;
     double *vari;
    vari=(double *) Calloc(N,double);
     double *sk;
    sk=(double *) Calloc(N,double);
    vari[0]=par[0]; vari[1]=par[1];
    par[0]=1;par[1]=1;
    sk[0]=nuis[2];sk[1]=nuis[3];
     //if(CheckCor(cormod,par)==-2){rho[0]=-2;return;}
    //if(par[0]<=0 || par[1]<=0 || par[2]<0 || par[3]<0){rho[0]=-2;return;}
           for(t=0;t<ntime[0];t++){
    for(i=0;i<ns[t];i++){
      for(v=t;v<ntime[0];v++){
      if(t==v){
         for(j=i;j<ns[t];j++){
            lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                cc=CorFct(cormod,lags,0,par,t,v);
   rho[h]=2*sk[t]*sk[v]*(sqrt(1-cc*cc)+cc*asin(cc)-1)/M_PI+sqrt(vari[t])*sqrt(vari[v])*cc;
            h++;}}
   else {  
        for(j=0;j<ns[v];j++){
        lags=dist(type[0],coordx[(i+NS[t])],coordx[(j+NS[v])],coordy[(i+NS[t])],coordy[(j+NS[v])],*REARTH);
                cc=CorFct(cormod,lags,0,par,t,v);
   rho[h]=2*sk[t]*sk[v]*(sqrt(1-cc*cc)+cc*asin(cc)-1)/M_PI+sqrt(vari[t])*sqrt(vari[v])*cc;
                h++;}}
    }}}
  return;
}



// Computation of the correlations for bivariate tapering:
void CorrelationMat_biv_tap(double *rho, double *coordx, double *coordy, double *coordt,int *cormod,
 double *nuis, double *par,double *radius, int *ns,int *NS)
{
  int i=0;
  //if(CheckCor(cormod,par)==-2){rho[0]=-2;return;}
      // if(par[0]<=0 || par[1]<=0 || par[2]<0 || par[3]<0){rho[0]=-2;return;}
    for(i=0;i<*npairs;i++) { rho[i]=CorFct(cormod,lags[i],0,par,first[i],second[i]);}
  return;
}

/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/****************** SPATIO-Temporal(BIVARIATE) CORRELATION Vector *******************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

//computation of correlation between a points and a vector (for kriging)
void Corr_c(double *cc,double *coordx, double *coordy, double *coordt, int *cormod, int *grid, double *locx,  double *locy,int *ncoord, int *nloc,int *tloc,
                 int *ns,int *NS,int *ntime, double *par, int *spt, int *biv, double *time,int *type, int *which,double *radius)
{
int i,j,h=0;
double dis=0.0;
    
if(!spt[0]&&!biv[0])	{   //spatial case     
     for(j=0;j<(*nloc);j++){
	     for(i=0;i<(*ncoord);i++){
           dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);
	      cc[h]=CorFct(cormod,dis,0,par,0,0);
	      h++;
	      }}}
else{    //spatio temporal  case or bivariate case
int t,v;double dit=0.0;
if(*spt) {   

        for(j=0;j<(*nloc);j++){
        for(v=0;v<(*tloc);v++){
           for(t=0;t<*ntime;t++){
                      dit=fabs(coordt[t]-time[v]);
                for(i=0;i<ns[t];i++){      
                   dis=dist(type[0],coordx[(i+NS[t])],locx[j],
                                    coordy[(i+NS[t])],locy[j],radius[0]);
       cc[h]=CorFct(cormod,dis,dit,par,0,0);
        h++;}}}}        
}

    if(*biv) {  
       for(j=0;j<(*nloc);j++){
            for(t=0;t<*ntime;t++){
                    for(i=0;i<ns[t];i++){
                      dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
                                cc[h]=CorFct(cormod,dis,0,par,which[0],t);
                                h++;}}}}
}
}

///compute the covariance btwen loc to predict and locaton sites for binomial and geometric RF
void Corr_c_poi(double *cc,double *coordx, double *coordy, double *coordt, int *cormod, int *grid, double *locx,  double *locy,int *ncoord, int *nloc,
                int *model,int *tloc,double *n, int *ns,int *NS,int *ntime, double *mean,double *nuis, double *par, int *spt, int *biv, double *time,int *type, int *which,double *radius)
{
    int i=0,j=0,h=0;
    int t=0,v=0;double dit=0.0;
    double dis=0.0,mui=0.0,muj=0.0,corr=0.0;
    if(!spt[0]&&!biv[0])  {   //spatial case
    
          for(j=0;j<(*nloc);j++){
                for(i=0;i<(*ncoord);i++){
                     dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);
                     corr=CorFct(cormod,dis,0,par,0,0);
                        mui=exp(mean[i]);
                        muj=exp(mean[j]);
                       cc[h]=sqrt(mui*muj)*corr_pois((1-nuis[0])*corr,mui, muj);  
                    h++;}}
           
      }
    else{    //spatio temporal  case or bivariate case

           if(*spt) {
        for(j=0;j<(*nloc);j++){
        for(v=0;v<(*tloc);v++){
           for(t=0;t<*ntime;t++){
                      dit=fabs(coordt[t]-time[v]);
               for(i=0;i<ns[t];i++){
                  //dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);   
                 ///dis=dist(type[0],coordx[(i+NS[t])],locx[(j+NS[v])],coordy[(i+NS[t])],locy[(j+NS[v])],radius[0]);
                   dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
                   corr=CorFct(cormod,dis,dit,par,0,0);
                             mui=exp(mean[j+*nloc * v]);      
                             muj=exp(mean[i+*nloc * t]);
                             cc[h]=sqrt(mui*muj)*corr_pois((1-nuis[0])*corr,mui, muj);
                    h++;}}}}
          }
      /* if(*biv) {
            
                // in case of an irregular grid of coordinates:
                for(j=0;j<(*nloc);j++){
                    for(i=0;i<*ncoord;i++){
                        dis=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                        for(t=0;t<*ntime;t++){
                            cc[h]=CorFct(cormod,dis,0,par,which[0],t);
                            h++;}}}
              }*/        
    }}

///compute the covariance btwen loc to predict and locaton sites for binomial and geometric RF
void Corr_c_bin(double *cc,double *coordx, double *coordy, double *coordt, int *cormod, int *grid, double *locx,  double *locy,int *ncoord, int *nloc,
                int *model,int *tloc,double *n, int *ns,int *NS,int *ntime, double *mean,double *nuis, double *par, int *spt, int *biv, double *time,int *type, int *which,double *radius)
{
    int i=0,j=0,h=0;
    int t=0,v=0;double dit=0.0;
    double dis=0.0,p1=0.0,p2=0.0,psj=0.0,ai=0.0,aj=0.0;
    if(!spt[0]&&!biv[0])  {   //spatial case
    
            for(j=0;j<(*nloc);j++){
                for(i=0;i<(*ncoord);i++){
                     dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);
                     //dis=dist(type[0],coordx[(i+NS[t])],locx[(j+NS[v])],coordy[(i+NS[t])],locy[(j+NS[v])],radius[0]);
                        ai=mean[i];aj=mean[j];
                        psj=pbnorm(cormod,dis,0,ai,aj,nuis[0],nuis[1],par,0);
                        p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
                       // compute the covariance!
                    if(*model==2||*model==11||*model==19) cc[h]=n[0]*(psj-p1*p2);//binomial
                   if(*model==14)            cc[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);  // geometric
                    h++;}}
           
      }
    else{    //spatio temporal  case or bivariate case

        if(*spt) {
        for(j=0;j<(*nloc);j++){
        for(v=0;v<(*tloc);v++){
           for(t=0;t<*ntime;t++){
                      dit=fabs(coordt[t]-time[v]);
               for(i=0;i<ns[t];i++){
                  //dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);   
                 ///dis=dist(type[0],coordx[(i+NS[t])],locx[(j+NS[v])],coordy[(i+NS[t])],locy[(j+NS[v])],radius[0]);
                   dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
                             ai=mean[j+*nloc * v];      
                             aj=mean[i+*nloc * t];
                                psj=pbnorm(cormod,dis,dit,ai,aj,nuis[0],nuis[1],par,0);
                                p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);
                               // Rprintf("%f %f %f %f %f %f %f %f \n",n[0],dis,dit,psj,p1,p2,nuis[0],nuis[1]);
                   if(*model==2||*model==11||*model==19) cc[h]=n[0]*(psj-p1*p2);//binomial
                   if(*model==14)            cc[h]=(psj-p1*p2)/((-psj+p1+p2)*p1*p2);  // geometric
                    h++;}}}}
          }
      /* if(*biv) {
            
                // in case of an irregular grid of coordinates:
                for(j=0;j<(*nloc);j++){
                    for(i=0;i<*ncoord;i++){
                        dis=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
                        for(t=0;t<*ntime;t++){
                            cc[h]=CorFct(cormod,dis,0,par,which[0],t);
                            h++;}}}
              }*/        
    }}




 void Corr_c_tap(double *cc,double *cc_tap,double *coordx, double *coordy, double *coordt, int *cormod, int *cormodtap, int *grid, double *locx,  double *locy,
                 double *mxd,double *mxt, int *ncoord, int *nloc, int *ns,int *NS,int*tloc,int *ntime, double *par, int *spt, int *biv, double *time,int *type,int *which,double *radius)
{
int i,j,h=0,*modtap;
double *partap,dis=0.0;

modtap=(int *) Calloc(1,int);*modtap=*cormodtap+1;
if(!spt[0]&&!biv[0])	{   //spatial case
    
partap=(double *) Calloc(1,double);
partap[0]=*mxd; // compact support
// in case of an irregular grid of coordinates:
    for(j=0;j<(*nloc);j++){
	  for(i=0;i<(*ncoord);i++){
	       dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);
	      cc[h]=CorFct(cormod,dis,0,par,0,0);
	      cc_tap[h]=cc[h]*CorFct(modtap,dis,0,partap,0,0);
	      h++;
	      }}
	 
    Free(partap);
}
else{    //spatio temporal  case or bivariate case
int t,v;double dit=0.0;
if(*spt) {  
    partap=(double *) Calloc(4,double);
    partap[0]=mxd[0]; // spatial compact support
    partap[2]=mxd[1]; // delta1 param for non separable taper
    partap[3]=mxd[2]; // delta1 param for non separable taper
    partap[1]=mxt[0]; // temporal compact support

 // in case of an irregular grid of coordinates:
	          for(j=0;j<(*nloc);j++){
        for(v=0;v<(*tloc);v++){
           for(t=0;t<*ntime;t++){
                      dit=fabs(coordt[t]-time[v]);
               for(i=0;i<ns[t];i++){
                  dis=dist(type[0],coordx[(i+NS[t])],locx[j],coordy[(i+NS[t])],locy[j],radius[0]);
		    cc[h]=CorFct(cormod,dis,dit,par,t,v);
	         cc_tap[h]=cc[h]*CorFct(modtap,dis,dit,partap,t,v);
		    h++;}}}}

    Free(partap);}
    if(*biv)  {
        
        partap=(double *) Calloc(3,double);
        partap[0]=mxd[0]; // spatial compact support
        partap[1]=mxd[1]; // temporal compact support
        partap[2]=mxd[2]; // temporal compact support
        partap[3]=mxd[3]; // colocated taper
            
             for(j=0;j<(*nloc);j++){
            for(t=0;t<*ntime;t++){
                    for(i=0;i<ns[t];i++){ //for(i=0;i<*ncoord;i++){
                        dis=dist(type[0],coordx[i],locx[j],coordy[i],locy[j],radius[0]);
                                cc[h]=CorFct(cormod,dis,0,par,which[0],t);
                                h++;}}}
      
        Free(partap);}
}
    Free(modtap);
}


// Derivatives with respect to R_power2 of the Cauchy correlation model:
double DCauchyPow(double R_power2, double scale, double rho)
{
  return -rho*log(R_pow(rho,-2/R_power2));
}
// Derivatives with respect to scale of the Cauchy correlation model:
double DCauchySc(double lag, double R_power2, double scale, double rho)
{
 return R_power2*rho*R_pow(rho, 2/R_power2)*R_pow(lag,2)/R_pow(scale,3);
}
// Derivatives with respect to scale of the Exponential correlation model:
double DExpoSc(double lag, double scale, double rho)
{
 return rho*lag/R_pow(scale,2);
}
// Derivatives with respect to scale of the Gaussian correlation model:
double DGaussSc(double lag, double scale, double rho)
{
  return 2*rho*R_pow(lag,2)/R_pow(scale,3);
}
// Derivatives with respect to R_power1 of the generalised Cauchy correlation model:
double DGenCauP1(double lag, double R_power1, double R_power2, double scale, double rho)
{
  if(lag)return R_power2*rho/R_power1*(log(1+R_pow(lag/scale,R_power1))/R_power1-
           R_pow(lag/scale,R_power1)*log(lag/scale)/
           (1+R_pow(lag/scale,R_power1)));
  else return 0.0;
}
// Derivatives with respect to R_power2 of the generalised Cauchy correlation model:
double DGenCauP2(double lag, double R_power1, double scale, double rho)
{
  return -rho*log(1+R_pow(lag/scale,R_power1))/R_power1;
}
// Derivatives with respect to scale of the generalised Cauchy correlation model:
double DGenCauSc(double lag, double R_power1, double R_power2, double scale, double rho)
{
  if(lag) return rho/(1+R_pow(lag/scale,2))*R_power2*R_pow(lag,R_power1)/R_pow(scale,R_power1+1);
  else return 0.0;
}
// Derivatives with respect to scale of the sferical correlation model:
double DSferiSc(double lag, double scale)
{
  if(lag<=scale) return 1.5*lag*(1-R_pow(lag/scale, 2))/R_pow(scale, 2);
  else return 0.0;
}

// Derivatives with respect to scale of the Wen1 correlation model:
double DWen1Sc(double lag, double scale, double smooth)
{
    
    if(lag<=scale) return (smooth+5)*(smooth+6)*lag*lag*R_pow(lag-scale,4)*R_pow(((scale-lag)/scale),smooth)/R_pow(scale,7);
    else return 0;
}

// Derivatives with respect to R_power of the Stable correlation model:
double DStabPow(double lag, double R_power, double scale, double rho)
{
  if(lag) return -rho*R_pow(lag/scale,R_power)*log(lag/scale);
  else return 0.0;
}
// Derivatives with respect to scale of the Stable correlation model:
double DStabSc(double lag, double R_power, double scale, double rho)
{
  if(lag) return rho*R_power*R_pow(lag/scale,R_power)/scale;
  else return 0.0;
}
// Derivatives with respect to scale of the Whittle-Matern correlation model:
double DWhMatSc(double eps, double lag, double scale, double smooth)
{
  if (lag){
    double pscale=0.0;
    pscale=(bessel_k(lag/(scale+eps),smooth,1)-
      bessel_k(lag/scale,smooth,1))/ eps;
    return R_pow(2,1-smooth)/gammafn(smooth)*R_pow(lag/scale,smooth)*
      (pscale-smooth*bessel_k(lag/scale,smooth,1)/scale);}
  else return 0;
}

// Derivatives with respect to scale of the wave model
double DWaveSc(double lag, double scale)
{
    if(lag==0) return 0.0;
    else return  sin(lag/scale)/lag-cos(lag/scale)/scale;
}

// Derivatives with respect to smooth of the Whittle-Matern correlation model:
double DWhMatSm(double eps, double lag, double scale, double smooth)
{
  if (lag){
    double psmooth=0.0;
    psmooth=(bessel_k(lag/scale,smooth+ eps,1)-
       bessel_k(lag/scale,smooth,1))/ eps;
    return R_pow(2,1-smooth)*R_pow(lag/scale,smooth)/gammafn(smooth)*
      ((log(lag/scale)-log(2)-digamma(smooth))*bessel_k(lag/scale,smooth,1)+psmooth);}
  else return 0;
}

// Derivatives with respect to smooth of the Wend1 correlation model:
double DWen1Sm(double lag, double scale, double smooth)
{
    if (lag<=scale) return R_pow(lag-scale,5)*R_pow((scale-lag)/scale,smooth)*( (log((scale-lag)/scale)*(smooth*lag+5*lag+scale))+lag)/R_pow(scale,6);
    else return 0;
}

double DMat_Cauchy_sc_t(double h,double u,double R_power2,double scale_s,double scale_t,double smooth)
{
  double arg=0,arg3=0;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else arg=1;
  return 2*R_pow(u,2)*R_power2*arg*arg3/(R_pow(scale_t,3)*(1+R_pow(u/scale_t, 2)));
}

double DMat_Cauchy_pw2(double h,double u,double R_power2,double scale_s,double scale_t,double smooth)
{
  double arg=0.0,arg2=0.0,arg3=0.0;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else  arg=1;
  if(1+R_pow(u/scale_t, 2)) arg2=arg*arg3*log(1+R_pow(u/scale_t, 2));
 return arg2;
}


double DMat_Cauchy_sc_s(double h,double u,double R_power2,double scale_s,double scale_t,double smooth)
{
  double arg1,arg2=0,arg3;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  if(h) {  arg2=2*smooth*scale_s*bessel_k(h/scale_s, smooth,1)-h*bessel_k(h/scale_s, smooth+1,1);}
  arg1=R_pow(2, 1-smooth)*R_pow(h/scale_s, smooth)*arg3;
  return -arg1*arg2/(gammafn(smooth)*R_pow(scale_s,2));
}

double DMat_Cauchy_sm(double h,double u,double eps, double R_power2, double scale_s,double scale_t,double smooth)
{
  double arg=0.0,arg2=0.0,arg3=0.0,psmooth=0.0;
  arg3=R_pow((1+R_pow(u/scale_t, 2)),-R_power2);
  psmooth=(bessel_k(h/scale_s,smooth+ eps,1)-bessel_k(h/scale_s,smooth,1))/ eps;
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
  else arg=1;
  if(h) arg2=-arg3*arg*(log(2)+digamma(smooth)-log(h/scale_s)-psmooth/bessel_k(h/scale_s, smooth, 1));
  else arg2=0;
  return arg2;
}

double DMat_Exp_sc_t(double h,double u,double scale_s,double scale_t,double smooth)
{
  double arg=0;
  if(h) arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
 else arg=1;
  return arg*u*exp(-u/scale_t)/R_pow(scale_t,2);
}

double DMat_Exp_sc_s(double h,double u,double scale_s,double scale_t,double smooth)
{
double arg1=0,arg2=0;

 if(h){arg1=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
   arg2=h*(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth+1, 1);}
 else {arg1=1;
   arg2=0;}

 return -arg1*2*smooth*exp(-u/scale_t)/scale_s + arg2*exp(-u/scale_t)/R_pow(scale_s,2);
}

double DMat_Exp_sm(double h,double u,double eps,double scale_s,double scale_t,double smooth)
{
  double arg=0.0,psmooth=0.0,arg2=0.0;
  psmooth=(bessel_k(h/scale_s,smooth+ eps,1)-bessel_k(h/scale_s,smooth,1))/ eps;
  if(h){arg=(R_pow(2, 1-smooth)/gammafn(smooth))*R_pow(h/scale_s, smooth)*bessel_k(h/scale_s, smooth, 1);
    arg2=-exp(-u/scale_t)*arg*(log(2)+digamma(smooth)-log(h/scale_s)-psmooth/bessel_k(h/scale_s,smooth,1));}
  else arg2=0;
  return arg2;
}

/***************************************************/
/* derivative of bivariate wendland2 model */
/***************************************************/
double DWen1sep_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                              rho= CorFunWend1_tap(h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))              rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)* CorFunWend1_tap(h,scale,smoo);
    return rho;
}
double DWen1sep_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)* CorFunWend1_tap(h,scale,smoo);
    if((c11==1)&&(c22==1))                             rho= CorFunWend1_tap(h,scale,smoo);
    return rho;
}
double DWen1sep_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1sep_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1sep_biv_scale(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWen1Sc(h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWen1Sc(h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWen1Sc(h, scale,smoo);
    return rho;
}
double DWen1sep_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)* CorFunWend1_tap(h,scale,smoo);
    return rho;
}
double DWen1sep_biv_smoo(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWen1Sm(h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWen1Sm(h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWen1Sm(h,scale,smoo);
    return rho;
}
/***************************************************/
/***************************************************/

/***************************************************/
/* derivative of bivariate matern separable  model */
/***************************************************/
double Dmatsep_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                  rho=CorFunWitMat(h, scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)*CorFunWitMat(h, scale,smoo);
    return rho;
}
double Dmatsep_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)*CorFunWitMat(h, scale,smoo);;
    if((c11==1)&&(c22==1))                             rho=CorFunWitMat(h, scale,smoo);
    return rho;
}
double Dmatsep_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double Dmatsep_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale ,double smoo,double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double Dmatsep_biv_scale(double h,double eps,double var11,double var22,double nug11,double nug22, double scale ,double smoo,double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWhMatSc(eps,h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWhMatSc(eps,h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWhMatSc(eps,h,scale,smoo);
    return rho;
}

double Dmatsep_biv_smo(double h,double eps,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWhMatSm(eps,h,scale,smoo);
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWhMatSm(eps,h,scale,smoo);
    if((c11==1)&&(c22==1))                   rho=var22*DWhMatSm(eps,h,scale,smoo);
    return rho;
}
double Dmatsep_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale,double smoo, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)*CorFunWitMat(h, scale,smoo);;
    return rho;
}
/***************************************************/
/***************************************************/
/***************************************************/
/* derivative of full bivariate matern model */
/***************************************************/
double DMat_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                     rho=CorFunWitMat(h, scale11,smoo11);
    if((c11==0&&c22==1)||(c11==1&&c22==0))     rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)*CorFunWitMat(h, scale12,smoo12);
    return rho;
}
double DMat_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)*CorFunWitMat(h, scale12,smoo12);
    if((c11==1)&&(c22==1))                             rho=CorFunWitMat(h, scale22,smoo22);
    return rho;
}
double DMat_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double DMat_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double DMat_biv_scale1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                 rho=var11*DWhMatSc(eps,h,scale11,smoo11);
    return rho;
}

double DMat_biv_scale2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))                 rho=var22*DWhMatSc(eps,h,scale22,smoo22);
    return rho;
}
double DMat_biv_scale1_contr(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWhMatSc(eps,h/0.5,scale11+scale22,0.5*(smoo11+smoo22));
    if((c11==0)&&(c22==0))                 rho=var11*DWhMatSc(eps,h,scale11,smoo11);
    return rho;
}
double DMat_biv_scale2_contr(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWhMatSc(eps,h/0.5,scale11+scale22,0.5*(smoo11+smoo22));
    if((c11==1)&&(c22==1))                 rho=var22*DWhMatSc(eps,h,scale22,smoo22);
    return rho;
}
double DMat_biv_scale12(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                        double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))    rho=col*sqrt(var11*var22)*DWhMatSc(eps,h,scale12,smoo12);
    return rho;
}
double DMat_biv_smo1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWhMatSm(eps,h,scale11,smoo11);
    return rho;
}
double DMat_biv_smo12(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                      double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWhMatSm(eps,h,scale12,smoo12);
    return rho;
}
double DMat_biv_smo2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))            rho=var22*DWhMatSm(eps,h,scale22,smoo22);
    return rho;
}
double DMat_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                    double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)*CorFunWitMat(h, scale12,smoo12);
    return rho;
}


/***************************************************/
/* derivative of full bivariate wend (contr) model */
/***************************************************/
double DWen1_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
  double rho=0.0;
    if((c11==0)&&(c22==0))                     rho=CorFunWend1_tap(h, scale11,smoo11);
    if((c11==0&&c22==1)||(c11==1&&c22==0))     rho=0.5*R_pow(var11,-0.5)*col*sqrt(var22)*CorFunWend1_tap(h, scale12,smoo12);
    return rho;
}
double DWen1_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
   double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))             rho=0.5*R_pow(var22,-0.5)*col*sqrt(var11)*CorFunWend1_tap(h, scale12,smoo12);
    if((c11==1)&&(c22==1))                             rho=CorFunWend1_tap(h, scale22,smoo22);
    return rho;
}
double DWen1_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;;
    if((c11==0)&&(c22==0))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==1)&&(c22==1))  {if(h==0)   rho=1;}
    return rho;
}
double DWen1_biv_scale1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                      double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                 rho=var11*DWen1Sc(h,scale11,smoo11);
    return rho;
}

double DWen1_biv_scale2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
   double rho=0.0;
    if((c11==1)&&(c22==1))                 rho=var22*DWen1Sc(h,scale22,smoo22);
    return rho;
    
}

double DWen1_biv_scale1_contr(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                 rho=var11*DWen1Sc(h,scale11,smoo11);
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWen1Sc(h,scale22+scale11,0.5*(smoo11+smoo12));
    return rho;
}


double DWen1_biv_scale2_contr(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0)) rho=col*sqrt(var11*var22)*DWen1Sc(h,scale22+scale11,0.5*(smoo11+smoo12));
    if((c11==1)&&(c22==1))                 rho=var22*DWen1Sc(h,scale22,smoo22);
    return rho;
}
double DWen1_biv_scale12(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))    rho=col*sqrt(var11*var22)*DWen1Sc(h,scale12,smoo12);
    return rho;
}
double DWen1_biv_smo1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0)&&(c22==0))                   rho=var11*DWen1Sm(h,scale11,smoo11);
    return rho;
}
double DWen1_biv_smo12(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))   rho=col*sqrt(var11*var22)*DWen1Sm(h,scale12,smoo12);
    return rho;
}
double DWen1_biv_smo2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
     double rho=0.0;
    if((c11==1)&&(c22==1))            rho=var22*DWen1Sm(h,scale22,smoo22);
    return rho;
}
double DWen1_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                    double smoo11, double smoo22,double smoo12, double col,int c11,int c22)
{
    double rho=0.0;
    if((c11==0&&c22==1)||(c11==1&&c22==0))  rho=sqrt(var11*var22)*CorFunWend1_tap(h,scale12,smoo12);
    return rho;
}

/***************************************************/
/* derivative of LMC  (contr) model */
/***************************************************/

double DLMC_contr_var1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo11;
    smoo11=CorFunStable(h, 1, scale11);
    if(h==0) {smoo11=smoo11+nug11;}
    if((c11==0)&&(c22==0))                  {rho=2*var11*smoo11;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=col*smoo11;}
    return rho;
}
double DLMC_contr_var2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo22;
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo22=smoo22+nug22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=2*var22*smoo22;}
    return rho;
}
double DLMC_contr_nug1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo11,smoo22;
    smoo11=CorFunStable(h, 1, scale11);
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo11=1;smoo22=0;}
    else   {smoo11=0;smoo22=0;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11+var22*col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(col,2)*smoo11+R_pow(var22,2)*smoo22;}
    return rho;
}
double DLMC_contr_nug2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0, smoo11,smoo22;
    smoo11=CorFunStable(h, 1, scale11);
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo11=0;smoo22=1;}
     else   {smoo11=0;smoo22=0;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11+R_pow(col,2)*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11+var22*col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(col,2)*smoo11+R_pow(var22,2)*smoo22;}
    return rho;
}

double DLMC_contr_col(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0, smoo11,smoo22;
    smoo11=CorFunStable(h, 1, scale11);
    smoo22=CorFunStable(h, 1, scale22);
    if(h==0) {smoo11=smoo11+nug11;smoo22=smoo22+nug22;}
    if((c11==0)&&(c22==0))                  {rho=2*col*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*smoo11+var22*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=2*col*smoo11;}
    return rho;
}
double DLMC_contr_scale11(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0, smoo11;
    smoo11=DStabSc(h, 1,scale11,CorFunStable(h, 1, scale11));
    if(h==0) {smoo11=smoo11+nug11;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(var11,2)*smoo11;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var11*col*smoo11;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(col,2)*smoo11;}
    return rho;
}
double DLMC_contr_scale22(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22)
{
    double rho=0.0,smoo22;
    smoo22=DStabSc(h, 1,scale22,CorFunStable(h, 1, scale22));
    if(h==0) {smoo22=smoo22+nug22;}
    if((c11==0)&&(c22==0))                  {rho=R_pow(col,2)*smoo22;}
    if((c11==0&&c22==1)||(c11==1&&c22==0))  {rho=var22*col*smoo22;}
    if((c11==1)&&(c22==1))                  {rho=R_pow(var22,2)*smoo22;}
    return rho;
}




void GradCorrFct(double rho, int *cormod, double eps, int *flag,
     double *grad, double h, double u, int c11, int c22,double *par)
{
  int i=0;
  double R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0, R_power_t=0;
  double scale=0.0, scale_s=0, scale_t=0, smooth=0.0, smooth_s=0.0,smooth_t=0.0;
  double var11=0.0, var22=0.0, nug11=0.0, nug22=0.0,scale11=0.0,scale12=0.0,scale22=0.0,col=0.0,smoo11=0.0,smoo12=0.0,smoo22=0.0;

  switch(*cormod)// Correlation functions are in alphabetical order
    {//spatial gradients of correlations:
    case 1:// Cauchy correlation function
      R_power2=par[0];
      scale=par[1];
      if(flag[0]==1){//R_power parameter
    grad[i]=DCauchyPow(R_power2,scale,rho);i++;}
      if(flag[1]==1)//scale parameter
  grad[i]=DCauchySc(h, R_power2, scale, rho);
      break;
    case 4:// Exponential correlation function
      scale=par[0];//scale parameter
      if(flag[0] == 1) grad[i]=DExpoSc(h, scale, rho);
      break;
    //case 6:// Gaussian correlation function
     // scale=par[0];//scale parameter
     // if(flag[0]==1) grad[i]=DGaussSc(h,scale,rho);
     // break;
    case 8:// Generalised Cuachy correlation function
      R_power1=par[0];
      R_power2=par[1];
      scale=par[2];
      if(flag[0]==1){//R_power1 parameter
  grad[i]=DGenCauP1(h,R_power1,R_power2,scale,rho);i++;}
      if(flag[1]==1){//R_power2 parameter
  grad[i]=DGenCauP2(h,R_power1,scale,rho);i++;}
      if(flag[2]==1)//scale parameter
  grad[i]=DGenCauSc(h,R_power1,R_power2,scale,rho);
      break;
    case 10:// Sferical correlation function
      scale=par[0];//scale parameter
      if(flag[0]==1) grad[i]=DSferiSc(h,scale);
      break;
   /* case 11:// wend2 correlation function
        scale=par[0];//scale parameter
        if(flag[0]==1) grad[i]=DwenSc(h,scale);
        break;*/
    case 12:// Stable correlation function
      R_power=par[0];
      scale=par[1];
      if(flag[0]==1){//R_power parameter
  grad[i]=DStabPow(h,R_power,scale,rho);i++;}
      //scale parameter
      if(flag[1]==1) grad[i]=DStabSc(h,R_power,scale,rho);
      break;
    case 14:// Whittle-Matern correlation function
      scale=par[0];
      smooth=par[1];
      if(flag[0]==1){//scale parameter
  grad[i]=DWhMatSc(eps,h,scale,smooth);i++;}
      //smooth parameter
      if(flag[1]==1) grad[i]=DWhMatSm(eps,h,scale,smooth);
      break;
   case 16:// wave  correlation function
      scale=par[0];//scale parameter
      if(flag[0]==1) grad[i]=DWaveSc(h,scale);
      break;
      //spatio-temproal gradients of correlations:

    case 84://Double Exponential
      scale_s=par[0];
      scale_t=par[1];
      if(flag[0]==1){//spatial-scale parameter
  grad[i]=DStabSc(h,1,scale_s,rho);i++;}
      //temporal-scale parameter
      if(flag[1]==1) grad[i]=DStabSc(u,1,scale_t,rho);
      break;
    case 86://Exponential-Gaussian
      scale_s=par[0];
      scale_t=par[1];
      smooth_s=par[2];
      smooth_t=par[3];
      // to do..../
      break;
    case 94://Stable-Stable
      R_power_s=par[0];
      R_power_t=par[1];
      scale_s=par[2];
      scale_t=par[3];
      if(flag[0]==1){//spatial-R_power parameter
  grad[i]=DStabPow(h,R_power_s,scale_s,rho);i++;}
      if(flag[1]==1){//temporal-R_power parameter
  grad[i]=DStabPow(u,R_power_t,scale_t,rho);i++;}
      if(flag[2]==1){//spatial-scale parameter///
  grad[i]=DStabSc(h,R_power_s,scale_s,rho);i++;}
      if(flag[3]==1){//temporal-scale parameter
  grad[i]=DStabSc(u,R_power_t,scale_t,rho);}
     break;
    case 110:
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale_s=par[5];
        smooth=par[6];
        if(flag[0]==1){//first variance
            grad[i]=DWen1sep_biv_var1(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DWen1sep_biv_var2(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DWen1sep_biv_nug1(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DWen1sep_biv_nug2(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DWen1sep_biv_col(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[5]==1){//scle
            grad[i]=DWen1sep_biv_scale(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);i++;}
        if(flag[6]==1){//smooh
            grad[i]=DWen1sep_biv_smoo(h,var11,var22,nug11,nug22,scale_s,col,smooth,c11,c22);}
        break;
        case 114:    /// bi wen1 full
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        if(flag[0]==1){//first variance
            grad[i]=DWen1_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DWen1_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DWen1_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DWen1_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DWen1_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DWen1_biv_scale1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[6]==1){//scale12
            grad[i]=DWen1_biv_scale12(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[7]==1){//scale2
            grad[i]=DWen1_biv_scale2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[8]==1){//smo1
            grad[i]=DWen1_biv_smo1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[9]==1){//smo12
            grad[i]=DWen1_biv_smo12(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[10]==1){//smo2
            grad[i]=DWen1_biv_smo2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);}
        break;
      case 120:  /// bi wen1 withcontr
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        if(flag[0]==1){//first variance
            grad[i]=DWen1_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DWen1_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DWen1_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DWen1_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DWen1_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DWen1_biv_scale1_contr(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[6]==1){//scale2
            grad[i]=DWen1_biv_scale2_contr(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[7]==1){//smo1
            grad[i]=DWen1_biv_smo1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[8]==1){//smo2
            grad[i]=DWen1_biv_smo2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);}
             break;
      case 122:  /// bi matern sep
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale=par[5];
        smooth=par[6];
        if(flag[0]==1){//first variance
            grad[i]=Dmatsep_biv_var1(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=Dmatsep_biv_var2(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=Dmatsep_biv_nug1(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=Dmatsep_biv_nug2(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=Dmatsep_biv_col(h,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[5]==1){//scle
            grad[i]=Dmatsep_biv_scale(h,eps,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);i++;}
        if(flag[6]==1){//smooth
            grad[i]=Dmatsep_biv_smo(h,eps,var11,var22,nug11,nug22,scale,smooth,col,c11,c22);}
        break;
      case 124:  /// LMC contr
        var11=par[0];
         col=par[1];
        var22=par[2];
        nug11=par[3];
        nug22=par[4];
        scale11=par[5];
        scale22=par[6];
        if(flag[0]==1){//first variance
            grad[i]=DLMC_contr_var1(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[2]==1){//second varuance
            grad[i]=DLMC_contr_var2(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[1]==1){//colocated
            grad[i]=DLMC_contr_col(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[3]==1){//first nugget
            grad[i]=DLMC_contr_nug1(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[4]==1){//second nugget
            grad[i]=DLMC_contr_nug2(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[5]==1){//scle1
            grad[i]=DLMC_contr_scale11(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);i++;}
        if(flag[6]==1){//scale2
            grad[i]=DLMC_contr_scale22(h,eps,var11,var22,nug11,nug22,scale11,scale22,col,c11,c22);}
      break;
        case 118:   /// bi matern contr
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale22=par[6];
        smoo11=par[7];
        smoo22=par[8];
        if(flag[0]==1){//first variance
            grad[i]=DMat_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DMat_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DMat_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DMat_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DMat_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DMat_biv_scale1_contr(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[6]==1){//scale2
            grad[i]=DMat_biv_scale2_contr(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
        if(flag[7]==1){//smo1
            grad[i]=DMat_biv_smo1(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);i++;}
          if(flag[8]==1){//smo2
            grad[i]=DMat_biv_smo2(h,eps,var11,var22,nug11,nug22,scale11,scale22,(scale11+scale22)/2,smoo11,smoo22,(smoo22+smoo11)/2,col,c11,c22);}
        break;
        case 128:    /// bi matern full
        var11=par[0];
        var22=par[1];
        nug11=par[2];
        nug22=par[3];
        col=par[4];
        scale11=par[5];
        scale12=par[6];
        scale22=par[7];
        smoo11=par[8];
        smoo12=par[9];
        smoo22=par[10];
        if(flag[0]==1){//first variance
            grad[i]=DMat_biv_var1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[1]==1){//second varuance
            grad[i]=DMat_biv_var2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[2]==1){//first nugget
            grad[i]=DMat_biv_nug1(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[3]==1){//second nugget
            grad[i]=DMat_biv_nug2(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[4]==1){//colocated
            grad[i]=DMat_biv_col(h,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[5]==1){//scale1
            grad[i]=DMat_biv_scale1(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[6]==1){//scale12
            grad[i]=DMat_biv_scale12(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[7]==1){//scale2
            grad[i]=DMat_biv_scale2(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[8]==1){//smo1
            grad[i]=DMat_biv_smo1(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[9]==1){//smo12
            grad[i]=DMat_biv_smo12(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);i++;}
        if(flag[10]==1){//smo2
            grad[i]=DMat_biv_smo2(h,eps,var11,var22,nug11,nug22,scale11,scale22,scale12,smoo11,smoo22,smoo12,col,c11,c22);}
        break;
    }

}





// compute the gradient matrix (numcoord...) for the spatial field:
void DCorrelationMat_biv_tap(int *cormod,double *coordx, double *coordy, double *coordt,double *drho,double *eps,int *flagcor,int *nparcor,double *parcor,double *rho)
{
    int i=0,m=0,k=0;
    double *gradcor,*derho;
    //Initializes gradients:
    gradcor=(double *) R_alloc(*nparcor, sizeof(double));
     derho= (double *) Calloc(*npairs * *nparcor,double);
    k=0;
    for(i=0;i<*npairs;i++){
        GradCorrFct(rho[i],cormod,eps[0],flagcor,gradcor,lags[i],0,first[i],second[i],parcor);
        for(m=0;m<*nparcor;m++){derho[k]=gradcor[m];k++;}}
    k=0;
    for(i=0;i<*nparcor;i++){
        for(m=0;m<*npairs;m++){
            drho[k]=derho[i+m* *nparcor];k++;}}
     Free(derho);
    return;
}










void DCorrelationMat_biv(int *cormod,double *coordx, double *coordy, double *coordt,double *drho,double *eps,int *flagcor,
      int *nparcor,double *parcor,double *rho)
{
 int i=0,j=0,h=0,k=0,s=0,t=0,v=0,st=0,npa=0;
 double *gradcor,*derho;

 st=ncoord[0] * *ntime;
 npa=0.5*st*(st-1)+st;
 gradcor=(double *) R_alloc(*nparcor, sizeof(double));
 derho= (double *) Calloc(npa * *nparcor,double);

for(i=0;i<ncoord[0];i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<ncoord[0];j++){
  if(i==j){
    for(v=t;v<*ntime;v++){
          GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,0,0,t,v,parcor);
          h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}
     else {
          for(v=0;v<*ntime;v++){
               GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,mlags[i][j],0,t,v,parcor);
               h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}}}}

 k=0;
 for(i=0;i<*nparcor;i++)
   for(j=0;j<npa;j++){
     drho[k]=derho[i+j* *nparcor];
     k++;}
    Free(derho);
 return;
}




void DCorrelationMat_biv2(int *cormod,double *coordx, double *coordy, double *coordt,double *drho,double *eps,int *flagcor,
      int *nparcor,double *parcor,double *rho)
{
 int i=0,j=0,h=0,k=0,s=0,t=0,v=0,st=0,npa=0;
 double *gradcor,*derho,lags=0.0;

 st=ncoord[0] * *ntime;
 npa=0.5*st*(st-1)+st;
 gradcor=(double *) R_alloc(*nparcor, sizeof(double));
 derho= (double *) Calloc(npa * *nparcor,double);

for(i=0;i<ncoord[0];i++){
    for(t=0;t<*ntime;t++){
      for(j=i;j<ncoord[0];j++){
  if(i==j){
    for(v=t;v<*ntime;v++){
          GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,0,0,t,v,parcor);
          h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}
     else {
          lags=dist(type[0],coordx[i],coordx[j],coordy[i],coordy[j],*REARTH);
          for(v=0;v<*ntime;v++){
               GradCorrFct(rho[h],cormod,eps[0],flagcor,gradcor,lags,0,t,v,parcor);
               h++;
   for(s=0;s<*nparcor;s++){
     derho[k]=gradcor[s];
     k++;}}}}}}

 k=0;
 for(i=0;i<*nparcor;i++)
   for(j=0;j<npa;j++){
     drho[k]=derho[i+j* *nparcor];
     k++;}
    Free(derho);
 return;
}




// Computes the spatio-temporal variogram:
double Variogram(int *cormod, double h, double u, double nugget, double var, double *par)
{
  double vario=0.0;
  //Computes the variogram
  vario=nugget+var*(1-CorFct(cormod,h,u,par,0,0));
  return vario;
}


// Vector of spatio-temporal correlations:
void VectCorrelation(double *rho, int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model, 
                     double *nuis,double *par, double *u)
{
  int i,j,t=0;
  double ai=0.0,aj=0.0,p1=0.0,p2=0.0,psj=0.0;
  for(j=0;j<*nlagt;j++)
    for(i=0;i<*nlags;i++){
      if(*model==1||*model==10||*model==12||*model==21||*model==30||
      *model==22||*model==24|*model==26||*model==27||*model==29||*model==34||*model==38||*model==39) rho[t]=CorFct(cormod, h[i], u[j], par,0,0);  // gaussian 
      //if(*model==12)                      rho[t]=R_pow(CorFct(cormod, h[i], u[j], par,0,0),2);  // chisq case
     // if(*model==10)  // skew gaussian case
      //{
      //    cc=CorFct(cormod,h[i], u[j],par,0,0);
      //    rho[t]=((2*R_pow(nuis[3],2)/M_PI)*(sqrt(1-cc*cc) + cc*asin(cc)-1) + cc*nuis[2])/(nuis[2]+R_pow(nuis[3],2)*(1-2/M_PI));
     // }
      /***************************************/
      if(*model==2||*model==11||*model==19||*model==14)   // binomial  type I or  II or geometric case
      {
          ai=mean[i];aj=mean[j];
          psj=pbnorm(cormod,h[i],u[j],ai,aj,nuis[0],nuis[1],par,0);
          p1=pnorm(ai,0,1,1,0); p2=pnorm(aj,0,1,1,0);

          if(*model==2||*model==11||*model==19) rho[t]=(psj-p1*p2)/sqrt(p1*p2*(1-p1)*(1-p2));  // binomyal type I II
          if(*model==14)                        rho[t]= ((psj-p1*p2)/((-psj+p1+p2)*p1*p2))/(sqrt((1-p1)*(1-p2))/(p1*p2));  //covar12/sqrt(var1var2)
      }
      /***************************************/
      t++;}
  return;
}



// Vector of spatio-temporal correlations:
void VectCorrelation_biv(double *rho, double *vario,int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model, 
                     double *nuis,double *par, double *u)
{
    int i,j,p,t=0;
    for(j=0;j<2;j++)
    for(p=0;p<2;p++)
    for(i=0;i<*nlags;i++){
        rho[t]=CorFct(cormod, h[i],0, par,j,p);
        vario[t]=CorFct(cormod, 0,0, par,j,p)-CorFct(cormod, h[i],0, par,j,p);
        t++;}
    return;
}



//  Bounds for the bivariate matern
double bi_matern_bounds(double scale11,double scale22,double scale12,double nu11,double nu22,double nu12,double t,int c){
  double scale,l,inf;
  scale=R_pow(scale11,2*nu11)*R_pow(scale22,2*nu22)/  R_pow(scale12,4*nu12);
  if(!c) inf=R_pow(R_pow(scale12,2)+t*t,(2*nu12+2))/(R_pow((R_pow(scale11,2)+t*t),(nu11+1))*R_pow((R_pow(scale22,2)+t*t),(nu22+1)));
  else  inf=1;
  l=(gammafn(nu11+1)*gammafn(nu22+1)*R_pow(gammafn(nu12),2)*scale*inf)/(gammafn(nu11)*gammafn(nu22)*R_pow(gammafn(nu12+1),2)); 
  return(l);
}

