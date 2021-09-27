#include "header.h"



void spectral_density(int *L,int *model,int *p, double *matrix ,double *matrix_out, 
                                double *C, double *a, double *nu1,double *Cg, double *ag, double *nu1g){
  int n_rows = (*L);
  int i,j,id=0,ig=n_rows,ih=0;
  double pi = acos(-1.0);
  int m = (*model);
  
  if(m==0){
    for(i=0;i<n_rows;i++){
      double norma2 = pow(matrix[id],2)+pow(matrix[ig],2);
      for(j=0;j<(*p);j++){
        if((a[j] > 0) && (nu1[j] > 0)){
        
          double rg1 = 2*log(2*pi*ag[0]);
          double rg2 = lgamma(nu1g[0]+1);
          double rg3 = - lgamma(nu1g[0]);
          double rg4 = - log(pi);
          double rg5 = - (nu1g[0]+1)*log(1+(pow(2*pi*ag[0],2))*norma2);
          double lnfg = rg1 + rg2 + rg3 + rg4 + rg5;
          
          double r1 = 2*log(2*pi*a[j]);
          double r2 = lgamma(nu1[j]+1);
          double r3 = - lgamma(nu1[j]);
          double r4 = - log(pi);
          double r5 = - (nu1[j]+1)*log(1+(pow(2*pi*a[j],2))*norma2);
          double lnf = r1 + r2 + r3 + r4 + r5;
          
          matrix_out[ih] = 2*(C[j]*exp(lnf)/ (Cg[0]*exp(lnfg)));
          ih = ih+1;
        }
        else {
          Rprintf("At least one parameter does not satisfy the model validity restrictions");
        }
      }
      id = id+1;
      ig = ig+1;
    }
  }
}

void matrix_temp(int *N ,double *matrix, double *l1 ,double *l2 ,double *v11 ,double *v21,double *v12,double *v22) {
  int i,j,id=0;
  for(i=0;i<(*N);i++){
    double a11 = l1[i]*pow(v11[i],2) + l2[i]*pow(v12[i],2);
    double a12 = l1[i]*v11[i]*v21[i] + l2[i]*v22[i]*v12[i] ;
    double a21 = l1[i]*v11[i]*v21[i] + l2[i]*v22[i]*v12[i];
    double a22 = l1[i]*pow(v12[i],2) + l2[i]*pow(v22[i],2);
    
    for(j=0;j<4;j++){
      if ( j==0 ) {                 
        matrix[id] = a11;
        id=id+1;
      }
      else if (j == 1) {
        matrix[id] = -a12;
        id=id+1;
      }
      else if (j == 2) {
        matrix[id] = -a21;
        id=id+1;
      }
      else if (j == 3) {
        matrix[id] = a22;
        id=id+1;
      }
    }  
  }
}


void vector_to_select(int *N, double *matrix) {
  int i,id=1;
  int a = 3;
  int b = 1;
  matrix[0]=1;
  for(i=1;i<=(*N);i++){
    int aux = id-1;
    int val = (i-1)%2;
    if ( val==0 ) {                 
      matrix[id] = matrix[aux]+a;
      id=id+1;
    }
    else if (val == 1) {
      matrix[id] = matrix[aux]+b;
      id=id+1;
    }
  }
}

void simu_on_coords(int *Ndim,int *Mcoords,int *Mu,double *coords,double *amatrix, 
                    double *matrix_phi,double *matrix_u,double *matrix_out){
  double pi = acos(-1.0);
  
  int n_rows_coords = (*Ndim);
  int n_rows_u = (*Mu);
  
  int i,j,ih=0,ip=0,row0 = 0;
  
  int ig=n_rows_coords,it=n_rows_u;
  int row1 = n_rows_coords;
  
  for(i=0;i<(*Ndim);i++){
    
    int amatrix_index0 = 0;
    int amatrix_index1 = 1;
    
    for(j=0;j<(*Mu);j++){
      
          int jaux = j;
          double val_phi = matrix_phi[jaux];
          double val_a0 = amatrix[amatrix_index0];
          double val_a1 = amatrix[amatrix_index1];
      
          double val1 = coords[ih];
          double val2 = coords[ig];
          double mul1 = val1*matrix_u[ip];
          double mul2 = val2*matrix_u[it];
        
          matrix_out[row0] = matrix_out[row0]+(val_a0*cos(2*pi*(mul1 + mul2)+val_phi));
          matrix_out[row1] = matrix_out[row1]+(val_a1*cos(2*pi*(mul1 + mul2)+val_phi));
          
          if((j+1)==(*Mu)){
            ip=0;
            it=n_rows_u;
          }
          else{
            ip=ip+1;
            it=it+1;
          }
        amatrix_index0 = amatrix_index0+2;
        amatrix_index1 = amatrix_index1+2;
    }
    row0=row0+1;
    row1=row1+1;
    ih=ih+1;
    ig=ig+1;
  }
}