#include "header.h"








void spectraldensityC(double u,int model,int d,int L,double *f,double *av,double *Cv,double *nu1v,double *nu2v){
    double norm_u = u;
    int h=0;
    for (int i=0;i<L*L;i++){if (av[i]>0){h = h+1;}}

    int positive_a;
    if (h==L*L){positive_a = 0;}
    
    h=0;
    int pos_nu1;
    for (int i=0;i<L*L;i++){if (av[i]>0){h = h+1;}}
    if (h==L*L){pos_nu1 = 0;}

    double *r = (double *) R_Calloc(L*L,double);
   
    if(model==14){   //Matern
        if((pos_nu1==0) &&(positive_a == 0 )){
            for (size_t i = 0; i < L*L; i++){
                double v = d/2;
                r[i] =  d*log(2*PI*av[i]) + lgamma(nu1v[i]+v) - lgamma(nu1v[i]) - v*log(PI) - (nu1v[i]+v)*log(1+pow(2*PI*av[i]*norm_u,2));
                r[i] =Cv[i]*exp(r[i]);
            }}else{Rprintf("Parameter does not satisfy the model validity restriction");}

    }
    //Wendland
    if (model==19){if(positive_a==0 && d==2){} }
    for (size_t i = 0; i < L*L; i++) {f[i] = r[i];}
    R_Free(r);

}



// preguntar si las coordenadas son enteras
void extraer(double *coord,int sequen1,double *sub_coord,int fila,int col, int d){
    int j=0;
    for (size_t i = sequen1; i < fila; i++){
        for (int h=0;h<d;h++){
            sub_coord[j+h*col] = coord[i+h*col]; 
             //Rprintf("%d %d %f %f --%d %d  %d\n",i,h,sub_coord[j+h*col] ,coord[i+h*col],d,fila,sequen1);
        }
        j=j+1;
   
    }
}


void rellenar_indice(int *index,int inicio, int final,int largo){
    for (size_t i=0;i<largo;i++){index[i] = inicio+i;}
}


void u_index_extraer(double *u_index,double *u, int *index,int largo,int d,int filas){
    for (int i=0;i<largo;i++){
        for (int h=0;h<d;h++){
            u_index[i+h*filas] = u[index[i]+h*filas]; }}
}


void mult_mat(double *z, double *x, int xrows, int xcols, double *y, int yrows, int ycols){ 
  char *notrans = "N";
  double one = 1.0, zero = 0.0, *tmp = NULL;
  tmp = (double *) R_Calloc(xrows * ycols, double);
  F77_CALL(dgemm)(notrans, notrans, &xrows, &ycols, &xcols, &one, x, &xrows, y, &yrows, &zero, tmp, &xrows FCONE FCONE);
  Memcpy(z, tmp, xrows * ycols);
  R_Free(tmp);
}

void tcrossprod(double *z, double *x,  int xrows, int xcols, double *y, int yrows, int ycols){ /* outer product of two given matrices. z <- x %*% t(y) */
  char *notrans = "N", *trans = "T";
  double one = 1.0, zero = 0.0, *tmp = NULL;
  /* use tmp so z can be either x or y */
  tmp = (double *) R_Calloc(xrows * yrows, double);
  F77_CALL(dgemm)(notrans, trans, &xrows, &yrows, &xcols, &one, x, &xrows, y, &yrows, &zero, tmp, &xrows FCONE FCONE);
  Memcpy(z, tmp, xrows * yrows);
  R_Free(tmp);
}


void mult_x_cons(double *x, double cte,int largo){
  int i; 
  for (i=0;i<largo;i++){x[i] = cte*x[i];}
}

void sumar_matrices(double *x0, double *x,double *y,int largo){
  int h; 
  for (h=0;h<largo;h++){x0[h] = x[h]+y[h];}
}

void restar_matrices(double *x0, double *x,double *y,int largo){
    int h; 
    for (h=0;h<largo;h++){x0[h] = x[h]-y[h];}
}


void cos_vec(double *x_cos,double *x,int largo){
  int h;
  for (h=0;h<largo;h++){x_cos[h] = cos(x[h]);}
}

void sen_vec(double *x_sen,double *x,int largo){
    int h;
    for (h=0;h<largo;h++){x_sen[h] = sin(x[h]);}
}

void llenar_simu(double *x,double *simu, int N,int *P, int m){

    int j=0;
    for (int i = N* *(P)*m;i<(N+1)* *(P)*m;i++){simu[i] = x[j];j = j+1;}

}
void extraer_col(int inicio, int final,double *x_original,double *x){
    int j=0;
    for(int i=inicio;i<=final;i++){x[j] = x_original[i];j=j+1;}
}

void llenar_simu1(double *simu1,double *simu,int *m,int *P,int *N,int lim, int i,double *L1){

    if (i==0){
        for (int j = 0; j < m[i] * *(P)* *(N)*lim; j++){simu1[j] = simu[j]/sqrt(*L1);}  
    }
    else{
        int k=0;
        int suma = m[i-1] * *(P)* *(N)*lim + m[i] * *(P)* *(N)*lim;
        for (size_t j = m[i] * *(P)* *(N)*lim ; j < suma ; j++){
            simu1[j]  = simu[k];
            k = k+1;}
    }
}


void C_tcrossprod(Rcomplex *z, Rcomplex *x,  int xrows, int xcols, Rcomplex *y, int yrows, int ycols){ /* outer product of two given matrices. z <- x %*% t(y) */
  char *notrans = "N", *trans = "T";
  Rcomplex one, zero, *tmp = NULL;
  one.r = 1;
  one.i = 0;
  zero.i =0;
  zero.r = 0;


      tmp = (Rcomplex *) R_Calloc(xrows * ycols, Rcomplex);
  F77_CALL(zgemm)(notrans, trans, &xrows, &yrows, &xcols, &one, x, &xrows, y, &yrows, &zero, tmp, &xrows FCONE FCONE);
  Memcpy(z, tmp, xrows * yrows);
  R_Free(tmp);
}


void C_mult_mat(Rcomplex *z, Rcomplex *x, int xrows, int xcols, Rcomplex *y, int yrows, int ycols){ 
  char *notrans = "N";
  Rcomplex one, zero, *tmp = NULL;
  one.r = 1;
  one.i = 0;
  zero.i =0;
  zero.r = 0;

   tmp = (Rcomplex *) R_Calloc(xrows * ycols, Rcomplex);
  F77_CALL(zgemm)(notrans, notrans, &xrows, &ycols, &xcols, &one, x, &xrows, y, &yrows, &zero, tmp, &xrows FCONE FCONE);
  Memcpy(z, tmp, xrows * ycols);
  R_Free(tmp);
}


void for_c(int *d_v,double *a_v,double *nu1_v,double *C_v,double *nu2_v,int *P, int *N, int *L,int *model,double *u,double *a0,double *nu0,double *A,double *B,
    int *sequen,int *largo_sequen,int *n,double *coord,double *phi, int *vtype,int *m1,double *simu1,double *L1){
    


/*Rprintf("  %f %f %f %f %f %d %d %d %f %f %f\n ", *u, *a0, *nu0, *A, *B,
     *sequen, *largo_sequen, *n,*simu1, *L1);*/


    // Iteraciones del for
    int prod = *P * *L * *N ;
    int d = *d_v;int p = *P;
    //variables utiles

    int model_1 = 14;

   double *f,*g,*u_1,*w,*trabajo_adicional;
   Rcomplex *autovalores,*autovectores;
   
       f = (double *) R_Calloc(p * p, double);
        g = (double *) R_Calloc(1, double);
         u_1 = (double *) R_Calloc(d, double);
          w = (double *) R_Calloc(p , double);

    autovalores = (Rcomplex *) R_Calloc(p, Rcomplex);
    autovectores = (Rcomplex *) R_Calloc(p*p,Rcomplex);
    trabajo_adicional = (double *) R_Calloc(3*p-1,double);



    int lwork = 3 * p - 1;
    int info;
    char *J = "V";
    char *U = "L";
    size_t N1= p;
    Rcomplex *vdm;
    vdm= (Rcomplex * )R_Calloc(p*p,Rcomplex);
    Rcomplex raizv[2];
    for (int l= 0;l<prod;l++){
        int h=0;
        for(int i=0;i<d;i++){u_1[i] = u[l+h];h = prod;}


        int g1 = 1;
        double norm_u = F77_CALL(dnrm2)(&d, u_1, &g1);

        spectraldensityC(norm_u,*model,d,p,f,a_v,C_v,nu1_v,nu2_v);
        double C0[] = {1};
        spectraldensityC(norm_u,model_1,2,1,g,a0,C0,nu0,nu2_v);
        for (int i =0;i<p*p;i++){f[i] = 2*f[i]/g[0];}
        if (p>=2){
            F77_CALL(dsyev)(J, U, &p, f, &p, w, trabajo_adicional, &lwork, &info,N1,N1);
            for (int i=0;i<p;i++){
                autovalores[i].r = w[p-i-1];
                autovalores[i].i = 0;
            }
            // modificar esta parte si se quiere para p >2
            autovectores[0].r = f[2];
            autovectores[1].r = f[3];
            autovectores[2].r = f[0];
            autovectores[3].r = f[1];

            autovectores[0].i = 0;
            autovectores[1].i = 0;
            autovectores[2].i = 0;
            autovectores[3].i = 0;
        }
        else{
            autovalores->r = *f;
            autovalores->i =0;
            autovectores->r = 1;
            autovectores->i = 0;
        }
        
        for (int i=0;i<p;i++){
            if (autovalores[i].r<0){
                autovalores[i].r =0;
            }
        }

        //Calculos previo para lo de A y B
        if(p==1){
            Rcomplex raiz;
            if (autovalores->r>=0){raiz.r =sqrt(autovalores->r);raiz.i =0;}
            else{raiz.i= sqrt(-1*autovalores->r);raiz.r =0;}

            A[l] = raiz.r;B[l] = raiz.i;
        }
        else{
            
            for (size_t i = 0; i < p; i++){
                if (autovalores[i].r>=0){raizv[i].r = sqrt(autovalores[i].r);raizv[i].i =0;}
                else{
                    raizv[i].i = sqrt(-1*autovalores[i].r);raizv[i].r =0;}
            }
            vdm[0].r = raizv[0].r;
            vdm[0].i = raizv[0].i;
            vdm[1].r = 0;
            vdm[2].r = 0;
            vdm[1].i = 0;
            vdm[2].i = 0;
            vdm[3].r = raizv[1].r;
            vdm[3].i = raizv[1].i;
            C_tcrossprod(vdm,vdm,2,2,autovectores,2,2);
            C_mult_mat(vdm,autovectores,2,2,vdm,2,2);
            double fo = (l+1-1)/p;
            int j = l+1 - p*floor(fo)-1;           //(h,l)
            double V; double imag_;
            for (int h=0;h<p;h++){
                if (j==0){ V = vdm[h].r;imag_ = -vdm[h].i;}
                else{V = vdm[2+h].r;imag_ = -vdm[2+h].i;}
                A[h+p*l] = V;B[h+p*l] = imag_;
            }
        } 
    }
    R_Free(f);
    R_Free(g);
    R_Free(u_1);
    R_Free(autovalores);
    R_Free(autovectores);
    R_Free(w);
    R_Free(trabajo_adicional);
    R_Free(vdm);
     //El otro for

    int lim = *(largo_sequen)-1;
    for (int i=0;i < lim;i++){
        int fila = sequen[i+1]-sequen[i];
        int col = *(n);
        double *sub_coord;
        sub_coord = (double *) R_Calloc(fila*d,double);  
        int sequen1 = sequen[i];
        extraer(coord,sequen1,sub_coord,fila,col,d);
        int m = fila;
        double *simu;
        simu = (double *) R_Calloc(m * *(P)* *(N),double);  
        for (size_t k = 0; k < *(N); k++){
            int *index;
            int largo = (k+1)**(L)* *(P) - k* *(L)* *(P);
            index =  (int *) R_Calloc(largo,int); 
            rellenar_indice(index,k* *(L)* *(P), (k+1)**(L)* *(P) ,largo);
            double *ui;
            ui =  (double *) R_Calloc(largo*d,double); 
            u_index_extraer(ui,u,index,largo,d,prod);
            double *x0;
            x0 = (double *) R_Calloc(fila*largo,double); 
            tcrossprod(x0,sub_coord,fila,d,ui,largo,d);
            double cte =2*PI;
            mult_x_cons(x0,cte,fila*largo);
            double *phi_cross;
            phi_cross= (double*) R_Calloc(largo,double); 
            u_index_extraer(phi_cross,phi,index,largo,1,prod);
            double *cbind;
            cbind = (double*) R_Calloc(m,double); 
            for (int h=0;h<m;h++){
                cbind[h] = 1;
            }
            double *x1;
            x1 = (double *) R_Calloc(fila*largo,double); 
            tcrossprod(x1,cbind,m,1,phi_cross,largo,1);
            double *x;
            x = (double *) R_Calloc(fila*largo,double);
            sumar_matrices(x,x1,x0,fila*largo);
            if (*vtype==0){
                if (*P==1){
                    double *x_cos; 
                    x_cos= (double *) R_Calloc(fila*largo,double);
                    cos_vec(x_cos,x,fila*largo);
                    double *x_sen;
                    x_sen = (double *) R_Calloc(fila*largo,double);
                    sen_vec(x_sen,x,fila*largo);
                    double *A_index;
                    A_index = (double *) R_Calloc(fila*largo,double);
                    double *B_index;
                    B_index = (double *) R_Calloc(fila*largo,double);
                    u_index_extraer(A_index,A,index,largo,1,prod);
                    u_index_extraer(B_index,B,index,largo,1,prod);

                    double *result1;
                    result1 = (double *) R_Calloc(m,double); 
                    mult_mat(result1,x_cos,m,largo,A_index,largo,1);

                    double *result2;
                    result2 = (double *) R_Calloc(m,double); 
                    mult_mat(result2,x_sen,m,largo,B_index,largo,1);
                    sumar_matrices(result2,result1,result2,m);
                    llenar_simu(result2,simu,k,P,m);
                    R_Free(x_cos);
                    R_Free(x_sen);
                    R_Free(result1);
                    R_Free(result2);
                    R_Free(A_index);
                    R_Free(B_index);
                }
                else{
                    double *x_cos;
                    x_cos = (double *) R_Calloc(fila*largo,double);
                    cos_vec(x_cos,x,fila*largo);
                    double *x_sen;
                    x_sen= (double *) R_Calloc(fila*largo,double);
                    sen_vec(x_sen,x,fila*largo);
                    double *A_col;
                    A_col = (double *) R_Calloc(*(P)*(index[largo-1]+1),double);
                    double *B_col;
                    B_col = (double *) R_Calloc(*(P)*(index[largo-1]+1),double);
                    int p = *P;
                    extraer_col(index[0],index[largo-1]*p+1,A,A_col);
                    extraer_col(index[0],index[largo-1]*p+1,B,B_col);
                    double *result1;
                    result1 = (double *) R_Calloc(m * *(P),double); 
                    tcrossprod(result1,x_cos,m,largo,A_col,p,largo);
                    double *result2; 
                    result2= (double *) R_Calloc(m * *(P),double); 
                    tcrossprod(result2,x_sen,m,largo,B_col,p,largo);

                    sumar_matrices(result2,result1,result2,m*p);
                    llenar_simu(result2,simu,k,P,m);

                    R_Free(result1);
                    R_Free(result2);
                    R_Free(A_col);
                    R_Free(B_col);
                    R_Free(x_cos);
                    R_Free(x_sen);
                }
            }
            else{
                double *cosphi;
                cosphi = (double *) R_Calloc(fila*largo,double);
                double *senphi;
                senphi = (double *) R_Calloc(fila*largo,double);
                cos_vec(cosphi,x1,fila*largo);
                sen_vec(senphi,x1,fila*largo);
                if (*P == 1){
                    Rprintf("d");
                }
                else{
                    double *x_cos;
                    x_cos = (double *) R_Calloc(fila*largo,double);
                    cos_vec(x_cos,x,fila*largo);
                    double *x_sen;
                    x_sen = (double *) R_Calloc(fila*largo,double);
                    sen_vec(x_sen,x,fila*largo);

                    double *A_col;
                    A_col = (double *) R_Calloc(*(P)*(index[largo-1]+1),double);
                    double *B_col; 
                    B_col= (double *) R_Calloc(*(P)*(index[largo-1]+1),double);
                    int p = *P;
                    extraer_col(index[0],index[largo-1]*p+1,A,A_col);
                    extraer_col(index[0],index[largo-1]*p+1,B,B_col);

                    restar_matrices(x_cos,x_cos,cosphi,fila*largo);
                    restar_matrices(x_sen,x_sen,senphi,fila*largo);
                    double *result1;
                    result1 = (double *) R_Calloc(m * *(P),double); 
                    tcrossprod(result1,x_cos,m,largo,A_col,p,largo);
                    double *result2; 
                    result2= (double *) R_Calloc(m * *(P),double); 
                    tcrossprod(result2,x_sen,m,largo,B_col,p,largo);
                    sumar_matrices(result2,result1,result2,m*p);
                    llenar_simu(result2,simu,k,P,m);
                    R_Free(result1);
                    R_Free(result2);
                    R_Free(A_col);
                    R_Free(B_col);
                    R_Free(x_sen);
                    R_Free(x_cos);
                }
                R_Free(cosphi);
                R_Free(senphi);
            }
            R_Free(phi_cross);
            R_Free(cbind);
            R_Free(x);
            R_Free(x1);
            R_Free(x0);
            R_Free(ui);
            R_Free(index);
        }
        //simu1 bla bla
        llenar_simu1(simu1,simu,m1,P,N,lim,i,L1);

        R_Free(simu);
        R_Free(sub_coord);
    }
}

