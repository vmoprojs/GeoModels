//#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <stdint.h>
#define LOW -1.0e15
#define MAXERR 1e-6
#define EPS DBL_EPSILON
/* for bivariate t distribution */
#define EPS1 1.0e-10
#define EPS2 1.0e-10
#define Inf 1.7976931348623157e+308
#define ETHRESH 1.0e-12
#define MACHEP   1.11022302462515654042E-16   /* 2**-53 */
#define MAXNUM   1.79769313486231570815E308    /* 2**1024*(1-MACHEP) */



//---------START GLOBAL VARIABLES-----------

double **dista;// 2x2 matrix of distance weight (for c)and commpact suppots (tap)
int *first;//vector of index in the bivariate case
int *second;//vector of index in the bivariate case
int *isbiv;//is bivariate?
int *ismem;//is with memoty allocation
int *isst;//is a spatio-temporal random field?
int *istap;//is tapering?
double *lags;// vector of spatial distances for tapering
double *lagt;// vector of temporal distance for tapering
double **mlags;// vector of spatial distances
double **mlagt;// vector of temporal distances
double *maxdist;// the threshould of the spatial distances
double *maxtime;// the threshould of the temporal distances below which the pairs are considered
double *maximdista;// the maximum spatial distance
double *maximtime;// the maximum temporal distance
double *minimdista; // the minimum spatial distance
double *minimtime;// the minimum temporal distance
int *ncoord;// number of total spatial coordinates
int *ncoordx;// number of the first spatial coordinates
int *ncoordy;// number of the second spatial coordinates
int *npairs;// effective number of pairs
int *nrep;// number of iid replicates of the random field
int *ntime;// number of times
double *REARTH; // radius of the sphere
double *tapsep; // parameter separability for space time quasi taper
double *tlags;
double *tlagt;
int *tfirst;
int *tsecond;
int *type;  //  type of distance
int *cdyn; // dynamic coords indicator
//char *fkernel; // path to kernel binary
//int *totpairs;// total number of pairs
//---------END GLOBAL VARIABLES-------------
//void indx(double *ind, int *n);

// fortran declaration for bivariate normal cdf:
double F77_NAME(chfm)(double *xr,double *xi,double *ar, double *ai,
  double *br, double *bi,double *r, double *ri,int *len, int *lnchf,int *ip);

double F77_NAME(bvnmvn)(double *lower, double *upper, int *infin, double *correl);

void mult_pmnorm( int *nvar , double *lower , double *upper , int *infin , double *corr , int *maxpts , double *abseps ,
     double *releps , double *esterror , double *result , int *fail );

void F77_NAME(sadmvn)( int* , double* , double*, int* , double* , int* , double* , double* , double* , double * , int*) ;

double F77_NAME(mvndst)(int *N,double *lower, double *upper, int *infin, double *correl);
// Internal function declaration:
// 1)
/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
Start
 ---------------------------------------------------------------*/
double Hyp_conf_laplace_approx(double a, double b, double z);
double asym_aprox_1F1(double a, double b, double z);
int fmax_int(int u,int v);
int fmin_int(int u,int v);
//double hypergeo(double aux, double *param);

double bi_matern_bounds(double scale11,double scale22,double scale12,double nu11,double nu22,double nu12,double t,int c);

double biv_binom (int NN, int u, int v, double p01,double p10,double p11);

double  biv_binom2(int NN_i,int NN_j, int k, int u, int v, double p01,double p10,double p11);

double biv_wrapped(double alfa,double u, double v, double mi, double mj, double nugget,double sill,double corr);

double biv_Weibull(double corr,double zi,double zj,double mui, double muj, double shape);
double  biv_Weibull2(double rho12,double zi,double zj,double mi,double mj, double shape1,double shape2);

double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape);
double biv_gamma2(double corr,double zi,double zj,double mui, double muj, double shape);
double biv_Kumara(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2);

//double log_biv_binom (int NN, double u, double v, double psm,double psj);
double biv_LogLogistic(double corr,double zi,double zj,double mui, double muj, double shape);
double biv_Logistic(double corr,double zi,double zj,double mui, double muj, double sill);

double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11);
double bin_aux(int a,int NN,int u,int v,double p1, double p2,double p11);
double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11);
double aux_biv_binomneg_simple(int NN, int u, double p01,double p10,double p11);
double aux_euv_binomneg (int N, double p1,double p2,double p11);
double corr_binomneg (int N, double p1,double p2,double p11);
double biv_poisbin (int NN,  int u, int v, double p01,double p10,double p11);
double biv_poisbinneg(int NN, int u, int v, double p01,double p10,double p11);


double biv_skew(double corr,double zi,double zj,double mi,double mj,double vari,double skew);

double biv_sinh(double corr,double zi,double zj,double mi,double mj,double skew,double tail,double vari);
//bivariate skew gaussian bivariate  distribution
double biv_skew2(double corr,double zi,double zj,double vari1,double vari2,double rho ,
         double skew1,double skew2);

double triv_skew(double x,double c_0i,double c_0j, double rho,double data_i,double data_j,double *nuis);


double LambertW(double z);
double inverse_lamb(double x,double tail);
double dbnorm(double x_i,double x_j,double mean_i,double mean_j,double sill,double corr);
double biv_tukey_h(double data_i, double data_j, double mean_i, double mean_j, double tail, double sill, double corr);


double cond_exp_skew(double c_0i,double c_0j,double rho, double data_i,double data_j,double *nuis);
double cond_exp_bin(int *cormod,double data_i,double data_j,double lags_i,double lags_j,double lags,double *nuis,double *par,double psm);


double marg_binom(int n,double x,double p);

double marg_binomneg(int n,int x,double p);

double marg_pois(int n,double x,double p);

double marg_geom(int x,double p);

double marg_p(double categ_0,double psm,int *model,int n);

double CheckCor(int *cormod, double *par);

void Comp_supp(double *c_supp,int *cormod, double h,double u, double *par);
double CorFct(int *cormod, double h, double u, double *par, int c11, int c22);
double CorFunCauchy(double lag, double power2, double scale);
double CorFunDagum(double lag, double power1, double power2, double scale);
double CorFunGenCauchy(double lag, double power1, double power2, double scale);
double CorFunGenCauchy2(double lag, double power1, double power2, double scale);
double CorFunWitMatCau(double h, double scale12,double smo12);
double Shkarofski(double lag, double a,double b, double k);
double CorFunSmoke(double h, double  scale,double  smooth);
double CorFunGenWitMatCau(double h, double scale,  double smoo,double beta);
double CorFunSferical(double lag, double scale);
double CorFunStable(double lag, double power, double scale);
double CorFunWave(double lag, double scale);
double CorFunWitMat(double lag, double scale, double smooth);
double CorFunWitMat1(double lag, double scale, double smooth);
double CorFunBohman(double lag,double scale);
double CorFunW0(double h,double scale,double power);
double CorFunW1(double h,double scale,double power);
double CorFunW2(double h,double scale,double power);
double CorFunW_gen(double h, double power1, double smooth, double scale);
double CorFunWend1(double lag,double scale);
double CorFunWend2(double lag,double scale);
double CorFunWend3(double lag,double scale);
double CorFunWend4(double lag,double scale);
double CorFunWend5(double lag,double scale);
double CorFunWendhole3(double lag,double scale);
double CorFunWendhole2(double lag,double scale);
double CorFunWendhole1(double lag,double scale);
double CorFunWendhole(double lag,double scale);
double CorFunWend0_tap(double lag,double scale,double smoo);
double CorFunWend1_tap(double lag,double scale,double smoo);
double CorFunWend2_tap(double lag,double scale,double smoo);

/****** derivative **************/
double DCauchyPow(double power2, double scale, double rho);
double DCauchySc(double lag, double power2, double scale, double rho);
double DExpoSc(double lag, double scale, double rho);
double DGaussSc(double lag, double scale, double rho);
double DGenCauP1(double lag, double power1, double power2, double scale, double rho);
double DGenCauP2(double lag, double power1, double scale, double rho);
double DGenCauSc(double lag, double power1, double power2, double scale, double rho);
double DSferiSc(double lag, double scale);
double DStabPow(double lag, double power, double scale, double rho);
double DStabSc(double lag, double power, double scale, double rho);
double DWaveSc(double lag, double scale);
double DWhMatSc(double eps, double lag, double scale, double smooth);
double DWhMatSm(double eps, double lag, double scale, double smooth);
double Dwen1Sc(double lag, double scale);
double DWen1Sc(double lag, double scale, double smooth);
double DWen1Sm(double lag, double scale, double smooth);
double DGneiting_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_32sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_32sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_32pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_32pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_32sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC2_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC2_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC2_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC2_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_GC2_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DIaco_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);
double DIaco_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);
double DIaco_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);
double DIaco_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);
double DIaco_pw2(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double power2);
double DPorcu_sc_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DPorcu_sc_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DPorcu_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DPorcu_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DPorcu_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DStein_sc_s(double h,double u, double power_t,double scale_s,double scale_t,double smooth);
double DStein_sc_t(double h,double u, double power_t,double scale_s,double scale_t,double smooth);
double DStein_pw_s(double h,double u, double power_t,double scale_s,double scale_t,double smooth);
double DStein_pw_t(double h,double u, double power_t,double scale_s,double scale_t,double smooth);
double DStein_sm(double h,double u,   double power_t,double scale_s,double scale_t,double smooth);
double DExp_Cauchy_sc_s(double h,double u,double scale_s,double scale_t,double power2);
double DExp_Cauchy_sc_t(double h,double u,double scale_s,double scale_t,double power2);
double DExp_Cauchy_pw2 (double h,double u,double scale_s,double scale_t,double power2);
double DExp_Exp_sc_s(double h,double u,double scale_s,double scale_t);
double DExp_Exp_sc_t(double h,double u,double scale_s,double scale_t);
double DMat_Exp_sc_s(double h,double u,double scale_s,double scale_t,double smooth);
double DMat_Exp_sc_t(double h,double u,double scale_s,double scale_t,double smooth);
double DMat_Exp_sm(double h,double u,double eps,double scale_s,double scale_t,double smooth);
double DMat_Cauchy_sc_s(double h,double u,double power2,double scale_s,double scale_t,double smooth);
double DMat_Cauchy_sc_t(double h, double u,double power2,double scale_s,double scale_t,double smooth);
double DMat_Cauchy_pw2(double h,double u,double power2,double scale_s,double scale_t,double smooth);
double DMat_Cauchy_sm(double h,double u,double eps, double power2,double scale_s,double scale_t,double smooth);

double DExpsep_biv_var11(double h,double var11,double var22,double nug11,double nug22, double scale, double col,int c11,int c22);
double DExpsep_biv_var22(double h,double var11,double var22,double nug11,double nug22, double scale, double col,int c11,int c22);
double DExpsep_biv_nug11(double h,double var11,double var22,double nug11,double nug22, double scale, double col,int c11,int c22);
double DExpsep_biv_nug22(double h,double var11,double var22,double nug11,double nug22, double scale, double col,int c11,int c22);
double DExpsep_biv_scale(double h,double var11,double var22,double nug11,double nug22, double scale, double col,int c11,int c22);
double DExpsep_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale, double col ,int c11,int c22);
double DMat_biv_scale1_contr(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_scale2_contr(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                             double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double Dmatsep_biv_var11(double h,double var11,double var22,double nug11,double nug22, double scale,double smooth, double col,int c11,int c22);
double Dmatsep_biv_var22(double h,double var11,double var22,double nug11,double nug22, double scale,double smooth, double col,int c11,int c22);
double Dmatsep_biv_nug11(double h,double var11,double var22,double nug11,double nug22, double scale,double smooth, double col,int c11,int c22);
double Dmatsep_biv_nug22(double h,double var11,double var22,double nug11,double nug22, double scale,double smooth, double col,int c11,int c22);
double Dmatsep_biv_scale(double h,double eps,double var11,double var22,double nug11,double nug22, double scale,double smooth, double col,int c11,int c22);
double Dmatsep_biv_smo(double h,double eps,double var11,double var22,double nug11,double nug22, double scale, double smooth,double col,int c11,int c22);
double Dmatsep_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale,double smooth, double col ,int c11,int c22);


double DMat_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_scale1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_scale12(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11,double eps, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_scale2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_smoo1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_smoo12(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DMat_biv_smoo2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22);



/***************************************************/
/* derivative of LMC  (contr) model */
/***************************************************/

double DLMC_contr_var1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22);
double DLMC_contr_var2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22);
double DLMC_contr_nug1(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22);
double DLMC_contr_nug2(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22);
double DLMC_contr_col(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22);
double DLMC_contr_scale11(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22);
double DLMC_contr_scale22(double h,double eps,double var11,double var22,double nug11,double nug22, double scale11,double scale22,double col,int c11,int c22);


/***************************************************/
/* derivative of bivariate wendland2 model */
/***************************************************/
double DWen1_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                     double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                    double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_scale1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_scale12(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                        double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_scale2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_smoo1(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                      double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_smoo12(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                       double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
double DWen1_biv_smoo2(double h,double var11,double var22,double nug11,double nug22, double scale11,double scale22, double scale12,
                        double smoo11, double smoo22,double smoo12, double col,int c11,int c22);
                      

/***************************************************/
/* derivative of bivariate separable wendland2 model */
/***************************************************/
double DWen1sep_biv_var1(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22);
double DWen1sep_biv_var2(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22);
double DWen1sep_biv_nug1(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22);
double DWen1sep_biv_nug2(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22);
double DWen1sep_biv_scale(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22);
double DWen1sep_biv_col(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22);
double DWen1sep_biv_smoo(double h,double var11,double var22,double nug11,double nug22, double scale, double col,double smoo,int c11,int c22);



void GradCorrFct(double rho, int *cormod, double eps, int *flag,
		 double *grad, double h, double u, int c11, int c22, double *par);
void GradVarioFct(double vario, int *cormod, double *eps, int *flag,
		  double *grad, double lag, double *par, double tsep);
void TapVectCorrelation(double *rho,int *cormod,double *tdists,int *ntdists,double *nuis,double *par);
void VectCorrelation(double *rho, int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model, 
                     double *nuis,double *par, double *u);
void VectCorrelation_biv(double *rho, double *vario,int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model, 
                     double *nuis,double *par, double *u);


double Variogram(int *cormod, double h, double u, double nugget, double var, double *par);

double VarioFct(int *cormod, double h, double *par, double u);

double VarioDobStable(double lag, double power_s, double power_t, double scale_s, double scale_t, double tsep);

double VarioGCauchy(double lag, double power1, double power2, double scale);


/*----------------------------------------------------------------
 Function for matrix manipulation
 ---------------------------------------------------------*/

void CoFactor(double **a,int n,double **b);
void Transpose(double **a,int n,double k);
void Matrix_prod(double **a,double **b,double **c,int n);
void Matrix_prod_vec(double **a,double *b,double *c,int n);
double Determinant(double **a,int n);
double QFORM(double **A,double *x,double *y,int n);
double QFORM2(double **A,double *x,double *y,int n, int m);
double Trace(double **A,int n);
void lubksb(double **a, int n, int *indx, double *b);
void ludcmp(double **a, int n, int *indx, double *dd);

/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
End
 ---------------------------------------------------------------*/

// 2)
/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
Start
 ---------------------------------------------------------------*/
void Comp_Pair_Kumaraswamy2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_Kumaraswamy_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_Weibull_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_SinhGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_2Gamma2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_Gamma2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_Gamma_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Cond_Gauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Cond_Gauss_st2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Cond_Gauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt, double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Cond_BinGauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data,int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Cond_BinGauss_st2( int *cormod, double *coordx, double *coordy, double *coordt, double *data,int *NN, double *par,int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Diff_Gauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Diff_Gauss_st2(int *cormod,double *coordx, double *coordy, double *coordt,  double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Diff_BinGauss2( int *cormod, double *coordx, double *coordy, double *coordt, double *data,int *NN,  double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Diff_BinGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Diff_Gauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_Gauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data,int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_WrapGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_Gauss_st2(int *cormod, double *coordx, double *coordy, double *coordt, double *data,int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int  *ns,int *NS,int *GPU, int *local);
void Comp_Pair_WrapGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_Gauss_biv2(int *cormod,  double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_WrapGauss_biv2(int *cormod,  double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Cond_Gauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_TapGauss(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_BinomGauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_Binom2Gauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_BinomnegGauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_GeomGauss2(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_PoisbinGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_PoisbinnegGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_PoisbinGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_PoisbinnegGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_Binom2Gauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_BinGauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_BinomGauss_biv2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_BinGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_BinomGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_BinomnegGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_GeomGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed,double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_SkewGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_SkewGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_T2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_TWOPIECET2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_TWOPIECET_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_TWOPIECEGauss2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                         int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_TWOPIECEGauss_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,
                        int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_T_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                        int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns,int *NS, int *GPU,int *local);
void Comp_Pair_LogLogistic_st2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                               int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *GPU, int *local);
void Comp_Pair_LogLogistic2(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                            int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *GPU, int *local);


/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
End
 ---------------------------------------------------------------*/



// 2.1)
/*----------------------------------------------------------------
 File name: CompositeLikelihood_OCL.c
 Description: functions for composite log-likelihood evaluation in OpenCL
 Start
 ---------------------------------------------------------------*/

void Comp_Cond_Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data,int *NN, double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS, int *local,int *GPU);
void Comp_Diff_Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data, int *NN,
                          double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis, int *ns, int *NS,int *local,int *GPU);
void Comp_Pair_WrapGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local,int *GPU);
void Comp_Pair_SinhGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_LogGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_Gauss_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                             double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_SkewGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                              double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_Gamma2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                          double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_PoisbinnegGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                    double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_BinomGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                               double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_BinomnegGauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_Binom2Gauss2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                                double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_Logistic2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                             int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);
void Comp_Pair_LogLogistic2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,
                                int *NN,  double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local, int *GPU);

void Comp_Pair_T2_OCL(int *cormod, double *coordx, double *coordy, double *coordt, double *data, int *NN,
                      double *par,  int *weigthed,double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,
                      int *local_wi, int *dev);


void Comp_Pair_T_st2_OCL(int *cormod, double *coordx, double *coordy, double *coordt,double *data,int *NN,
                         double *par, int *weigthed, double *res,double *mean,double *mean2,double *nuis,int *ns, int *NS,int *local_wi, int *dev);

void wendintegral_call(double *x, double *param, double *res);
/*----------------------------------------------------------------
 File name: CompositeLikelihood_OCL.c
 Description: functions for composite log-likelihood evaluation in OpenCL
 End
 ---------------------------------------------------------------*/



// 3)
/*----------------------------------------------------------------
File name: Distributions.c
Description: procedures for the computation of useful distributions
Start
 ---------------------------------------------------------------*/
double biv_two_pieceTukeyh(double rho,double zi,double zj,double sill,double eta,double tail,double p11,double mui,double muj);
double cdf_norm(double lim1,double lim2,double a11,double a12);
double cdf_norm2(double lim1,double lim2,double a11,double a12, double a22);

double log_biv_Norm(double corr,double zi,double zj,double mi,double mj,double vari, double nugget);

double d2norm(double x, double y, double rho);

double d2lognorm(double x, double y, double sill,double nugget, double mux,double muy,double rho);

double dNnorm(int N,double **M, double *dat);

double int_pt(double x, double df);

double int_gen(double x,double mu, double alpha,double lag,double supp);
double int_hyp(double x,double a, double b,double c,double z);
double int_gen_hyp(double x,double a, double b,double z,double c);
void integr_pt(double *x, int n, void *ex);
void integr_hyp(double *x, int n, void *ex);
void integr_gen_hyp(double *x, int n, void *ex);
double HyperG_integral(double x, double *param);
double IntHyp(double x, double *param);
void integr_gen(double *x, int n, void *ex);

double int_gen_skew(double x,double data_i, double data_j,double c_0i,double c_0j,double rho,double *nuis);
void integr_gen_skew(double *x, int n, void *ex);
double cond_exp_skew(double c_0i,double c_0j,double rho, double data_i,double data_j,double *nuis);

double  wendintegral(double x, double *param);

double pbnorm(int *cormod, double h, double u, double lim1, double lim2, double nugget, double var,double *par, double thr);
double phalf_gauss (double z);
double pbnorm22(double lim1,double lim2,double corr,double nugget);
double ptnorm(int which,int *cormod, double h0,double h1,double h2, double u0, double u1,double u2, 
                 double *nuis, double *par, double thr);
double pbhalf_gauss (double zi,double zj,double rho,double nugget);
double pnorm_two_piece(double x, double eta);
double pbnorm_two_piece( int *cormod, double h, double u, 
    double xi, double xj, double nugget, double var,double eta,double *par);
void vpbnorm(int *cormod, double *h, double *u, int *nlags, int *nlagt,
	     double *nuis, double *par, double *rho, double *thr);


/*----------------------------------------------------------------
File name: Distributions.c
Description: procedures for the computation of useful distributions
End
 ---------------------------------------------------------------*/

// 4)
/*----------------------------------------------------------------
File name: Godambe.c
Description: procedures for the computation of the Godambe matrix
Start
 ---------------------------------------------------------------*/


void GodambeMat(double *betas,int *biv,double *coordx, double *coordy, double *coordt, int *cormod, double *data, int *dst,
          double *eps,int *flagcor, int *flagnuis, int *grid, int *like, double *mean,int *model,double *NN, int *nbetas,
          int *npar, int *nparc,int *nparcT, double *parcor, double *nuis, double *score,
          double *sensmat, int *spt,  int *type_lik, double *varimat,
          int *vartype, double *winc, double *winstp,double *winct,double *winstp_t,int *weigthed, double *X,int *ns,int *NS);


void Sens_Cond_Gauss(double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
		     int *npar, int *nparc, double *parcor, double *score, double *sensmat,int *weigthed);

void Sens_Cond_Gauss_st(double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
			int *npar, int *nparc, double *parcor, double *score, double *sensmat,int *weigthed);

void Sens_Cond_Gauss_biv(double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis, double *nuis, int *np,
                         int *npar, int *nparc, double *parcor, double *score, double *sensmat,int *weigthed);


void Sens_Cond_Gauss_ij(double rho, int *flag, double *gradcor, int *npar,
			int *nparc, double *par, double *sensmat);

void Sens_Cond_Gauss_biv_ij(double *gradcorttii,double *gradcorvvii,double *gradcorvtii,double *gradcortvii ,
                            double *gradcorttij ,double *gradcorvvij, double *gradcorvtij,double *gradcortvij,
                            double **inverse,int *flag,int *npar, int * nparc, double *par,int N,double *sens);

void Sens_Diff_Gauss(double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
		     double *nuis, int *np,int *npar, int *nparc, double *parcor,double *score, double *sensmat,int *weigthed);

void Sens_Diff_Gauss_st(double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
			double *nuis, int *np,int *npar, int *nparc, double *parcor, double *score, double *sensmat,int *weigthed);

void Sens_Diff_Gauss_ij(double *gradient, int *npar, double *sensmat,double weigths);


void Sens_Pair(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,double *NN,
         double *nuis, int *np, int *nbetas,int *npar,  int *nparc,int *nparcT,double *mean, int *model, double *parcor,double *score,  double *sensmat,int *weigthed,double *Z);

void Sens_Pair_st(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, 
         double *eps, int *flagcor, int *flagnuis, double *NN,double *nuis, int *np,
      int *nbetas, int *npar,int *nparc,int *nparcT, double *mean, int *model,
            double *parcor, double *score, double *sensmat,  int *weigthed,double *Z, int *ns, int *NS);

void Sens_Pair_biv(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, double *eps, int *flagcor, int *flagnuis,
                         double *NN, double *nuis, int *np,int *npar, int *nparc,int *nparcT, double *mean,int *mode,
                         double *parcor, double *score, double *sens,int *weigthed,double *Z, int *ns, int *NS);

void Sens_Pair_Gauss_ij(double rho, int *flag, double *gradcor, int *npar, int nbetas,
      int *nparc, double *par, double *sensmat, double *Xi, double *Xj);

void Sens_Pair_Gauss_biv_ij(double rhott,double rhotv,double rhovt,double rhovv,double *gradcortt,double *gradcortv,double *gradcorvt,double *gradcorvv,
                            int *flag, int *npar,int *nparc, double *par, double *sensmat);

void Sensitivity(double *betas,int *biv,double *coordx,double *coordy,double *coordt,int *cormod,  double *data, double *eps, int *flagcor, int *flagnuis, int *like,
     double *mean, int *model, double *NN,int *nbetas, int *npar, int *nparc,int *nparcT, double *parcor, double *nuis, int *np,double *score,
     double *sensmat, int *spt, int *type_lik,int *weigthed,double *Z,int *ns,int *NS);
 

void Vari_SubSamp(double *betas,double *coordx, double *coordy, double *coordt,int *cormod, double *data,
                  int *dst,double *eps, int *flagcor, int *flagnuis, int *grid,
                  int *like, double *mean,int *model, double *NN, int *nbetas,int *npar, int *nparc,int *nparcT, double *nuis, int *np,
                  double *parcor,  int *type_lik, double *varimat,
                  double *winc, double *winstp,int *weigthed,double *Z);
void Vari_SubSamp_biv(double *betas,double *coordx, double *coordy, double *coordt,int *cormod, double *data,
                  int *dst,double *eps, int *flagcor, int *flagnuis, int *grid,
                  int *like, int *model, double *NN,int *npar, int *nparc,int *nparcT, double *nuis, int *np,
                  double *parcor,  int *type_lik, double *varimat,
                  double *winc, double *winstp,int *weigthed,double *Z,int *ns,int *NS);

void Vari_SubSamp_st2(double *betas,double *coordx, double *coordy, double *coordt, int *cormod, double *data, int *dst, double *eps, int *flagcor,
                   int *flagnuis,int *like,double *mean,int *model,double *NN,int *nbetas,int *npar, int *nparc,int *nparcT, double *nuis,int *np, double *parcor,
                   int *type_lik,double *varimat, double *winc, double *winstp,double *winc_t,double *winstp_t,int *weigthed,double *Z,int *ns,int *NS);
void cumvec(int *ns,int *res,int len);
void Rep(double *coordt,int *ns, double *res);
/*----------------------------------------------------------------
File name: Godambe.c
Description: procedures for the computation of the Godambe matrix
End
 ---------------------------------------------------------------*/

// 5)
/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of the gredients
Start
 ---------------------------------------------------------------*/

void Grad_Cond_Bin(double rho,double pij, double p,int *flag, double *gradcor,
		   double *grad, int *npar, double *nuis,  double u, double v);

void Grad_Cond_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v);

void Grad_Diff_Bin(double rho,double pij, double p,int *flag, double *gradcor,
		   double *grad, int *npar, double *nuis,  double u, double v);

void Grad_Diff_Gauss(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par, double u, double v);

void Grad_Diff_Vario(double rho, int *flag, double *gradcor, double *grad,
		     int *npar, double *par);


void Grad_Pair_Binom(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc,int *nparcT, int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);
void Grad_Pair_Wrapped(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc,int *nparcT, int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_Binomneg(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas);

void Grad_Pair_Gamma(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc,int *nparcT, int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_LogGauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_LogLogistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_Logistic(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_StudenT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_Tukeyh(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas);

void Grad_Pair_Weibull(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc, int *nparcT,int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas);

void Grad_Pair_Sinh(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_Skewgauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_Twopiecegauss(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas);

void Grad_Pair_TwopieceT(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,
  double NN,int *npar,int *nparc,int *nparcT, int nbetas, double *nuis, double *par, double u, double v,
       double ai, double aj,double *Xl,double *Xm,double **sX,int l,int m,double *betas);

void Grad_Pair_Gauss2(double rho,int *cormod,int *flag,int *flagcor, double *gradcor, double *grad, double lag, double lagt,double NN,
       int *npar,int *nparc, int *nparcT,int nbetas,double *nuis, double *par,  double u, double v,double ai,double aj,double *Xl, double *Xm
       ,double **sX,int l,int m,double *betas);

void Grad_Pair_Gauss(double rho, int *flag,int *flagcor,  double *gradcor, double *grad,
                     int *npar,int *nparc, int nbetas, double *par, double u, double v, double *Xl,double *Xm);

void Grad_Pair_Gauss_biv(double rhott, double rhotv, double rhovt,double rhovv,int *flag,
                         double *gradcortt, double *gradcortv, double *gradcorvt, double *gradcorvv,
                         double *grad,int *npar, double *par, double u, double v);


void Grad_Cond_Gauss_biv(double *gradcorttii,double *gradcorvvii,double *gradcorvtii,double *gradcortvii ,
                        double *gradcorttij ,double *gradcorvvij, double *gradcorvtij,double *gradcortvij,
                        double **inverse,int *flag,double *grad,int *npar, double *par,double *dat,int N);



double mij(double qij, double w, double pij, double p);

double nij(double qij, double w, double pij, double p);

/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of the gradients
End
 ---------------------------------------------------------------*/

// 5)
/*----------------------------------------------------------------
File name: Likelihood.c
Description: procedures for computation of the full likelihood
Start
 ---------------------------------------------------------------*/

// Insert here

/*----------------------------------------------------------------
File name: Likelihood.c
Description: procedures for computation of the full likelihood
End
 ---------------------------------------------------------------*/


// 6)
/*----------------------------------------------------------------
File name: WeightedLeastSquare.c
Description: procedures for the estimation of the model parameters
via the weighted least square method.
Start
 ---------------------------------------------------------------*/



void Binned_Variogram2(double *bins, double *coordx, double *coordy, double *coordt,double *data, int *lbins, double *moms, int *nbins);

void Binned_Variogram_22(double *bins,double *coordx, double *coordy, double *coordt, double *data, int *lbins, double *moms, int *nbins);

void Binned_Variogram_st2(double *bins, double *bint,double *coordx, double *coordy, double *coordt, double *data, int *lbins,
                         int *lbinst, int *lbint, double *moms, double *momst,
                         double *momt, int *nbins, int *nbint, int *ns,int *NS);


void Binned_Variogram_biv2(double *bins,double *coordx, double *coordy, double *coordt, double *data, 
  int *cross_lbins, double *cross_moms, int *nbins,int *marg_lbins, double *marg_moms,int *ns, int *NS);


void Cloud_Variogram2(double *bins,double *coordx, double *coordy, double *coordt, double *data, int *lbins, double *moms, int *nbins);

void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);



void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
		      int *nbins, int *nbint, double *nuis, double *par, double *res);



/*----------------------------------------------------------------
File name: WeightedLeastSquare.c
Description: procedures for the estimation of the model parameters
via the weighted least square method.
End
 ---------------------------------------------------------------*/

// 7)
/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
Start
 ---------------------------------------------------------------*/

void DeleteGlobalVar();

void RangeDist(double *max, double *min);

double fac(int x,int j);

double logfac(int n);

/*for bivariate t distributions*/
double Poch(double q,double n);
double biv_T(double rho,double zi,double zj,double nuu);
double biv_two_pieceT(double rho,double zi,double zj,double sill,double nuu,double eta,double p11,
  double mui,double muj);
/**********************************************/
double biv_half_Gauss(double rho,double zi,double zj);
double biv_two_pieceGaussian(double rho,double zi,double zj,double sill,double eta,double p11,
  double mui,double muj);



double stirling(double x);

double dist(int type_dist,double coordx,double locx,double coordy,double locy, double radius); 

double Dist_geodesic(double loni, double lati, double lonj, double latj,double radius);

double Dist_chordal(double loni, double lati, double lonj, double latj,double radius);


int is_equal(double val1, double val2);

double Maxima(double *x, int *size);
double Maxima_i(int *x, int size);

double Minima(double *x, int *size);

void Maxima_Minima_dist(double *res,double *coordx,double *coordy,int *nsize,int *type_dist,double *radius);

void Maxima_Minima_time(double *res,double *coordt,int *nsize);

void Range(double *x, double *ran, int *size);

void Seq(double *x, int len, double *res);

void SeqStep(double *x, int len, double step, double *res);

void SetSampling(double *coordx, double *coordy, double *data, int n, int *npts, int nbetas,
     double *scoordx, double *scoordy, double *sdata, double xmax,
     double xmin, double ymax, double ymin, double **sX,double **X);


void SetSampling_biv(double *coordx, double *coordy, double *data, int n, int *npts,
                     double *scoordx, double *scoordy, double *sdata, double xmax,
                     double xmin, double ymax, double ymin);

void SetSampling_st(double *data,double *sdata,int *ncoord,int *ntime, int nbetas,
        int wint,int k, double **sX,double **X);

void SetSampling_s(double *coordx, double *coordy, double *data, int *npts, int nbetas,
                   double *scoordx, double *scoordy, double *sdata, double xmax,
                   double xmin, double ymax, double ymin,double **sX,double **X,int *ns,int *NS, int *NS_sub, double *res_sub, double *coordt,int *ns_sub);

void SetSampling_t(double *data,double *sdata, int nbetas,int npts,
                   int nt,int wint,int k,double **sX,double **X,int *ns_sub,int *NS_sub,int nsub_t, int *ntimeS, double *s2cx, double *s2cy, double *scoordx, double *scoordy);



void SetGlobalVar(int *biv,double *coordx,double *coordy,double *coordt,
      int *grid,int *ia,
		  int *idx,int *ismal,int *ja,int *mem, int *nsite,int *nsitex,int *nsitey,
		  int *npair, double *radius, int *replic,double *srange, double *sep,int *st, int *times,double *trange,
		  int *tap,int *tapmodel,int *tp,int *weighted, int *dyn);

void Space_Dist(double *coordx,double *coordy,int grid,int *ia,int *idx,
		int *ismal,int *ja,double thres);

void SpaceTime_Dist(int biv,double *coordx,double *coordy,double *coordt,int *grid,int *ia,int *idx,int *ismal,int *ja,
                    int *tapmodel,double *thres,double *thret);


/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for the computation of useful quantities.
End
 ---------------------------------------------------------------*/





/*----------------------------------------------------------------
 File name: Host.c
 Description: procedures for OCL computation.
 Start
 ---------------------------------------------------------------*/
#define CL_SILENCE_DEPRECATION
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif

#define LOW -1.0e15
#define MAXERR 1e-6
#define EPS DBL_EPSILON

#pragma once

#include <string.h>

#define MAX_PLATFORMS     8
#define MAX_DEVICES      16
#define MAX_INFO_STRING 256



#ifdef __cplusplus
#include <cstdio>
#endif





const char *err_code (cl_int err_in);
void check_error(cl_int err, const char *operation, char *filename, int line);
#define checkError(E, S) check_error(E,S,__FILE__,__LINE__)
unsigned getDeviceList(cl_device_id devices[MAX_DEVICES]);
void getDeviceName(cl_device_id device, char name[MAX_INFO_STRING]);
int parseUInt(const char *str, cl_uint *output);
double sum_total(double *arr, int ngrid);
void param_OCL(int *cormod,int *NN,double *par,int *weigthed,double *nuis,int *int_par, double *dou_par);
void exec_kernel(double *h_x, double *h_y, double *h_mean, double *h_data, int *int_par,double *dou_par,
                 int *local_wi, int *dev, double *res, char *f_name);
void create_binary_kernel(int *dev, char **fname);
int DeviceInfo();
void exec_kernel_st(double *h_x, double *h_y,double *h_t, double *h_mean, double *h_data, int *int_par,double *dou_par,
                    int *local_wi, int *dev, double *res, char *f_name, int *ns, int *NS);
void param_st_OCL(int *cormod,int *NN,double *par,int *weigthed,double *nuis,int *int_par, double *dou_par);
void CBessel(double *xxx, double *nuu, int *expscale,double *res, int *tipo);
void exec_kernel_source(double *h_x, double *h_y, double *h_mean, double *h_data, int *int_par,double *dou_par,
                        int *local_wi, int *dev, double *res, char *f_name);
void exec_kernel_st_dyn(double *h_x, double *h_y,double *h_t, double *h_mean, double *h_data, int *int_par,double *dou_par,
                        int *local_wi, int *dev, double *res, char *f_name,int *ns, int *NS);
void cdf_norm_call(double *lim1,double *lim2,double *a11,double *a12, double *res);


/*----------------------------------------------------------------
 File name: Host.c
 Description: procedures for OCL computation.
 End
 ---------------------------------------------------------------*/





/*for bivariate t distributions*/
double hyt2f1( double a, double b, double c, double x, double *loss );
double hys2f1( double a,double b,double c,double x,double *loss );
double hyp2f1( double a,double b,double c,double x);
double hypergeo(double a,double b,double c,double x);
void hypergeo_call(double *a,double *b,double *c,double *x, double *res);
