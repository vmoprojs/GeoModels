
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#define extern
#include "header.h"
#undef extern



/********************** utility C calls ****************************************************/

extern void hyperg_U_e_call( double *a,  double *b,  double *x, double *val);
extern void SetGlobalVar2 (int *nsite, int *times,
                           double *h,int *nn,  double *maxh,
                           double *u,int *tt,  double *maxu,
                           int *st,int *biv,int *one, int *two);

extern void SetGlobalVar(int *biv,double *coordx,double *coordy,double *coordz,double *coordt,int *grid,int *ia,
      int *idx,int *ismal,int *ja,int *mem, int *nsite,int *nsitex,int *nsitey,int *nsitez,
      int *npair,double *radius,double *srange, double *sep,int *st, int *times,double *trange,
      int *tap,int *tapmodel,int *tp,int *weighted, int *colidx,int *rowidx, int *ns, int *NS, int *dyn);

extern void DeleteGlobalVar(void);
extern void DeleteGlobalVar2(void);

extern void Maxima_Minima_dist(double *res,double *coordx,double *coordy,double *coordz,int *nsize,int *type_dist,double *radius);

/* for Turning band */
extern void TBD1d(double *ux, double *uy, double *sx, double *sy, double *phi, int *L, int *N, double *result);
extern void spectraldensityC(double u,int model,int d,int L,double *f,double *av,double *Cv,double *nu1v,double *nu2v, double *params_other);
extern void extraer(double *coord,int sequen1,double *sub_coord,int fila,int col, int d);
extern void rellenar_indice(int *index,int inicio, int final,int largo);
extern void u_index_extraer(double *u_index,double *u, int *index,int largo,int d,int filas);
extern void mult_mat(double *z, double *x, int xrows, int xcols, double *y, int yrows, int ycols);
extern void tcrossprod(double *z, double *x,  int xrows, int xcols, double *y, int yrows, int ycols);
extern void mult_x_cons(double *x, double cte,int largo);
extern void sumar_matrices(double *x0, double *x,double *y,int largo);
extern void restar_matrices(double *x0, double *x,double *y,int largo);
extern void cos_vec(double *x_cos,double *x,int largo);
extern void sen_vec(double *x_sen,double *x,int largo);
extern void llenar_simu(double *x,double *simu, int N,int *P, int m);
extern void extraer_col(int inicio, int final,double *x_original,double *x);
extern void llenar_simu1(double *simu1,double *simu,int *m,int *P,int *N,int lim, int i,double *L1);
extern void C_tcrossprod(Rcomplex *z, Rcomplex *x,  int xrows, int xcols, Rcomplex *y, int yrows, int ycols);
extern void C_mult_mat(Rcomplex *z, Rcomplex *x, int xrows, int xcols, Rcomplex *y, int yrows, int ycols);
extern void for_c(int *d_v, double *a_v, double *nu1_v, double *C_v, double *nu2_v, 
           int *P, int *N, int *L, int *model, double *u,
           double *a0, double *nu0, double *A, double *B,
           int *sequen, int *largo_sequen, int *n,
           double *coord, double *phi, int *vtype, int *m1, double *simu1, double *L1, double *params_other);
extern void spectral_density_1d(double *norm_u, int *N, double *av, double *params_other, double *nu1v, int *model, double *result);
/*****/

extern void pairs(int *ncoords,double *data,double *coordx,double *coordy,double *coordz,double *numbins,
                       double *bins,double *v0,double *v1,double *v2,double *maxdist);


/********************** for variogrms computations  ****************************************************/

extern void Binned_Variogram_biv2(double *bins,double *coordx, double *coordy,double *coordz, double *coordt, double *data,
  int *cross_lbins, double *cross_moms, int *nbins,int *marg_lbins, double *marg_moms,int *ns, int *NS);

extern void Binned_Variogram_st2(double *bins, double *bint, double *coordx, double *coordy,double *coordz, double *coordt,double *data, int *lbins, int *lbinst,
       int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint, int *ns,int *NS);

extern void Binned_Variogram2(double *bins, double *coordx, double *coordy, double *coordz, double *coordt,double *data, int *lbins, double *moms, int *nbins);

extern void Binned_Variogram2new(double *bins, int *np,double *data1, double *data2, double *vdist, int *lbins, double *moms, int *nbins,double *mm);


extern void Binned_Variogram_biv2new(double *bins, int *np,double *data1, double *data2,  double *vdist, double *mm,
     double *moms00,double *moms10,double *moms11,int *lbins00,int *lbins10,int *lbins11,int *nbins,int *first, int *second);


extern void Binned_Variogram_st2_dyn(double *bins, double *bint, double *coordx, double *coordy,double *coordz, double *coordt,double *data, int *lbins, int *lbinst,
                              int *lbint, double *moms,double *momst, double *momt, int *nbins, int *nbint, int *ns,int *NS);

extern void Cloud_Variogram2(double *bins,double *coordx, double *coordy,double *coordz, double *coordt, double *data, int *lbins, double *moms, int *nbins);
extern void LeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
          int *nbins, int *nbint, double *nuis, double *par, double *res);
extern void WLeastSquare_G(double *bins, double *bint, int *cormod, double *lbins, double *moms,
           int *nbins, int *nbint, double *nuis, double *par, double *res);

extern void Binned_Variogram_22(double *bins, double *coordx, double *coordy,double *coordz, double *coordt,double *data, int *lbins, double *moms, int *nbins);

/********************** for distribution computations  ****************************************************/

extern void biv_unif_CopulaClayton_call(double *x,double *y,double *rho, double *nu, double *res);
extern void biv_unif_CopulaGauss_call(double *x,double *y,double *rho, double *res);
extern void corr_kuma_vec(double *rho,double *eta,double *gam,double *res, int *n);
extern void vpbnorm(int *cormod, double *h, double *u, int *nlags, int *nlagt,
                    double *nuis, double *par, double *rho, double *thr);

/*********************** for correlation computations ***************************************************/

extern void Corr_c(double *cc,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod, int *grid, double *locx,double *locy,double *locz,
        int *ncoord, int *nloc,int *tloc,int *ns,int *NS,int *ntime, double *par, int *spt, int *biv, double *time,int *type, int *which,double *radius);

extern void Corr_c_bin(double *cc,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod, int *grid, double *locx,double *locy,double *locz,
            int *ncoord, int *nloc,int *model,int *tloc,int *nn,int *n, int *ns,int *NS,int *ntime, double *mean,
            double *nuis, double *par, int *spt, int *biv, 
                      double *time,int *type, int *which,double *radius);


extern void Corr_c_tap(double *cc,double *cc_tap,double *coordx, double *coordy, double *coordz, double *coordt, int *cormod, int *cormodtap, int *grid,
 double *locx,  double *locy,double *locz,
                double *mxd,double *mxt, int *ncoord, int *nloc, int *ns,int *NS,int *tloc,int *ntime, double *par, int *spt, int *biv, double *time,
                int *type,int *which,double *radius);

extern void VectCorrelation(double *rho, int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model,
                            double *nuis,double *par, double *u,int *N);

extern void VectCorrelation_biv(double *rho, double *vario,int *cormod, double *h, int *nlags, int *nlagt,double *mean,int *model,
                                double *nuis,double *par, double *u,int *N);

extern void CorrelationMat2(double *rho,double *coordx, double *coordy, double *coordz, double *coordt,  int *cormod,
 double *nuis, double *par,double *radius,int *ns, int *NS);

extern void CorrelationMat_dis2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt,  int *cormod, double *mean,
                         int *nn,double *nuis, double *par,double *radius, int *ns, int *NS,int *model);

extern void CorrelationMat_st_dyn2(double *rho, double *coordx, double *coordy,double *coordz, double *coordt,int *cormod,
  double *nuis, double *par,double *radius, int *ns,int *NS);
extern void CorrelationMat_st_dyn_dis2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt,  int *cormod,  double *mean,int *n,
                                double *nuis, double *par,double *radius, int *ns, int *NS, int *model);
extern void CorrelationMat_biv_dyn2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod,  double *nuis,
  double *par,double *radius, int *ns,int *NS);

extern void CorrelationMat_biv_tap(double *rho, double *coordx, double *coordy,double *coordz, double *coordt,int *cormod,
 double *nuis, double *par,double *radius, int *ns,int *NS);

extern void DCorrelationMat_biv(int *cormod,double *coordx, double *coordy, double *coordz,double *coordt,double *drho,double *eps,int *flagcor,
      int *nparcor,double *parcor,double *rho);

extern void DCorrelationMat_biv2(int *cormod,double *coordx, double *coordy, double *coordz,double *coordt,double *drho,double *eps,int *flagcor,
      int *nparcor,double *parcor,double *rho);

extern void DCorrelationMat_biv_tap(int *cormod,double *coordx, double *coordy, double *coordz,double *coordt,double *drho,
            double *eps,int *flagcor,int *nparcor,double *parcor,double *rho);


extern void CorrelationMat_biv_skew_dyn2(double *rho,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod,
 double *nuis, double *par,double *radius, int *ns,int *NS);
extern void CorrelationMat_tap(double *rho,double *coordx, double *coordy, double *coordz,double *coordt, int *cormod,  double *nuis, double *par,double *radius,
  int *ns, int *NS);
extern void CorrelationMat_dis_tap(double *rho,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod,  double *nuis, double *par,double *radius,
           int *ns, int *NS, int *n1,int *n2, double *mu1,double *mu2,int *model);
extern void CorrelationMat_st_tap(double *rho,double *coordx, double *coordy, double *coordz,double *coordt, int *cormod, double *nuis,
  double *par,double *radius, int *ns, int *NS);

/*********************** for bivariate composite likelihood ***************************************************/
extern void Comp_Pair_Gauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);

extern void Comp_Diff_Gauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);

extern void Comp_Pair_WrapGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_SinhGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_SkewGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Weibull2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy22mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Beta2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_LogGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisbinnegGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisbinGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegGaussZINB2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNGauss_misp2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_LogLogistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Logistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Pois2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisGamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisGammaZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_PoisGamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_Pois2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_SkewT2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_T2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Tukeyhh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Tukeyh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_Tukeygh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_T2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECETukeyh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECET2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECEGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECEBIMODAL2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_WrapGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_T_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Pois_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisGamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_PoisGamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_Pois_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Tukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Tukeyhh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void  Comp_Pair_TWOPIECEGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void  Comp_Pair_TWOPIECETukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void  Comp_Pair_TWOPIECET_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisbinGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisbinnegGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegGaussZINB_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Diff_Gauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_SkewGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_SinhGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy2_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Beta_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Weibull_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_LogGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_LogLogistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Logistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECEBIMODAL_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_biv2mem(int *cormod, double *data1,double *data2,int *NN,
    double *par, int *weigthed,double *res,double *mean1,double *mean2,
    double *nuis, int *local,int *GPU);
extern void Comp_Pair_SkewGauss_biv2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
    double *par, int *weigthed,double *res,double *mean1,double *mean2,
    double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_T_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_GaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);

extern void Comp_Pair_PoisCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern  void Comp_Pair_LogisticCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond);
extern void Comp_Pair_BetaCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Beta2Cop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_GammaCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond);
extern void Comp_Pair_WeibullCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond);
extern void Comp_Pair_LogGaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond);
extern void Comp_Pair_KumaraswamyCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy2Cop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond);
extern void Comp_Pair_BinomNNGaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond);
extern void Comp_Pair_BinomnegGaussCop2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop,int *cond);
extern  void Comp_Pair_BinomnegBinary2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);

/*********************** for conditional composite likelihood ***************************************************/
extern void Comp_Cond_Gauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                         double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Tukeyh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                          double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Tukeyhh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_WrapGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);

extern void Comp_Cond_SkewGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                             double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_T2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                     double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_T2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_SkewT2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                    double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_Tukeygh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_SinhGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                             double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                         double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Weibull2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_LogGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                            double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Beta2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                        double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Kumaraswamy2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Kumaraswamy22mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_Pois2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNGauss_misp2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                     double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_PoisGamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                        double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_PoisGamma2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                             double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_PoisGammaZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                             double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Pois2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                        double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                              double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                             double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern  void Comp_Cond_BinomnegBinary2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegLogi2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECETukeyh2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                  double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECET2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                             double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECEGauss2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECEBIMODAL2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegGaussZINB2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                     double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_PoisZIP2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_LogLogistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Logistic2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                            double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                            double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_T_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                        double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_T_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Tukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                             double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Tukeyhh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                              double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_SkewGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_SinhGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                            double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Weibull_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                              double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_LogGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Beta_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Kumaraswamy_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                  double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Kumaraswamy2_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_Pois_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_PoisGamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_PoisGamma_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Pois_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                  double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                    double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegLogi_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECETukeyh_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                     double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECET_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECEGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                    double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECEBIMODAL_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegGaussZINB_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                        double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                              double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_PoisZIP_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                         double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_LogLogistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                  double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Logistic_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_WrapGauss_st2mem(int *cormod, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);



/********************** spatial anisotropic  marginal************************************/
extern void Comp_Pair_Gauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Diff_Gauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_WrapGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_SinhGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_SkewGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gamma2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                               double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Weibull2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                     double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy22mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Beta2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                              double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_LogGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                  double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisbinnegGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                         double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisbinGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                       double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegLogi2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomnegGaussZINB2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                    double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomLogi2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNGauss_misp2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BinomNNLogi2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                     double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_LogLogistic2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                     double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Logistic2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                  double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Pois2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                              double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisGamma2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_PoisGamma2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                              double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisZIP2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_Pois2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                         double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_PoisZIP2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                            double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_SkewT2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                          double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_T2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                      double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Tukeyhh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Tukeyh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Gauss_misp_Tukeygh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                            double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_T2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                           double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECETukeyh2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                        double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECET2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECEGauss2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                       double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TWOPIECEBIMODAL2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                         double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_GaussCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_PoisCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_TCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_BetaCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                  double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Beta2Cop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                   double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_KumaraswamyCop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                         double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Pair_Kumaraswamy2Cop2mem_aniso(int *cormod, double *coord1,double *coord2, double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
                                          double *nuis, int *local,int *GPU,int *type_cop, int *cond);


/********************** spatial anisotropic conditional ************************************/
extern void Comp_Cond_Gauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Tukeyhh2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_SkewGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_T2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_T2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_SkewT2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_Tukeygh2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_SinhGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gamma2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Weibull2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_LogGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Beta2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Kumaraswamy2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Kumaraswamy22mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_Pois2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNGauss_misp2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_PoisGamma2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_PoisGamma2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Pois2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomLogi2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomNNLogi2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegLogi2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECETukeyh2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECET2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECEGauss2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_TWOPIECEBIMODAL2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_BinomnegGaussZINB2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_PoisZIP2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Gauss_misp_PoisZIP2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_LogLogistic2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
extern void Comp_Cond_Logistic2mem_aniso(int *cormod, double *coord1, double *coord2,double *data1,double *data2,int *N1,int *N2,
 double *par, int *weigthed, double *res,double *mean1,double *mean2,
 double *nuis, int *local,int *GPU,int *type_cop, int *cond);
/*********************************************************************************************/

extern void CorrelationMat_st_dis_tap(double *rho,double *coordx, double *coordy,double *coordz, double *coordt, int *cormod,  double *nuis, double *par,double *radius,
                               int *ns, int *NS, int *n1,int *n2, double *mu1,double *mu2,int  *model);

extern void hypergeo_call(double *a,double *b,double *c,double *x, double *res);
extern void hyperg_1F2_e_call( double *a,  double *b,double *c,  double *x, double *val);

extern void hyperg_call(double *a,double *b,double *x,double *res);
extern void biv_pois_call(double *corr,int *r, int *t, double *mean_i, double *mean_j,double *res);
extern void appellF4_call(double *a,double *b,double *c,double *d,double *x,double *y, double *res);
extern void biv_gamma_call(double *corr,double *zi,double *zj,double *mui, double *muj, double *shape, double *res);

extern void biv_binomneg_call(int *NN,int *u, int *v, double *p01, double *p10,double *p11,double *res);

extern void biv_binom_call(int *NN,int *u, int *v, double *p01, double *p10,double *p11,double *res);

/*extern void matrix_temp(int *N ,double *matrix, double *l1 ,double *l2 ,double *v11 ,double *v21,double *v12,double *v22);*/

/*extern void vector_to_select(int *N, double *matrix);*/


extern void lgnd (int *lmax,double *x, double *p);
extern void qnorm55_call(double *p, double *mu, double *sigma, int *lower_tail, int *log_p, double *res);


extern void integr_kuma(double *x, int n, void *ex);


static const R_CMethodDef CEntries[] = {
    {"hyperg_U_e_call",             (DL_FUNC) &hyperg_U_e_call,          4},
    {"integr_kuma",                 (DL_FUNC) &integr_kuma,              3},
    {"lgnd",                        (DL_FUNC) &lgnd,                     3},
    {"qnorm55_call",                (DL_FUNC) &qnorm55_call,             6},
    {"biv_binomneg_call",           (DL_FUNC) &biv_binomneg_call,        7},
    {"biv_binom_call",              (DL_FUNC) &biv_binom_call,           7},
    {"biv_gamma_call",              (DL_FUNC) &biv_gamma_call,           7},
    {"hyperg_call",                 (DL_FUNC) &hyperg_call,              4},
    {"biv_pois_call",               (DL_FUNC) &biv_pois_call,            6},
    {"appellF4_call",               (DL_FUNC) &appellF4_call,            7},
    {"hypergeo_call",               (DL_FUNC) &hypergeo_call,            5},
    {"hyperg_1F2_e_call",           (DL_FUNC) &hyperg_1F2_e_call,        5},
    {"CorrelationMat_st_dis_tap",   (DL_FUNC) &CorrelationMat_st_dis_tap,16},
    {"CorrelationMat2",             (DL_FUNC) &CorrelationMat2,             11},
    {"CorrelationMat_dis2",         (DL_FUNC) &CorrelationMat_dis2,         14},
    {"CorrelationMat_st_dyn2",      (DL_FUNC) &CorrelationMat_st_dyn2,      11},
    {"CorrelationMat_biv_dyn2",     (DL_FUNC) &CorrelationMat_biv_dyn2,     11},
    {"CorrelationMat_st_dyn_dis2",  (DL_FUNC) &CorrelationMat_st_dyn_dis2,  14},
    {"CorrelationMat_biv_tap",      (DL_FUNC) &CorrelationMat_biv_tap,      10},
    {"DCorrelationMat_biv",         (DL_FUNC) &DCorrelationMat_biv,         11},
    {"DCorrelationMat_biv2",        (DL_FUNC) &DCorrelationMat_biv2,        11},
    {"DCorrelationMat_biv_tap",     (DL_FUNC) &DCorrelationMat_biv_tap,     11},
    {"CorrelationMat_biv_skew_dyn2",(DL_FUNC) &CorrelationMat_biv_skew_dyn2,11},
    {"CorrelationMat_tap",          (DL_FUNC) &CorrelationMat_tap,          10},
    {"CorrelationMat_dis_tap",      (DL_FUNC) &CorrelationMat_dis_tap,      16},
    {"CorrelationMat_st_tap",       (DL_FUNC) &CorrelationMat_st_tap,       11},
    {"Binned_Variogram2",           (DL_FUNC) &Binned_Variogram2,            9},
    {"Binned_Variogram2new",        (DL_FUNC) &Binned_Variogram2new,         9},
     {"Binned_Variogram_biv2new",   (DL_FUNC) &Binned_Variogram_biv2new,    15},
    {"Binned_Variogram_biv2",       (DL_FUNC) &Binned_Variogram_biv2,       13},
    {"Binned_Variogram_st2",        (DL_FUNC) &Binned_Variogram_st2,        17},
    {"Binned_Variogram_st2_dyn",    (DL_FUNC) &Binned_Variogram_st2_dyn,    17},
    {"Cloud_Variogram2",            (DL_FUNC) &Cloud_Variogram2,             9},
    {"LeastSquare_G",               (DL_FUNC) &LeastSquare_G,               10},
    {"WLeastSquare_G",              (DL_FUNC) &WLeastSquare_G,              10},
    {"Binned_Variogram_22",         (DL_FUNC) &Binned_Variogram_22,          9},
    {"biv_unif_CopulaClayton_call", (DL_FUNC) &biv_unif_CopulaClayton_call,  5},
    {"biv_unif_CopulaGauss_call",   (DL_FUNC) &biv_unif_CopulaGauss_call,    4},
    {"Corr_c",                      (DL_FUNC) &Corr_c,                      23},
    {"Corr_c_bin",                  (DL_FUNC) &Corr_c_bin,                  28},
    {"Corr_c_tap",                  (DL_FUNC) &Corr_c_tap,                  27},
    {"corr_kuma_vec",               (DL_FUNC) &corr_kuma_vec,                5},
    {"DeleteGlobalVar",             (DL_FUNC) &DeleteGlobalVar,              0},
    {"DeleteGlobalVar2",            (DL_FUNC) &DeleteGlobalVar2,             0},
    {"Maxima_Minima_dist",          (DL_FUNC) &Maxima_Minima_dist,           7},
    {"pairs",                       (DL_FUNC) &pairs,                       11},
    {"SetGlobalVar2",               (DL_FUNC) &SetGlobalVar2,               12},
    {"SetGlobalVar",               (DL_FUNC) &SetGlobalVar,                 31},
   /* {"simu_on_coords",              (DL_FUNC) &simu_on_coords,               8},*/

/* for Turning band */

    {"TBD1d",                       (DL_FUNC) &TBD1d,                      8},
    {"spectraldensityC",            (DL_FUNC) &spectraldensityC,          10},
    {"spectral_density_1d",         (DL_FUNC) &spectral_density_1d,        7},
    {"extraer",                     (DL_FUNC) &extraer,                    6},
    {"rellenar_indice",             (DL_FUNC) &rellenar_indice,            4},
    {"u_index_extraer",             (DL_FUNC) &u_index_extraer,            6},
    {"mult_mat",                    (DL_FUNC) &mult_mat,                   7},
    {"tcrossprod",                  (DL_FUNC) &tcrossprod,                 7},
    {"mult_x_cons",                 (DL_FUNC) &mult_x_cons,                3},
    {"sumar_matrices",              (DL_FUNC) &sumar_matrices,             4},
    {"restar_matrices",             (DL_FUNC) &restar_matrices,            4},
    {"cos_vec",                     (DL_FUNC) &cos_vec,                    3},
    {"sen_vec",                     (DL_FUNC) &sen_vec,                    3},
    {"llenar_simu",                 (DL_FUNC) &llenar_simu,                5},
    {"extraer_col",                 (DL_FUNC) &extraer_col,                4},
    {"llenar_simu1",                (DL_FUNC) &llenar_simu1,               8},
    {"C_tcrossprod",                (DL_FUNC) &C_tcrossprod,               7},
    {"C_mult_mat",                  (DL_FUNC) &C_mult_mat,                 7},
    {"for_c",                       (DL_FUNC) &for_c,                     24},
    {"VectCorrelation",             (DL_FUNC) &VectCorrelation,             11},
    {"VectCorrelation_biv",         (DL_FUNC) &VectCorrelation_biv,         12},
    {"vpbnorm",                     (DL_FUNC) &vpbnorm,                      9},

  /*********************** for pairwise composite likelihood ***************************************************/
    {"Comp_Pair_Gauss2mem",         (DL_FUNC) &Comp_Pair_Gauss2mem,         15},
    {"Comp_Pair_Tukeyh2mem",         (DL_FUNC) &Comp_Pair_Tukeyh2mem,         15},
    {"Comp_Pair_Tukeyhh2mem",         (DL_FUNC) &Comp_Pair_Tukeyhh2mem,         15},
    {"Comp_Pair_SkewGauss2mem",         (DL_FUNC) &Comp_Pair_SkewGauss2mem,         15},
    {"Comp_Pair_WrapGauss2mem",         (DL_FUNC) &Comp_Pair_WrapGauss2mem,         15},
    {"Comp_Pair_T2mem",         (DL_FUNC) &Comp_Pair_T2mem,         15},
    {"Comp_Pair_Gauss_misp_T2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_T2mem,         15},
    {"Comp_Pair_Gauss_misp_SkewT2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_SkewT2mem,         15},
    {"Comp_Pair_Gauss_misp_Tukeygh2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_Tukeygh2mem,         15},
    {"Comp_Pair_SinhGauss2mem",         (DL_FUNC) &Comp_Pair_SinhGauss2mem,         15},
    {"Comp_Pair_Gamma2mem",         (DL_FUNC) &Comp_Pair_Gamma2mem,         15},
    {"Comp_Pair_Weibull2mem",         (DL_FUNC) &Comp_Pair_Weibull2mem,         15},
    {"Comp_Pair_LogGauss2mem",         (DL_FUNC) &Comp_Pair_LogGauss2mem,         15},
    {"Comp_Pair_Beta2mem",         (DL_FUNC) &Comp_Pair_Beta2mem,         15},
    {"Comp_Pair_Kumaraswamy2mem",         (DL_FUNC) &Comp_Pair_Kumaraswamy2mem,         15},
    {"Comp_Pair_Kumaraswamy22mem",         (DL_FUNC) &Comp_Pair_Kumaraswamy22mem,         15},
    {"Comp_Pair_Gauss_misp_Pois2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_Pois2mem,         15},
    {"Comp_Pair_BinomNNGauss_misp2mem",         (DL_FUNC) &Comp_Pair_BinomNNGauss_misp2mem,         15},
    {"Comp_Pair_Gauss_misp_PoisGamma2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_PoisGamma2mem,         15},
    {"Comp_Pair_PoisGamma2mem",         (DL_FUNC) &Comp_Pair_PoisGamma2mem,         15},
    {"Comp_Pair_PoisGammaZIP2mem",         (DL_FUNC) &Comp_Pair_PoisGammaZIP2mem,         15},
    {"Comp_Pair_Pois2mem",         (DL_FUNC) &Comp_Pair_Pois2mem,         15},
    {"Comp_Pair_BinomGauss2mem",         (DL_FUNC) &Comp_Pair_BinomGauss2mem,         15},
    {"Comp_Pair_BinomLogi2mem",         (DL_FUNC) &Comp_Pair_BinomLogi2mem,         15},
    {"Comp_Pair_BinomNNGauss2mem",         (DL_FUNC) &Comp_Pair_BinomNNGauss2mem,         15},
    {"Comp_Pair_BinomNNLogi2mem",         (DL_FUNC) &Comp_Pair_BinomNNLogi2mem,         15},
    {"Comp_Pair_BinomnegGauss2mem",         (DL_FUNC) &Comp_Pair_BinomnegGauss2mem,         15},
    {"Comp_Pair_TWOPIECETukeyh2mem",         (DL_FUNC) &Comp_Pair_TWOPIECETukeyh2mem,         15},
    {"Comp_Pair_TWOPIECET2mem",         (DL_FUNC) &Comp_Pair_TWOPIECET2mem,         15},
    {"Comp_Pair_TWOPIECEGauss2mem",         (DL_FUNC) &Comp_Pair_TWOPIECEGauss2mem,         15},
    {"Comp_Pair_TWOPIECEBIMODAL2mem",         (DL_FUNC) &Comp_Pair_TWOPIECEBIMODAL2mem,         15},
    {"Comp_Pair_BinomnegGaussZINB2mem",         (DL_FUNC) &Comp_Pair_BinomnegGaussZINB2mem,         15},
    {"Comp_Pair_PoisZIP2mem",         (DL_FUNC) &Comp_Pair_PoisZIP2mem,         15},
    {"Comp_Pair_Gauss_misp_PoisZIP2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_PoisZIP2mem,         15},
    {"Comp_Pair_LogLogistic2mem",         (DL_FUNC) &Comp_Pair_LogLogistic2mem,         15},
    {"Comp_Pair_Logistic2mem",         (DL_FUNC) &Comp_Pair_Logistic2mem,         15},
    {"Comp_Pair_Gauss_st2mem",         (DL_FUNC) &Comp_Pair_Gauss_st2mem,         15},
    {"Comp_Pair_T_st2mem",         (DL_FUNC) &Comp_Pair_T_st2mem,         15},
    {"Comp_Pair_Gauss_misp_T_st2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_T_st2mem,         15},
    {"Comp_Pair_Tukeyh_st2mem",         (DL_FUNC) &Comp_Pair_Tukeyh_st2mem,         15},
    {"Comp_Pair_Tukeyhh_st2mem",         (DL_FUNC) &Comp_Pair_Tukeyhh_st2mem,         15},
    {"Comp_Pair_SkewGauss_st2mem",         (DL_FUNC) &Comp_Pair_SkewGauss_st2mem,         15},
    {"Comp_Pair_SinhGauss_st2mem",         (DL_FUNC) &Comp_Pair_SinhGauss_st2mem,         15},
    {"Comp_Pair_Gamma_st2mem",         (DL_FUNC) &Comp_Pair_Gamma_st2mem,         15},
    {"Comp_Pair_Weibull_st2mem",         (DL_FUNC) &Comp_Pair_Weibull_st2mem,         15},
    {"Comp_Pair_LogGauss_st2mem",         (DL_FUNC) &Comp_Pair_LogGauss_st2mem,         15},
    {"Comp_Pair_Beta_st2mem",         (DL_FUNC) &Comp_Pair_Beta_st2mem,         15},
    {"Comp_Pair_Kumaraswamy_st2mem",         (DL_FUNC) &Comp_Pair_Kumaraswamy_st2mem,         15},//85
    {"Comp_Pair_Kumaraswamy2_st2mem",         (DL_FUNC) &Comp_Pair_Kumaraswamy2_st2mem,         15},
    {"Comp_Pair_Gauss_misp_Pois_st2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_Pois_st2mem,         15},
    {"Comp_Pair_Gauss_misp_PoisGamma_st2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_PoisGamma_st2mem,         15},
    {"Comp_Pair_PoisGamma_st2mem",         (DL_FUNC) &Comp_Pair_PoisGamma_st2mem,         15},
    {"Comp_Pair_Pois_st2mem",         (DL_FUNC) &Comp_Pair_Pois_st2mem,         15},
    {"Comp_Pair_BinomGauss_st2mem",         (DL_FUNC) &Comp_Pair_BinomGauss_st2mem,         15},
    {"Comp_Pair_BinomLogi_st2mem",         (DL_FUNC) &Comp_Pair_BinomLogi_st2mem,         15},
    {"Comp_Pair_BinomNNGauss_st2mem",         (DL_FUNC) &Comp_Pair_BinomNNGauss_st2mem,         15},
    {"Comp_Pair_BinomNNLogi_st2mem",         (DL_FUNC) &Comp_Pair_BinomNNLogi_st2mem,         15},
    {"Comp_Pair_BinomnegBinary2mem",         (DL_FUNC) &Comp_Pair_BinomnegBinary2mem,         15},
    {"Comp_Pair_BinomnegGauss_st2mem",         (DL_FUNC) &Comp_Pair_BinomnegGauss_st2mem,         15},
    {"Comp_Pair_BinomnegLogi_st2mem",         (DL_FUNC) &Comp_Pair_BinomnegLogi_st2mem,         15},//96
    {"Comp_Pair_TWOPIECETukeyh_st2mem",         (DL_FUNC) &Comp_Pair_TWOPIECETukeyh_st2mem,         15},
    {"Comp_Pair_TWOPIECET_st2mem",         (DL_FUNC) &Comp_Pair_TWOPIECET_st2mem,         15},
    {"Comp_Pair_TWOPIECEGauss_st2mem",         (DL_FUNC) &Comp_Pair_TWOPIECEGauss_st2mem,         15},
    {"Comp_Pair_TWOPIECEBIMODAL_st2mem",         (DL_FUNC) &Comp_Pair_TWOPIECEBIMODAL_st2mem,         15},
    {"Comp_Pair_BinomnegGaussZINB_st2mem",         (DL_FUNC) &Comp_Pair_BinomnegGaussZINB_st2mem,         15},
    {"Comp_Pair_PoisZIP_st2mem",         (DL_FUNC) &Comp_Pair_PoisZIP_st2mem,         15},
    {"Comp_Pair_Gauss_misp_PoisZIP_st2mem",         (DL_FUNC) &Comp_Pair_Gauss_misp_PoisZIP_st2mem,         15},
    {"Comp_Pair_LogLogistic_st2mem",         (DL_FUNC) &Comp_Pair_LogLogistic_st2mem,         15},
    {"Comp_Pair_Logistic_st2mem",         (DL_FUNC) &Comp_Pair_Logistic_st2mem,         15},
    {"Comp_Pair_WrapGauss_st2mem",         (DL_FUNC) &Comp_Pair_WrapGauss_st2mem,         15},
  /*********************** for copula ***************************************************/
    {"Comp_Pair_GaussCop2mem",         (DL_FUNC) &Comp_Pair_GaussCop2mem,         14},
    {"Comp_Pair_LogGaussCop2mem",         (DL_FUNC) &Comp_Pair_LogGaussCop2mem,         14},
    {"Comp_Pair_LogisticCop2mem",         (DL_FUNC) &Comp_Pair_LogisticCop2mem,         14},
    {"Comp_Pair_PoisCop2mem",         (DL_FUNC) &Comp_Pair_PoisCop2mem,         14},
    {"Comp_Pair_WeibullCop2mem",         (DL_FUNC) &Comp_Pair_WeibullCop2mem,         14},
    {"Comp_Pair_GammaCop2mem",         (DL_FUNC) &Comp_Pair_GammaCop2mem,         14},
    {"Comp_Pair_BinomNNGaussCop2mem",         (DL_FUNC) &Comp_Pair_BinomNNGaussCop2mem,         14}, 
    {"Comp_Pair_BinomnegGaussCop2mem",         (DL_FUNC) &Comp_Pair_BinomnegGaussCop2mem,         14},    
    {"Comp_Pair_TCop2mem",         (DL_FUNC) &Comp_Pair_TCop2mem,         14},
    {"Comp_Pair_BetaCop2mem",         (DL_FUNC) &Comp_Pair_BetaCop2mem,         14},
    {"Comp_Pair_Beta2Cop2mem",         (DL_FUNC) &Comp_Pair_Beta2Cop2mem,         14},
    {"Comp_Pair_KumaraswamyCop2mem",         (DL_FUNC) &Comp_Pair_KumaraswamyCop2mem,         14},
    {"Comp_Pair_Kumaraswamy2Cop2mem",         (DL_FUNC) &Comp_Pair_Kumaraswamy2Cop2mem,         14},//90
    /*********************** for conditional composite likelihood ***************************************************/
    {"Comp_Cond_Gauss2mem",         (DL_FUNC) &Comp_Cond_Gauss2mem,         15},
    {"Comp_Cond_Tukeyh2mem",         (DL_FUNC) &Comp_Cond_Tukeyh2mem,         15},
    {"Comp_Cond_Tukeyhh2mem",         (DL_FUNC) &Comp_Cond_Tukeyhh2mem,         15},
    {"Comp_Cond_SkewGauss2mem",         (DL_FUNC) &Comp_Cond_SkewGauss2mem,         15},
    {"Comp_Cond_WrapGauss2mem",         (DL_FUNC) &Comp_Cond_WrapGauss2mem,         15},
    {"Comp_Cond_T2mem",         (DL_FUNC) &Comp_Cond_T2mem,         15},
    {"Comp_Cond_Gauss_misp_T2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_T2mem,         15},
    {"Comp_Cond_Gauss_misp_SkewT2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_SkewT2mem,         15},
    {"Comp_Cond_Gauss_misp_Tukeygh2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_Tukeygh2mem,         15},
    {"Comp_Cond_SinhGauss2mem",         (DL_FUNC) &Comp_Cond_SinhGauss2mem,         15},
    {"Comp_Cond_Gamma2mem",         (DL_FUNC) &Comp_Cond_Gamma2mem,         15},
    {"Comp_Cond_Weibull2mem",         (DL_FUNC) &Comp_Cond_Weibull2mem,         15},
    {"Comp_Cond_LogGauss2mem",         (DL_FUNC) &Comp_Cond_LogGauss2mem,         15},
    {"Comp_Cond_Beta2mem",         (DL_FUNC) &Comp_Cond_Beta2mem,         15},
    {"Comp_Cond_Kumaraswamy2mem",         (DL_FUNC) &Comp_Cond_Kumaraswamy2mem,         15},
    {"Comp_Cond_Kumaraswamy22mem",         (DL_FUNC) &Comp_Cond_Kumaraswamy22mem,         15},
    {"Comp_Cond_Gauss_misp_Pois2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_Pois2mem,         15},
    {"Comp_Cond_BinomNNGauss_misp2mem",         (DL_FUNC) &Comp_Cond_BinomNNGauss_misp2mem,         15},
    {"Comp_Cond_Gauss_misp_PoisGamma2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_PoisGamma2mem,         15},
    {"Comp_Cond_PoisGamma2mem",         (DL_FUNC) &Comp_Cond_PoisGamma2mem,         15},
    {"Comp_Cond_PoisGammaZIP2mem",         (DL_FUNC) &Comp_Cond_PoisGammaZIP2mem,         15},
    {"Comp_Cond_Pois2mem",         (DL_FUNC) &Comp_Cond_Pois2mem,         15},
    {"Comp_Cond_BinomGauss2mem",         (DL_FUNC) &Comp_Cond_BinomGauss2mem,         15},
    {"Comp_Cond_BinomLogi2mem",         (DL_FUNC) &Comp_Cond_BinomLogi2mem,         15},
    {"Comp_Cond_BinomNNGauss2mem",         (DL_FUNC) &Comp_Cond_BinomNNGauss2mem,         15},
    {"Comp_Cond_BinomnegBinary2mem",         (DL_FUNC) &Comp_Cond_BinomnegBinary2mem,         15},
    {"Comp_Cond_BinomNNLogi2mem",         (DL_FUNC) &Comp_Cond_BinomNNLogi2mem,         15},
    {"Comp_Cond_BinomnegGauss2mem",         (DL_FUNC) &Comp_Cond_BinomnegGauss2mem,         15},
    {"Comp_Cond_TWOPIECETukeyh2mem",         (DL_FUNC) &Comp_Cond_TWOPIECETukeyh2mem,         15},
    {"Comp_Cond_TWOPIECET2mem",         (DL_FUNC) &Comp_Cond_TWOPIECET2mem,         15},
    {"Comp_Cond_TWOPIECEGauss2mem",         (DL_FUNC) &Comp_Cond_TWOPIECEGauss2mem,         15},
    {"Comp_Cond_TWOPIECEBIMODAL2mem",         (DL_FUNC) &Comp_Cond_TWOPIECEBIMODAL2mem,         15},
    {"Comp_Cond_BinomnegGaussZINB2mem",         (DL_FUNC) &Comp_Cond_BinomnegGaussZINB2mem,         15},
    {"Comp_Cond_PoisZIP2mem",         (DL_FUNC) &Comp_Cond_PoisZIP2mem,         15},
    {"Comp_Cond_Gauss_misp_PoisZIP2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_PoisZIP2mem,         15},
    {"Comp_Cond_LogLogistic2mem",         (DL_FUNC) &Comp_Cond_LogLogistic2mem,         15},
    {"Comp_Cond_Logistic2mem",         (DL_FUNC) &Comp_Cond_Logistic2mem,         15},
    {"Comp_Cond_Gauss_st2mem",         (DL_FUNC) &Comp_Cond_Gauss_st2mem,         15},
    {"Comp_Cond_T_st2mem",         (DL_FUNC) &Comp_Cond_T_st2mem,         15},
    {"Comp_Cond_Gauss_misp_T_st2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_T_st2mem,         15},
    {"Comp_Cond_Tukeyh_st2mem",         (DL_FUNC) &Comp_Cond_Tukeyh_st2mem,         15},
    {"Comp_Cond_Tukeyhh_st2mem",         (DL_FUNC) &Comp_Cond_Tukeyhh_st2mem,         15},
    {"Comp_Cond_SkewGauss_st2mem",         (DL_FUNC) &Comp_Cond_SkewGauss_st2mem,         15},
    {"Comp_Cond_SinhGauss_st2mem",         (DL_FUNC) &Comp_Cond_SinhGauss_st2mem,         15},
    {"Comp_Cond_Gamma_st2mem",         (DL_FUNC) &Comp_Cond_Gamma_st2mem,         15},
    {"Comp_Cond_Weibull_st2mem",         (DL_FUNC) &Comp_Cond_Weibull_st2mem,         15},
    {"Comp_Cond_LogGauss_st2mem",         (DL_FUNC) &Comp_Cond_LogGauss_st2mem,         15},
    {"Comp_Cond_Beta_st2mem",         (DL_FUNC) &Comp_Cond_Beta_st2mem,         15},
    {"Comp_Cond_Kumaraswamy_st2mem",         (DL_FUNC) &Comp_Cond_Kumaraswamy_st2mem,         15},//85
    {"Comp_Cond_Kumaraswamy2_st2mem",         (DL_FUNC) &Comp_Cond_Kumaraswamy2_st2mem,         15},
    {"Comp_Cond_Gauss_misp_Pois_st2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_Pois_st2mem,         15},
    {"Comp_Cond_Gauss_misp_PoisGamma_st2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_PoisGamma_st2mem,         15},
    {"Comp_Cond_PoisGamma_st2mem",         (DL_FUNC) &Comp_Cond_PoisGamma_st2mem,         15},
    {"Comp_Cond_Pois_st2mem",         (DL_FUNC) &Comp_Cond_Pois_st2mem,         15},
    {"Comp_Cond_BinomGauss_st2mem",         (DL_FUNC) &Comp_Cond_BinomGauss_st2mem,         15},
    {"Comp_Cond_BinomLogi_st2mem",         (DL_FUNC) &Comp_Cond_BinomLogi_st2mem,         15},
    {"Comp_Cond_BinomNNGauss_st2mem",         (DL_FUNC) &Comp_Cond_BinomNNGauss_st2mem,         15},
    {"Comp_Cond_BinomNNLogi_st2mem",         (DL_FUNC) &Comp_Cond_BinomNNLogi_st2mem,         15},
    {"Comp_Cond_BinomnegGauss_st2mem",         (DL_FUNC) &Comp_Cond_BinomnegGauss_st2mem,         15},
    {"Comp_Cond_BinomnegLogi_st2mem",         (DL_FUNC) &Comp_Cond_BinomnegLogi_st2mem,         15},//96
    {"Comp_Cond_TWOPIECETukeyh_st2mem",         (DL_FUNC) &Comp_Cond_TWOPIECETukeyh_st2mem,         15},
    {"Comp_Cond_TWOPIECET_st2mem",         (DL_FUNC) &Comp_Cond_TWOPIECET_st2mem,         15},
    {"Comp_Cond_TWOPIECEGauss_st2mem",         (DL_FUNC) &Comp_Cond_TWOPIECEGauss_st2mem,         15},
    {"Comp_Cond_TWOPIECEBIMODAL_st2mem",         (DL_FUNC) &Comp_Cond_TWOPIECEBIMODAL_st2mem,         15},
    {"Comp_Cond_BinomnegGaussZINB_st2mem",         (DL_FUNC) &Comp_Cond_BinomnegGaussZINB_st2mem,         15},
    {"Comp_Cond_PoisZIP_st2mem",         (DL_FUNC) &Comp_Cond_PoisZIP_st2mem,         15},
    {"Comp_Cond_Gauss_misp_PoisZIP_st2mem",         (DL_FUNC) &Comp_Cond_Gauss_misp_PoisZIP_st2mem,         15},
    {"Comp_Cond_LogLogistic_st2mem",         (DL_FUNC) &Comp_Cond_LogLogistic_st2mem,         15},
    {"Comp_Cond_Logistic_st2mem",         (DL_FUNC) &Comp_Cond_Logistic_st2mem,         15},
    {"Comp_Cond_WrapGauss_st2mem",         (DL_FUNC) &Comp_Cond_WrapGauss_st2mem,         15},

    /********************** spatial anisotropic +************************************/
    {"Comp_Pair_Gauss2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss2mem_aniso,         17},
    {"Comp_Pair_Pois2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss2mem_aniso,         17},
    {"Comp_Diff_Gauss2mem_aniso",         (DL_FUNC) &Comp_Diff_Gauss2mem_aniso,         17},
    {"Comp_Pair_WrapGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_WrapGauss2mem_aniso,         17},
    {"Comp_Pair_SinhGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_SinhGauss2mem_aniso,         17},
    {"Comp_Pair_SkewGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_SkewGauss2mem_aniso,         17},
    {"Comp_Pair_Gamma2mem_aniso",         (DL_FUNC) &Comp_Pair_Gamma2mem_aniso,         17},
    {"Comp_Pair_Weibull2mem_aniso",         (DL_FUNC) &Comp_Pair_Weibull2mem_aniso,         17},
    {"Comp_Pair_Kumaraswamy22mem_aniso",         (DL_FUNC) &Comp_Pair_Kumaraswamy22mem_aniso,         17},
    {"Comp_Pair_Beta2mem_aniso",         (DL_FUNC) &Comp_Pair_Beta2mem_aniso,         17},
    {"Comp_Pair_LogGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_LogGauss2mem_aniso,         17},
    {"Comp_Pair_PoisbinnegGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_PoisbinnegGauss2mem_aniso,         17},
    {"Comp_Pair_PoisbinGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_PoisbinGauss2mem_aniso,         17},
    {"Comp_Pair_BinomnegGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomnegGauss2mem_aniso,         17},//109
    {"Comp_Pair_BinomnegLogi2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomnegLogi2mem_aniso,         17},
    {"Comp_Pair_BinomnegGaussZINB2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomnegGaussZINB2mem_aniso,         17},
    {"Comp_Pair_BinomGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomGauss2mem_aniso,         17},//112
    {"Comp_Pair_BinomLogi2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomLogi2mem_aniso,         17},
    {"Comp_Pair_BinomNNGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomNNGauss2mem_aniso,         17},
    {"Comp_Pair_BinomNNGauss_misp2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomNNGauss_misp2mem_aniso,         17},
    {"Comp_Pair_BinomNNLogi2mem_aniso",         (DL_FUNC) &Comp_Pair_BinomNNLogi2mem_aniso,         17},
    {"Comp_Pair_LogLogistic2mem_aniso",         (DL_FUNC) &Comp_Pair_LogLogistic2mem_aniso,         17},
    {"Comp_Pair_Logistic2mem_aniso",         (DL_FUNC) &Comp_Pair_Logistic2mem_aniso,         17},
    {"Comp_Pair_PoisGamma2mem_aniso",         (DL_FUNC) &Comp_Pair_PoisGamma2mem_aniso,         17},
    {"Comp_Pair_Gauss_misp_PoisGamma2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss_misp_PoisGamma2mem_aniso,         17},
    {"Comp_Pair_PoisZIP2mem_aniso",         (DL_FUNC) &Comp_Pair_PoisZIP2mem_aniso,         17},
    {"Comp_Pair_Gauss_misp_Pois2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss_misp_Pois2mem_aniso,         17},
    {"Comp_Pair_Gauss_misp_PoisZIP2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss_misp_PoisZIP2mem_aniso,         17},
    {"Comp_Pair_Gauss_misp_SkewT2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss_misp_SkewT2mem_aniso,         17},
    {"Comp_Pair_Gauss_misp_T2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss_misp_T2mem_aniso,         17},
    {"Comp_Pair_Tukeyhh2mem_aniso",         (DL_FUNC) &Comp_Pair_Tukeyhh2mem_aniso,         17},
    {"Comp_Pair_Tukeyh2mem_aniso",         (DL_FUNC) &Comp_Pair_Tukeyh2mem_aniso,         17},
    {"Comp_Pair_Gauss_misp_Tukeygh2mem_aniso",         (DL_FUNC) &Comp_Pair_Gauss_misp_Tukeygh2mem_aniso,         17},
    {"Comp_Pair_T2mem_aniso",         (DL_FUNC) &Comp_Pair_T2mem_aniso,         17},
    {"Comp_Pair_TWOPIECETukeyh2mem_aniso",         (DL_FUNC) &Comp_Pair_TWOPIECETukeyh2mem_aniso,         17},
    {"Comp_Pair_TWOPIECET2mem_aniso",         (DL_FUNC) &Comp_Pair_TWOPIECET2mem_aniso,         17},
    {"Comp_Pair_TWOPIECEGauss2mem_aniso",         (DL_FUNC) &Comp_Pair_TWOPIECEGauss2mem_aniso,         17},
    {"Comp_Pair_TWOPIECEBIMODAL2mem_aniso",         (DL_FUNC) &Comp_Pair_TWOPIECEBIMODAL2mem_aniso,         17},
    {"Comp_Pair_GaussCop2mem_aniso",         (DL_FUNC) &Comp_Pair_GaussCop2mem_aniso,         17},
    {"Comp_Pair_TCop2mem_aniso",         (DL_FUNC) &Comp_Pair_TCop2mem_aniso,         17},
    {"Comp_Pair_BetaCop2mem_aniso",         (DL_FUNC) &Comp_Pair_BetaCop2mem_aniso,         17},
    {"Comp_Pair_Beta2Cop2mem_aniso",         (DL_FUNC) &Comp_Pair_Beta2Cop2mem_aniso,         17},
    {"Comp_Pair_KumaraswamyCop2mem_aniso",         (DL_FUNC) &Comp_Pair_KumaraswamyCop2mem_aniso,         17},
    {"Comp_Pair_Kumaraswamy2Cop2mem_aniso",         (DL_FUNC) &Comp_Pair_Kumaraswamy2Cop2mem_aniso,         17},
    {"Comp_Cond_Gauss2mem_aniso",         (DL_FUNC) &Comp_Cond_Gauss2mem_aniso,         17},
    {"Comp_Cond_Tukeyhh2mem_aniso",         (DL_FUNC) &Comp_Cond_Tukeyhh2mem_aniso,         17},
    {"Comp_Cond_SkewGauss2mem_aniso",         (DL_FUNC) &Comp_Cond_SkewGauss2mem_aniso,         17},
    {"Comp_Cond_T2mem_aniso",         (DL_FUNC) &Comp_Cond_T2mem_aniso,         17},
    {"Comp_Cond_Gauss_misp_T2mem_aniso",         (DL_FUNC) &Comp_Cond_Gauss_misp_T2mem_aniso,         17},
    {"Comp_Cond_Gauss_misp_SkewT2mem_aniso",         (DL_FUNC) &Comp_Cond_Gauss_misp_SkewT2mem_aniso,         17},
    {"Comp_Cond_SinhGauss2mem_aniso",         (DL_FUNC) &Comp_Cond_SinhGauss2mem_aniso,         17},
    {"Comp_Cond_Gamma2mem_aniso",         (DL_FUNC) & Comp_Cond_Gamma2mem_aniso,         17},
    {"Comp_Cond_Weibull2mem_aniso",         (DL_FUNC) &Comp_Cond_Weibull2mem_aniso,         17},
    {"Comp_Cond_LogGauss2mem_aniso",         (DL_FUNC) &Comp_Cond_LogGauss2mem_aniso,         17},
    {"Comp_Cond_Beta2mem_aniso",         (DL_FUNC) &Comp_Cond_Beta2mem_aniso,         17},
    {"Comp_Cond_Kumaraswamy2mem_aniso",         (DL_FUNC) &Comp_Cond_Kumaraswamy2mem_aniso,         17},
    {"Comp_Cond_Gauss_misp_Pois2mem_aniso",         (DL_FUNC) &Comp_Cond_Gauss_misp_Pois2mem_aniso,         17},
    {"Comp_Cond_BinomNNGauss_misp2mem_aniso",         (DL_FUNC) &Comp_Cond_BinomNNGauss_misp2mem_aniso,         17},
    {"Comp_Cond_Gauss_misp_PoisGamma2mem_aniso",         (DL_FUNC) &Comp_Cond_Gauss_misp_PoisGamma2mem_aniso,         17},
    {"Comp_Cond_PoisGamma2mem_aniso",         (DL_FUNC) &Comp_Cond_PoisGamma2mem_aniso,         17},
    {"Comp_Cond_Pois2mem_aniso",         (DL_FUNC) &Comp_Cond_Pois2mem_aniso,         17},
    {"Comp_Cond_BinomGauss2mem_aniso",         (DL_FUNC) &Comp_Cond_BinomGauss2mem_aniso,         17},
    {"Comp_Cond_BinomLogi2mem_aniso",         (DL_FUNC) &Comp_Cond_BinomLogi2mem_aniso,         17},
    {"Comp_Cond_BinomNNGauss2mem_aniso",         (DL_FUNC) &Comp_Cond_BinomNNGauss2mem_aniso,         17},
    {"Comp_Cond_BinomNNLogi2mem_aniso",         (DL_FUNC) & Comp_Cond_BinomNNLogi2mem_aniso,         17},
    {"Comp_Cond_BinomnegGauss2mem_aniso",         (DL_FUNC) &Comp_Cond_BinomnegGauss2mem_aniso,         17},
    {"Comp_Cond_TWOPIECETukeyh2mem_aniso",         (DL_FUNC) &Comp_Cond_TWOPIECETukeyh2mem_aniso,         17},
    {"Comp_Cond_TWOPIECET2mem_aniso",         (DL_FUNC) &Comp_Cond_TWOPIECET2mem_aniso,         17},
    {"Comp_Cond_TWOPIECEGauss2mem_aniso",         (DL_FUNC) &Comp_Cond_TWOPIECEGauss2mem_aniso,         17},
    {"Comp_Cond_TWOPIECEBIMODAL2mem_aniso",         (DL_FUNC) &Comp_Cond_TWOPIECEBIMODAL2mem_aniso,         17},
    {"Comp_Cond_BinomnegGaussZINB2mem_aniso",         (DL_FUNC) &Comp_Cond_BinomnegGaussZINB2mem_aniso,         17},
    {"Comp_Cond_PoisZIP2mem_aniso",         (DL_FUNC) &Comp_Cond_PoisZIP2mem_aniso,         17},
    {"Comp_Cond_LogLogistic2mem_aniso",         (DL_FUNC) &Comp_Cond_LogLogistic2mem_aniso,         17},
    {"Comp_Cond_Logistic2mem_aniso",         (DL_FUNC) &Comp_Cond_Logistic2mem_aniso,         17},
/*************************/
    {"Comp_Pair_Gauss_biv2mem",         (DL_FUNC) &Comp_Pair_Gauss_biv2mem,         12},
    {NULL, NULL, 0}
};

void R_init_GeoModels(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
void R_unload_GeoModels(DllInfo *info) {
  // just to avoid warning from compiler
}
