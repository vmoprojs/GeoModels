#define LOW -1.0e15
#define MAXERR 1e-6
//#define EPS DBL_EPSILON
double CorFunWitMat(double lag, double scale, double smooth);
double Dist_chordal(double loni, double lati, double lonj, double latj,double radius);
double Dist_geodesic(double loni, double lati, double lonj, double latj,double radius);
double dist(int type_dist,double coordx,double locx,double coordy,double locy,double radius);
double CorFunCauchy(double lag, double R_power2, double scale);
double CorFunStable(double lag, double R_power, double scale);
double CorFunDagum(double lag, double R_power1, double R_power2, double scale);
double CorFunGenCauchy(double lag, double R_power1, double R_power2, double scale);
double CorFunSferical(double lag, double scale);
double CorFunW0(double lag,double scale,double smoo);
double CorFunW1(double lag,double scale,double smoo);
double CorFunW2(double lag,double scale,double smoo);
double CorFunWave(double lag, double scale);
double CorFunWitMat(double lag, double scale, double smooth);
double CorFct(int cormod, double h, double u, double par0,double par1,double par2,double par3, int c11, int c22);
double Variogram(int cormod, double h, double u, double nugget, double var, double par0,double par1,double par2,double par3);
double CorFunBohman(double lag,double scale);
double biv_sinh(double corr,double zi,double zj,double skew,double tail,double vari);
double d2lognorm(double x, double y, double sill, double mux,double muy,double rho);
double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape);

double pnorm_OCL(double x, double mu, double sd);
double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11);
double pbnorm(int cormod, double h, double u, double mean1, double mean2, double nugget, double var,double par0,double par1,double par2,double par3, double thr);
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11);
double biv_binom(int NN, int u, int v, double p01,double p10,double p11);
double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11);
double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11);

// ===================================== Bessel  =====================================//
/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2014  R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* Constants und Documentation that apply to several of the
 * ./bessel_[ijky].c  files */

/* *******************************************************************
 Explanation of machine-dependent constants
 beta	  = Radix for the floating-point system
 minexp = Smallest representable power of beta
 maxexp = Smallest power of beta that overflows
 it = p = Number of bits (base-beta digits) in the mantissa
 (significand) of a working precision (floating-point) variable
 NSIG	  = Decimal significance desired.  Should be set to
 INT(LOG10(2)*it+1).	 Setting NSIG lower will result
 in decreased accuracy while setting NSIG higher will
 increase CPU time without increasing accuracy.  The
 truncation error is limited to a relative error of
 T=.5*10^(-NSIG).
 ENTEN  = 10 ^ K, where K is the largest int such that
 ENTEN is machine-representable in working precision
 ENSIG  = 10 ^ NSIG
 RTNSIG = 10 ^ (-K) for the smallest int K such that K >= NSIG/4
 ENMTEN = Smallest ABS(X) such that X/4 does not underflow
 XINF	  = Largest positive machine number; approximately beta ^ maxexp
 == DBL_MAX (defined in  #include <float.h>)
 SQXMIN = Square root of beta ^ minexp = sqrt(DBL_MIN)
 EPS	  = The smallest positive floating-point number such that 1.0+EPS > 1.0
 = beta ^ (-p)	 == DBL_EPSILON
 For I :
 EXPARG = Largest working precision argument that the library
 EXP routine can handle and upper limit on the
 magnitude of X when IZE=1; approximately LOG(beta ^ maxexp)
 For I and J :
 xlrg_IJ = xlrg_BESS_IJ (was = XLARGE). Upper limit on the magnitude of X
 (when IZE=2 for I()).  Bear in mind that if floor(abs(x)) =: N, then
 at least N iterations of the backward recursion will be executed.
 The value of 10 ^ 4 was used till Feb.2009, when it was increased
 to 10 ^ 5 (= 1e5).
 For j :
 XMIN_J  = Smallest acceptable argument for RBESY; approximately
 max(2*beta ^ minexp, 2/XINF), rounded up
 For Y :
 xlrg_Y =  (was = XLARGE). Upper bound on X;
 approximately 1/DEL, because the sine and cosine functions
 have lost about half of their precision at that point.
 EPS_SINC = Machine number below which sin(x)/x = 1; approximately SQRT(EPS).
 THRESH = Lower bound for use of the asymptotic form;
 approximately AINT(-LOG10(EPS/2.0))+1.0
 For K :
 xmax_k =  (was = XMAX). Upper limit on the magnitude of X when ize = 1;
 i.e. maximal x for UNscaled answer.
 Solution to equation:
 W(X) * (1 -1/8 X + 9/128 X^2) = beta ^ minexp
 where  W(X) = EXP(-X)*SQRT(PI/2X)
 --------------------------------------------------------------------
 Approximate values for some important machines are:
 beta minexp maxexp it NSIG ENTEN ENSIG RTNSIG ENMTEN	 EXPARG
 IEEE (IBM/XT,
 SUN, etc.) (S.P.)  2	  -126	128  24	  8  1e38   1e8	  1e-2	4.70e-38     88
 IEEE	(...) (D.P.)  2	 -1022 1024  53	 16  1e308  1e16  1e-4	8.90e-308   709
 CRAY-1	      (S.P.)  2	 -8193 8191  48	 15  1e2465 1e15  1e-4	1.84e-2466 5677
 Cyber 180/855
 under NOS  (S.P.)  2	  -975 1070  48	 15  1e322  1e15  1e-4	1.25e-293   741
 IBM 3033     (D.P.) 16	   -65	 63  14	  5  1e75   1e5	  1e-2	2.16e-78    174
 VAX	      (S.P.)  2	  -128	127  24	  8  1e38   1e8	  1e-2	1.17e-38     88
 VAX D-Format (D.P.)  2	  -128	127  56	 17  1e38   1e17  1e-5	1.17e-38     88
 VAX G-Format (D.P.)  2	 -1024 1023  53	 16  1e307  1e16  1e-4	2.22e-308   709
 And routine specific :
 xlrg_IJ xlrg_Y xmax_k EPS_SINC XMIN_J    XINF   THRESH
 IEEE (IBM/XT,
 SUN, etc.) (S.P.)	1e4  1e4   85.337  1e-4	 2.36e-38   3.40e38	8.
 IEEE	(...) (D.P.)	1e4  1e8  705.342  1e-8	 4.46e-308  1.79e308   16.
 CRAY-1	      (S.P.)	1e4  2e7 5674.858  5e-8	 3.67e-2466 5.45e2465  15.
 Cyber 180/855
 under NOS  (S.P.)	1e4  2e7  672.788  5e-8	 6.28e-294  1.26e322   15.
 IBM 3033     (D.P.)	1e4  1e8  177.852  1e-8	 2.77e-76   7.23e75    17.
 VAX	      (S.P.)	1e4  1e4   86.715  1e-4	 1.18e-38   1.70e38	8.
 VAX e-Format (D.P.)	1e4  1e9   86.715  1e-9	 1.18e-38   1.70e38    17.
 VAX G-Format (D.P.)	1e4  1e8  706.728  1e-8	 2.23e-308  8.98e307   16.
 */
#define nsig_BESS	16
#define ensig_BESS	1e16
#define rtnsig_BESS	1e-4
#define enmten_BESS	8.9e-308
#define enten_BESS	1e308

#define exparg_BESS	709.
#define xlrg_BESS_IJ	1e5
#define xlrg_BESS_Y	1e8
#define thresh_BESS_Y	16.

#define xmax_BESS_K	705.342/* maximal x for UNscaled answer */
#define sqxmin_BESS_K	1.49e-154


#define M_eps_sinc	2.149e-8



#define DBL_MAX 1.7976931348623157e+308
#define DBL_EPSILON 2.2204460492503131e-16
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#define DBL_MIN 2.2250738585072014e-308
#define ML_NAN 0.0/0.0
#define POSINF 1.0 /0.0
#define ML_NEGINF -1.0 /0.0
#define HSQRT 1.414213562373095048801688724209698078569671

#define min0(x, y) (((x) <= (y)) ? (x) : (y))
#define max0(x, y) (((x) <= (y)) ? (y) : (x))

typedef void integr_fn(double *x, int n, void *ex);
typedef int Rboolean;
enum { FALSE, TRUE };



double bessel_ii(double x, double alpha, double expo);
double bessel_kk(double x, double alpha, double expo);
static void II_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bi, int *ncalc);
static void KK_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bk, int *ncalc);

double bessel_ii(double x, double alpha, double expo)
{
    int nb, ncalc, ize;
    double na, *bi;
    /*#ifndef MATHLIB_STANDALONE
     const void *vmax;
     #endif*/
    
    /*#ifdef IEEE_754
     // NaNs propagated correctly
     if (isnan(x) || isnan(alpha)) return x + alpha;
     #endif*/
    if (x < 0) {
        printf("ME_RANGE bessel_i");
        return ML_NAN;
    }
    ize = (int)expo;
    
    na = floor(alpha);
    if (alpha < 0) {
        /* Using Abramowitz & Stegun  9.6.2 & 9.6.6
         * this may not be quite optimal (CPU and accuracy wise) */
        return(bessel_ii(x, -alpha, expo) +
               ((alpha == na) ? /* sin(pi * alpha) = 0 */ 0 :
                bessel_kk(x, -alpha, expo) *
                ((ize == 1)? 2. : 2.*exp(-2.*x))/M_PI_F * sin(-alpha*M_PI_F)));
    }
    nb = 1 + (int)na;/* nb-1 <= alpha < nb */
    alpha -= (double)(nb-1);
    
    /*#ifdef MATHLIB_STANDALONE
     bi = (double *) calloc(nb, sizeof(double));
     if (!bi) printf("bessel_i allocation error\n");
     #else*/
    //vmax = vmaxget();
    bi = (double *) calloc( nb, sizeof(double));
    //#endif
    //printf("bessel_i(x=%g): ncalc (=%d) != nb (=%d); alpha=%g. ANTES ize: %d\n",
    //     x, ncalc, nb, alpha, ize);
    II_bessel(&x, &alpha, &nb, &ize, bi, &ncalc);
    if(ncalc != nb) {/* error input */
        if(ncalc < 0)
            printf("bessel_i(%g): ncalc (=%d) != nb (=%d); alpha=%g. Arg. out of range?\n",
                   x, ncalc, nb, alpha);
        else
            printf("bessel_i(%g,nu=%g): precision lost in result\n",
                   x, alpha+(double)nb-1);
    }
    x = bi[nb-1];
    /*#ifdef MATHLIB_STANDALONE
     free(bi);
     #else
     //vmaxset(vmax);
     #endif*/
    return x;
}





double bessel_kk(double x, double alpha, double expo)
{
    int nb, ncalc, ize;
    double *bk;
    /*#ifndef MATHLIB_STANDALONE
     const void *vmax;
     #endif*/
    
    /*#ifdef IEEE_754
     // NaNs propagated correctly
     if (ISNAN(x) || ISNAN(alpha)) return x + alpha;
     #endif*/
    if (x < 0) {
        printf("ME_RANGE bessel_k");
        return ML_NAN;
    }
    ize = (int)expo;
    if(alpha < 0)
        alpha = -alpha;
    nb = 1+ (int)floor(alpha);/* nb-1 <= |alpha| < nb */
    alpha -= (double)(nb-1);
    /*#ifdef MATHLIB_STANDALONE
     bk = (double *) calloc(nb, sizeof(double));
     if (!bk) MATHLIB_ERROR("%s", _("bessel_k allocation error"));
     #else*/
    //vmax = vmaxget();
    bk = (double *) calloc( nb, sizeof(double));
    //#endif
    KK_bessel(&x, &alpha, &nb, &ize, bk, &ncalc);
    if(ncalc != nb) {/* error input */
        if(ncalc < 0)
            printf("bessel_k(%g): ncalc (=%d) != nb (=%d); alpha=%g. Arg. out of range?\n",
                   x, ncalc, nb, alpha);
        else
            printf("bessel_k(%g,nu=%g): precision lost in result\n",
                   x, alpha+(double)nb-1);
    }
    x = bk[nb-1];
    /*#ifdef MATHLIB_STANDALONE
     free(bk);
     #else
     vmaxset(vmax);
     #endif*/
    return x;
}


static void II_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bi, int *ncalc)
{
    
    /* -------------------------------------------------------------------
     This routine calculates Bessel functions I_(N+ALPHA) (X)
     for non-negative argument X, and non-negative order N+ALPHA,
     with or without exponential scaling.
     Explanation of variables in the calling sequence
     X     - Non-negative argument for which
     I's or exponentially scaled I's (I*EXP(-X))
     are to be calculated.	If I's are to be calculated,
     X must be less than exparg_BESS (IZE=1) or xlrg_BESS_IJ (IZE=2),
     (see bessel.h).
     ALPHA - Fractional part of order for which
     I's or exponentially scaled I's (I*EXP(-X)) are
     to be calculated.  0 <= ALPHA < 1.0.
     NB    - Number of functions to be calculated, NB > 0.
     The first function calculated is of order ALPHA, and the
     last is of order (NB - 1 + ALPHA).
     IZE   - Type.	IZE = 1 if unscaled I's are to be calculated,
     = 2 if exponentially scaled I's are to be calculated.
     BI    - Output vector of length NB.	If the routine
     terminates normally (NCALC=NB), the vector BI contains the
     functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
     corresponding exponentially scaled functions.
     NCALC - Output variable indicating possible errors.
     Before using the vector BI, the user should check that
     NCALC=NB, i.e., all orders have been calculated to
     the desired accuracy.	See error returns below.
     *******************************************************************
     *******************************************************************
     Error returns
     In case of an error,	NCALC != NB, and not all I's are
     calculated to the desired accuracy.
     NCALC < 0:  An argument is out of range. For example,
     NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >= EXPARG_BESS.
     In this case, the BI-vector is not calculated, and NCALC is
     set to MIN0(NB,0)-1 so that NCALC != NB.
     NB > NCALC > 0: Not all requested function values could
     be calculated accurately.	This usually occurs because NB is
     much larger than ABS(X).  In this case, BI[N] is calculated
     to the desired accuracy for N <= NCALC, but precision
     is lost for NCALC < N <= NB.  If BI[N] does not vanish
     for N > NCALC (because it is too small to be represented),
     and BI[N]/BI[NCALC] = 10**(-K), then only the first NSIG-K
     significant figures of BI[N] can be trusted.
     Intrinsic functions required are:
     DBLE, EXP, gamma_cody, INT, MAX, MIN, REAL, SQRT
     Acknowledgement
     This program is based on a program written by David J.
     Sookne (2) that computes values of the Bessel functions J or
     I of float argument and long order.  Modifications include
     the restriction of the computation to the I Bessel function
     of non-negative float argument, the extension of the computation
     to arbitrary positive order, the inclusion of optional
     exponential scaling, and the elimination of most underflow.
     An earlier version was published in (3).
     References: "A Note on Backward Recurrence Algorithms," Olver,
     F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
     pp 941-947.
     "Bessel Functions of Real Argument and Integer Order,"
     Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
     125-132.
     "ALGORITHM 597, Sequence of Modified Bessel Functions
     of the First Kind," Cody, W. J., Trans. Math. Soft.,
     1983, pp. 242-245.
     Latest modification: May 30, 1989
     Modified by: W. J. Cody and L. Stoltz
     Applied Mathematics Division
     Argonne National Laboratory
     Argonne, IL  60439
     */
    
    /*-------------------------------------------------------------------
     Mathematical constants
     -------------------------------------------------------------------*/
    double const__ = 1.585;
    
    /* Local variables */
    int nend, intx, nbmx, k, l, n, nstart;
    double pold, test,	p, em, en, empal, emp2al, halfx,
    aa, bb, cc, psave, plast, tover, psavel, sum, nu, twonu;
    
    /*Parameter adjustments */
    --bi;
    nu = *alpha;
    twonu = nu + nu;
    
    /*-------------------------------------------------------------------
     Check for X, NB, OR IZE out of range.
     ------------------------------------------------------------------- */
    
    if (*nb > 0 && *x >= 0. &&	(0. <= nu && nu < 1.) &&
        (1 <= *ize && *ize <= 2) ) {
        
        *ncalc = *nb;
        if(*ize == 1 && *x > exparg_BESS) {
            for(k=1; k <= *nb; k++)
                bi[k]=POSINF; /* the limit *is* = Inf */
            return;
            printf("bessel_i(%g): ncalc (=%d) != nb (=%d); alpha=%g. PILAS!!!\n",
                   *x, *ncalc, *nb, *alpha);
        }
        if(*ize == 2 && *x > xlrg_BESS_IJ) {
            for(k=1; k <= *nb; k++)
                bi[k]= 0.; /* The limit exp(-x) * I_nu(x) --> 0 : */
            return;
        }
        intx = (int) (*x);/* fine, since *x <= xlrg_BESS_IJ <<< LONG_MAX */
        if (*x >= rtnsig_BESS) { /* "non-small" x ( >= 1e-4 ) */
            
            /* -------------------------------------------------------------------
             Initialize the forward sweep, the P-sequence of Olver
             ------------------------------------------------------------------- */
            nbmx = *nb - intx;
            n = intx + 1;
            en = (double) (n + n) + twonu;
            plast = 1.;
            p = en / *x;
            /* ------------------------------------------------
             Calculate general significance test
             ------------------------------------------------ */
            test = ensig_BESS + ensig_BESS;
            if (intx << 1 > nsig_BESS * 5) {
                test = sqrt(test * p);
            } else {
                test /= pow(const__, intx);
            }
            if (nbmx >= 3) {
                /* --------------------------------------------------
                 Calculate P-sequence until N = NB-1
                 Check for possible overflow.
                 ------------------------------------------------ */
                tover = enten_BESS / ensig_BESS;
                nstart = intx + 2;
                nend = *nb - 1;
                for (k = nstart; k <= nend; ++k)
                {
                    n = k;
                    en += 2.;
                    pold = plast;
                    plast = p;
                    p = en * plast / *x + pold;
                    if (p > tover) {
                        /* ------------------------------------------------
                         To avoid overflow, divide P-sequence by TOVER.
                         Calculate P-sequence until ABS(P) > 1.
                         ---------------------------------------------- */
                        tover = enten_BESS;
                        p /= tover;
                        plast /= tover;
                        psave = p;
                        psavel = plast;
                        nstart = n + 1;
                        do {
                            ++n;
                            en += 2.;
                            pold = plast;
                            plast = p;
                            p = en * plast / *x + pold;
                        }
                        while (p <= 1.);
                        
                        bb = en / *x;
                        /* ------------------------------------------------
                         Calculate backward test, and find NCALC,
                         the highest N such that the test is passed.
                         ------------------------------------------------ */
                        test = pold * plast / ensig_BESS;
                        test *= .5 - .5 / (bb * bb);
                        p = plast * tover;
                        --n;
                        en -= 2.;
                        nend = min0(*nb,n);
                        for (l = nstart; l <= nend; ++l) {
                            *ncalc = l;
                            pold = psavel;
                            psavel = psave;
                            psave = en * psavel / *x + pold;
                            if (psave * psavel > test) {
                                goto L90;
                            }
                        }
                        *ncalc = nend + 1;
                    L90:
                        --(*ncalc);
                        goto L120;
                    }
                }
                n = nend;
                en = (double)(n + n) + twonu;
                /*---------------------------------------------------
                 Calculate special significance test for NBMX > 2.
                 --------------------------------------------------- */
                test = fmax(test,sqrt(plast * ensig_BESS) * sqrt(p + p));
            }
            /* --------------------------------------------------------
             Calculate P-sequence until significance test passed.
             -------------------------------------------------------- */
            do {
                ++n;
                en += 2.;
                pold = plast;
                plast = p;
                p = en * plast / *x + pold;
            } while (p < test);
            
        L120:
            /* -------------------------------------------------------------------
             Initialize the backward recursion and the normalization sum.
             ------------------------------------------------------------------- */
            ++n;
            en += 2.;
            bb = 0.;
            aa = 1. / p;
            em = (double) n - 1.;
            empal = em + nu;
            emp2al = em - 1. + twonu;
            sum = aa * empal * emp2al / em;
            nend = n - *nb;
            if (nend < 0) {
                /* -----------------------------------------------------
                 N < NB, so store BI[N] and set higher orders to 0..
                 ----------------------------------------------------- */
                bi[n] = aa;
                nend = -nend;
                for (l = 1; l <= nend; ++l) {
                    bi[n + l] = 0.;
                }
            } else {
                if (nend > 0) {
                    /* -----------------------------------------------------
                     Recur backward via difference equation,
                     calculating (but not storing) BI[N], until N = NB.
                     --------------------------------------------------- */
                    
                    for (l = 1; l <= nend; ++l) {
                        --n;
                        en -= 2.;
                        cc = bb;
                        bb = aa;
                        /* for x ~= 1500,  sum would overflow to 'inf' here,
                         * and the final bi[] /= sum would give 0 wrongly;
                         * RE-normalize (aa, sum) here -- no need to undo */
                        if(nend > 100 && aa > 1e200) {
                            /* multiply by  2^-900 = 1.18e-271 */
                            cc	= ldexp(cc, -900);
                            bb	= ldexp(bb, -900);
                            sum = ldexp(sum,-900);
                        }
                        aa = en * bb / *x + cc;
                        em -= 1.;
                        emp2al -= 1.;
                        if (n == 1) {
                            break;
                        }
                        if (n == 2) {
                            emp2al = 1.;
                        }
                        empal -= 1.;
                        sum = (sum + aa * empal) * emp2al / em;
                    }
                }
                /* ---------------------------------------------------
                 Store BI[NB]
                 --------------------------------------------------- */
                bi[n] = aa;
                if (*nb <= 1) {
                    sum = sum + sum + aa;
                    goto L230;
                }
                /* -------------------------------------------------
                 Calculate and Store BI[NB-1]
                 ------------------------------------------------- */
                --n;
                en -= 2.;
                bi[n] = en * aa / *x + bb;
                if (n == 1) {
                    goto L220;
                }
                em -= 1.;
                if (n == 2)
                    emp2al = 1.;
                else
                    emp2al -= 1.;
                
                empal -= 1.;
                sum = (sum + bi[n] * empal) * emp2al / em;
            }
            nend = n - 2;
            if (nend > 0) {
                /* --------------------------------------------
                 Calculate via difference equation
                 and store BI[N], until N = 2.
                 ------------------------------------------ */
                for (l = 1; l <= nend; ++l) {
                    --n;
                    en -= 2.;
                    bi[n] = en * bi[n + 1] / *x + bi[n + 2];
                    em -= 1.;
                    if (n == 2)
                        emp2al = 1.;
                    else
                        emp2al -= 1.;
                    empal -= 1.;
                    sum = (sum + bi[n] * empal) * emp2al / em;
                }
            }
            /* ----------------------------------------------
             Calculate BI[1]
             -------------------------------------------- */
            bi[1] = 2. * empal * bi[2] / *x + bi[3];
        L220:
            sum = sum + sum + bi[1];
            
        L230:
            /* ---------------------------------------------------------
             Normalize.  Divide all BI[N] by sum.
             --------------------------------------------------------- */
            if (nu != 0.)
                sum *= (tgamma(1. + nu) * pow(*x * .5, -nu));
            if (*ize == 1)
                sum *= exp(-(*x));
            aa = enmten_BESS;
            if (sum > 1.)
                aa *= sum;
            for (n = 1; n <= *nb; ++n) {
                if (bi[n] < aa)
                    bi[n] = 0.;
                else
                    bi[n] /= sum;
            }
            return;
        } else { /* small x  < 1e-4 */
            
            /* -----------------------------------------------------------
             Two-term ascending series for small X.
             -----------------------------------------------------------*/
            aa = 1.;
            empal = 1. + nu;
            //#ifdef IEEE_754
            /* No need to check for underflow */
            //            halfx = .5 * *x;
            //#else
            if (*x > enmten_BESS)
                halfx = .5 * *x;
            else
                halfx = 0.;
            //#endif
            
            if (nu != 0.)
                aa = pow(halfx, nu) / tgamma(empal);
            if (*ize == 2)
                aa *= exp(-(*x));
            bb = halfx * halfx;
            bi[1] = aa + aa * bb / empal;
            if (*x != 0. && bi[1] == 0.)
                *ncalc = 0;
            if (*nb > 1) {
                if (*x == 0.) {
                    for (n = 2; n <= *nb; ++n)
                        bi[n] = 0.;
                } else {
                    /* -------------------------------------------------
                     Calculate higher-order functions.
                     ------------------------------------------------- */
                    cc = halfx;
                    tover = (enmten_BESS + enmten_BESS) / *x;
                    if (bb != 0.)
                        tover = enmten_BESS / bb;
                    for (n = 2; n <= *nb; ++n) {
                        aa /= empal;
                        empal += 1.;
                        aa *= cc;
                        if (aa <= tover * empal)
                            bi[n] = aa = 0.;
                        else
                            bi[n] = aa + aa * bb / empal;
                        if (bi[n] == 0. && *ncalc > n)
                            *ncalc = n - 1;
                    }
                }
            }
        }
    } else { /* argument out of range */
        *ncalc = min0(*nb,0) - 1;
        
    }
}
static void KK_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bk, int *ncalc)
{
    /*-------------------------------------------------------------------
     This routine calculates modified Bessel functions
     of the third kind, K_(N+ALPHA) (X), for non-negative
     argument X, and non-negative order N+ALPHA, with or without
     exponential scaling.
     Explanation of variables in the calling sequence
     X     - Non-negative argument for which
     K's or exponentially scaled K's (K*EXP(X))
     are to be calculated.	If K's are to be calculated,
     X must not be greater than XMAX_BESS_K.
     ALPHA - Fractional part of order for which
     K's or exponentially scaled K's (K*EXP(X)) are
     to be calculated.  0 <= ALPHA < 1.0.
     NB    - Number of functions to be calculated, NB > 0.
     The first function calculated is of order ALPHA, and the
     last is of order (NB - 1 + ALPHA).
     IZE   - Type.	IZE = 1 if unscaled K's are to be calculated,
     = 2 if exponentially scaled K's are to be calculated.
     BK    - Output vector of length NB.	If the
     routine terminates normally (NCALC=NB), the vector BK
     contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
     or the corresponding exponentially scaled functions.
     If (0 < NCALC < NB), BK(I) contains correct function
     values for I <= NCALC, and contains the ratios
     K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
     NCALC - Output variable indicating possible errors.
     Before using the vector BK, the user should check that
     NCALC=NB, i.e., all orders have been calculated to
     the desired accuracy.	See error returns below.
     *******************************************************************
     Error returns
     In case of an error, NCALC != NB, and not all K's are
     calculated to the desired accuracy.
     NCALC < -1:  An argument is out of range. For example,
     NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >= XMAX_BESS_K.
     In this case, the B-vector is not calculated,
     and NCALC is set to MIN0(NB,0)-2	 so that NCALC != NB.
     NCALC = -1:  Either  K(ALPHA,X) >= XINF  or
     K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) >= XINF.	 In this case,
     the B-vector is not calculated.	Note that again
     NCALC != NB.
     0 < NCALC < NB: Not all requested function values could
     be calculated accurately.  BK(I) contains correct function
     values for I <= NCALC, and contains the ratios
     K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
     Intrinsic functions required are:
     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT
     Acknowledgement
     This program is based on a program written by J. B. Campbell
     (2) that computes values of the Bessel functions K of float
     argument and float order.  Modifications include the addition
     of non-scaled functions, parameterization of machine
     dependencies, and the use of more accurate approximations
     for SINH and SIN.
     References: "On Temme's Algorithm for the Modified Bessel
     Functions of the Third Kind," Campbell, J. B.,
     TOMS 6(4), Dec. 1980, pp. 581-586.
     "A FORTRAN IV Subroutine for the Modified Bessel
     Functions of the Third Kind of Real Order and Real
     Argument," Campbell, J. B., Report NRC/ERB-925,
     National Research Council, Canada.
     Latest modification: May 30, 1989
     Modified by: W. J. Cody and L. Stoltz
     Applied Mathematics Division
     Argonne National Laboratory
     Argonne, IL  60439
     -------------------------------------------------------------------
     */
    /*---------------------------------------------------------------------
     * Mathematical constants
     *	A = LOG(2) - Euler's constant
     *	D = SQRT(2/PI)
     ---------------------------------------------------------------------*/
    double a = .11593151565841244881;
    
    /*---------------------------------------------------------------------
     P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA + Euler's constant
     Coefficients converted from hex to decimal and modified
     by W. J. Cody, 2/26/82 */
    double p[8] = { .805629875690432845,20.4045500205365151,
        157.705605106676174,536.671116469207504,900.382759291288778,
        730.923886650660393,229.299301509425145,.822467033424113231 };
    double q[7] = { 29.4601986247850434,277.577868510221208,
        1206.70325591027438,2762.91444159791519,3443.74050506564618,
        2210.63190113378647,572.267338359892221 };
    /* R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA) */
    double r[5] = { -.48672575865218401848,13.079485869097804016,
        -101.96490580880537526,347.65409106507813131,
        3.495898124521934782e-4 };
    double s[4] = { -25.579105509976461286,212.57260432226544008,
        -610.69018684944109624,422.69668805777760407 };
    /* T    - Approximation for SINH(Y)/Y */
    double t[6] = { 1.6125990452916363814e-10,
        2.5051878502858255354e-8,2.7557319615147964774e-6,
        1.9841269840928373686e-4,.0083333333333334751799,
        .16666666666666666446 };
    /*---------------------------------------------------------------------*/
    double estm[6] = { 52.0583,5.7607,2.7782,14.4303,185.3004, 9.3715 };
    double estf[7] = { 41.8341,7.1075,6.4306,42.511,1.35633,84.5096,20.};
    
    /* Local variables */
    int iend, i, j, k, m, ii, mplus1;
    double x2by4, twox, c, blpha, ratio, wminf;
    double d1, d2, d3, f0, f1, f2, p0, q0, t1, t2, twonu;
    double dm, ex, bk1, bk2, nu;
    
    ii = 0; /* -Wall */
    
    ex = *x;
    nu = *alpha;
    *ncalc = min0(*nb,0) - 2;
    if (*nb > 0 && (0. <= nu && nu < 1.) && (1 <= *ize && *ize <= 2)) {
        if(ex <= 0 || (*ize == 1 && ex > xmax_BESS_K)) {
            if(ex <= 0) {
                if(ex < 0) printf("ME_RANGE, K_bessel\n");
                for(i=0; i < *nb; i++)
                    bk[i] = POSINF;
            } else /* would only have underflow */
                for(i=0; i < *nb; i++)
                    bk[i] = 0.;
            *ncalc = *nb;
            return;
        }
        k = 0;
        if (nu < sqxmin_BESS_K) {
            nu = 0.;
        } else if (nu > .5) {
            k = 1;
            nu -= 1.;
        }
        twonu = nu + nu;
        iend = *nb + k - 1;
        c = nu * nu;
        d3 = -c;
        if (ex <= 1.) {
            /* ------------------------------------------------------------
             Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
             Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
             ------------------------------------------------------------ */
            d1 = 0.; d2 = p[0];
            t1 = 1.; t2 = q[0];
            for (i = 2; i <= 7; i += 2) {
                d1 = c * d1 + p[i - 1];
                d2 = c * d2 + p[i];
                t1 = c * t1 + q[i - 1];
                t2 = c * t2 + q[i];
            }
            d1 = nu * d1;
            t1 = nu * t1;
            f1 = log(ex);
            f0 = a + nu * (p[7] - nu * (d1 + d2) / (t1 + t2)) - f1;
            q0 = exp(-nu * (a - nu * (p[7] + nu * (d1-d2) / (t1-t2)) - f1));
            f1 = nu * f0;
            p0 = exp(f1);
            /* -----------------------------------------------------------
             Calculation of F0 =
             ----------------------------------------------------------- */
            d1 = r[4];
            t1 = 1.;
            for (i = 0; i < 4; ++i) {
                d1 = c * d1 + r[i];
                t1 = c * t1 + s[i];
            }
            /* d2 := sinh(f1)/ nu = sinh(f1)/(f1/f0)
             *	   = f0 * sinh(f1)/f1 */
            if (fabs(f1) <= .5) {
                f1 *= f1;
                d2 = 0.;
                for (i = 0; i < 6; ++i) {
                    d2 = f1 * d2 + t[i];
                }
                d2 = f0 + f0 * f1 * d2;
            } else {
                d2 = sinh(f1) / nu;
            }
            f0 = d2 - nu * d1 / (t1 * p0);
            if (ex <= 1e-10) {
                /* ---------------------------------------------------------
                 X <= 1.0E-10
                 Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
                 --------------------------------------------------------- */
                bk[0] = f0 + ex * f0;
                if (*ize == 1) {
                    bk[0] -= ex * bk[0];
                }
                ratio = p0 / f0;
                c = ex * DBL_MAX;
                if (k != 0) {
                    /* ---------------------------------------------------
                     Calculation of K(ALPHA,X)
                     and  X*K(ALPHA+1,X)/K(ALPHA,X),	ALPHA >= 1/2
                     --------------------------------------------------- */
                    *ncalc = -1;
                    if (bk[0] >= c / ratio) {
                        return;
                    }
                    bk[0] = ratio * bk[0] / ex;
                    twonu += 2.;
                    ratio = twonu;
                }
                *ncalc = 1;
                if (*nb == 1)
                    return;
                
                /* -----------------------------------------------------
                 Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),
                 L = 1, 2, ... , NB-1
                 ----------------------------------------------------- */
                *ncalc = -1;
                for (i = 1; i < *nb; ++i) {
                    if (ratio >= c)
                        return;
                    
                    bk[i] = ratio / ex;
                    twonu += 2.;
                    ratio = twonu;
                }
                *ncalc = 1;
                goto L420;
            } else {
                /* ------------------------------------------------------
                 10^-10 < X <= 1.0
                 ------------------------------------------------------ */
                c = 1.;
                x2by4 = ex * ex / 4.;
                p0 = .5 * p0;
                q0 = .5 * q0;
                d1 = -1.;
                d2 = 0.;
                bk1 = 0.;
                bk2 = 0.;
                f1 = f0;
                f2 = p0;
                do {
                    d1 += 2.;
                    d2 += 1.;
                    d3 = d1 + d3;
                    c = x2by4 * c / d2;
                    f0 = (d2 * f0 + p0 + q0) / d3;
                    p0 /= d2 - nu;
                    q0 /= d2 + nu;
                    t1 = c * f0;
                    t2 = c * (p0 - d2 * f0);
                    bk1 += t1;
                    bk2 += t2;
                } while (fabs(t1 / (f1 + bk1)) > DBL_EPSILON ||
                         fabs(t2 / (f2 + bk2)) > DBL_EPSILON);
                bk1 = f1 + bk1;
                bk2 = 2. * (f2 + bk2) / ex;
                if (*ize == 2) {
                    d1 = exp(ex);
                    bk1 *= d1;
                    bk2 *= d1;
                }
                wminf = estf[0] * ex + estf[1];
            }
        } else if (DBL_EPSILON * ex > 1.) {
            /* -------------------------------------------------
             X > 1./EPS
             ------------------------------------------------- */
            *ncalc = *nb;
            bk1 = 1. / (M_SQRT_2dPI * sqrt(ex));
            for (i = 0; i < *nb; ++i)
                bk[i] = bk1;
            return;
            
        } else {
            /* -------------------------------------------------------
             X > 1.0
             ------------------------------------------------------- */
            twox = ex + ex;
            blpha = 0.;
            ratio = 0.;
            if (ex <= 4.) {
                /* ----------------------------------------------------------
                 Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0
                 ----------------------------------------------------------*/
                d2 = trunc(estm[0] / ex + estm[1]);
                m = (int) d2;
                d1 = d2 + d2;
                d2 -= .5;
                d2 *= d2;
                for (i = 2; i <= m; ++i) {
                    d1 -= 2.;
                    d2 -= d1;
                    ratio = (d3 + d2) / (twox + d1 - ratio);
                }
                /* -----------------------------------------------------------
                 Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
                 recurrence and K(ALPHA,X) from the wronskian
                 -----------------------------------------------------------*/
                d2 = trunc(estm[2] * ex + estm[3]);
                m = (int) d2;
                c = fabs(nu);
                d3 = c + c;
                d1 = d3 - 1.;
                f1 = DBL_MIN;
                f0 = (2. * (c + d2) / ex + .5 * ex / (c + d2 + 1.)) * DBL_MIN;
                for (i = 3; i <= m; ++i) {
                    d2 -= 1.;
                    f2 = (d3 + d2 + d2) * f0;
                    blpha = (1. + d1 / d2) * (f2 + blpha);
                    f2 = f2 / ex + f1;
                    f1 = f0;
                    f0 = f2;
                }
                f1 = (d3 + 2.) * f0 / ex + f1;
                d1 = 0.;
                t1 = 1.;
                for (i = 1; i <= 7; ++i) {
                    d1 = c * d1 + p[i - 1];
                    t1 = c * t1 + q[i - 1];
                }
                p0 = exp(c * (a + c * (p[7] - c * d1 / t1) - log(ex))) / ex;
                f2 = (c + .5 - ratio) * f1 / ex;
                bk1 = p0 + (d3 * f0 - f2 + f0 + blpha) / (f2 + f1 + f0) * p0;
                if (*ize == 1) {
                    bk1 *= exp(-ex);
                }
                wminf = estf[2] * ex + estf[3];
            } else {
                /* ---------------------------------------------------------
                 Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by
                 backward recurrence, for  X > 4.0
                 ----------------------------------------------------------*/
                dm = trunc(estm[4] / ex + estm[5]);
                m = (int) dm;
                d2 = dm - .5;
                d2 *= d2;
                d1 = dm + dm;
                for (i = 2; i <= m; ++i) {
                    dm -= 1.;
                    d1 -= 2.;
                    d2 -= d1;
                    ratio = (d3 + d2) / (twox + d1 - ratio);
                    blpha = (ratio + ratio * blpha) / dm;
                }
                bk1 = 1. / ((M_SQRT_2dPI + M_SQRT_2dPI * blpha) * sqrt(ex));
                if (*ize == 1)
                    bk1 *= exp(-ex);
                wminf = estf[4] * (ex - fabs(ex - estf[6])) + estf[5];
            }
            /* ---------------------------------------------------------
             Calculation of K(ALPHA+1,X)
             from K(ALPHA,X) and  K(ALPHA+1,X)/K(ALPHA,X)
             --------------------------------------------------------- */
            bk2 = bk1 + bk1 * (nu + .5 - ratio) / ex;
        }
        /*--------------------------------------------------------------------
         Calculation of 'NCALC', K(ALPHA+I,X),	I  =  0, 1, ... , NCALC-1,
         &	  K(ALPHA+I,X)/K(ALPHA+I-1,X),	I = NCALC, NCALC+1, ... , NB-1
         -------------------------------------------------------------------*/
        *ncalc = *nb;
        bk[0] = bk1;
        if (iend == 0)
            return;
        
        j = 1 - k;
        if (j >= 0)
            bk[j] = bk2;
        
        if (iend == 1)
            return;
        
        m = min0((int) (wminf - nu),iend);
        for (i = 2; i <= m; ++i) {
            t1 = bk1;
            bk1 = bk2;
            twonu += 2.;
            if (ex < 1.) {
                if (bk1 >= DBL_MAX / twonu * ex)
                    break;
            } else {
                if (bk1 / ex >= DBL_MAX / twonu)
                    break;
            }
            bk2 = twonu / ex * bk1 + t1;
            ii = i;
            ++j;
            if (j >= 0) {
                bk[j] = bk2;
            }
        }
        
        m = ii;
        if (m == iend) {
            return;
        }
        ratio = bk2 / bk1;
        mplus1 = m + 1;
        *ncalc = -1;
        for (i = mplus1; i <= iend; ++i) {
            twonu += 2.;
            ratio = twonu / ex + 1./ratio;
            ++j;
            if (j >= 1) {
                bk[j] = ratio;
            } else {
                if (bk2 >= DBL_MAX / ratio)
                    return;
                
                bk2 *= ratio;
            }
        }
        *ncalc = max0(1, mplus1 - k);
        if (*ncalc == 1)
            bk[0] = bk2;
        if (*nb == 1)
            return;
        
    L420:
        for (i = *ncalc; i < *nb; ++i) { /* i == *ncalc */
            //#ifndef IEEE_754
            //            if (bk[i-1] >= DBL_MAX / bk[i])
            //                return;
            //#endif
            bk[i] *= bk[i-1];
            (*ncalc)++;
        }
    }
}

// ===================================== END Bessel  =====================================//

// ===================================== START Integrate  =====================================//
// https://github.com/wch/r-source/blob/trunk/src/appl/integrate.c

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2014  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  Most of this file is C translations of Fortran routines in
 *  QUADPACK: the latter is part of SLATEC 'and therefore in the public
 *  domain' (https://en.wikipedia.org/wiki/QUADPACK).
 *
 *
 */

//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif

//#include <math.h>
//#include <float.h>
//#include <Rmath.h> /* for fmax2, fmin2, imin2 */
//#include <R_ext/Applic.h> /* exporting the API , particularly */
/*--- typedef void integr_fn(double *x, int n, void *ex) ---
 * vectorizing function   f(x[1:n], ...) -> x[]  {overwriting x[]}.
 * Vectorization can be used to speed up the integrand
 * instead of calling it  n  times.
 */

/* f2c-ed translations + modifications of QUADPACK functions */



void Rdqags(integr_fn f, void *ex, double *a, double *b,
                double *epsabs, double *epsrel,
                double *result, double *abserr, int *neval, int *ier,
                int *limit, int *lenw, int *last, int *iwork, double *work);


static void rdqagse(integr_fn f, void *ex, double *, double *,
                        double *, double *, int *, double *, double *,
                        int *, int *, double *, double *, double *,
                        double *, int *, int *);

static void rdqk21(integr_fn f, void *ex,
                   double *, double *, double *, double *, double *, double *);
static void rdqelg(int *, double *, double *, double *, double *, int *);


static void rdqpsrt(int *, int *, int *, double *, double *, int *, int *);


/* Table of constant values */

//static double c_b6 = 0.; // QUITAR STATIC PARA OPENCL
//static double c_b7 = 1.; // QUITAR STATIC PARA OPENCL



void Rdqags(integr_fn f, void *ex, double *a, double *b,
                double *epsabs, double *epsrel,
                double *result, double *abserr, int *neval, int *ier,
                int *limit, int *lenw, int *last, int *iwork, double *work)
{
    int l1, l2, l3;
    
    //         check validity of limit and lenw.
    
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*limit < 1 || *lenw < *limit *4) return;
    
    //         prepare call for dqagse.
    
    l1 = *limit;
    l2 = *limit + l1;
    l3 = *limit + l2;
    
    rdqagse(f, ex, a, b, epsabs, epsrel, limit, result, abserr, neval, ier,
                work, &work[l1], &work[l2], &work[l3], iwork, last);
    
    return;
}

static void rdqagse(integr_fn f, void *ex, double *a, double *b, double *
                        epsabs, double *epsrel, int *limit, double *result,
                        double *abserr, int *neval, int *ier, double *alist,
                        double *blist, double *rlist, double *elist, int *
                        iord, int *last)
{
    // Local variables
    Rboolean noext, extrap;
    int k,ksgn, nres;
    int ierro;
    int ktmin, nrmax;
    int iroff1, iroff2, iroff3;
    int id;
    int numrl2;
    int jupbnd;
    int maxerr;
    double res3la[3];
    double rlist2[52];
    double abseps, area, area1, area2, area12, dres, epmach;
    double a1, a2, b1, b2, defabs, defab1, defab2, oflow, uflow, resabs, reseps;
    double error1, error2, erro12, errbnd, erlast, errmax, errsum;
    
    double correc = 0.0, erlarg = 0.0, ertest = 0.0, small = 0.0;
    
    
    // ===first executable statement  dqagse
    // Parameter adjustments
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist;
    
    // Function Body
    epmach = DBL_EPSILON;
    
    //            test on validity of parameters
    //            ------------------------------
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    alist[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.;
    elist[1] = 0.;
    if (*epsabs <= 0. && *epsrel < max0(epmach * 50., 5e-29)) {
        *ier = 6;
        return;
    }
    
    //           first approximation to the integral
    //           -----------------------------------
    
    uflow = DBL_MIN;
    oflow = DBL_MAX;
    ierro = 0;
    rdqk21(f, ex, a, b, result, abserr, &defabs, &resabs);
    
    //           test on accuracy.
    
    dres = fabs(*result);
    errbnd = max0(*epsabs, *epsrel * dres);
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    if (*abserr <= epmach * 100. * defabs && *abserr > errbnd)
        *ier = 2;
    if (*limit == 1)
        *ier = 1;
    if (*ier != 0 || (*abserr <= errbnd && *abserr != resabs)
        || *abserr == 0.) goto L140;
    
    //           initialization
    //           --------------
    
    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    numrl2 = 2;
    ktmin = 0;
    extrap = FALSE;
    noext = FALSE;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (1. - epmach * 50.) * defabs) {
        ksgn = 1;
    }
    
    //           main do-loop
    //           ------------
    
    for (*last = 2; *last <= *limit; ++(*last)) {
        
        //           bisect the subinterval with the nrmax-th largest error estimate.
        
        a1 = alist[maxerr];
        b1 = (alist[maxerr] + blist[maxerr]) * .5;
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        rdqk21(f, ex, &a1, &b1, &area1, &error1, &resabs, &defab1);
        rdqk21(f, ex, &a2, &b2, &area2, &error2, &resabs, &defab2);
        
        //           improve previous approximations to integral
         //and error and test for accuracy.
        
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if (!(defab1 == error1 || defab2 == error2)) {
            
            if (fabs(rlist[maxerr] - area12) <= fabs(area12) * 1e-5 &&
                erro12 >= errmax * .99) {
                if (extrap)
                    ++iroff2;
                else //if(! extrap)
                    ++iroff1;
            }
            if (*last > 10 && erro12 > errmax)
                ++iroff3;
        }
        rlist[maxerr] = area1;
        rlist[*last] = area2;
        errbnd = max0(*epsabs, *epsrel * fabs(area));
        
        //          test for roundoff error and eventually set error flag.
        
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
            *ier = 2;
        if (iroff2 >= 5)
            ierro = 3;
        
        //set error flag in the case that the number of subintervals equals limit.
        if (*last == *limit)
            *ier = 1;
        
        //          set error flag in the case of bad integrand behaviour
         //at a point of the integration range.
        
        if (max0(fabs(a1), fabs(b2)) <=
            (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) {
            *ier = 4;
        }
        
        //           append the newly-created intervals to the list.
        
        if (error2 > error1) {
            alist[maxerr] = a2;
            alist[*last] = a1;
            blist[*last] = b1;
            rlist[maxerr] = area2;
            rlist[*last] = area1;
            elist[maxerr] = error2;
            elist[*last] = error1;
        } else {
            alist[*last] = a2;
            blist[maxerr] = b1;
            blist[*last] = b2;
            elist[maxerr] = error1;
            elist[*last] = error2;
        }
        
        //           call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
         //with nrmax-th largest error estimate (to be bisected next).
        
        //L30:
        rdqpsrt(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
        
        if (errsum <= errbnd)   goto L115;// ===jump out of do-loop
        if (*ier != 0)		break;
        if (*last == 2)	{ // L80:
            small = fabs(*b - *a) * .375;
            erlarg = errsum;
            ertest = errbnd;
            rlist2[1] = area;	continue;
        }
        if (noext)		continue;
        
        erlarg -= erlast;
        if (fabs(b1 - a1) > small) {
            erlarg += erro12;
        }
        if (!extrap) {
            
            //         test whether the interval to be bisected next is the
            // smallest interval.
            
            if (fabs(blist[maxerr] - alist[maxerr]) > small) {
                continue;
            }
            extrap = TRUE;
            nrmax = 2;
        }
        
        if (ierro != 3 && erlarg > ertest) {
            
            //           the smallest interval has the largest error.
            // before bisecting decrease the sum of the errors over the
            // larger intervals (erlarg) and perform extrapolation.
            
            id = nrmax;
            jupbnd = *last;
            if (*last > *limit / 2 + 2) {
                jupbnd = *limit + 3 - *last;
            }
            for (k = id; k <= jupbnd; ++k) {
                maxerr = iord[nrmax];
                errmax = elist[maxerr];
                if (fabs(blist[maxerr] - alist[maxerr]) > small) {
                    goto L90;
                }
                ++nrmax;
                // L50:
            }
        }
        //           perform extrapolation.  L60:
        
        ++numrl2;
        rlist2[numrl2 - 1] = area;
        rdqelg(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
        ++ktmin;
        if (ktmin > 5 && *abserr < errsum * .001) {
            *ier = 5;
        }
        if (abseps < *abserr) {
            ktmin = 0;
            *abserr = abseps;
            *result = reseps;
            correc = erlarg;
            ertest = max0(*epsabs, *epsrel * fabs(reseps));
            if (*abserr <= ertest) {
                break;
            }
        }
        
        //           prepare bisection of the smallest interval.  L70:
        
        if (numrl2 == 1) {
            noext = TRUE;
        }
        if (*ier == 5) {
            break;
        }
        maxerr = iord[1];
        errmax = elist[maxerr];
        nrmax = 1;
        extrap = FALSE;
        small *= .5;
        erlarg = errsum;
    L90:
        ;
    }
    
    
    // L100:	set final result and error estimate.
    //		------------------------------------
    
    if (*abserr == oflow) 	goto L115;
    if (*ier + ierro == 0) 	goto L110;
    if (ierro == 3)
        *abserr += correc;
    if (*ier == 0)
        *ier = 3;
    if (*result == 0. || area == 0.) {
        if (*abserr > errsum) 	goto L115;
        if (area == 0.) 	goto L130;
    }
    else { // L105:
        if (*abserr / fabs(*result) > errsum / fabs(area))
            goto L115;
    }
    
L110://		test on divergence.
    if (ksgn == -1 && max0(fabs(*result), fabs(area)) <= defabs * .01) {
        goto L130;
    }
    if (.01 > *result / area || *result / area > 100. || errsum > fabs(area)) {
        *ier = 5;
    }
    goto L130;
    
L115://		compute global integral sum.
    *result = 0.;
    for (k = 1; k <= *last; ++k)
        *result += rlist[k];
    *abserr = errsum;
L130:
    if (*ier > 2)
        L140:
        *neval = *last * 42 - 21;
    return;
}


static void rdqelg(int *n, double *epstab, double *
                   result, double *abserr, double *res3la, int *nres)
{
    // Local variables
    int i__, indx, ib, ib2, ie, k1, k2, k3, num, newelm, limexp;
    double delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epmach, epsinf;
    double oflow, ss, res;
    double errA, err1, err2, err3, tol1, tol2, tol3;
    
   
    
    // ===first executable statement  dqelg
    // Parameter adjustments
    --res3la;
    --epstab;
    
    // Function Body
    epmach = DBL_EPSILON;
    oflow = DBL_MAX;
    ++(*nres);
    *abserr = oflow;
    *result = epstab[*n];
    if (*n < 3) {
        goto L100;
    }
    limexp = 50;
    epstab[*n + 2] = epstab[*n];
    newelm = (*n - 1) / 2;
    epstab[*n] = oflow;
    num = *n;
    k1 = *n;
    for (i__ = 1; i__ <= newelm; ++i__) {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab[k1 + 2];
        e0 = epstab[k3];
        e1 = epstab[k2];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2 - e1;
        err2 = fabs(delta2);
        tol2 = max0(fabs(e2), e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = max0(e1abs, fabs(e0)) * epmach;
        if (err2 <= tol2 && err3 <= tol3) {
            //           if e0, e1 and e2 are equal to within machine
            // accuracy, convergence is assumed.
            *result = res;//		result = e2
            *abserr = err2 + err3;//	abserr = fabs(e1-e0)+fabs(e2-e1)
            
            goto L100;	// ===jump out of do-loop
        }
        
        e3 = epstab[k1];
        epstab[k1] = e1;
        delta1 = e1 - e3;
        err1 = fabs(delta1);
        tol1 = max0(e1abs, fabs(e3)) * epmach;
        
        //          if two elements are very close to each other, omit
        // a part of the table by adjusting the value of n
        
        if (err1 > tol1 && err2 > tol2 && err3 > tol3) {
            ss = 1. / delta1 + 1. / delta2 - 1. / delta3;
            epsinf = fabs(ss * e1);
            
            //           test to detect irregular behaviour in the table, and
            // eventually omit a part of the table adjusting the value of n.
            
            if (epsinf > 1e-4) {
                goto L30;
            }
        }
        
        *n = i__ + i__ - 1;
        goto L50;// ===jump out of do-loop
        
        
    L30:// compute a new element and eventually adjust the value of result.
        
        res = e1 + 1. / ss;
        epstab[k1] = res;
        k1 += -2;
        errA = err2 + fabs(res - e2) + err3;
        if (errA <= *abserr) {
            *abserr = errA;
            *result = res;
        }
    }
    
    //           shift the table.
    
L50:
    if (*n == limexp) {
        *n = (limexp / 2 << 1) - 1;
    }
    
    if (num / 2 << 1 == num) ib = 2; else ib = 1;
    ie = newelm + 1;
    for (i__ = 1; i__ <= ie; ++i__) {
        ib2 = ib + 2;
        epstab[ib] = epstab[ib2];
        ib = ib2;
    }
    if (num != *n) {
        indx = num - *n + 1;
        for (i__ = 1; i__ <= *n; ++i__) {
            epstab[i__] = epstab[indx];
            ++indx;
        }
    }
    //L80:
    if (*nres >= 4) {
        // L90:
        *abserr = fabs(*result - res3la[3]) +
        fabs(*result - res3la[2]) +
        fabs(*result - res3la[1]);
        res3la[1] = res3la[2];
        res3la[2] = res3la[3];
        res3la[3] = *result;
    } else {
        res3la[*nres] = *result;
        *abserr = oflow;
    }
    
L100:// compute error estimate
    *abserr = max0(*abserr, epmach * 5. * fabs(*result));
    return;
}



static void  rdqk21(integr_fn f, void *ex, double *a, double *b, double *result,
                    double *abserr, double *resabs, double *resasc)
{
    // Initialized data
    
     double wg[5] = { .066671344308688137593568809893332,
        .149451349150580593145776339657697,
        .219086362515982043995534934228163,
        .269266719309996355091226921569469,
        .295524224714752870173892994651338 };
     double xgk[11] = { .995657163025808080735527280689003,
        .973906528517171720077964012084452,
        .930157491355708226001207180059508,
        .865063366688984510732096688423493,
        .780817726586416897063717578345042,
        .679409568299024406234327365114874,
        .562757134668604683339000099272694,
        .433395394129247190799265943165784,
        .294392862701460198131126603103866,
        .14887433898163121088482600112972,0. };
     double wgk[11] = { .011694638867371874278064396062192,
        .03255816230796472747881897245939,
        .05475589657435199603138130024458,
        .07503967481091995276704314091619,
        .093125454583697605535065465083366,
        .109387158802297641899210590325805,
        .123491976262065851077958109831074,
        .134709217311473325928054001771707,
        .142775938577060080797094273138717,
        .147739104901338491374841515972068,
        .149445554002916905664936468389821 };
    
    
    // Local variables
    double fv1[10], fv2[10], vec[21];
    double absc, resg, resk, fsum, fval1, fval2;
    double hlgth, centr, reskh, uflow;
    double fc, epmach, dhlgth;
    int j, jtw, jtwm1;
    
    
    
    // ===first executable statement  dqk21
    epmach = DBL_EPSILON;
    uflow = DBL_MIN;
    
    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    dhlgth = fabs(hlgth);
    
    //           compute the 21-point kronrod approximation to
    // the integral, and estimate the absolute error.
    
    resg = 0.;
    vec[0] = centr;
    for (j = 1; j <= 5; ++j) {
        jtw = j << 1;
        absc = hlgth * xgk[jtw - 1];
        vec[(j << 1) - 1] = centr - absc;
        // L5:
        vec[j * 2] = centr + absc;
    }
    for (j = 1; j <= 5; ++j) {
        jtwm1 = (j << 1) - 1;
        absc = hlgth * xgk[jtwm1 - 1];
        vec[(j << 1) + 9] = centr - absc;
        vec[(j << 1) + 10] = centr + absc;
    }
    f(vec, 21, ex);
    fc = vec[0];
    resk = wgk[10] * fc;
    *resabs = fabs(resk);
    for (j = 1; j <= 5; ++j) {
        jtw = j << 1;
        absc = hlgth * xgk[jtw - 1];
        fval1 = vec[(j << 1) - 1];
        fval2 = vec[j * 2];
        fv1[jtw - 1] = fval1;
        fv2[jtw - 1] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j - 1] * fsum;
        resk += wgk[jtw - 1] * fsum;
        *resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
        // L10:
    }
    for (j = 1; j <= 5; ++j) {
        jtwm1 = (j << 1) - 1;
        absc = hlgth * xgk[jtwm1 - 1];
        fval1 = vec[(j << 1) + 9];
        fval2 = vec[(j << 1) + 10];
        fv1[jtwm1 - 1] = fval1;
        fv2[jtwm1 - 1] = fval2;
        fsum = fval1 + fval2;
        resk += wgk[jtwm1 - 1] * fsum;
        *resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
        // L15:
    }
    reskh = resk * .5;
    *resasc = wgk[10] * fabs(fc - reskh);
    for (j = 1; j <= 10; ++j) {
        *resasc += wgk[j - 1] * (fabs(fv1[j - 1] - reskh) +
                                 fabs(fv2[j - 1] - reskh));
        // L20:
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs((resk - resg) * hlgth);
    if (*resasc != 0. && *abserr != 0.) {
        *abserr = *resasc * min0(1., pow(*abserr * 200. / *resasc, 1.5));
    }
    if (*resabs > uflow / (epmach * 50.)) {
        *abserr = max0(epmach * 50. * *resabs, *abserr);
    }
    return;
}



static void rdqpsrt(int *limit, int *last, int *maxerr,
                        double *ermax, double *elist, int *iord, int *nrmax)
{
    // Local variables
    int i, j, k, ido, jbnd, isucc, jupbn;
    double errmin, errmax;

    
    // Parameter adjustments
    --iord;
    --elist;
    
    // Function Body
    
    //           check whether the list contains more than
    // two error estimates.
    if (*last <= 2) {
        iord[1] = 1;
        iord[2] = 2;
        goto Last;
    }
    //           this part of the routine is only executed if, due to a
    // difficult integrand, subdivision increased the error
    // estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.
    
    errmax = elist[*maxerr];
    if (*nrmax > 1) {
        ido = *nrmax - 1;
        for (i = 1; i <= ido; ++i) {
            isucc = iord[*nrmax - 1];
            if (errmax <= elist[isucc])
                break; // out of for-loop
            iord[*nrmax] = isucc;
            --(*nrmax);
            // L20:
        }
    }
    
    //L30:       compute the number of elements in the list to be maintained
    // in descending order. this number depends on the number of
    // subdivisions still allowed.
    if (*last > *limit / 2 + 2)
        jupbn = *limit + 3 - *last;
    else
        jupbn = *last;
    
    errmin = elist[*last];
    
    //          insert errmax by traversing the list top-down,
    // starting comparison from the element elist(iord(nrmax+1)).
    
    jbnd = jupbn - 1;
    for (i = *nrmax + 1; i <= jbnd; ++i) {
        isucc = iord[i];
        if (errmax >= elist[isucc]) {// ===jump out of do-loop
            // L60: insert errmin by traversing the list bottom-up.
            iord[i - 1] = *maxerr;
            for (j = i, k = jbnd; j <= jbnd; j++, k--) {
                isucc = iord[k];
                if (errmin < elist[isucc]) {
                    // goto L80; ===jump out of do-loop
                    iord[k + 1] = *last;
                    goto Last;
                }
                iord[k + 1] = isucc;
            }
            iord[i] = *last;
            goto Last;
        }
        iord[i - 1] = isucc;
    }
    
    iord[jbnd] = *maxerr;
    iord[jupbn] = *last;
    
Last:// set maxerr and ermax.
    
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
    return;
}


double beta (double x, double y)
{
    return( (  tgamma(x)*tgamma(y)  )/ tgamma(x+y)  );
}

// integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp)
{
    double res=0.0,y;
    y=lag/supp;
    res=pow(1-x,mu-1)*pow(x*x-y*y,alpha)/beta(2*alpha+1,mu);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}

// function generalized wendland  to integrate

void integr_gen(double *x, int n, void *ex)
{
    int i;double mu,alpha,beta,y;
    mu =    ((double*)ex)[0];  //mu
    alpha = ((double*)ex)[1];  //alpha
    beta =     ((double*)ex)[2];  //csupp
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_gen(x[i],mu,alpha,y,beta);}
    return;
}

// function computing generalized wendland
double wendintegral(double x, double *param)
{
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = pow(DBL_EPSILON, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;		     // as instructed in WRE
    iwork =   (int *) calloc(subdiv, sizeof(int));  // idem
    work = (double *) calloc(lenw, sizeof(double)); // idem
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=x/param[2];
    upper=1;
    // Compute the integral
    if(x<=param[2]) {
        Rdqags(integr_gen, (void *) &ex,
                   &lower, &upper, &epsabs, &epsrel, &result,
                   &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
        
    }else   {result=0;}
    free(iwork);free(work);
    return(result);
}

// ===================================== END Integrate  =====================================//

// ===================================== Bivariate Normal  =====================================//

double Phi(double x);
double Phi2diag( double x, double a, double px, double pxs );
double Phi2help( double x, double y, double rho );
double Phi2( double x, double y, double rho );
double cdf_norm_OCL(double lim1,double lim2,double a11,double a12);


// https://www.jstatsoft.org/article/view/v052i10/v52i10.pdf

double pnorm_OCL(double x, double mu, double sd)
{
    double z = (x-mu)/(sd);
    return  (  (2-erfc(z/HSQRT))/2  )  ;
}

double Phi(double x)
{
    double val =(1+     (1-erfc(x/HSQRT) )    )/2;
    //printf("val: %f\n",val);
    return ( val );
}


double Phi2diag( double x, double a, double px, double pxs )
{
    if( a <= 0.0 ) return px;
    if( a >= 1.0 ) return px * px;
    double b = 2.0 - a, sqrt_ab = sqrt( a * b );
    double c1 = 6.36619772367581343e-001;
    double c2 = 1.25331413731550025;
    double c3 = 1.57079632679489662;
    double c4 = 1.591549430918953358e-001;
    double asr = ( a > 0.1 ? asin( 1.0 - a ) : acos( sqrt_ab ) );
    double comp = px * pxs;
    if( comp * ( 1.0 - a - c1 * asr ) < 5e-17 )
        return b * comp;
    double tmp = c2 * x;
    double alpha = a * x * x / b;
    double a_even = -tmp * a;
    double a_odd = -sqrt_ab * alpha;
    double beta = x * x;
    double b_even = tmp * sqrt_ab;
    double b_odd = sqrt_ab * beta;
    double delta = 2.0 * x * x / b;
    double d_even = ( 1.0 - a ) * c3 - asr;
    double d_odd = tmp * ( sqrt_ab - a );
    double res = 0.0, res_new = d_even + d_odd;
    int k = 2;
    while( res != res_new )
    {
        d_even = ( a_odd + b_odd + delta * d_even ) / k;
        a_even *= alpha / k;
        b_even *= beta / k;
        k++;
        a_odd *= alpha / k;
        b_odd *= beta / k;
        d_odd = ( a_even + b_even + delta * d_odd ) / k;
        k++;
        res = res_new;
        res_new += d_even + d_odd;
    }
    res *= exp( -x * x / b ) * c4;
    return max0( ( 1.0 + c1 * asr ) * comp, b * comp - max0( 0.0, res ) );
}


double Phi2help( double x, double y, double rho )
{
    double s = sqrt( ( 1.0 - rho ) * ( 1.0 + rho ) );
    double a = 0.0, b1 = -fabs( x ), b2 = 0.0;
    if( rho > 0.99 )
    {
        double tmp = sqrt( ( 1.0 - rho ) / ( 1.0 + rho ) );
        b2 = -fabs( ( x - y ) / s - x * tmp );
        a = pow( ( x - y ) / x / s - tmp,2 );
    }
    else if( rho < -0.99 )
    {
        double tmp = sqrt( ( 1.0 + rho ) / ( 1.0 - rho ) );
        b2 = -fabs( ( x + y ) / s - x * tmp );
        a = pow( ( x + y ) / x / s - tmp,2 );
    }
    else
    {
        b2 = -fabs( rho * x - y ) / s;
        a = pow( b2 / x ,2);
    }
    
    double p1 = Phi( b1 ), p2 = Phi( b2 ), q = 0.0;
    if( a <= 1.0 )
        q = 0.5 * Phi2diag( b1, 2.0 * a / ( 1.0 + a ), p1, p2 );
    else
        q = p1 * p2 - 0.5 * Phi2diag( b2, 2.0 / ( 1.0 + a ), p2, p1 );
    int c1 = ( y / x >= rho ), c2 = ( x < 0.0 ), c3 = c2 && ( y >= 0.0 );
    return ( c1 && c3 ? q - 0.5
            : c1 && c2 ? q
            : c1 ? 0.5 - p1 + q
            : c3 ? p1 - q - 0.5
            : c2 ? p1 - q
            : 0.5 - q );
}

double Phi2( double x, double y, double rho )
{
    if( ( 1.0 - rho ) * ( 1.0 + rho ) <= 0.0 )
    {
        if( rho > 0.0 )
        {
            return (Phi( min0( x, y ) ));
        }
        else
        {
            return (max0( 0.0, min0( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ));
        }
    }
    if( x == 0.0 && y == 0.0 )
    {
        if( rho > 0.0 )
        {
            return (Phi2diag( 0.0, 1.0 - rho, 0.5, 0.5 ));
        }
        else
        {
            return (0.5 - Phi2diag( 0.0, 1.0 + rho, 0.5, 0.5 ));
        }
    }
    
    return (max0( 0.0,
                 min0( 1.0,
                      Phi2help( x, y, rho ) + Phi2help( y, x, rho ) ) ));
}


double cdf_norm_OCL(double lim1,double lim2,double a11,double a12)
{
    double res=0;
    double  uppe[2]={lim1/sqrt(a11),lim2/sqrt(a11)}, corre[1] ={a12/a11};
    //int infin[2]={0,0};
    double auxil=1-pow(corre[0],2);
    double det=pow(a11,2)-pow(a12,2) ;
    //res= a11* sqrt(auxil/det) *  F77_CALL(bvnmvn)(lowe,uppe,infin,corre);
    res= a11* sqrt(auxil/det) *  Phi2(uppe[0],uppe[1],corre[0]);
    return(res);
}



// ===================================== Distance Functions  =====================================//

// Utility.c
double Dist_chordal(double loni, double lati, double lonj, double latj,double radius)
{
    double ai, bi, aj, bj, val=0.0;
    if (loni == lonj && lati == latj) return val;
    ai = (lati)*M_PI_F/180;
    bi = (loni)*M_PI_F/180;
    aj = (latj)*M_PI_F/180;
    bj = (lonj)*M_PI_F/180;
    val=radius  *sqrt(pow(cos(ai) * cos(bi)-cos(aj)  *cos(bj) ,2) +
                      pow(cos(ai) * sin(bi)-cos(aj) * sin(bj) ,2)+
                      pow(sin(ai)-sin(aj) ,2));
    return(val);
}

// Computes the Geodesic distance between to coordinates:
double Dist_geodesic(double loni, double lati, double lonj, double latj,double radius)
{
    double ai, bi, aj, bj, val=0.0,val2=0.0;
    if (loni == lonj && lati == latj) return val;
    ai = (lati)*M_PI_F/180;
    bi = (loni)*M_PI_F/180;
    aj = (latj)*M_PI_F/180;
    bj = (lonj)*M_PI_F/180;
    val = sin(ai) * sin(aj) + cos(ai) * cos(aj) * cos(bi - bj);
    if(val<= -1)  val2=M_PI_F*radius;
    if(val>=1) val2=0;
    val2 = acos(val)*radius;
    return(val2);
}
/*
 void GeoDist(double *coordx, double *coordy, int *ncoord, double *res,int *type_dist,double radius)
 {
 int i=0, h=0, j=0;
 
 for(i=0;i<(*ncoord-1);i++)
 for(j=(i+1);j<ncoord[0];j++){
 if(*type_dist==1) res[h]=Dist_chordal(coordx[i],coordy[i],coordx[j],coordy[j],radius);
 if(*type_dist==2) res[h]=Dist_geodesic(coordx[i],coordy[i],coordx[j],coordy[j,radius]);
 h++;}
 
 return;
 }
 */

double dist(int type_dist,double coordx,double locx,double coordy,double locy,double radius)
{
    double lags=0.0;
    
    if(type_dist==0) lags=hypot(coordx-locx,coordy-locy);                        //euclidean
    if(type_dist==2) lags=Dist_geodesic(coordx,coordy,locx,locy,radius);           //great circle
    if(type_dist==1) lags=Dist_chordal(coordx,coordy,locx,locy,radius);      //chordal
    
    return(lags);
}

// ===================================== CorrelationFunction.c  ==================================//

// Cauhcy class of correlation models:
double CorFunCauchy(double lag, double R_power2, double scale)
{
    double rho=0.0;
    // Computes the correlation:
    rho=pow((1+pow(lag/scale,2)),-R_power2/2);
    return rho;
}
// Stable class of correlation models:
double CorFunStable(double lag, double R_power, double scale)
{
    double rho=0.0;
    // Computes the correlation:
    rho=exp(-pow(lag/scale,R_power));
    return rho;
}
// Dagum:
double CorFunDagum(double lag, double R_power1, double R_power2, double scale)
{
    double rho=0.0;
    rho=1-pow(pow(lag/scale,R_power1)/(1+pow(lag/scale,R_power1)), R_power2/R_power1);
    return rho;
}
// Generalised Cauhcy class of correlation models:
double CorFunGenCauchy(double lag, double R_power1, double R_power2, double scale)
{
    double rho=0.0;
    rho=pow((1+pow(lag/scale,R_power1)), -R_power2/R_power1);
    return rho;
}

// Sferical class of correlation models:
double CorFunSferical(double lag, double scale)
{
    double rho=0.0;
    if(lag<=scale) rho=1-1.5*lag/scale+0.5*pow(lag/scale, 3);
    else rho=0;
    return rho;
}
/* wendland function alpha=0*/
double CorFunW0(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1) rho=  pow(1-x,smoo);
    else rho=0;
    return rho;
}
/* wendland function alpha=1*/
double CorFunW1(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)  rho=pow(1-x,smoo+1)*(1+(smoo+1)*x);
    else rho=0;
    return rho;
}
/* wendland function alpha=2*/
double CorFunW2(double lag,double scale,double smoo)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)  rho=pow(1-x,smoo+2)*(3+x*(3*smoo+6)+x*x*(pow(smoo,2)+4*smoo+3))/3;
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
// Whittle=matern class of correlation models:
double CorFunWitMat(double lag, double scale, double smooth)
{
    double rho=0.0;
    // Computes the correlation:
    if(lag==0) rho=1;
    else  rho=(pow(lag/scale,smooth)*bessel_kk(lag/scale,smooth,1))/(pow(2,smooth-1)*tgamma(smooth));
    return rho;
}

/* generalized wendland function*/
double CorFunW_gen(double lag,double R_power1,double smooth,double scale)  // mu alpha beta
{
    double rho=0.0,x=0;
    
    /* case alpha=0--Askey funcyion*/
    if(smooth==0) {
        x=lag/scale;
        if(x<=1) rho=pow(1-x,R_power1);
        else rho=0;
    }
    /* case alpha>0*/
    if(smooth>0) {
        x=lag;
        double *param;
        param=(double *) calloc(3,sizeof(double));
        param[0]=R_power1;param[1]=smooth;param[2]=scale;  //mu,alpha //beta
        rho=wendintegral(x,param);
        free(param);
    }
    return rho;
}


double CorFct(int cormod, double h, double u, double par0,double par1,double par2,double par3, int c11, int c22)
{
    double arg=0.0, col=0.0,R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0.0, R_power_t=0.0, var11=0.0, var22=0.0;
    double rho=0.0, sep=0, scale=0.0, smooth=0.0,smooth_s=0.0,smooth_t=0.0, scale_s=0.0, scale_t=0, x=0, nug11=0.0, nug22=0.0;
    double scale11=0.0, scale22=0.0, scale12=0.0, smoo11=0.0, smoo22=0.0, smoo12=0.0,R_power11=0.0, R_power22=0.0, R_power12=0.0;
    switch(cormod) // Correlation functions are in alphabetical order
    {
        case 1:// Cauchy correlation function
            R_power1=2;
            R_power2=par0;
            scale=par1;
            rho=CorFunCauchy(h, R_power2, scale);
            break;
        case 4:// Exponential correlation function
            R_power=1;
            scale=par0;
            rho=CorFunStable(h, R_power, scale);
            break;
        case 5: // Dagum
            R_power1=par0;
            R_power2=par1;
            scale=par2;
            rho=CorFunDagum(h, R_power1, R_power2, scale);
            break;
        case 6:// Gaussian correlation function
            R_power=2;
            scale=par0;
            rho=CorFunStable(h, R_power, scale);
            break;
        case 8: // Generalised Cuachy correlation function
            R_power1=par0;
            R_power2=par1;
            scale=par2;
            rho=CorFunGenCauchy(h, R_power1, R_power2, scale);
            break;
        case 10:// Sferical correlation function
            scale=par0;
            rho=CorFunSferical(h, scale);
            break;
        case 11://wen0
            R_power=par0;
            scale=par1;
            rho=CorFunW0(h,scale,R_power);
            break;
        case 13://wen1
            R_power=par0;
            scale=par1;
            rho=CorFunW1(h,scale,R_power);
            break;
        case 12:// Stable correlation function
            R_power=par0;
            scale=par1;
            rho=CorFunStable(h, R_power, scale);
            break;
        case 14://  Whittle-Matern correlation function
            scale=par0;
            smooth=par1;
            rho=CorFunWitMat(h, scale, smooth);
            break;
        case 15://wen2
            R_power=par0;
            scale=par1;
            rho=CorFunW2(h,scale,R_power);
            break;
        case 16: //wave
            scale=par0;
            rho=CorFunWave(h,scale);
            break;
        case 17://  multiquadric correlation function valid on sphere
            R_power=par0;
            scale=par1;
            rho=pow(1-R_power/2,2*scale)/pow(1+pow(R_power/2,2)-R_power*cos(h),scale);
            break;
        case 18://  sinsphere correlation function valid on sphere
            R_power=par0;
            rho=1-pow(sin(h/(2)),R_power);
            break;
        case 19: // Generalised wend correlation function
            R_power1=par0;
            scale=par1;
            smooth=par2;
            rho=CorFunW_gen(h, R_power1, smooth, scale);
            break;
            
        case 84:// Double exp:
            scale_s=par0;
            scale_t=par1;
            rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
            break;
        

    }
    return rho;
}

// Computes the spatio-temporal variogram:
double Variogram(int cormod, double h, double u, double nugget, double var, double par0,double par1,double par2,double par3)
{
    double vario=0.0;
    //Computes the variogram
    vario=nugget+var*(1-CorFct(cormod,h,u,par0,par1,par2,par3,0,0));
    return vario;
}
double CorFunBohman(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)
    {
        if (x>0) rho=(1-x)*(sin(2*M_PI_F*x)/(2*M_PI_F*x))+(1-cos(2*M_PI_F*x))/(2*M_PI_F*M_PI_F*x);
        else   rho=1;
    }
    else rho=0;
    return rho;
}


// ===================================== START: Distributions.c  ==================================//
double biv_sinh(double corr,double zi,double zj,double skew,double tail,double vari)
{
    double b1=0.0,b2=0.0,A=0.0,B=0.0,k=0.0,res=0.0,Z1,Z2;
    b1=tail * asinh(zi)-skew;
    b2=tail * asinh(zj)-skew;
    k=1-pow(corr,2);
    A=pow(2 * M_PI_F * pow(k,0.5) * vari,-1) * cosh(b1) * cosh(b2) * pow(tail,2)/sqrt((pow(zi,2)+1) * (pow(zj,2)+1));
    Z1=sinh(b1);Z2=sinh(b2);
    B=exp(- (Z1*Z1 + Z2*Z2 - 2*corr*Z1*Z2)/(2*k)  );
    res=A*B;
    return(res);
}
// compute  bivariate log-normal pdf:
double d2lognorm(double x, double y, double sill, double mux,double muy,double rho)
{
    double res=0.0, q=0.0, omr=1-pow(rho,2);
    q=(pow((log(x)-mux)/sqrt(sill),2)+
       pow((log(y)-muy)/sqrt(sill),2)
       -2*rho*(log(x)-mux)/sqrt(sill)*
       (log(y)-muy)/sqrt(sill))/omr;
    res=exp(-q/2)/(2*x*y*M_PI_F*sill*sqrt(omr));
    //Rprintf("%f\n",res);
    return(res);
}
//bivariate skew gaussian distribution
double biv_skew(double corr,double zi,double zj,double vari,double skew)
{
    double aux1=0.0,aux2=0.0,pdf1,pdf2,quadr;
    double dens,det,fact1,fact2, fact3,c1,lim1,lim2,cdf1,cdf2,a11,a12;
    double nu2=pow(skew,2.0);
    double tau2 =vari;
    if(-LOW<skew<LOW)
    { det=pow(tau2,2)-pow(tau2*corr,2);
        quadr=(1.0/det)*( aux1*(pow(zj,2.0)+pow(zi,2.0)) - 2.0*corr*zi*zj  ) ;
        dens=(0.5/M_PI_F) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;}
    else{
        
        c1= 1.0 /  (1-pow(corr,2)) ;
        aux1 = tau2 + nu2 ;
        aux2 =  tau2*corr + nu2*corr ;
        det =   pow(aux1,2) - pow(aux2,2) ;
        quadr= (1.0/det)*( aux1*(pow(zj,2.0)+pow(zi,2.0)) - 2.0*aux2*zi*zj  ) ;
        pdf1=  (0.5/M_PI_F) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
        //////
        aux1 = (nu2/tau2)*c1 + c1 ;
        aux2 = (nu2/tau2)*c1*corr + c1*corr;
        det=  pow(aux1,2) - pow(aux2,2) ;
        a11  =  (1/det)* aux1   ;
        a12  =  (1/det)* aux2   ;
        fact3 =( c1*skew)/ ( tau2 * det ) ;
        fact1 = fact3 * (aux1 - corr*aux2) ;
        fact2 = fact3 * (aux2 - corr*aux1) ;
        lim1 =   fact1*zi +  fact2*zj   ;
        lim2 =   fact2*zi +  fact1*zj   ;
        cdf1 = cdf_norm_OCL(lim1,lim2,a11,a12)  ;//cdf1 = cdf_norm_OCL(lim1,lim2,a11,a12) ;
        //printf("lim1: %f\tlim2: %f\ta11: %f\ta12: %f\tcdf1: %f\n",lim1,lim2,a11,a12,cdf1);
        //////////////////////////////////////////////////////////////////////////////////
        aux1 = tau2 + nu2 ;
        aux2 =  tau2*corr - nu2*corr ;
        det =   pow(aux1,2) - pow(aux2,2) ;
        quadr= (1.0/det)*( aux1*(pow(zj,2.0)+pow(zi,2.0)) - 2.0*aux2*zi*zj  ) ;
        pdf2=  (0.5/M_PI_F) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
        //////
        aux1 = (nu2/tau2)*c1 + c1 ;
        aux2 = (nu2/tau2)*c1*corr - c1*corr;
        det=  pow(aux1,2) - pow(aux2,2) ;
        a11  =  (1/det)* aux1   ;
        a12  =  (1/det)* aux2   ;
        fact3 =( c1*skew)/ ( tau2 * det  ) ;
        fact1 = fact3 * (aux1 - corr*aux2) ;
        fact2 = fact3 * (aux2 - corr*aux1) ;
        lim1 =   fact1*zi +  fact2*zj   ;
        lim2 =   fact2*zi +  fact1*zj   ;
        cdf2 = cdf_norm_OCL(lim1,lim2,a11,a12) ; //cdf2 = cdf_norm_OCL(lim1,lim2,a11,a12)
        //printf("lim1: %f\tlim2: %f\ta11: %f\ta12: %f\tcdf2: %f\n",lim1,lim2,a11,a12,cdf2);
        
        dens=2*(pdf1*cdf1+pdf2*cdf2);
    }
    
    //printf("dens: %f\n",dens);
    return(dens);
}

double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double a=0.0,A=0.0,z=0.0,res=0.0,B=0.0,C=0.0;
    double ci=exp(mui);
    double cj=exp(muj);
    if(corr)   {
        a=1-pow(corr,2);
        z=shape*corr*sqrt(zi*zj)/(sqrt(ci*cj)*a);
        
        A=pow(zi*zj,shape/2-1) * exp(-shape*((zi/ci)+(zj/cj))/(2*a));
        C=pow(z/2,1-shape/2);
        B=tgamma(shape/2)*pow(a,shape/2)*pow(2,shape)*pow(shape,-shape)*pow(ci*cj,shape/2);
        res=A*C*bessel_ii(z,shape/2-1,1)/B;//A*C*bessel_i(z,shape/2-1,1)/B
    }
    else
    {
        B=(pow((shape/(2*ci)),shape/2)*pow(zi,shape/2-1)*exp(-(shape*zi/(2*ci))) )/tgamma(shape/2);
        C=(pow((shape/(2*cj)),shape/2)*pow(zj,shape/2-1)*exp(-(shape*zj/(2*cj))) )/tgamma(shape/2);
        res=B*C;
    }
    return(res);
    
}

double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0, kk=0.0;
    a_u=min0(u,v);//min0(u,v)
    dens=0;
    pp=p11*p01*p10;
    bb=(1-p11)/p11;
    for(a=0;a<=a_u;a++){
        double aux1 = a+1;
        double aux2 = u-a+1;
        double aux3 = v-a+1;
        kk=exp(-(lgamma(aux1)+lgamma(aux2)+lgamma(aux3)));//exp(-(lgamma(a+1)+lgamma(u-a+1)+lgamma(v-a+1)))
        dens+=kk*(pow(NN*((p11*(p01+p10)-p01*p10*(p11+1))/pp),a)*pow(NN*((p10-p11)/(p10*p11)),u-a)*pow(NN*bb,v-a));
    }
    return(exp(-NN*bb)*dens);
}

// compute the bivariate normal cdf for the bernoulli RF:

double pbnorm(int cormod, double h, double u, double mean1, double mean2, double nugget, double var,double par0,double par1,double par2,double par3, double thr)
{
    double res=0;
    //double lim_inf[2]={0,0};//lower bound for the integration
    double lim_sup[2]={mean1,mean2};
    //int infin[2]={0,0};//set the bounds for the integration
    double corr[1]={var*CorFct(cormod,h,u,par0,par1,par2,par3,0,0)};
    //res=F77_CALL(bvnmvn)(lim_inf,lim_sup,infin,corr);
    res = Phi2(lim_sup[0],lim_sup[1],corr[0]);
    return(res);
}

// bivariate pois-binomial
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0, kk=0.0;
    a_u=min0(u,v);
    dens=0;
    for(a=0;a<=a_u;a++){
        double aux1 = a+1;
        double aux2 = u-a+1;
        double aux3 = v-a+1;
        kk=exp(-(lgamma(aux1)+lgamma(aux2)+lgamma(aux3)));//exp(-(lgamma(aux1)+lgamma(aux2)+lgamma(aux3)))
        dens+=kk*(pow(NN*p11,a)*pow(NN*(p01-p11),u-a)*pow(NN*(p10-p11),v-a));
    }
    return(exp(-NN*(p01+p10-p11))*dens);
}

double biv_binom(int NN, int u, int v, double p01,double p10,double p11)
{
    
    int a;
    double kk=0.0,dens=0.0;
    for(a=max0(0,u+v-NN);a<=min0(u,v);a++)
    {
        double aux1 = a+1;
        double aux2 = u-a+1;
        double aux3 = v-a+1;
        double aux4 = NN+1;
        double aux5 = NN-u-v+a+1;
        //kk=exp(logfac(NN)-(logfac(a)+logfac(u-a)+logfac(v-a)+logfac(NN-u-v+a)));
        kk=exp(lgamma(aux4)-(lgamma(aux1)+lgamma(aux2)+lgamma(aux3)+lgamma(aux5)));
        dens+=kk*(pow(p11,a)*pow(p01-p11,u-a)*pow(p10-p11,v-a)*pow(1+p11-(p01+p10),NN-u-v+a));
    }
    return(dens);
}


double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11)
{
    int a=0,i=0;
    double kk1=0.0,kk2=0.0,dens1=0.0,dens2=0.0;
    
    for(a=max0(0,u-v+NN-1);a<=NN-2;a++){
        for(i=max0(0,a-u);i<=min0(a,NN-1);i++){
            double aux1 = NN+u;
            double aux2 = i+1;
            double aux3 = NN-i;
            double aux4 = a-i+1;
            double aux5 = u-a+i+1;
            double aux6 = v-u;
            double aux7 = v+a-NN-u+2;
            double aux8 = NN-a-1;
            kk1=exp(lgamma(aux1)-(lgamma(aux2)+lgamma(aux3)+lgamma(aux4)+lgamma(aux5)));
            kk2=exp(lgamma(aux6)-(lgamma(aux7)+lgamma(aux8)));
            dens1+=kk1*kk2*pow(p11,i+1)*pow(1+p11-(x+y),u-a+i)*
            pow(x-p11,NN-i-1)*pow(y-p11,a-i)*pow(1-y,v-u-NN+a+1)*pow(y,NN-a-1);
        }}
    
    for(a=max0(0,u-v+NN);a<=NN-1;a++){
        for(i=max0(0,a-u);i<=min0(a,NN-1);i++){
            double aux1 = NN+u;
            double aux2 = i+1;
            double aux3 = NN-i;
            double aux4 = a-i+1;
            double aux5 = u-a+i+1;
            double aux6 = v-u;
            double aux7 = v+a-NN-u+1;
            double aux8 = NN-a;
            kk1=exp(lgamma(aux1)-(lgamma(aux2)+lgamma(aux3)+lgamma(aux4)+lgamma(aux5)));
            kk2=exp(lgamma(aux6)-(lgamma(aux7)+lgamma(aux8)));
            dens2+=kk1*kk2*pow(p11,i)*pow(1+p11-(x+y),u-a+i)*
            pow(x-p11,NN-i)*pow(y-p11,a-i)*pow(1-y,v-u-NN+a)*pow(y,NN-a);
        }}
    return(dens1+dens2);
}

// bivariate negative  binomial
double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
    double kk1=0.0,dens=0.0;int i=0;
    if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);
    
    if(u==v)            {
        for(i=max0(0,NN-u-1);i<=NN-1;i++){
            double aux1 = NN+u;
            double aux2 = i+1;
            double aux3 = NN-i;
            double aux4 = u-NN+2+i;
            kk1=exp(lgamma(aux1)-(lgamma(aux2)+lgamma(aux3)+lgamma(aux3)+lgamma(aux4)));
            dens+=kk1*pow(p11,i+1)*pow(1+p11-(p01+p10),u-NN+1+i)*pow(p01-p11,NN-1-i)*pow(p10-p11,NN-1-i); }
    }
    
    if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
    return(dens);
}

double  biv_binom2 (int NN_i,int NN_j, int k, int u, int v, double p01,double p10,double p11)
{
    
    int a,i,j;
    double const1=0.0,const2=0.0,const3=0.0,dens=0.0,dens1;
    double P10=p01-p11;
    double P01=p10-p11;
    double P00=1+p11-(p01+p10);
    double P11=p11;
    
    for(i=0;i<=min0(NN_i-k,u);i++){
        for(j=0;j<=min0(NN_j-k,v);j++){
            
            for(a=max0(0,u+v-k-i-j);a<=min0(u-i,v-j);a++){
                double aux1 = k+1;
                double aux2 = a+1;
                double aux3 = u-i-a+1;
                double aux4 = v-j-a+1;
                double aux5 = k-u-v+i+j+a+1;
                const1=exp(lgamma(aux1)-(lgamma(aux2)+lgamma(aux3)+lgamma(aux4)+lgamma(aux5)));
                dens1=const1*pow(P11,a)*pow(P00,k-u-v+i+j+a)*
                pow(P10,u-i-a)*pow(P01,v-j-a);
                double aux6 =NN_i-k+1;
                double aux7 =NN_i-k-i+1;
                double aux8 = i+1;
                double aux9 =NN_j-k+1;
                double aux10 = NN_j-k-j+1;
                double aux11 = j+1;
                const2=exp(lgamma(aux6)-(lgamma(aux7)+lgamma(aux8)));
                const3=exp(lgamma(aux9)-(lgamma(aux10)+lgamma(aux11)));
                
                //Rprintf("%f %f %f\n",const2,const3,dens1);
                dens+=dens1*const2*const3*
                pow(P11+P10,i) * pow(P11+P01,j) *
                pow(P00+P01,NN_i-k-i) * pow(P00+P10,NN_j-k-j);
            }}}
    return(dens);
}

// ===================================== END Distributions.c  ==================================//



__kernel void Comp_Cond_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par) // main test kernel
{
    
    int j, gid = get_global_id(0);
    double s1=0.0, s12=0.0, lags=0.0,weights=1.0, sum=0.0;
    double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    //double mean = dou_par[2];
    
    //double2 par = (double2)(dou_par[0],dou_par[1]);
    
    int ncoord  = int_par[1];
    int cormod  = int_par[0];
    int type    = int_par[3];
    
    
    
    //if(nuis1<0 || nuis0<0) {res[gid]=LOW;}
    s1=nuis0+nuis1;//set nugget + sill
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   ) //((gid+j)!= j) && ((gid+j) < n)
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            //printf("%f\t%d\t%d\t%f\t%f\t\n",lags,gid+j,j,data[gid+j],data[j]);
            if(lags<=maxdist){
                //s12=par0*exp(-lags/par1); // par = varianza, parametro de escala y media
                s12=nuis1*CorFct(cormod, lags, 0, par0,par1,par2,par3,0,0); //sill * corr
                det=pow(s1,2)-pow(s12,2);
                //printf("GPU\t%f\t%f\t%f\t%f\t%f\t%d\t%f\n",s1,s12,par0,par1,nuis1,cormod);
                //printf("GPU\t%f\t%d\t%f\t%f\t%f\t%f\n",nuis1,cormod,lags,par0,par1,s12);
                //printf("CPU:%f\t%f\t%f\t%f\n",par0,par1,par2,par3);
                u=data[gid+j]-mean[gid+j]; //data[si] - mean // nuis es la media
                v=data[j]-mean[j]; //data[sj] - mean
                
                if(!isnan(u)&&!isnan(v) )
                {
                    u2=pow(u,2);
                    v2=pow(v,2);
                    sum+= (-log(2*M_PI_F)-log(det)+log(s1)+
                            (u2+v2)*(0.5/s1-s1/det)+2*s12*u*v/det)*weights;
                }
                
            }
            //printf("IF= GLOBAL %d LOCAL %d SUM %f\n",gid,local_id,sum);
        }
        
        else
            continue;
    }
    //barrier(CLK_GLOBAL_MEM_FENCE);
    res[gid] = sum;

}


__kernel void Comp_Pair_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par) // main test kernel
{
    
    int j, gid = get_global_id(0);
    //res[gid] = coordx[gid]+coordy[gid];
    double s1=0.0, s12=0.0, lags=0.0,weights=1.0, sum=0.0;
    double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    //double2 par = (double2)(dou_par[0],dou_par[1]);
    int cormod  = int_par[0];
    int ncoord  = int_par[1];
    int type    = int_par[3];
    
    //if(nuis1<0 || nuis0<0) {res[gid]=LOW;}
    s1=nuis0+nuis1;//set nugget + sill
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   ) //((gid+j)!= j) && ((gid+j) < n)
        {
            //double dx = (coordx[j]-coordx[gid+j]);
            //double dy = (coordy[j]-coordy[gid+j]);
            //lags = sqrt( dx * dx + dy * dy);
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            //printf("lags:%f\n",lags);
            if(lags<=maxdist){
                //s12=par0*exp(-lags/par1); // par = varianza, parametro de escala y media
                s12=nuis1*CorFct(cormod, lags, 0, par0,par1,par2,par3,0,0); //sill * corr
                det=pow(s1,2)-pow(s12,2);
                
                u=data[gid+j]-mean[gid+j]; //data[si] - mean // nuis es la media
                v=data[j]-mean[j]; //data[sj] - mean
                
                if(!isnan(u)&&!isnan(v) )
                {
                    u2=pow(u,2);
                    v2=pow(v,2);
                    sum+= -0.5*(2*log(2*M_PI_F)+log(det)+
                                (s1*(u2+v2)-2*s12*u*v)/det)*weights;
                }
                
            }
            //printf("IF= GLOBAL %d LOCAL %d SUM %f\n",gid,local_id,sum);
            //printf("SUM %f\n",sum);
        }
        
        else
           continue;
    }
    //barrier(CLK_GLOBAL_MEM_FENCE);
    res[gid] = sum;
    
}

// Composite marginal (difference) log-likelihood for the spatial Gaussian model:
__kernel void Comp_Diff_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double lags=0.0,weights=1.0, sum=0.0,vario=0.0;
    double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    

    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                vario=Variogram(cormod,lags,0,nuis0,nuis1,par0,par1,par2,par3);
                u=data[gid+j];
                v=data[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    if(weigthed) weights=CorFunBohman(lags,maxdist);
                    sum+=  -0.5*(log(2*M_PI_F)+log(vario)+
                                 pow(u-v,2)/(2*vario))*weights;
                }
                
            }
        }
        
        else
            continue;
    }
    //barrier(CLK_GLOBAL_MEM_FENCE);
    res[gid] = sum;
    
}



__kernel void Comp_Pair_WrapGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double s1=0.0, s12=0.0,sum=0.0;
    double det=0.0, u=0.0,v=0.0,weights=1.0;
    double quadr=0.0,wrap_gauss;
    double lags=0.0,alfa=2.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    s1=nuis0+nuis1;//set nugget + sill
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                s12=nuis1*CorFct(cormod,lags,0,par0,par1,par2,par3,0,0); //sill * corr
                det=pow(s1,2)-pow(s12,2);
                
                u=data[gid+j]-2*atan(mean[gid+j])-M_PI_F;; //data[si] - mean
                v=data[j]-2*atan(mean[j])-M_PI_F;; //data[sj] - mean
                //printf("u,v: %f\t%f\n",u,v);
                if(!isnan(u)&&!isnan(v) )
                {
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    double k1=-alfa,k2=-alfa;
                    wrap_gauss = 5;
                    while(k1<=alfa){
                        while(k2<=alfa){
                            
                            quadr = -0.5*(1.0/det)*(s1*pow((v+2*k1*M_PI_F),2.0)+s1*pow((u+2*k2*M_PI_F),2.0)
                                                    -2.0*s12*(u+2*k2*M_PI_F)*(v+2*k1*M_PI_F));
                            wrap_gauss +=  (1/2.0*M_PI_F)*(1/sqrt(det)*exp(quadr)) ;
                            k2 = k2+1;}
                        k1 = k1+1;k2 = -alfa;}
                    sum+=  log(wrap_gauss)*weights ;
                    //printf("sum: %f\n",sum);
                }
                
            }
        }
        
        else
            continue;
    }
    //barrier(CLK_GLOBAL_MEM_FENCE);
    res[gid] = sum;
    
}


__kernel void Comp_Pair_SinhGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,bb=0.0,weights=1.0,sum=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
   
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                zi=(data[gid+j]-mean[gid+j])/sqrt(nuis1);
                zj=(data[j]-mean[j])/sqrt(nuis1);
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    bb=log(biv_sinh(corr,zi,zj,nuis2,nuis3,nuis1));
                    sum+=  weights*bb;
                    //printf("sum: %f\n",sum);
                }
                
            }
        }
        
        else
            continue;
    }
    res[gid] = sum;
    
}


__kernel void Comp_Pair_LogGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double corr,zi,zj,lags,bb=0.0,weights=1.0,sum=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                zi=data[gid+j];
                zj=data[j];
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    bb=log(d2lognorm(zi,zj,nuis1, mean[gid+j], mean[j],corr));
                    sum+=  weights*bb;
                    //printf("sum: %f\n",sum);
                }
                
            }
        }
        
        else
            continue;
    }
    res[gid] = sum;
    
}



__kernel void Comp_Pair_SkewGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double corr,zi,zj,lags,bb=0.0,weights=1.0,sum=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                zi=data[gid+j]-mean[gid+j];
                zj=data[j]-mean[j];
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    sum+=  weights*log(biv_skew((1-nuis0)*corr,zi,zj,nuis1,nuis2));
                }
                
            }
        }
        
        else
            continue;
    }
    /*if((isfinite(sum)))
    {
        printf("SI ES FINITO: %f\n",sum);
    }*/
    res[gid] = sum;
    
}


__kernel void Comp_Pair_Gamma2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double corr,zi,zj,lags,bb=0.0,weights=1.0,sum=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                zi=data[gid+j];
                zj=data[j];
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    sum+=  weights*log(biv_gamma(corr,zi,zj,mean[gid+j],mean[j],nuis2)); //
                }
                
            }
        }
        
        else
            continue;
    }
 
    res[gid] = sum;
    
}

__kernel void Comp_Pair_PoisbinnegGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par) // PILAS!!!!! FALTA CONFIRMAR!!!
{
    
    int j, gid = get_global_id(0),uu=0,vv=0;
    double u,v,dens=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0, sum=0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success

    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[5];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                ai=mean[gid+j];
                aj=mean[j];
                
                psj=pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);//pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);
                p1=pnorm_OCL(ai,0,1);//pnorm_OCL(ai,0,1)
                p2=pnorm_OCL(aj,0,1);//pnorm_OCL(aj,0,1)
                
                u=data[gid+j];
                v=data[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    //corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    uu=(int) u;
                    vv=(int) v;
                    dens=biv_poisbinneg(NN,uu,vv,p1,p2,psj);//biv_poisbinneg(NN,uu,vv,p1,p2,psj)
                    sum+=  log(dens)*weights; //
                    //printf("sum: %f\n",sum);
                    //if(!(isfinite(sum))) sum = LOW;
                }
                
            }
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}

__kernel void Comp_Pair_PoisbinGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par) // PILAS!!!!! FALTA CONFIRMAR!!!
{
    
    int j, gid = get_global_id(0),uu=0,vv=0;
    double u,v,dens=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0, sum=0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[5];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                ai=mean[gid+j];
                aj=mean[j];
                
                psj=pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);//pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);
                p1=pnorm_OCL(ai,0,1);//pnorm_OCL(ai,0,1)
                p2=pnorm_OCL(aj,0,1);//pnorm_OCL(aj,0,1)
                
                u=data[gid+j];
                v=data[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    //corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    uu=(int) u;
                    vv=(int) v;
                    dens=biv_poisbin(NN,uu,vv,p1,p2,psj);//biv_poisbinneg(NN,uu,vv,p1,p2,psj)
                    sum+=  log(dens)*weights; //
                    //printf("sum: %f\n",sum);
                    //if(!(isfinite(sum))) sum = LOW;
                }
                
            }
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}


__kernel void Comp_Pair_BinomGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0),uu=0,vv=0;
    double u,v,dens=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0, sum=0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[5];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                ai=mean[gid+j];
                aj=mean[j];
                
                psj=pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);//pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);
                p1=pnorm_OCL(ai,0,1);//pnorm_OCL(ai,0,1)
                p2=pnorm_OCL(aj,0,1);//pnorm_OCL(aj,0,1)
                
                u=data[gid+j];
                v=data[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    uu=(int) u;
                    vv=(int) v;
                    dens=biv_binom(NN,uu,vv,p1,p2,psj);//biv_binom(NN,uu,vv,p1,p2,psj)
                    sum+=  log(dens)*weights; //
                    
                }
                
            }
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}

__kernel void Comp_Pair_BinomnegGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0),uu=0,vv=0;
    double u,v,dens=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0, sum=0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[5];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                ai=mean[gid+j];
                aj=mean[j];
                
                psj=pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);//pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);
                p1=pnorm_OCL(ai,0,1);//pnorm_OCL(ai,0,1)
                p2=pnorm_OCL(aj,0,1);//pnorm_OCL(aj,0,1)
                
                u=data[gid+j];
                v=data[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    uu=(int) u;
                    vv=(int) v;
                    dens=biv_binomneg(NN,uu,vv,p1,p2,psj);//
                    sum+=  log(dens)*weights; //
                    
                }
                
            }
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}

__kernel void Comp_Pair_Binom2Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0),uu=0,vv=0;
    double u,v,dens=0.0,lags=0.0,weights=1.0,ai=0.0,aj=0.0, sum=0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[5];
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                
                ai=mean[gid+j];
                aj=mean[j];
                
                psj=pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);//pbnorm(cormod,lags,0,ai,aj,nuis0,nuis1,par0,par1,par2,par3,0);
                p1=pnorm_OCL(ai,0,1);//pnorm_OCL(ai,0,1)
                p2=pnorm_OCL(aj,0,1);//pnorm_OCL(aj,0,1)
                
                u=data[gid+j];
                v=data[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    uu=(int) u;
                    vv=(int) v;
                    dens=biv_binomneg(NN,uu,vv,p1,p2,psj);//
                    sum+=  log(dens)*weights; //
                    
                }
                
            }
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}
/******************************************************************************************/
/********************* SPACE TIME CASE *****************************************************/
/******************************************************************************************/


__kernel void Comp_Pair_Gauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    double maxdist = dou_par[6];
    double maxtime	=	dou_par[11];
    double nuis0 = dou_par[4];
    double nuis1 = dou_par[5];
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    //int ls = get_global_size(0);
    //int ms = get_global_size(1);
    
    int m=0,v =0;
    double s1=0.0, s12=0.0, lags=0.0, lagt=0.0,weights=1.0,sum=0.0;
    double det=0.0, u=0.0, u2=0.0, w=0.0, w2=0.0;
    
    //int gid = (npts*t+l);
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    //int j = (ntime*gidx+gidy);
    
    
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    if(l >= ncoord) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        //res[i]=0;
        //mom_cond0[i] = 0;mom_cond1[i] = 0;mom_cond2[i] = 0;mom_cond3[i] = 0;
        s1=nuis0+nuis1;
        for(m = l;m<ncoord;m++)
        {
            if(l==m)
            {
                for(v = t+1;v<ntime;v++)
                {
                    lagt=fabs(coordt[t]-coordt[v]);
                    
                    if(lagt<=maxtime)
                    {
                        //printf("%d\t%d\n",(t+ntime*l),(v+ntime*l));
                        s12=nuis1*CorFct(cormod,0, lagt,par0,par1,par2,par3,t,v);
                        det=pow(s1,2)-pow(s12,2);
                        
                        u=data[(t+ntime*l)]-mean[(t+ntime*l)];
                        w=data[(v+ntime*l)]-mean[(v+ntime*l)];
                        if(!isnan(u)&&!isnan(w) ){
                            u2=pow(u,2);
                            w2=pow(w,2);
                            if(weigthed) {weights=CorFunBohman(lagt,maxtime);}
                            sum+= (-0.5*(2*log(2*M_PI_F)+log(det)+
                                          (s1*(u2+w2)-2*s12*u*w)/det))*weights;}
                        
                    }
                }
            }
            else
            {
                lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                for(v=0;v<ntime;v++)
                {
                    lagt=fabs(coordt[t]-coordt[v]);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        //printf("%d\t%d\n",(t+ntime*l),(v+ntime*m));
                        //printf("GPU: %f\t%f\n",lags,lagt);
                        s12=nuis1*CorFct(cormod,lags, lagt,par0,par1,par2,par3,t,v);
                        det=pow(s1,2)-pow(s12,2);
                        
                        u=data[(t+ntime*l)]-mean[(t+ntime*l)];
                        w=data[(v+ntime*m)]-mean[(v+ntime*m)];
                        if(!isnan(u)&&!isnan(w) ){
                            u2=pow(u,2);
                            w2=pow(w,2);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= (-0.5*(2*log(2*M_PI_F)+log(det)+
                                          (s1*(u2+w2)-2*s12*u*w)/det))*weights;}
                    }
                }
            }
        }
        res[i] = sum;
    }
    //printf("GPU: %f\t%d\n",sum,i);
    
}

