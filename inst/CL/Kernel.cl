double biv_sinh(double corr,double zi,double zj,double mi,double mj,double skew,double tail,double vari);
double d2lognorm(double x, double y, double sill,double nugget, double mux,double muy,double rho);
double biv_skew(double corr,double zi,double zj,double mi,double mj,double vari,double skew);
double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape);
double biv_binom(int NN, int u, int v, double p01,double p10,double p11);
double pbnorm(int cormod, double h, double u, double mean1, double mean2, double nugget, double var,double par0,double par1,double par2,double par3, double thr);
double hyp2f1( double a,double b,double c,double x);
double hyt2f1( double a,double b,double c,double x,double *loss );
double hys2f1( double a,double b,double c,double x,double *loss );
double hypergeo(double a,double b,double c,double x);
long double digammal(long double x);
double biv_T(double rho,double zi,double zj,double nuu);
double appellF4(double a,double b,double c,double d,double x,double y);
#define LOW -1.0e15
// ===================================== Bessel  =====================================//

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





#define DBL_MAX 1.7976931348623157e+308
#define DBL_EPSILON 2.2204460492503131e-16
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#define DBL_MIN 2.2250738585072014e-308
#define HSQRT 1.414213562373095048801688724209698078569671

// For biv_T
#define EPS1 1.0e-10
#define ETHRESH 1.0e-12
#define MACHEP   1.11022302462515654042E-16   /* 2**-53 */
#define MAXNUM   1.79769313486231570815E308    /* 2**1024*(1-MACHEP) */


double bessel_ii(double x, double alpha, double expo);
double bessel_kk(double x, double alpha, double expo);
static void II_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bi, int *ncalc);
static void KK_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bk, int *ncalc);
double bessel_ii_minus(double x, double alpha, double expo);




/*************************************
 An ANSI-C implementation of the digamma-function for real arguments based
 on the Chebyshev expansion proposed in appendix E of
 http://arXiv.org/abs/math.CA/0403344 . This is identical to the implementation
 by Jet Wimp, Math. Comp. vol 15 no 74 (1961) pp 174 (see Table 1).
 For other implementations see
 the GSL implementation for Psi(Digamma) in
 http://www.gnu.org/software/gsl/manual/html_node/Psi-_0028Digamma_0029-Function.html
 
 Richard J. Mathar, 2005-11-24
 **************************************/
#ifndef M_PIl
/** The constant Pi in high precision */
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
/** the natural logarithm of 2 in high precision */
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

/** The digamma function in long double precision.
 * @param x the real value of the argument
 * @return the value of the digamma (psi) function at that point
 * @author Richard J. Mathar
 * @since 2005-11-24
 */
long double digammal(long double x)
{
    /* force into the interval 1..3 */
    if( x < 0.0L )
        return digammal(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
    else if( x < 1.0L )
        return digammal(1.0L+x)-1.0L/x ;
    else if ( x == 1.0L)
        return -M_GAMMAl ;
    else if ( x == 2.0L)
        return 1.0L-M_GAMMAl ;
    else if ( x == 3.0L)
        return 1.5L-M_GAMMAl ;
    else if ( x > 3.0L)
    /* duplication formula */
        return 0.5L*(digammal(x/2.0L)+digammal((x+1.0L)/2.0L))+M_LN2l ;
    else
    {
        /* Just for your information, the following lines contain
         * the Maple source code to re-generate the table that is
         * eventually becoming the Kncoe[] array below
         * interface(prettyprint=0) :
         * Digits := 63 :
         * r := 0 :
         *
         * for l from 1 to 60 do
         * 	d := binomial(-1/2,l) :
         * 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
         * 	evalf(r) ;
         * 	print(%,evalf(1+Psi(1)-r)) ;
         *o d :
         *
         * for N from 1 to 28 do
         * 	r := 0 :
         * 	n := N-1 :
         *
         *	for l from iquo(n+3,2) to 70 do
         *		d := 0 :
         *		for s from 0 to n+1 do
         *		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
         *		od :
         *		if 2*l-n > 1 then
         *		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
         *		fi :
         *	od :
         *	print(evalf((-1)^n*2*r)) ;
         *od :
         *quit :
         */
        long double Kncoe[] = { .30459198558715155634315638246624251L,
            .72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
            .27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
            .17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
            .11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
            .83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
            .59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
            .42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
            .304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
            .21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
            .15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
            .11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
            .80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
            .58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
            .41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;
        
        long double Tn_1 = 1.0L ;	/* T_{n-1}(x), started at n=1 */
        long double Tn = x-2.0L ;	/* T_{n}(x) , started at n=1 */
        long double resul = Kncoe[0] + Kncoe[1]*Tn ;
        
        x -= 2.0L ;
        
        for(int n = 2 ; n < sizeof(Kncoe)/sizeof(long double) ;n++)
        {
            const long double Tn1 = 2.0L * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
            resul += Kncoe[n]*Tn1 ;
            Tn_1 = Tn ;
            Tn = Tn1 ;
        }
        return resul ;
    }
}



// https://github.com/wch/r-source/blob/trunk/src/nmath/bessel_i.c
double bessel_ii(double x, double alpha, double expo)
{
    int nb, ncalc, ize;
    double na;
    
    if (x < 0 ) {
        //printf("A: %f\t%f\n",x,alpha);
        x=NAN;
        return x;
    }
    if ( isinf(x)) {
        //printf("B: %f\t%f\n",x,alpha);
        x=INFINITY;
        return x;
    }
    ize = (int)expo;
    
    na = floor(alpha);
   
    double cond = fabs(alpha-na);
    double B;
    /*if(alpha < 0)
    {
        //printf("C: %f\t%f\n",x,alpha);
        if(cond<DBL_EPSILON)//alpha==na cond<DBL_EPSILON
        {
            B =0;
        }
        else
        {
            B = bessel_kk(x, -alpha, expo) * 2.;
        }
        x = bessel_ii(x, -alpha, expo) + B;
        return(x);
    }*/
    
    if (alpha < 0) {
     
        return(bessel_ii(x, -alpha, expo) +
               ((alpha == na) ?  0 :
                bessel_kk(x, -alpha, expo) *
                ((ize == 1)? 2. : 2.*exp(-2.*x))/M_PI * sin(-alpha*M_PI)));
    }
    
    
    nb = 1 + (int)na;/* nb-1 <= alpha < nb */
    //printf("-----ALPHA:%f\n",alpha);
    alpha -= (double)(nb-1);
    
    //printf("NB:%d\n",nb);
    
    double bi[nb];
    //bi = (double *) calloc( nb, sizeof(double));
    //printf("BEFORE: x:%f\talpha:%f\tnb:%d\tize:%d\tncalc:%d\n\n",x,alpha,nb,ize,ncalc);
    II_bessel(&x, &alpha, &nb, &ize, bi, &ncalc);
    //printf("AFTER: x:%f\talpha:%f\tnb:%d\tize:%d\tncalc:%d\tbi:%f\n\n",x,alpha,nb,ize,ncalc,bi[nb-1]);
    /*if(ncalc != nb) {
        printf("E: %f\t%f\n",x,alpha);
        //x=NAN;
        //return x;
    }*/
    x = bi[nb-1];
    //bi[nb];
    return x;
}


double bessel_ii_minus(double x, double alpha, double expo)
{
    
    double na;
    
    if (x < 0) {
        x=NAN;
        return x;
    }
    na = floor(alpha);
    
    double cond = fabs(alpha-na);
    double B;
    
    if(cond<=DBL_EPSILON)
     {
     B =0;
     }
     else
     {
     B = bessel_kk(x, -alpha, expo) * 2.;
     }
    //printf("x0:%f alpha:%.3f expo:%.3f\n",x,-alpha,expo);
    //x = bessel_ii(x, -alpha, expo);
    //printf("x1:%.3f alpha:%.3f expo:%.3f\n",x,-alpha,expo);
    x = bessel_ii(x, -alpha, expo) + B;
    //x =  B;
    
    return(x);
}

double bessel_kk(double x, double alpha, double expo)
{
    int nb, ncalc, ize;
    
    
    if (x < 0) {
        //printf("ME_RANGE bessel_k");
        x=NAN;
        return x;
    }
    ize = (int)expo;
    if(alpha < 0) alpha = -alpha;
    nb = 1+ (int)floor(alpha);
    
    alpha -= (double)(nb-1);
    double bk[nb];
    //bk = (double *) calloc( nb, sizeof(double));
    //printf("NB:%d\n",nb);
    //bk[2];
    //bk = (double *) malloc( nb);
    
    KK_bessel(&x, &alpha, &nb, &ize, bk, &ncalc);
    /*if(ncalc != nb) {
        x=NAN;
        return x;
    }*/
    x = bk[nb-1];
    //x = alpha-expo;
    return x;
}


static void II_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bi, int *ncalc)
{

    
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
                bi[k]=INFINITY; /* the limit *is* = Inf */
            return;
            //printf("bessel_i(%g): ncalc (=%d) != nb (=%d); alpha=%g. PILAS!!!\n",
            //       *x, *ncalc, *nb, *alpha);
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
                        nend = min(*nb,n);
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
                test = max(test,sqrt(plast * ensig_BESS) * sqrt(p + p));
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
                            cc	= ldexp(cc, -900);;//ldexp(cc, -900);
                            bb	= ldexp(bb, -900);;//ldexp(bb, -900);
                            sum = ldexp(sum,-900);;//ldexp(sum,-900);
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
        *ncalc = min(*nb,0) - 1;
        
    }
}



static void KK_bessel(double *x, double *alpha, int *nb,
                      int *ize, double *bk, int *ncalc)
{
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
    *ncalc = min(*nb,0) - 2;
    if (*nb > 0 && (0. <= nu && nu < 1.) && (1 <= *ize && *ize <= 2)) {
        if(ex <= 0 || (*ize == 1 && ex > xmax_BESS_K)) {
            if(ex <= 0) {
                //if(ex < 0) printf("ME_RANGE, K_bessel\n");
                for(i=0; i < *nb; i++)
                    bk[i] = INFINITY;
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
        
        m = min((int) (wminf - nu),iend);
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
        *ncalc = max(1, mplus1 - k);
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

typedef void integr_fn(double *x, int n, void *ex);
typedef int Rboolean;
enum { FALSE, TRUE };
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
    if (*epsabs <= 0. && *epsrel < max(epmach * 50., 5e-29)) {
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
    errbnd = max(*epsabs, *epsrel * dres);
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
        errbnd = max(*epsabs, *epsrel * fabs(area));
        
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
        
        if (max(fabs(a1), fabs(b2)) <=
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
            ertest = max(*epsabs, *epsrel * fabs(reseps));
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
    if (ksgn == -1 && max(fabs(*result), fabs(area)) <= defabs * .01) {
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
        tol2 = max(fabs(e2), e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = max(e1abs, fabs(e0)) * epmach;
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
        tol1 = max(e1abs, fabs(e3)) * epmach;
        
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
    *abserr = max(*abserr, epmach * 5. * fabs(*result));
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
        *abserr = *resasc * min(1., pow(*abserr * 200. / *resasc, 1.5));
    }
    if (*resabs > uflow / (epmach * 50.)) {
        *abserr = max(epmach * 50. * *resabs, *abserr);
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
    return (res);
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
    //free(iwork);free(work);
    return(result);
}

// ===================================== END Integrate  =====================================//




// ============= START qnorm



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

// ============= END: qnorm


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
    double sol=NAN;
    if( a <= 0.0 ) sol = px;
    if( a >= 1.0 ) sol =  px * px;
    double b = 2.0 - a, sqrt_ab = sqrt( a * b );
    double c1 = 6.36619772367581343e-001;
    double c2 = 1.25331413731550025;
    double c3 = 1.57079632679489662;
    double c4 = 1.591549430918953358e-001;
    //double asr = ( a > 0.1 ? asin( 1.0 - a ) : acos( sqrt_ab ) );
    double asr;
    if(a > 0.1)
    {
        asr = asin( 1.0 - a );
    }else
    {
        asr = acos( sqrt_ab );
    }
    
    double comp = px * pxs;
    if( comp * ( 1.0 - a - c1 * asr ) < 5e-17 ) sol =  b * comp;
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
    /*while( res != res_new )
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
    }*/
    double cond = fabs(res-res_new);
    while( cond>DBL_EPSILON )
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
         cond = fabs(res-res_new);
     }
    double sol1;
    if(isnan(sol))
    {
        res *= exp( -x * x / b ) * c4;
        sol1 =  max( ( 1.0 + c1 * asr ) * comp, b * comp - max( 0.0, res ) );
    }else{
        sol1 = sol;
    }
    return sol1;
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
    
    bool c13 = (c1 && c3);
    bool c12 = (c1 && c2);
    double sol;
    if(c13) {sol = (q - 0.5);}
    else if(c12) {sol = (q);}
    else if(c1 ) {sol = (0.5 - p1 + q);}
    else if(c3 ) {sol = (p1 - q - 0.5);}
    else if(c2 ) {sol = (p1 - q);}
    else {sol = (0.5 - q);}
    
    //return (0.5 - p1 + q );
    return (sol );
    //return ( c1 && c3 ? q - 0.5
    //        : c1 && c2 ? q
    //        : c1 ? 0.5 - p1 + q
    //        : c3 ? p1 - q - 0.5
    //        : c2 ? p1 - q
    //        : 0.5 - q );
}

double Phi2( double x, double y, double rho )
{
    double sol = NAN;
    if( ( 1.0 - rho ) * ( 1.0 + rho ) <= 0.0 )
    {
        if( rho > 0.0 )
        {
            //return (Phi( min( x, y ) ));
            sol = (Phi( min( x, y ) ));
        }
        else
        {
            //return (max( 0.0, min( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ));
            sol = (max( 0.0, min( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ));
        }
    }
    if( x == 0.0 && y == 0.0 )
    {
        if( rho > 0.0 )
        {
            //return (Phi2diag( 0.0, 1.0 - rho, 0.5, 0.5 ));
            sol = (Phi2diag( 0.0, 1.0 - rho, 0.5, 0.5 ));
        }
        else
        {
            //return (0.5 - Phi2diag( 0.0, 1.0 + rho, 0.5, 0.5 ));
            sol = (0.5 - Phi2diag( 0.0, 1.0 + rho, 0.5, 0.5 ));
        }
    }
    else
    {
        sol = (max( 0.0,
                   min( 1.0,
                       Phi2help( x, y, rho ) + Phi2help( y, x, rho ) ) ));
    }
    
    return (sol);
}


double cdf_norm_OCL(double lim1,double lim2,double a11,double a12)
{
    double res=0;
    //double  uppe[2]={lim1/sqrt(a11),lim2/sqrt(a11)}, corre[1] ={a12/a11};
    //int infin[2]={0,0};
    double auxil=1-pow(a12/a11,2);
    double det=pow(a11,2)-pow(a12,2) ;
    //res= a11* sqrt(auxil/det) *  Phi2(uppe[0],uppe[1],corre[0]);
    res= a11* sqrt(auxil/det) *  Phi2(lim1/sqrt(a11),lim2/sqrt(a11),a12/a11);
    //res= a11* sqrt(auxil/det) *  Phi(a12/a11);
    return(res);
}


// ===================================== START Distance Functions  =====================================//

// Utility.c
double Dist_chordal(double loni, double lati, double lonj, double latj,double radius)
{
    double ai, bi, aj, bj, val=0.0;
    if (loni == lonj && lati == latj) return val;
    ai = (lati)*M_PI/180;
    bi = (loni)*M_PI/180;
    aj = (latj)*M_PI/180;
    bj = (lonj)*M_PI/180;
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
    ai = (lati)*M_PI/180;
    bi = (loni)*M_PI/180;
    aj = (latj)*M_PI/180;
    bj = (lonj)*M_PI/180;
    val = sin(ai) * sin(aj) + cos(ai) * cos(aj) * cos(bi - bj);
    if(val<= -1)  val2=M_PI*radius;
    if(val>=1) val2=0;
    val2 = acos(val)*radius;
    return(val2);
}

double dist(int type_dist,double coordx,double locx,double coordy,double locy,double radius)
{
    double lags=0.0;
    
    if(type_dist==0) lags=hypot(coordx-locx,coordy-locy);                        //euclidean
    if(type_dist==2) lags=Dist_geodesic(coordx,coordy,locx,locy,radius);           //great circle
    if(type_dist==1) lags=Dist_chordal(coordx,coordy,locx,locy,radius);      //chordal
    
    return(lags);
}

// ===================================== END Distance Functions  =====================================//

// ===================================== START CorrelationFunction.c  ==================================//
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
    if(x<=1)
    {rho=pow(1-x,smoo+1)*(1+(smoo+1)*x);}
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
        //free(param);
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
            // ========================   SPACE
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
        case 12:// Stable correlation function
            R_power=par0;
            scale=par1;
            rho=CorFunStable(h, R_power, scale);
            break;
        case 13://wen1
            R_power=par0;
            scale=par1;
            rho=CorFunW1(h,scale,R_power);
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
        case 18://  (sinpower) sinsphere correlation function valid on sphere
            R_power=par0;
            rho=1-pow(sin(h/(2)),R_power);
            break;
        case 19: // Generalised wend correlation function
            R_power1=par0;
            scale=par1;
            smooth=par2;
            rho=CorFunW_gen(h, R_power1, smooth, scale);
            break;
    }
    return rho;
}


double CorFct_st(int cormod, double h, double u, double par0,double par1,double par2,double par3,double par4,double par5,double par6, int c11, int c22)
{
    double arg=0.0, col=0.0,R_power=0.0, R_power1=0.0, R_power2=0.0, R_power_s=0.0, R_power_t=0.0, var11=0.0, var22=0.0;
    double rho=0.0, sep=0, scale=0.0, smooth=0.0,smooth_s=0.0,smooth_t=0.0, scale_s=0.0, scale_t=0, x=0, nug11=0.0, nug22=0.0;
    double scale11=0.0, scale22=0.0, scale12=0.0, smoo11=0.0, smoo22=0.0, smoo12=0.0,R_power11=0.0, R_power22=0.0, R_power12=0.0;
    switch(cormod) // Correlation functions are in alphabetical order
    {
            // ========================   SPACE TIME
        case 42:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            sep=1/(1+exp(-par4));par4;
            arg=1+pow(u/scale_t, R_power_t);
            rho=exp(-(pow(h/scale_s, R_power_s))*pow(arg, -0.5*sep*R_power_s))/arg;
            break;
        case 44:// Iaco-Cesare model as in (14) of Gneitint (2006): note that R_power parameters are in [0,1]
            R_power2=par0;
            R_power_s=par1;
            R_power_t=par2;
            scale_s=par3;
            scale_t=par4;
            rho=pow(1+pow(h/scale_s, R_power_s)+pow(u/scale_t, R_power_t),-R_power2);
            break;
        case 46://Porcu model as in (4) of Porcu Bevilacqua et al (2010), with beta_1=beta_2=1
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            sep=par4;
            if(sep>0) rho=pow(0.5*pow(1+pow(h/scale_s, R_power_s),sep)+0.5*pow(1+pow(u/scale_t, R_power_t),sep),-1/sep);
            else rho=pow((1+pow(h/scale_s, R_power_s))*(1+pow(u/scale_t,R_power_t)),-1);
            break;
        case 50:   //Gneiting correlation with prac ranges "giusto"
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            sep=par4;
            arg=1+pow(u, R_power_t)/scale_t;
            rho=exp(-(pow(h, R_power_s)/scale_s)*pow(arg, 0.5*sep*R_power_s))/pow(arg,1.5);
            break;
        case 52:// Gneiting correlation model valid on the sphere (equation 8)
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            sep=par4;
            arg=1+pow(h/scale_s, 2*R_power_s);
            rho=exp(-pow(u/scale_t, 2*R_power_t)/(pow(arg, sep*R_power_t)))/arg;
            break;
        case 54:// Gneiting correlation model valid on the sphere  (equazione 9)
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            sep=par4;
            // arg=1+R_pow(h/scale_s, 2*R_power_s);
            // rho=R_pow(1+R_pow(u/scale_t, 2*R_power_t)/R_pow(arg,R_power_t*sep),-1)/R_pow(arg,1);
            arg=1+pow(u/scale_t, 2*R_power_t);
            rho=exp(-pow(h/scale_s, R_power_s)*pow(arg,R_power_s*sep))/pow(arg,1);
            break;
        case 58:  //st multiquaderic
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            arg=pow(1+pow(u/scale_t,R_power_t),-1);
            rho= pow(pow(1-R_power_s/2,2)/(1+pow(R_power_s/2,2)-R_power_s*arg*cos(h)),scale_s);   // model B2 in the paper  (eq 9 right part)
            break;
        case 61:  //no sep gneiting  with temporal matern margin
            R_power_s=par0;
            R_power=par1;
            scale_s=par2;
            scale_t=par3;
            sep=par4;
            smooth=par5;
            arg=1+pow(h/scale_s, R_power_s);
            if(u>0) rho=pow(arg,-R_power)*CorFunWitMat(u,scale_t*pow(arg,sep/2),smooth);
            else  rho=pow(arg,-R_power);
            
            break;
        case 62:  //no sep gneiting  with spatial matern margin
            
            R_power_t=par0;
            R_power=par1;
            scale_s=par2;
            scale_t=par3;
            sep=par4;
            smooth=par5;
            arg=1+pow(u/scale_t, R_power_t);
            if(h>0)  rho=pow(arg,-R_power)*CorFunWitMat(h,scale_s*pow(arg,sep/2),smooth);
            else  rho=pow(arg,-R_power);
            
            break;
        case 63:  //
            R_power_t=par0;
            R_power_s=par1;
            R_power=par2;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            arg=pow(1+pow(u/scale_t,R_power_t/2),-1/(R_power_t/2));
            rho=pow(arg,R_power)*CorFunW0(h,scale_s*pow(arg,sep),R_power_s);
            break;
        case 64:
            R_power_s=par0;
            R_power=par1;
            R_power_t=par2;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            arg=pow(1+pow(h/scale_s,R_power_s/2),-1/(R_power_s/2));
            rho=pow(arg,R_power)*CorFunW0(u,scale_t*pow(arg,sep),R_power_t);  //2.5+2*0
            break;
        case 65:
            R_power_t=par0;
            R_power_s=par1;
            R_power=par2;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            arg=pow(1+pow(u/scale_t,R_power_t/2),-1/(R_power_t/2));
            rho=pow(arg,R_power)*CorFunW1(h,scale_s*pow(arg,sep),R_power_s);
            break;
        case 66:
            R_power_s=par0;
            R_power=par1;
            R_power_t=par2;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            arg=pow(1+pow(h/scale_s,R_power_s/2),-1/(R_power_s/2));
            rho=pow(arg,R_power)*CorFunW1(u,scale_t*pow(arg,sep),R_power_t); //2.5+2*1
            break;
        case 67:  //
            R_power_t=par0;
            R_power_s=par1;
            R_power=par2;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            arg=pow(1+pow(u/scale_t,R_power_t/2),-1/(R_power_t/2));
            rho=pow(arg,R_power)*CorFunW2(h,scale_s*pow(arg,sep),R_power_s);
            break;
        case 68:
            R_power_s=par0;
            R_power_t=par2;
            R_power=par1;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            arg=pow(1+pow(h/scale_s,R_power_s/2),-1/(R_power_s/2));
            rho=pow(arg,R_power)*CorFunW2(u,scale_t*pow(arg,sep),R_power_t); ////2.5+2*2
            break;
        case 88:
            R_power_s=par0;
            R_power=par1;
            R_power_t=par2;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            smooth=par6;
            arg=pow(1+pow(h/scale_s,R_power_s/2),-1/(R_power_s/2));
            rho=pow(arg,R_power)*CorFunW_gen(u,R_power_t,smooth,scale_t*pow(arg,sep));
            break;
        case 87:
            R_power_t=par0;
            R_power_s=par1;
            R_power=par2;
            scale_s=par3;
            scale_t=par4;
            sep=par5;
            smooth=par6;
            arg=pow(1+pow(u/scale_t,R_power_t/2),-1/(R_power_t/2));
            rho=pow(arg,R_power)*CorFunW_gen(h,R_power_s,smooth,scale_s*pow(arg,sep));
            break;
        case 69:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW0(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
            break;
        case 70:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW0(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
            break;
        case 71:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW0(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
            break;
        case 72:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW1(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
            break;
        case 73:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW1(h,scale_s,R_power_s)*CorFunW1(u,scale_t,R_power_t);
            break;
        case 74:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW1(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
            break;
        case 75:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW2(h,scale_s,R_power_s)*CorFunW0(u,scale_t,R_power_t);
            break;
        case 77:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunW2(h,scale_s,R_power_s)*CorFunW2(u,scale_t,R_power_t);
            break;
            // END non-separable correlation functions
            // START separable correlation functions:
        case 84:// Double exp:
            scale_s=par0;
            scale_t=par1;
            rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
            break;
        case 94:// Stable-stab:
            R_power_s=par0;
            R_power_t=par1;
            scale_s=par2;
            scale_t=par3;
            rho=CorFunStable(h,R_power_s,scale_s)*CorFunStable(u,R_power_t,scale_t);
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
double Variogram_st(int cormod, double h, double u, double nugget, double var, double par0,double par1,double par2,double par3,double par4,double par5,double par6)
{
    double vario=0.0;
    //Computes the variogram
    vario=nugget+var*(1-CorFct_st(cormod,h,u,par0,par1,par2,par3,par4,par5,par6,0,0));
    return vario;
}


double CorFunBohman(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)
    {
        if (x>0) rho=(1-x)*(sin(2*M_PI*x)/(2*M_PI*x))+(1-cos(2*M_PI*x))/(2*M_PI*M_PI*x);
        else   rho=1;
    }
    else rho=0;
    return rho;
}


// ===================================== END CorrelationFunction.c  ==================================//




/**************** for bivaraite T distribution */////////////////////////
double hyp2f1( double a,double b,double c,double x)
{
    double d, d1, d2, e;
    double p, q, r, s, y, ax;
    double ia, ib, ic, id, err;
    int flag, i, aid;
    
    err = 0.0;
    ax = fabs(x);
    s = 1.0 - x;
    flag = 0;
    ia = round(a); /* nearest integer to a */
    ib = round(b);
    
    if( a <= 0 )
    {
        if( fabs(a-ia) < EPS1 )    /* a is a negative integer */
            flag |= 1;
    }
    
    if( b <= 0 )
    {
        if( fabs(b-ib) < EPS1 )    /* b is a negative integer */
            flag |= 2;
    }
    
    if( ax < 1.0 )
    {
        if( fabs(b-c) < EPS1 )   /* b = c */
        {
            y = pow( s, -a ); /* s to the -a power */
            goto hypdon;
        }
        if( fabs(a-c) < EPS1 )   /* a = c */
        {
            y = pow( s, -b ); /* s to the -b power */
            goto hypdon;
        }
    }
    
    
    
    if( c <= 0.0 )
    {
        ic = round(c);  /* nearest integer to c */
        if( fabs(c-ic) < EPS1 )    /* c is a negative integer */
        {
            /* check if termination before explosion */
            if( (flag & 1) && (ia > ic) )
                goto hypok;
            if( (flag & 2) && (ib > ic) )
                goto hypok;
            goto hypdiv;
        }
    }
    
    if( flag )      /* function is a polynomial */
        goto hypok;
    
    if( ax > 1.0 )      /* series diverges  */
        goto hypdiv;
    
    p = c - a;
    ia = round(p); /* nearest integer to c-a */
    if( (ia <= 0.0) && (fabs(p-ia) < EPS1) ) /* negative int c - a */
        flag |= 4;
    
    r = c - b;
    ib = round(r); /* nearest integer to c-b */
    if( (ib <= 0.0) && (fabs(r-ib) < EPS1) ) /* negative int c - b */
        flag |= 8;
    
    d = c - a - b;
    id = round(d); /* nearest integer to d */
    q = fabs(d-id);
    
    /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
     * for reporting a bug here.  */
    if( fabs(ax-1.0) < EPS1 )      /* |x| == 1.0 */
    {
        if( x > 0.0 )
        {
            if( flag & 12 ) /* negative int c-a or c-b */
            {
                //Rprintf("c-a or c-b negative");
                if( d >= 0.0 )
                    goto hypf;
                else
                    goto hypdiv;
            }
            if( d <= 0.0 )
            {
                goto hypdiv;}
            //y = exp(lgammafn(c)+lgammafn(d) -(lgammafn(p) + lgammafn(r)));
            y = tgamma(c)*tgamma(d)/(tgamma(p)*tgamma(r));
            goto hypdon;
        }
        
        if( d <= -1.0 )
            goto hypdiv;
        
    }
    
    /* Conditionally make d > 0 by recurrence on c
     * AMS55 #15.2.27
     */
    if( d < 0.0 )
    {
        /* Try the power series first */
        y = hyt2f1( a, b, c, x, &err );
        if( err < ETHRESH )
            goto hypdon;
        /* Apply the recurrence if power series fails */
        err = 0.0;
        aid = 2 - id;
        e = c + aid;
        d2 = hyp2f1(a,b,e,x);
        d1 = hyp2f1(a,b,e+1.0,x);
        q = a + b + 1.0;
        for( i=0; i<aid; i++ )
        {
            r = e - 1.0;
            y = (e*(r-(2.0*e-q)*x)*d2 + (e-a)*(e-b)*x*d1)/(e*r*s);
            e = r;
            d1 = d2;
            d2 = y;
        }
        goto hypdon;
    }
    
    
    if( flag & 12 )
        goto hypf; /* negative integer c-a or c-b */
    
hypok:
    y = hyt2f1( a, b, c, x, &err );
    
    
hypdon:
    if( err > ETHRESH )
    {
        //Rprintf( "hyp2f1, PLOSS\n");
    }
    return(y);
    
    /* The transformation for c-a or c-b negative integer
     * AMS55 #15.3.3
     */
hypf:
    y = pow( s, d ) * hys2f1( c-a, c-b, c, x, &err );
    goto hypdon;
    
    /* The alarm exit */
hypdiv:
    //Rprintf( "hyp2f1, OVERFLOW\n");
    return( MAXNUM );
}


/* Apply transformations for |x| near 1
 * then call the power series
 */
double hyt2f1( double a,double b,double c,double x,double *loss )
{
    double p, q, r, s, t, y, d, err, err1;
    double ax, id, d1, d2, e, y1;
    int i, aid;
    
    err = 0.0;
    s = 1.0 - x;
    if( x < -0.5 )
    {
        if( b > a )
            y = pow( s, -a ) * hys2f1( a, c-b, c, -x/s, &err );
        
        else
            y = pow( s, -b ) * hys2f1( c-a, b, c, -x/s, &err );
        
        goto done;
    }
    
    d = c - a - b;
    id = round(d);  /* nearest integer to d */
    
    // Rprintf("%lf  %lf %lf %lf\n" , d, a, b, c);
    if( x > 0.9 )
    {
        if( fabs(d-id) > EPS1 ) /* test for integer c-a-b */
        {
            //  Rprintf("integer case\n");
            /* Try the power series first */
            y = hys2f1( a, b, c, x, &err );
            if( err < ETHRESH ) {
                //Rprintf("Power series failed");
                goto done;
            }
            /* If power series fails, then apply AMS55 #15.3.6 */
            q = hys2f1( a, b, 1.0-d, s, &err );
            if (d < 0)
                q *= tgamma(d) /(tgamma(c-a) * tgamma(c-b));
            else q *= exp(lgamma(d)  - (lgamma(c-a) + lgamma(c-b)));
            r = pow(s,d) * hys2f1( c-a, c-b, d+1.0, s, &err1 );
            if ( d > 0)
                r *= tgamma(-d) /tgamma(a) * tgamma(b);
            else  r *= exp(lgamma(-d) - (lgamma(a) + lgamma(b)));
            y = q + r;
            //  Rprintf("\n %lf %lf %lf \n", y, q, r);
            q = fabs(q); /* estimate cancellation error */
            r = fabs(r);
            if( q > r )
                r = q;
            err += err1 + (MACHEP*r)/y;
            
            y *= tgamma(c);
            goto done;
        }
        else
        {
            //  Rprintf("non-int case R2 %lf\n", x);
            //  Rprintf("a = %lf, b=%lf, c=%lf\n ", a, b, c);
            /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12 */
            if( id >= 0.0 )
            {
                //    Rprintf("id >= 0\n");
                e = d;
                d1 = d;
                d2 = 0.0;
                aid = id;
            }
            else
            {
                e = -d;
                d1 = 0.0;
                d2 = d;
                aid = -id;
            }
            
            ax = log(s);
            
            /* sum for t = 0 */
            y = digammal(1.0) + digammal(1.0+e) - digammal(a+d1) - digammal(b+d1) - ax;
            y /= tgamma(e+1.0);
            
            p = (a+d1) * (b+d1) * s / tgamma(e+2.0); /* Poch for t=1 */
            t = 1.0;
            do
            {
                r = digammal(1.0+t) + digammal(1.0+t+e) - digammal(a+t+d1)
                - digammal(b+t+d1) - ax;
                q = p * r;
                y += q;
                p *= s * (a+t+d1) / (t+1.0);
                p *= (b+t+d1) / (t+1.0+e);
                t += 1.0;
            }
            while( fabs(q/y) > EPS1 );
            
            
            if( id == 0.0 )
            {
                y *= tgamma(c)/(tgamma(a)*tgamma(b));
                goto psidon;
            }
            
            y1 = 1.0;
            
            if( aid == 1 )
                goto nosum;
            
            //  Rprintf("sum case b=%lf\n", b);
            t = 0.0;
            p = 1.0;
            for( i=1; i<aid; i++ )
            {
                r = 1.0-e+t;
                p *= s * (a+t+d2) * (b+t+d2) / r;
                t += 1.0;
                p /= t;
                y1 += p;
            }
        nosum:
            p = tgamma(c);
            //  Rprintf("e = %lf, a=%lf, b=%lf, d1= %lf, d2=%lf\n ", e, a, b, d1, d2);
            //  y1 *= gammafn(e) * p / (gammafn(a+d1) * gammafn(b+d1));
            y1 *= exp(lgamma(e) + log(p) - (lgamma(a+d1) - lgamma(b+d1)));
            //  y *= p / (gammafn(a+d2) * gammafn(b+d2));
            y *= exp(log(p) - lgamma(a+d2+ .00001))/tgamma(b+d2 + .00001);
            if( (aid & 1) != 0 )
                y = -y;
            
            q = pow( s, id ); /* s to the id power */
            if( id > 0.0 )
                y *= q;
            else
                y1 *= q;
            
            y += y1;
        psidon:
            goto done;
        }
        
    }
    
    /* Use defining power series if no special cases */
    y = hys2f1( a, b, c, x, &err );
    
done:
    *loss = err;
    return(y);
}


/* Defining power series expansion of Gauss hypergeometric function */

double hys2f1( double a,double b,double c,double x,double *loss )
//double *loss; /* estimates loss of significance */
{
    double f, g, h, k, m, s, u, umax;
    int i;
    
    i = 0;
    umax = 0.0;
    f = a;
    g = b;
    h = c;
    s = 1.0;
    u = 1.0;
    k = 0.0;
    do
    {
        if( fabs(h) < EPS1 )
        {
            *loss = 1.0;
            return( MAXNUM );
        }
        m = k + 1.0;
        u = u * ((f+k) * (g+k) * x / ((h+k) * m));
        s += u;
        k = fabs(u);  /* remember largest term summed */
        if( k > umax )
            umax = k;
        k = m;
        if( ++i > 10000 ) /* should never happen */
        {
            *loss = 1.0;
            return(s);
        }
    }
    while( fabs(u/s) > MACHEP );
    
    /* return estimated relative error */
    *loss = (MACHEP*umax)/fabs(s) + (MACHEP*i);
    
    return(s);
}





double hypergeo(double a,double b,double c,double x)
{
    double sol;
    sol =(hyp2f1( a, b, c, x ));
    return(sol);
}



// END aux fun for biv_T

// ===================================== START: Distributions.c  ==================================//
double biv_sinh(double corr,double zi,double zj,double mi,double mj,double skew,double tail,double vari)
{
    double b1=0.0,b2=0.0,A=0.0,B=0.0,k=0.0,res=0.0,Z1,Z2;
    double xi=(zi-mi)/sqrt(vari);
    double xj=(zj-mj)/sqrt(vari);
    b1=tail * asinh(xi)-skew;
    b2=tail * asinh(xj)-skew;
    k=1-pow(corr,2);
    A=pow(2 * M_PI * pow(k,0.5) * vari,-1) * cosh(b1) * cosh(b2) * pow(tail,2)/sqrt((pow(xi,2)+1) * (pow(xj,2)+1));
    Z1=sinh(b1);Z2=sinh(b2);
    B=exp(- (Z1*Z1 + Z2*Z2 - 2*corr*Z1*Z2)/(2*k)  );
    res=A*B;
    return(res);
}
// compute  bivariate log-normal pdf:

double d2lognorm(double x, double y, double sill,double nugget, double mux,double muy,double rho)
{
    double res=0.0, q=0.0, b2=sill+pow(nugget,1), omr=pow(b2,2)-pow(rho*sill,2);
    
    q=(b2*pow((log(x)-mux),2)+
       b2*pow((log(y)-muy),2)
       -2*rho*sill*(log(x)-mux)*(log(y)-muy))/omr;
    res=exp(-q/2)/(2*x*y*M_PI*sqrt(omr));
    return(res);
}
//bivariate skew gaussian distribution
double biv_skew(double corr,double zi,double zj,double mi,double mj,double vari,double skew)
{
    double aux1=0.0,aux2=0.0,pdf1,pdf2,quadr;
    zi=zi-mi;
    zj=zj-mj;
    double dens,det,fact1,fact2, fact3,c1,lim1,lim2,cdf1,cdf2,a11,a12;
    double nu2=pow(skew,2.0);
    double tau2 =vari;
    if(-LOW<skew<LOW)
    { det=pow(tau2,2)-pow(tau2*corr,2);
        quadr=(1.0/det)*( aux1*(pow(zj,2.0)+pow(zi,2.0)) - 2.0*corr*zi*zj  ) ;
        dens=(0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;}
    else{
        
        c1= 1.0 /  (1-pow(corr,2)) ;
        aux1 = tau2 + nu2 ;
        aux2 =  tau2*corr + nu2*corr ;
        det =   pow(aux1,2) - pow(aux2,2) ;
        quadr= (1.0/det)*( aux1*(pow(zj,2.0)+pow(zi,2.0)) - 2.0*aux2*zi*zj  ) ;
        pdf1=  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
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
        pdf2=  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
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
    double a=0.0,A=0.0,z=0.0,res=0.0,B=0.0,C=0.0, D=0.0;
    double ci=exp(mui);
    double cj=exp(muj);
    double gam = tgamma(shape/2);
    if(corr)   {
        a=1-pow(corr,2);
        z=shape*sqrt(zi*zj*corr*corr)/(sqrt(ci*cj)*a);
        A=pow(zi*zj,shape/2-1) * exp(-shape*((zi/ci)+(zj/cj))/(2*a));
        C=pow(z/2,1-shape/2);
        //B=tgamma(shape/2)*pow(a,shape/2)*pow(2,shape)*pow(shape,-shape)*pow(ci*cj,shape/2);
        B=gam*pow(a,shape/2)*pow(2,shape)*pow(shape,-shape)*pow(ci*cj,shape/2);
        /*if(isinf(z))
        {
            D =1.7976931348623157e+308;
            //printf("MIERDAAAAA\n");
        }else
        {
            D = bessel_ii(z,shape/2-1,1);
        }*/
        //D =isinf(z)?1.7976931348623157e+308:bessel_ii(z,shape/2-1,1);
        D = bessel_ii(z,shape/2-1,1);//bessel_ii(z,shape/2-1,1)
        
        //if(isinf(D)) D = 1.7976931348623157e+308;
        //if(isinf(C)) C = 1.7976931348623157e+308;
        //if(B<DBL_EPSILON) B = 4.9406564584124654E-324;
        res=A*C*D/B;
        //if((shape/2-1)<0) {printf("%f\t%f\tB:%f\t%f\t%f\tsh:%f\n\n",res,A,B,C,D,shape/2-1);}
        //printf("res:%f\tA:%f\tB:%f\tC:%f\tD:%f\tz:%f\tsha:%f\n\n",res,A,B,C,D,z,shape/2-1);
        //if(isnan(res)) res =4.9406564584124654E-324;
        //printf("D:%f\tz:%f\tshape:%f\tB:%f\tres:%f\n\n",D,z,shape/2-1,B,res);
        //printf("res:%f\tA:%f\tB:%f\tC:%f\tD:%f\tz:%f\tsha:%f\n\n",res,A,B,C,D,z,shape/2-1);
    }
    else
    {
        B=(pow((shape/(2*ci)),shape/2)*pow(zi,shape/2-1)*exp(-(shape*zi/(2*ci))) )/gam;
        C=(pow((shape/(2*cj)),shape/2)*pow(zj,shape/2-1)*exp(-(shape*zj/(2*cj))) )/gam;
        res=B*C;
        
    }
    
    //if(res < DBL_EPSILON || isnan(res)) res=4.9406564584124654E-324;
    //bool c1 = ( isnan(res) || res<DBL_EPSILON);
    //if(c1 ) res=4.9406564584124654E-324;
    //if(isinf(res)) res=DBL_MAX;//1.7976931348623157e+308;
    
    return(res);
    
}

double biv_binom(int NN, int u, int v, double p01,double p10,double p11)
{
    
    int a;
    double kk=0.0,dens=0.0;
    for(a=max(0,u+v-NN);a<=min(u,v);a++)
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

double biv_wrapped (double alfa,double u, double v,double mi,double mj,double nugget,double sill,double corr)
{
    double x,y,s1=0.0,s12=0.0,quadr=0.0,det=0.0,wrap_gauss=0.0; // 5???
    //2*atan(mean[i])-M_PI
    x=u-2*atan(mi)-M_PI; y=v-2*atan(mj)-M_PI;
    s1=nugget+sill;
    s12=sill*corr; //sill * corr
    det=pow(s1,2)-pow(s12,2);
    double k1=-alfa,k2=-alfa;
    while(k1<=alfa){
        while(k2<=alfa){
            quadr = -0.5*(1.0/det)*(s1*pow((y+2*k1*M_PI),2.0)+s1*pow((x+2*k2*M_PI),2.0)
                                    -2.0*s12*(x+2*k2*M_PI)*(y+2*k1*M_PI));
            wrap_gauss +=  (1/2.0*M_PI)*(1/sqrt(det)*exp(quadr)) ;
            k2 = k2+1;}
        k1 = k1+1;k2 = -alfa;}
    return(wrap_gauss);
}


// compute the bivariate normal cdf for the bernoulli RF:

double pbnorm(int cormod, double h, double u, double mean1, double mean2, double nugget, double var,double par0,double par1,double par2,double par3, double thr)
{
    double res=0;
    //double lim_inf[2]={0,0};//lower bound for the integration
    double lim_sup[2]={mean1,mean2};
    //int infin[2]={0,0};//set the bounds for the integration
    double corr[1]={(1-nugget)*CorFct(cormod,h,u,par0,par1,par2,par3,0,0)};
    //res=F77_CALL(bvnmvn)(lim_inf,lim_sup,infin,corr);
    res = Phi2(lim_sup[0],lim_sup[1],corr[0]);
    return(res);
}


double pbnorm_st(int cormod, double h, double u, double mean1, double mean2, double nugget, double var,double par0,double par1,double par2,double par3,double par4,double par5,double par6, double thr)
{
    double res=0;
    //double lim_inf[2]={0,0};//lower bound for the integration
    double lim_sup[2]={mean1,mean2};
    //int infin[2]={0,0};//set the bounds for the integration
    double corr[1]={var*CorFct_st(cormod,h,u,par0,par1,par2,par3,par4,par5,par6,0,0)};
    //res=F77_CALL(bvnmvn)(lim_inf,lim_sup,infin,corr);
    res = Phi2(lim_sup[0],lim_sup[1],corr[0]);
    return(res);
}


double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11)
{
    int a=0,i=0;
    double kk1=0.0,kk2=0.0,dens1=0.0,dens2=0.0;
    
    for(a=max(0,u-v+NN-1);a<=NN-2;a++){
        for(i=max(0,a-u);i<=min(a,NN-1);i++){
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
    
    for(a=max(0,u-v+NN);a<=NN-1;a++){
        for(i=max(0,a-u);i<=min(a,NN-1);i++){
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
        for(i=max(0,NN-u-1);i<=NN-1;i++){
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

// bivariate pois-binomial
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0, kk=0.0;
    a_u=min(u,v);
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

// bivariate pois-bineg
double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0, kk=0.0;
    a_u=min(u,v);
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

/**************/
double biv_Logistic(double corr,double zi,double zj,double mui, double muj, double beta)
{
    double a=0.0,A=0.0,D=0.0,res=0.0,B=0.0,C=0.0;double ci=exp(mui);double cj=exp(muj);
    if(corr)   {
        a=1-pow(corr,2);
        
        A=(exp((zi-ci)/beta)*exp((zj-cj)/beta))/(pow(a,-3)*pow(beta,2));
        B=(exp((zi-ci)/beta)+a)*(exp((zj-cj)/beta)+a);
        C=(pow(corr,2)*exp((zi-ci)/beta)*exp((zj-cj)/beta))/B;
        D=(C+1)/pow(1-C,3);
        res=A*D/pow(B,2);
    }
    else
    {
        B=exp((zi-ci)/beta)*pow((exp((zi-ci)/beta)+1),-2)/beta;
        C=exp((zi-ci)/beta)*pow((exp((zi-ci)/beta)+1),-2)/beta;
        res=B*C;
    }
    //printf("%f\n",res);
    return(res);
    
}

/**************/
double biv_LogLogistic(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double a=0.0,A=0.0,D=0.0,res=0.0,B=0.0,C=0.0;double ci=exp(mui);double cj=exp(muj);
    if(corr)   {
        a=1-pow(corr,2);
        
        A=pow(shape/sqrt(ci*cj),2)*pow((zi*zj)/(ci*cj),shape-1)/pow(a,-3);
        B=(pow((zi/ci),shape)+a)*(pow((zj/cj),shape)+a);
        C=(pow(corr,2)*pow(zi*zj,shape))/(pow(ci*cj,shape)*B);
        D=(C+1)/pow(1-C,3);
        res=A*D/(pow(B,2));
    }
    else
    {
        B=(shape/ci)*pow((zi/ci),shape-1)*pow((pow((zi/ci),shape)+1),-2);
        C=(shape/cj)*pow((zj/cj),shape-1)*pow((pow((zj/cj),shape)+1),-2);
        res=B*C;
    }
    return(res);
    
}



/*******************************Weibull**/
double biv_Weibull(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double z=0.0,a=0.0,A=0.0,k=0.0,res=0.0,B=0.0,C=0.0;double ci=exp(mui);double cj=exp(muj);;
    k=tgamma(1+1/shape);
    if(corr)   {
        a=1-pow(corr,2);
        z=2*fabs(corr)*pow(zi*zj,shape/2)*pow(k,shape)/(pow(ci*cj,shape/2)*a);
        A=pow(shape,2)*pow(k,2*shape)*pow(zi*zj,shape-1);
        B= exp(-pow(k,shape)*(pow(zi/ci,shape)+pow(zj/cj,shape))/a);
        C=a*pow(ci*cj,shape);
        res=A*B*bessel_ii(z,0,2)/(exp(-z)*C);
        
    }
    else
    {
        B=shape*pow(k,shape)*pow(zi,shape-1)* exp(-pow(k,shape)*pow(zi/ci,shape))/pow(ci,shape);
        C=shape*pow(k,shape)*pow(zj,shape-1)* exp(-pow(k,shape)*pow(zj/cj,shape))/pow(cj,shape);
        res=B*C;
    }
    return(res);
    
}


/*********** bivariate T distribution********************/

double biv_T(double rho,double zi,double zj,double nuu)
{
    double nu=1/nuu;
    int k=0;
    double B=0.0,C=0.0,res0=0.0,RR=0.0,pp1=0.0,pp2=0.0;
    double bb1,bb2;
    double x=zi;double y=zj;
    double cc=(nu+1)/2; double nu2=nu/2;
    double x1=(x*x+nu); double y1=(y*y+nu);
    double rho2=pow(1-rho*rho,-cc);
    double b1 = (pow(nu,nu))*pow(x1*y1,-cc)*pow(tgamma(cc),2);
    double c1 = M_PI*pow(tgamma(nu/2),2)*rho2;
    double b2 = rho*x*y*pow(nu,nu+2)*pow(x1*y1,-nu2-1);
    double c2 = 2*M_PI*rho2;
    double a1 = 0; double a2 = 0;
    double aux  = pow(rho*x*y,2)/(x1*y1);
    double aux1 = pow(rho*nu,2)/(x1*y1);
    /// indipendent t distributions
    if(fabs(rho)<=EPS1)
    {
    C = lgamma(cc)+log(pow((1+x*x/nu),-cc))-log(sqrt(M_PI*nu))-lgamma(nu/2);
    B = lgamma(cc)+log(pow((1+y*y/nu),-cc))-log(sqrt(M_PI*nu))-lgamma(nu/2);
    return(exp(B)*exp(C));
    }
    while( k<=6000 )
    {
        //pp1=log(hypergeo(cc+k,cc+k,0.5,aux));
        pp1=(0.5-2*(cc+k))*log(1-aux)+log(hypergeo(0.5-(cc+k),0.5-(cc+k),0.5,aux)); //euler
        bb1=pp1+k*log(aux1)+2*(lgamma(cc+k)-lgamma(cc))-lgamma(k+1.0)-lgamma(nu2+k)+lgamma(nu2);
        a1 = a1 + exp(bb1);
        //pp2=log(hypergeo(nu2+1+k,nu2+1+k,1.5,aux));
        pp2=(1.5-2*(nu2+1+k))*log(1-aux)+log(hypergeo(1.5-(nu2+1+k),1.5-(nu2+1+k),1.5,aux));//euler
        bb2=pp2+k*log(aux1)+2*log((1+k/nu2))+lgamma(nu2+k)-lgamma(k+1.0)-lgamma(nu2);
        a2 = a2 + exp(bb2);
        RR=(b1/c1)*a1+(b2/c2)*a2;
        if(RR>DBL_MAX) return(res0);
        if((fabs(RR-res0)<1e-30) ) {break;}
        else {res0=RR;}
        k++;
    }
        return(RR);
}

double appellF4(double a,double b,double c,double d,double x,double y)

{

    double RR=0.0,bb=0.0,res0=0.0;
    int k=0;
    while( k<=5000 )
    {
        bb=k*log(y)+(lgamma(a+k)+lgamma(b+k)+lgamma(d))-(lgamma(a)+lgamma(b)+
                                                         lgamma(d+k)+lgamma(k+1.0))+
        (c-(a+k)-(b+k))*log(1-x)+log(hypergeo(c-a-k,c-b-k,c,x)); //euler
        RR=RR+exp(bb);
           if(RR>DBL_MAX) return(res0);
        if((fabs(RR-res0)<1e-30)  ) {break;}
        else {res0=RR;}
        k++;
    }
    return(RR);
}


double appellF4_mod(double nu,double rho2,double x,double y)
{
    double xx,yy,x2,y2,arg,arg1,pp1,pp2,app;
    xx=x*x;yy=y*y;
    x2=xx+nu;
    y2=yy+nu;
    arg=(nu+1)/2;
    arg1=nu/2;
    pp1=pow(nu,nu)*pow(x2*y2,-arg)*pow(tgamma(arg),2);
    pp2=M_PI*pow(tgamma(arg1),2)*pow(1-rho2,-arg);
    app=appellF4(arg,arg,0.5,arg1,rho2*xx*yy/(x2*y2), nu*nu*rho2/(x2*y2));
    return(4*pp1*app/pp2);
}




/*********** bivariate two piece-T distribution********************/
double biv_two_pieceT(double rho,double zi,double zj,double sill,double nuu,double eta,
                      double p11,double mui,double muj)
{
    double res;
    double nu=1/nuu;
    double etamas=1+eta;
    double etamos=1-eta;
    double rho2=rho*rho;
    double zistd=(zi-mui)/sqrt(sill);
    double zjstd=(zj-muj)/sqrt(sill);
    if(zi>=mui&&zj>=muj)
    {res=          (p11/pow(etamos,2))*appellF4_mod(nu,rho2,zistd/etamos,zjstd/etamos);}
    if(zi>=mui&&zj<muj)
    {res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho2,zistd/etamos,zjstd/etamas);}
    if(zi<mui&&zj>=muj)
    {res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho2,zistd/etamas,zjstd/etamos);}
    if(zi<mui&&zj<muj)
    {res=    ((p11+eta)/pow(etamas,2))*appellF4_mod(nu,rho2,zistd/etamas,zjstd/etamas);}
    return(res/sill);
}

// Start: biv_two_pieceGaussian:
double biv_half_Gauss(double rho,double zi,double zj)
{
    double kk=0, dens=0,a=0,b=0,rho2=rho*rho;
    kk=(M_PI)*sqrt(1-rho2);
    a=exp(- (1/(2*(1-rho2)))*(pow(zi,2)+pow(zj,2)-2*rho*zi*zj));
    b=exp(- (1/(2*(1-rho2)))*(pow(zi,2)+pow(zj,2)+2*rho*zi*zj));
    dens=(a+b)/kk;
    return(dens);
}

double biv_two_pieceGaussian(double rho,double zi,double zj,double sill,double eta,
                             double p11,double mui,double muj)
{
    double res;
    double etamas=1+eta;
    double etamos=1-eta;
    double zistd=(zi-mui)/sqrt(sill);
    double zjstd=(zj-muj)/sqrt(sill);
    if(zi>=mui&&zj>=muj)
    {res=          (p11/pow(etamos,2))*biv_half_Gauss(rho,zistd/etamos,zjstd/etamos);}
    if(zi>=mui&&zj<muj)
    {res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Gauss(rho,zistd/etamos,zjstd/etamas);}
    if(zi<mui&&zj>=muj)
    {res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Gauss(rho,zistd/etamas,zjstd/etamos);}
    if(zi<mui&&zj<muj)
    {res=    ((p11+eta)/pow(etamas,2))*biv_half_Gauss(rho,zistd/etamas,zjstd/etamas);}
    return(res/sill);
}
// End: biv_two_pieceGaussian:

double log_biv_Norm(double corr,double zi,double zj,double mi,double mj,double vari, double nugget)
{
    double u,v,u2,v2,det,s1,s12,dens;
    u=zi-mi;
    v=zj-mj;
    u2=pow(u,2);v2=pow(v,2);
    s1=vari+nugget;s12=vari*corr;
    det=pow(s1,2)-pow(s12,2);
    dens=(-0.5*(2*log(2*M_PI)+log(det)+(s1*(u2+v2)-2*s12*u*v)/det));
    return(dens);
}
// ===================================== END: Distributions.c  ==================================//

/******************************************************************************************/
/********************* SPATIAL CASE *****************************************************/
/******************************************************************************************/
__kernel void Comp_Cond_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
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
    
    
    int ncoord  = int_par[1];
    int cormod  = int_par[0];
    int type    = int_par[3];
    
    
    
    
    s1=nuis0+nuis1;
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            
            if(lags<=maxdist){
            
                s12=nuis1*CorFct(cormod, lags, 0, par0,par1,par2,par3,0,0);
                det=pow(s1,2)-pow(s12,2);
            
                u=data[gid+j]-mean[gid+j];
                v=data[j]-mean[j];
                
                if(!isnan(u)&&!isnan(v) )
                {
                    u2=pow(u,2);
                    v2=pow(v,2);
                    sum+= (-log(2*M_PI)-log(det)+log(s1)+
                           (u2+v2)*(0.5/s1-s1/det)+2*s12*u*v/det)*weights;
                }
                
            }
            
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}


__kernel void Comp_Pair_Gauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    
    double  corr=0.0, lags=0.0,weights=1.0, sum=0.0;
    //double det=0.0, u=0.0, u2=0.0, v=0.0, v2=0.0;
    
    double maxdist = dou_par[6];
    double nuis0 = dou_par[4];//nugget
    double nuis1 = dou_par[5]; //sill
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    
    
    int cormod  = int_par[0];
    int ncoord  = int_par[1];
    int type    = int_par[3];
    
    //s1=nuis0+nuis1;
    
    
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            if(lags<=maxdist){
                corr=CorFct(cormod, lags, 0, par0,par1,par2,par3,0,0);
                if(!isnan(data[gid+j])&&!isnan(data[j]) )
                {
                    sum+=log_biv_Norm(corr,data[gid+j],data[j],mean[gid+j],mean[j],nuis1,nuis0)*weights;
                }}}
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}


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
                    sum+=  -0.5*(log(2*M_PI)+log(vario)+
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
    double u=0.0,v=0.0,weights=1.0;
    double lags=0.0, corr=0.0,sum=0.0,wrap_gauss,alfa=2.0;
    
    
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
                
                u=data[gid+j];//-2*atan(mean[gid+j])-M_PI;
                v=data[j];//-2*atan(mean[j])-M_PI;;
                //printf("u,v: %f\t%f\n",u,v);
                if(!isnan(u)&&!isnan(v) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    wrap_gauss=biv_wrapped(alfa,u,v,mean[gid+j],mean[j],nuis0,nuis1,corr);
                    
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    sum+=  log(wrap_gauss)*weights ;
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
                zi=(data[gid+j]);
                zj=(data[j]);
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    bb=log(biv_sinh(corr,zi,zj,mean[gid+j],mean[j],nuis2,nuis3,nuis1));
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
                    bb=log(d2lognorm(zi,zj,nuis1,nuis0, mean[gid+j], mean[j],corr));
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
                zi=data[gid+j];
                zj=data[j];
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    sum+=  weights*log(biv_skew((1-nuis0)*corr,zi,zj,mean[gid+j],mean[j],nuis1,nuis2));
                }
                
            }
        }
        
        else
            continue;
    }
    res[gid] = sum;
    
}



__kernel void Comp_Pair_Gamma2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    
    
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
    double corr,zi,zj,lags,bb=0.0,weights=1.0,sum=0.0,sill=1-nuis0;;
    
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
                    
                    sum+=  weights*log(biv_gamma(sill*corr,zi,zj,mean[gid+j],mean[j],nuis2));
                    //printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n\n\n",sum,corr,zj,zi,mean[gid+j],mean[j],nuis2/2-1);
                    
                }
                
            }
        }
        
        else
            continue;
    }
    if(isnan(sum)) sum=LOW;
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
                    
                } } }
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

__kernel void Comp_Pair_PoisbinGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
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
                    dens=biv_poisbin(NN,uu,vv,p1,p2,psj);//
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


__kernel void Comp_Pair_PoisbinnegGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
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
                    //printf("GPU lags:%f\n",lags);
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    uu=(int) u;
                    vv=(int) v;
                    //printf("dens:%f\n",dens);
                    dens=biv_poisbinneg(NN,vv,uu,p2,p1,psj);// EL ORDEN IMPORTA
                    //printf("GPU dens %f\t%d\t%d\t%f\t%f\t%f\n",dens,uu,vv,p1,p2,psj);
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





__kernel void Comp_Pair_Logistic2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double corr,zi,zj,lags,weights=1.0,sum=0.0;
    
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
                    sum+=  weights*log(biv_Logistic(corr,zi,zj,mean[gid+j],mean[j],nuis1));
                    //sum+=  weights*0.5;
                }
                
            }
        }
        
        else
            continue;
    }
    res[gid] = sum;
    
}



__kernel void Comp_Pair_LogLogistic2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    double corr,zi,zj,lags,weights=1.0,sum=0.0;
    
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
                    sum+=  weights*log(biv_LogLogistic(corr,zi,zj,mean[gid+j],mean[j],nuis2));
                    //sum+=  weights*0.5;
                }
                
            }
        }
        
        else
            continue;
    }
    res[gid] = sum;
    
}



__kernel void Comp_Pair_Weibull2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    
    
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
    
    double corr,zi,zj,lags,bb=0.0,weights=1.0,sum=0.0,sill=1-nuis0;
    
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
                    
                    sum+=  weights*log(biv_Weibull(sill*corr,zi,zj,mean[gid+j],mean[j],nuis2));
                    //printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n\n\n",sum,corr,zj,zi,mean[gid+j],mean[j],nuis2/2-1);
                    
                }
                
            }
        }
        
        else
            continue;
    }
    //if(isnan(sum)) sum=LOW;
    res[gid] = sum;
    
}



__kernel void Comp_Pair_T2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    
    double lags,weights=1.0, sum=0.0;
    double zi, zj, bl,corr;
    
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
            //lags=0.5;
            //printf("lags:%f \n",lags);
            if(lags<=maxdist){
                
                zi=data[j];
                zj=data[gid+j];
                //printf("zi:%f zj:%f \n",zi,zj);
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    //corr=CorFct(cormod,lags,0,par,0,0);
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    bl=biv_T(corr*(1-nuis1),(zi-mean[j])/sqrt(nuis2),(zj-mean[gid+j])/sqrt(nuis2),nuis0)/nuis2;
                    //printf("corr:%f bl:%f\n",corr,bl);
                    sum+= weights*log(bl);
                    //printf("corr:%f bl:%f sum:%f\n",corr,bl,sum);
                    
                }
                
            }
            
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}



__kernel void Comp_Pair_TWOPIECET2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    
    double lags,weights=1.0, sum=0.0;
    double zi, zj, bl,corr,p11,qq;
    
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
    
    
    
    qq=qnorm55((1-nuis3)/2,0,1,1,0);
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            //lags=0.5;
            //printf("lags:%f \n",lags);
            if(lags<=maxdist){
                
                zi=data[j];
                zj=data[gid+j];
                //printf("zi:%f zj:%f \n",zi,zj);
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    //corr=CorFct(cormod,lags,0,par,0,0);
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    p11=pbnorm(cormod,lags,0,qq,qq,0,1,par0,par1,par2,par3,0);
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    bl=biv_two_pieceT(corr,zi,zj,nuis2,nuis0,nuis3,p11,mean[j],mean[gid+j]);
                    //printf("corr:%f bl:%f\n",corr,bl);
                    sum+= weights*log(bl);
                    //printf("corr:%f bl:%f sum:%f\n",corr,bl,sum);
                    
                }
                
            }
            
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}

__kernel void Comp_Pair_TWOPIECEGauss2_OCL(__global const double *coordx,__global const double *coordy,__global const double *mean, __global const double *data, __global double *res,__global const int *int_par,__global const double *dou_par)
{
    
    int j, gid = get_global_id(0);
    
    double lags,weights=1.0, sum=0.0;
    double zi, zj, bl,corr,p11,qq;
    
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
    
    
    
    qq=qnorm55((1-nuis2)/2,0,1,1,0);
    for (j = 0; j < ncoord; j++) {
        if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
        {
            
            lags = dist(type,coordx[j],coordx[gid+j],coordy[j],coordy[gid+j],REARTH);
            //lags=0.5;
            //printf("lags:%f \n",lags);
            if(lags<=maxdist){
                
                zi=data[j];
                zj=data[gid+j];
                //printf("zi:%f zj:%f \n",zi,zj);
                
                if(!isnan(zi)&&!isnan(zj) )
                {
                    //corr=CorFct(cormod,lags,0,par,0,0);
                    corr=CorFct(cormod,lags,0,par0,par1,par2,par3,0,0);
                    p11=pbnorm(cormod,lags,0,qq,qq,0,1,par0,par1,par2,par3,0);
                    if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                    bl=biv_two_pieceGaussian(corr,zi,zj,nuis1,nuis2,p11,mean[j],mean[gid+j]);
                    //printf("corr:%f bl:%f\n",corr,bl);
                    sum+= weights*log(bl);
                    //printf("corr:%f bl:%f sum:%f\n",corr,bl,sum);
                    
                }
                
            }
            
        }
        
        else
            continue;
    }
    
    res[gid] = sum;
    
}

/******************************************************************************************/
/********************* SPACE TIME CASE ***************************************************/
/******************************************************************************************/


__kernel void Comp_Pair_Gauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
{
    
    double maxdist = dou_par[6];
    double maxtime	=	dou_par[11];
    double nuis0 = dou_par[4];//nugget
    double nuis1 = dou_par[5];//sill
    double nuis2 = dou_par[9];
    double nuis3 = dou_par[10];
    double par0 = dou_par[0];
    double par1 = dou_par[1];
    double par2 = dou_par[2];
    double par3 = dou_par[3];
    double REARTH = dou_par[8];
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
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
    double  lags=0.0, lagt=0.0,weights=1.0,sum=0.0, corr=0.0;
    double u=0.0,  w=0.0;
    
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
    
    //for (j = 0; j < ncoord; j++) {
    //    if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    //if(l >= ncoord) isValid = false;
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        corr=CorFct_st(cormod,lags,0,par0,par1,par2,par3,par4,par5,par6,t,v);
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= log_biv_Norm(corr,u,w,mean[(l+NS[t])],mean[(m+NS[v])],nuis1,nuis0)*weights;}
                    }
                }
            }
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        
                        corr=CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,t,v);
                        
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= log_biv_Norm(corr,u,w,mean[(l+NS[t])],mean[(m+NS[v])],nuis1,nuis0)*weights;}
                }
                }
            }
        }
      res[i] = sum;
    }
}




__kernel void Comp_Pair_WrapGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
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
    double lags=0.0, lagt=0.0,weights=1.0,sum=0.0,corr=0.0;
    double u=0.0, u2=0.0, w=0.0;
    double wrap_gauss,alfa=2.0;
    
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
    
    //for (j = 0; j < ncoord; j++) {
    //    if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    //if(l >= ncoord) isValid = false;
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            
                            corr=CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,t,v);
                            wrap_gauss=biv_wrapped(alfa,u,w,mean[(l+NS[t])],mean[(m+NS[v])],nuis0,nuis1,corr);
                            
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= log(wrap_gauss)*weights;
                        }
                        //printf("GPU: %d\t%d\t%d\t%d\t%f\n",l,t,v,m,sum);
                    }}}
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            
                            corr=CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,t,v);
                            wrap_gauss=biv_wrapped(alfa,u,w,mean[(l+NS[t])],mean[(m+NS[v])],nuis0,nuis1,corr);
                            
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= log(wrap_gauss)*weights;
                        }
                    }
                }}}
        res[i] = sum;
    }
}


__kernel void Comp_Pair_PoisbinGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[6];
    
    
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    //int ls = get_global_size(0);
    //int ms = get_global_size(1);
    
    int m=0,v =0;
    int uu=0,ww=0;
    double dens=0.0,lags=0.0,lagt=0.0,weights=1.0,u,w, sum=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
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
    
    //for (j = 0; j < ncoord; j++) {
    //    if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    //if(l >= ncoord) isValid = false;
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        psj=pbnorm_st(cormod,lags,0,mean[(l+NS[t])],mean[(m+NS[v])],nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(mean[(l+NS[t])],0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(mean[(m+NS[v])],0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                         
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            dens=biv_poisbin(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                        //printf("GPU: %d\t%d\t%d\t%d\t%f\n",l,t,v,m,sum);
                    }}}
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        psj=pbnorm_st(cormod,lags,lagt,mean[(l+NS[t])],mean[(m+NS[v])],nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(mean[(l+NS[t])],0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(mean[(m+NS[v])],0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        
                        if(!isnan(u)&&!isnan(w) ){
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            dens=biv_poisbin(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                    }}}}
        res[i] = sum;
    }
}


__kernel void Comp_Pair_PoisbinnegGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[6];
    
    
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    //int ls = get_global_size(0);
    //int ms = get_global_size(1);
    
    int m=0,v =0;
    int uu=0,ww=0;
    double dens=0.0,lags=0.0,lagt=0.0,weights=1.0,u,w, sum=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
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
    
    //for (j = 0; j < ncoord; j++) {
    //    if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    //if(l >= ncoord) isValid = false;
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        psj=pbnorm_st(cormod,lags,0,mean[(l+NS[t])],mean[(m+NS[v])],nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(mean[(l+NS[t])],0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(mean[(m+NS[v])],0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            dens=biv_poisbinneg(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                        //printf("GPU: %d\t%d\t%d\t%d\t%f\n",l,t,v,m,sum);
                    }}}
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        psj=pbnorm_st(cormod,lags,lagt,mean[(l+NS[t])],mean[(m+NS[v])],nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(mean[(l+NS[t])],0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(mean[(m+NS[v])],0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        
                        if(!isnan(u)&&!isnan(w) ){
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            dens=biv_poisbinneg(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                    }}}}
        res[i] = sum;
    }
}





__kernel void Comp_Cond_Gauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
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
    
    //for (j = 0; j < ncoord; j++) {
    //    if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    //if(l >= ncoord) isValid = false;
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        s1=nuis0+nuis1;
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                
                for(m=l+1;m<ns[t];m++)
                {
                    
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        s12=nuis1*CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,t,v);
                        det=pow(s1,2)-pow(s12,2);
                        u=data[(l+NS[t])]-mean[(l+NS[t])];
                        w=data[(m+NS[v])]-mean[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            u2=pow(u,2);
                            w2=pow(w,2);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= (-log(2*M_PI)-log(det)+log(s1)+
                                   (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det)*weights;}
                    }}
            }
            else
            {
                
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        
                        s12=nuis1*CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,t,v);
                        det=pow(s1,2)-pow(s12,2);
                        
                        u=data[(l+NS[t])]-mean[(l+NS[t])];
                        w=data[(m+NS[v])]-mean[(m+NS[v])];
    
                        if(!isnan(u)&&!isnan(w) ){
                            u2=pow(u,2);
                            w2=pow(w,2);
                            
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= (-log(2*M_PI)-log(det)+log(s1)+
                                   (u2+w2)*(0.5/s1-s1/det)+2*s12*u*w/det)*weights;}
                    }
                }
            }
        }
        res[i] = sum;
    }
}



__kernel void Comp_Diff_Gauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
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
    double  lags=0.0, lagt=0.0,weights=1.0,sum=0.0;
    double vario=0.0,u=0.0, w=0.0;
    
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
    
    //for (j = 0; j < ncoord; j++) {
    //    if (   ((gid+j)!= j) && ((gid+j) < ncoord)   )
    
    bool isValid = true;
    //printf("%d\t%d\n",l,t);
    
    //if(l >= ncoord) isValid = false;
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        vario=Variogram_st(cormod,lags,0,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6);
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= (-0.5*(log(2*M_PI)+log(vario)+
                                         pow(u-w,2)/(2*vario)))*weights;}
                    }}
            }
            else
            {
                
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        vario=Variogram_st(cormod,lags,lagt,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6);
                        u=data[(l+NS[t])];
                        w=data[(m+NS[v])];
                        
                        if(!isnan(u)&&!isnan(w) ){
                            
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= (-0.5*(log(2*M_PI)+log(vario)+
                                         pow(u-w,2)/(2*vario)))*weights;}
                    }
                }
            }
        }
        res[i] = sum;
    }
}



__kernel void Comp_Pair_SkewGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
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
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                
                for(m=l+1;m<ns[t];m++)
                {
                    
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,t,v);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= weights*log(biv_skew((1-nuis0)*corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis1,nuis2));
                        }
                    }}
            }
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,t,v);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= weights*log(biv_skew((1-nuis0)*corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis1,nuis2));}
                    }
                }
            }
        }
        res[i] = sum;
    }
}

__kernel void Comp_Pair_SinhGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);

    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                
                for(m=l+1;m<ns[t];m++)
                {
                    
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        
                        zi=(data[(l+NS[t])]);
                        zj=(data[(m+NS[v])]);
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= weights*log(biv_sinh(corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2,nuis3,nuis1));
                        }
                    }}
            }
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= weights*log(biv_sinh(corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2,nuis3,nuis1));}
                    }
                }
            }
        }
        res[i] = sum;
    }
}


__kernel void Comp_Pair_Gamma_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0,sill=1-nuis0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= weights*log(biv_gamma(sill*corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2));
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= weights*log(biv_gamma(sill*corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2));}
                    }}}}
        res[i] = sum;
    }
}




__kernel void Comp_Pair_LogGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= weights*log(d2lognorm(zi,zj,nuis1,nuis0, mean[(l+NS[t])], mean[(m+NS[v])],corr));
                            
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= weights*log(d2lognorm(zi,zj,nuis1,nuis0, mean[(l+NS[t])], mean[(m+NS[v])],corr));}
                    }}}}
        res[i] = sum;
    }
}


__kernel void Comp_Pair_BinomGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[6];
    
    
    
    int l = get_global_id(0);
    int t = get_global_id(1);

    
    int m=0,v =0;
    int uu=0,ww=0;
    double dens=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0, a=0.0,b=0.0,sum=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);

    
    bool isValid = true;

    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        a = mean[(l+NS[t])];
                        b = mean[(m+NS[v])];
                        
                        psj=pbnorm_st(cormod,lags,0,a,b,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(a,0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(b,0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            dens=biv_binom(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                        //printf("GPU: %d\t%d\t%d\t%d\t%f\n",l,t,v,m,sum);
                    }}}
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        a = mean[(l+NS[t])];
                        b = mean[(m+NS[v])];
                        
                        psj=pbnorm_st(cormod,lags,lagt,a,b,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(mean[(l+NS[t])],0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(mean[(m+NS[v])],0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        
                        if(!isnan(u)&&!isnan(w) ){
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            dens=biv_binom(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                    }}}}
        res[i] = sum;
    }
}


__kernel void Comp_Pair_BinomnegGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    int NN        = int_par[6];
    
    
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    int uu=0,ww=0;
    double dens=0.0,lags=0.0,lagt=0.0,weights=1.0,u=0.0,w=0.0, a=0.0,b=0.0,sum=0.0;
    double p1=0.0,p2=0.0;//probability of marginal success
    double psj=0.0;//probability of joint success
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
    {
        for(v = t;v<ntime;v++)
        {
            if(t==v)
            {
                for(m=l+1;m<ns[t];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist)
                    {
                        a = mean[(l+NS[t])];
                        b = mean[(m+NS[v])];
                        
                        psj=pbnorm_st(cormod,lags,0,a,b,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(a,0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(b,0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        if(!isnan(u)&&!isnan(w) ){
                            
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            dens=biv_binomneg(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                        //printf("GPU: %d\t%d\t%d\t%d\t%f\n",l,t,v,m,sum);
                    }}}
            else
            {
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++)
                {
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        a = mean[(l+NS[t])];
                        b = mean[(m+NS[v])];
                        
                        psj=pbnorm_st(cormod,lags,lagt,a,b,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                        p1=pnorm_OCL(mean[(l+NS[t])],0,1);//pnorm_OCL(ai,0,1)
                        p2=pnorm_OCL(mean[(m+NS[v])],0,1);
                        u = data[(l+NS[t])];
                        w = data[(m+NS[v])];
                        
                        if(!isnan(u)&&!isnan(w) ){
                            uu=(int) u;
                            ww=(int) w;
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            dens=biv_binomneg(NN,uu,ww,p1,p2,psj);
                            sum+= log(dens)*weights;
                        }
                    }}}}
        res[i] = sum;
    }
}


__kernel void Comp_Pair_LogLogistic_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= weights*log(biv_LogLogistic(corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2));
                            
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= weights*log(biv_LogLogistic(corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2));
                        }}}}}
        res[i] = sum;
    }
}



__kernel void Comp_Pair_Logistic_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= weights*log(biv_Logistic(corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis1));
                            
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= weights*log(biv_Logistic(corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis1));
                        }}}}}
        res[i] = sum;
    }
}




__kernel void Comp_Pair_Weibull_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0,sill=1-nuis0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            sum+= weights*log(biv_Weibull(sill*corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2));
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            sum+= weights*log(biv_Weibull(sill*corr,zi,zj,mean[(l+NS[t])],mean[(m+NS[v])],nuis2));}
                    }}}}
        res[i] = sum;
    }
}





__kernel void Comp_Pair_T_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0,bl=0.0;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    
    if(isValid)
        
    {
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        //printf("a:%f b:%f i:%d j:%d\n",data[(l+NS[t])],data[(m+NS[v])],(l+NS[t]),(m+NS[v]));
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            //printf("a:%f b:%f\n",(mean[(l+NS[t])]),(mean[(m+NS[v])]));
                            bl=biv_T(corr*(1-nuis1),(zi-mean[(l+NS[t])])/sqrt(nuis2), (zj-mean[(m+NS[v])])/sqrt(nuis2),nuis0)/nuis2;
                            sum+= weights*log(bl);
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            bl=biv_T(corr*(1-nuis1),(zi-mean[(l+NS[t])])/sqrt(nuis2),(zj-mean[(m+NS[v])])/sqrt(nuis2),nuis0)/nuis2;
                            sum+= weights*log(bl);
                        }
                    }}}
        }
        res[i] = sum;
    }
}





__kernel void Comp_Pair_TWOPIECET_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0,bl=0.0,qq,p11;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    
    if(isValid)
        
    {
        qq=qnorm55((1-nuis3)/2,0,1,1,0);
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        //printf("a:%f b:%f i:%d j:%d\n",data[(l+NS[t])],data[(m+NS[v])],(l+NS[t]),(m+NS[v]));
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            //printf("a:%f b:%f\n",(mean[(l+NS[t])]),(mean[(m+NS[v])]));
                            
                            //p11=pbnorm_st(cormod,lags,0,qq,qq,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                            p11=pbnorm(cormod,lags,0,qq,qq,0,1,par0,par1,par2,par3,0);
                            bl=biv_two_pieceT(corr,zi,zj,nuis2,nuis0,nuis3,p11,mean[(l+NS[t])],mean[(m+NS[v])]);
                            
                            
                            sum+= weights*log(bl);
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            
                            //p11=pbnorm_st(cormod,lags,0,qq,qq,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                            p11=pbnorm(cormod,lags,0,qq,qq,0,1,par0,par1,par2,par3,0);
                            bl=biv_two_pieceT(corr,zi,zj,nuis2,nuis0,nuis3,p11,mean[(l+NS[t])],mean[(m+NS[v])]);
                            
                            sum+= weights*log(bl);
                        }
                    }}}
        }
        res[i] = sum;
    }
}





__kernel void Comp_Pair_TWOPIECEGauss_st2_OCL(__global const double *coordt,__global const double *coordx,__global const double *coordy,__global const double *data,__global const double *mean,  __global double *res,__global const int *int_par,__global const double *dou_par,__global const int *ns,__global const int *NS)
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
    double par4 = dou_par[12];
    double par5 = dou_par[13];
    double par6 = dou_par[14];
    
    int cormod      = int_par[0];
    int ncoord      = int_par[1];
    int ntime      = int_par[5];
    int weigthed    = int_par[2];
    int type        = int_par[3];
    
    int l = get_global_id(0);
    int t = get_global_id(1);
    
    
    int m=0,v =0;
    double corr=0.0,zi=0.0,zj=0.0,lags=0.0,lagt=0.0,weights=1.0, sum=0.0,bl=0.0,qq,p11;
    
    int m1 = get_local_id(0);
    int v1 = get_local_id(1);
    
    int lsize_m = get_local_size(0);
    int lsize_v = get_local_size(1);
    
    int wx = (l-m1)/lsize_m;
    int wy = (t-v1)/lsize_v;
    
    int gidx = (wx*lsize_m+m1);
    int gidy = (wy*lsize_v+v1);
    
    int i = (ncoord*gidy+gidx);
    
    bool isValid = true;
    
    if(l >= ns[t]) isValid = false;
    
    if(t >= ntime) isValid = false;
    
    
    if(isValid)
        
    {
        qq=qnorm55((1-nuis2)/2,0,1,1,0);
        for(v = t;v<ntime;v++){
            if(t==v){
                for(m=l+1;m<ns[t];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lags<=maxdist){
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        //printf("a:%f b:%f i:%d j:%d\n",data[(l+NS[t])],data[(m+NS[v])],(l+NS[t]),(m+NS[v]));
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, 0,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist);}
                            //printf("a:%f b:%f\n",(mean[(l+NS[t])]),(mean[(m+NS[v])]));
                            
                            //p11=pbnorm_st(cormod,lags,0,qq,qq,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                            p11=pbnorm(cormod,lags,0,qq,qq,0,1,par0,par1,par2,par3,0);
                            bl=biv_two_pieceGaussian(corr,zi,zj,nuis1,nuis2,p11,mean[(l+NS[t])],mean[(m+NS[v])]);
                            
                            
                            sum+= weights*log(bl);
                        }}}}
            else{
                lagt=fabs(coordt[t]-coordt[v]);
                for(m=0;m<ns[v];m++){
                    lags=dist(type,coordx[l],coordx[m],coordy[l],coordy[m],REARTH);
                    if(lagt<=maxtime && lags<=maxdist)
                    {
                        zi=data[(l+NS[t])];
                        zj=data[(m+NS[v])];
                        if(!isnan(zi)&&!isnan(zj) ){
                            corr =CorFct_st(cormod,lags, lagt,par0,par1,par2,par3,par4,par5,par6,0,0);
                            if(weigthed) {weights=CorFunBohman(lags,maxdist)*CorFunBohman(lagt,maxtime);}
                            
                            //p11=pbnorm_st(cormod,lags,0,qq,qq,nuis0,nuis1,par0,par1,par2,par3,par4,par5,par6,0);
                            p11=pbnorm(cormod,lags,0,qq,qq,0,1,par0,par1,par2,par3,0);
                            bl=biv_two_pieceGaussian(corr,zi,zj,nuis1,nuis2,p11,mean[(l+NS[t])],mean[(m+NS[v])]);
                            
                            sum+= weights*log(bl);
                        }
                    }}}
        }
        res[i] = sum;
    }
}
