#include "header.h"

//************************************* START hyperg.c*****************************************
// FUNCTION: 1F1
//Source: https://github.com/scipy/scipy/blob/master/scipy/special/cephes/hyperg.c


void  hyperg_call(double *a,double *b,double *x,double *res)
{
    *res = hyperg(*a,*b,*x);
}

double hyperg(a, b, x)
double a, b, x;
{
    double asum, psum, acanc, pcanc, temp;

    /* See if a Kummer transformation will help */
    temp = b - a;
    if (fabs(temp) < 0.001 * fabs(a))
    return (exp(x) * hyperg(temp, b, -x));


    /* Try power & asymptotic series, starting from the one that is likely OK */
    if (fabs(x) < 10 + fabs(a) + fabs(b)) {
    psum = hy1f1p(a, b, x, &pcanc);
    if (pcanc < 1.0e-15)
        goto done;
    asum = hy1f1a(a, b, x, &acanc);
    }
    else {
    psum = hy1f1a(a, b, x, &pcanc);
    if (pcanc < 1.0e-15)
        goto done;
    asum = hy1f1p(a, b, x, &acanc);
    }

    /* Pick the result with less estimated error */

    if (acanc < pcanc) {
    pcanc = acanc;
    psum = asum;
    }

  done:
    if (pcanc > 1.0e-12)
    //sf_error("hyperg", SF_ERROR_LOSS, NULL);
        printf("hyperg SF_ERROR_LOSS\n");

    return (psum);
}




/* Power series summation for confluent hypergeometric function                */


double hy1f1p(double a, double b, double x, double *err)
{
    double n, a0, sum, t, u, temp, maxn;
    double an, bn, maxt;
    double y, c, sumc;

    /* set up for power series summation */
    an = a;
    bn = b;
    a0 = 1.0;
    sum = 1.0;
    c = 0.0;
    n = 1.0;
    t = 1.0;
    maxt = 0.0;
    *err = 1.0;

    maxn = 200.0 + 2 * fabs(a) + 2 * fabs(b);

    while (t > MACHEP) {
    if (bn == 0) {        /* check bn first since if both   */
        //sf_error("hyperg", SF_ERROR_SINGULAR, NULL);
        printf("hyperg SF_ERROR_SINGULAR\n");
        return (NPY_INFINITY);    /* an and bn are zero it is     */
    }
    if (an == 0)        /* a singularity            */
        return (sum);
    if (n > maxn) {
        /* too many terms; take the last one as error estimate */
        c = fabs(c) + fabs(t) * 50.0;
        goto pdone;
    }
    u = x * (an / (bn * n));

    /* check for blowup */
    temp = fabs(u);
    if ((temp > 1.0) && (maxt > (DBL_MAX / temp))) {
        *err = 1.0;        /* blowup: estimate 100% error */
        return sum;
    }

    a0 *= u;

    y = a0 - c;
    sumc = sum + y;
    c = (sumc - sum) - y;
    sum = sumc;

    t = fabs(a0);

    an += 1.0;
    bn += 1.0;
    n += 1.0;
    }

  pdone:

    /* estimate error due to roundoff and cancellation */
    if (sum != 0.0) {
    *err = fabs(c / sum);
    }
    else {
    *err = fabs(c);
    }

    if (*err != *err) {
    /* nan */
    *err = 1.0;
    }

    return (sum);
}


/*                                                     hy1f1a()        */
/* asymptotic formula for hypergeometric function:
 *
 *        (    -a
 *  --    ( |z|
 * |  (b) ( -------- 2f0( a, 1+a-b, -1/x )
 *        (  --
 *        ( |  (b-a)
 *
 *
 *                                x    a-b                     )
 *                               e  |x|                        )
 *                             + -------- 2f0( b-a, 1-a, 1/x ) )
 *                                --                           )
 *                               |  (a)                        )
 */

double hy1f1a(double a, double b, double x, double *err)

{
    double h1, h2, t, u, temp, acanc, asum, err1, err2;
    if (x == 0) {
    acanc = 1.0;
    asum = NPY_INFINITY;
    goto adone;
    }
    temp = log(fabs(x));
    t = x + temp * (a - b);
    u = -temp * a;

    if (b > 0) {
    temp = lgam(b);
    t += temp;
    u += temp;
    }

    h1 = hyp2f0(a, a - b + 1, -1.0 / x, 1, &err1);

    temp = exp(u) / gamma(b - a);
    h1 *= temp;
    err1 *= temp;

    h2 = hyp2f0(b - a, 1.0 - a, 1.0 / x, 2, &err2);

    if (a < 0)
    temp = exp(t) / gamma(a);
    else
    temp = exp(t - lgam(a));

    h2 *= temp;
    err2 *= temp;

    if (x < 0.0)
    asum = h1;
    else
    asum = h2;

    acanc = fabs(err1) + fabs(err2);

    if (b < 0) {
    temp = gamma(b);
    asum *= temp;
    acanc *= fabs(temp);
    }


    if (asum != 0.0)
    acanc /= fabs(asum);

    if (acanc != acanc)
    /* nan */
    acanc = 1.0;

    if (asum == NPY_INFINITY || asum == -NPY_INFINITY)
    /* infinity */
    acanc = 0;

    acanc *= 30.0;        /* fudge factor, since error of asymptotic formula
                 * often seems this much larger than advertised */

  adone:


    *err = acanc;
    return (asum);
}

/*                                                     hyp2f0()        */

double hyp2f0(double a, double b, double x, int type, double *err)

{
    double a0, alast, t, tlast, maxt;
    double n, an, bn, u, sum, temp;

    an = a;
    bn = b;
    a0 = 1.0e0;
    alast = 1.0e0;
    sum = 0.0;
    n = 1.0e0;
    t = 1.0e0;
    tlast = 1.0e9;
    maxt = 0.0;

    do {
    if (an == 0)
        goto pdone;
    if (bn == 0)
        goto pdone;

    u = an * (bn * x / n);

    /* check for blowup */
    temp = fabs(u);
    if ((temp > 1.0) && (maxt > (DBL_MAX / temp)))
        goto error;

    a0 *= u;
    t = fabs(a0);

    /* terminating condition for asymptotic series:
     * the series is divergent (if a or b is not a negative integer),
     * but its leading part can be used as an asymptotic expansion
     */
    if (t > tlast)
        goto ndone;

    tlast = t;
    sum += alast;        /* the sum is one term behind */
    alast = a0;

    if (n > 200)
        goto ndone;

    an += 1.0e0;
    bn += 1.0e0;
    n += 1.0e0;
    if (t > maxt)
        maxt = t;
    }
    while (t > MACHEP);


  pdone:            /* series converged! */

    /* estimate error due to roundoff and cancellation */
    *err = fabs(MACHEP * (n + maxt));

    alast = a0;
    goto done;

  ndone:            /* series did not converge */

    /* The following "Converging factors" are supposed to improve accuracy,
     * but do not actually seem to accomplish very much. */

    n -= 1.0;
    x = 1.0 / x;

    switch (type) {        /* "type" given as subroutine argument */
    case 1:
    alast *=
        (0.5 + (0.125 + 0.25 * b - 0.5 * a + 0.25 * x - 0.25 * n) / x);
    break;

    case 2:
    alast *= 2.0 / 3.0 - b + 2.0 * a + x - n;
    break;

    default:
    ;
    }

    /* estimate error due to roundoff, cancellation, and nonconvergence */
    *err = MACHEP * (n + maxt) + fabs(a0);

  done:
    sum += alast;
    return (sum);

    /* series blew up: */
  error:
    *err = NPY_INFINITY;
    //sf_error("hyperg", SF_ERROR_NO_RESULT, NULL);
    printf("hyperg SF_ERROR_NO_RESULT\n");
    return (sum);
}



//************************************* END hyperg.c*****************************************




// ===================================== START: Bivariate Normal  =====================================//


//#define HSQRT 1.414213562373095048801688724209698078569671
#define HSQRT 1.4142

// https://www.jstatsoft.org/article/view/v052i10/v52i10.pdf

double Phi(double x)
{
    double val =(1+     (1-erfc(x/HSQRT) )    )/2;
    //double val =(1+     (1-0.1 )    )/2;
    
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
    while( cond>DEPSILON )
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
        sol1 =  fmax( ( 1.0 + c1 * asr ) * comp, b * comp - fmax( 0.0, res ) );
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
        a = R_pow( ( x - y ) / x / s - tmp,2 );
    }
    else if( rho < -0.99 )
    {
        double tmp = sqrt( ( 1.0 + rho ) / ( 1.0 - rho ) );
        b2 = -fabs( ( x + y ) / s - x * tmp );
        a = R_pow( ( x + y ) / x / s - tmp,2 );
    }
    else
    {
        b2 = -fabs( rho * x - y ) / s;
        a = R_pow( b2 / x ,2);
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
            sol = (Phi( fmin( x, y ) ));
        }
        else
        {
            //return (max( 0.0, min( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ));
            sol = (fmax( 0.0, fmin( 1.0, Phi( x ) + Phi( y ) - 1.0 ) ));
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
        sol = (fmax( 0.0,
                   fmin( 1.0,
                       Phi2help( x, y, rho ) + Phi2help( y, x, rho ) ) ));
    }
    
    return (sol);
}



/*for bivariate t distributions*/
double A[] = {
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2
};

double B[] = {
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5
};

double C[] = {
    /* 1.00000000000000000000E0, */
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6
};
double polevl(double x, const double coef[], int N)
{
    double ans;
    int i;
    const double *p;
    
    p = coef;
    ans = *p++;
    i = N;
    
    do
        ans = ans * x + *p++;
    while (--i);
    
    return (ans);
}

/*                                                     p1evl() */
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

double p1evl(double x, const double coef[], int N)
{
    double ans;
    const double *p;
    int i;
    
    p = coef;
    ans = x + *p++;
    i = N - 1;
    
    do
        ans = ans * x + *p++;
    while (--i);
    
    return (ans);
}

// Integrand for hypergeometric computation:



/* integrand  in  hypergeometric function*/
double int_gen_hyp(double x,double a, double b,double z,double c)
{
    double res=0.0;
    res=R_pow(x,b-1)*R_pow(1-x,c-b-1)* R_pow(1-z*x,-a);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}
void integr_gen_hyp(double *x, int n, void *ex){
    int i;double a,b,c,y;
    a =    ((double*)ex)[0];  //a
    b = ((double*)ex)[1];  //b
    c=     ((double*)ex)[2];  //c
    y =     ((double*)ex)[3];  //y
    for (i=0;i<n;i++) {x[i]=int_gen_hyp(x[i],a,b,y,c);}
    return;
}
// function computing generalized wendland
double HyperG_integral(double x, double *param) {
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;         /* as instructed in WRE */
    iwork =   (int *) Calloc(subdiv, int);  /* idem */
    work = (double *) Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=0;
    upper=1;
    // Compute the integral
    Rdqags(integr_gen_hyp, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);

    Free(iwork);Free(work);
    return(result);
}




//************************************** ST igam.c*****************************************



void igam_call(double *a,double *x,double *res)
{
    *res = igam(*a,*x);
}

double igam(a, x)
double a, x;
{
    double absxma_a;
    
    /* Check zero integration limit first */
    if (x == 0)
        return (0.0);
    
    if ((x < 0) || (a <= 0)) {
        //sf_error("gammainc", SF_ERROR_DOMAIN, NULL);
       // printf("gammainc  SF_ERROR_DOMAIN\n");
        return (NAN);
    }
    
    /* Asymptotic regime where a ~ x; see [2]. */
    absxma_a = fabs(x - a) / a;
    if ((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) {
        return asymptotic_series(a, x, IGAM);
    } else if ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a))) {
        return asymptotic_series(a, x, IGAM);
    }
    
    if ((x > 1.0) && (x > a)) {
        return (1.0 - igamc(a, x));
    }
    
    return igam_series(a, x);
}


double igamc(double a, double x)
{
    double absxma_a;
    
    if ((x < 0) || (a <= 0)) {
        //sf_error("gammaincc", SF_ERROR_DOMAIN, NULL);
        //printf("gammainc  SF_ERROR_DOMAIN\n");
        return (NAN);
    } else if (x == 0) {
        return 1;
    } else if (isinf(x)) {
        return 0.0;
    }
    
    /* Asymptotic regime where a ~ x; see [2]. */
    absxma_a = fabs(x - a) / a;
    if ((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) {
        return asymptotic_series(a, x, IGAMC);
    } else if ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a))) {
        return asymptotic_series(a, x, IGAMC);
    }
    
    /* Everywhere else; see [2]. */
    if (x > 1.1) {
        if (x < a) {
            return 1.0 - igam_series(a, x);
        } else {
            return igamc_continued_fraction(a, x);
        }
    } else if (x <= 0.5) {
        if (-0.4 / log(x) < a) {
            return 1.0 - igam_series(a, x);
        } else {
            return igamc_series(a, x);
        }
    } else {
        if (x * 1.1 < a) {
            return 1.0 - igam_series(a, x);
        } else {
            return igamc_series(a, x);
        }
    }
}




double igam_fac(double a, double x)
{
    double ax, fac, res, num;
    
    if (fabs(a - x) > 0.4 * fabs(a)) {
        ax = a * log(x) - x - lgam(a);
        if (ax < -MAXLOG) {
            //sf_error("igam", SF_ERROR_UNDERFLOW, NULL);
            //printf("gammainc  SF_ERROR_DOMAIN\n");
            return 0.0;
        }
        return exp(ax);
    }
    
    fac = a + lanczos_g - 0.5;
    res = sqrt(fac / exp(1)) / lanczos_sum_expg_scaled(a);
    
    if ((a < 200) && (x < 200)) {
        res *= exp(a - x) * pow(x / fac, a);
    } else {
        num = x - a - lanczos_g + 0.5;
        res *= exp(a * log1pmx(num / fac) + x * (0.5 - lanczos_g) / fac);
    }
    
    return res;
}


/* Compute igamc using DLMF 8.9.2. */
double igamc_continued_fraction(double a, double x)
{
    int i;
    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;
    
    ax = igam_fac(a, x);
    if (ax == 0.0) {
        return 0.0;
    }
    
    /* continued fraction */
    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1 / qkm1;
    
    for (i = 0; i < MAXITER; i++) {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        yc = y * c;
        pk = pkm1 * z - pkm2 * yc;
        qk = qkm1 * z - qkm2 * yc;
        if (qk != 0) {
            r = pk / qk;
            t = fabs((ans - r) / r);
            ans = r;
        }
        else
            t = 1.0;
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if (fabs(pk) > big) {
            pkm2 *= biginv;
            pkm1 *= biginv;
            qkm2 *= biginv;
            qkm1 *= biginv;
        }
        if (t <= MACHEP) {
            break;
        }
    }
    
    return (ans * ax);
}


/* Compute igam using DLMF 8.11.4. */
double igam_series(double a, double x)
{
    int i;
    double ans, ax, c, r;
    
    ax = igam_fac(a, x);
    if (ax == 0.0) {
        return 0.0;
    }
    
    /* power series */
    r = a;
    c = 1.0;
    ans = 1.0;
    
    for (i = 0; i < MAXITER; i++) {
        r += 1.0;
        c *= x / r;
        ans += c;
        if (c <= MACHEP * ans) {
            break;
        }
    }
    
    return (ans * ax / a);
}


/* Compute igamc using DLMF 8.7.3. This is related to the series in
 * igam_series but extra care is taken to avoid cancellation.
 */
double igamc_series(double a, double x)
{
    int n;
    double fac = 1;
    double sum = 0;
    double term, logx;
    
    for (n = 1; n < MAXITER; n++) {
        fac *= -x / n;
        term = fac / (a + n);
        sum += term;
        if (fabs(term) <= MACHEP * fabs(sum)) {
            break;
        }
    }
    
    logx = log(x);
    term = -expm1(a * logx - lgam1p(a));
    return term - exp(a * logx - lgam(a)) * sum;
}


/* Compute igam/igamc using DLMF 8.12.3/8.12.4. */
double asymptotic_series(double a, double x, int func)
{
    int k, n, sgn;
    int maxpow = 0;
    double lambda = x / a;
    double sigma = (x - a) / a;
    double eta, res, ck, ckterm, term, absterm;
    double absoldterm = NPY_INFINITY;
    double etapow[NIC] = {1};
    double sum = 0;
    double afac = 1;
    
    if (func == IGAM) {
        sgn = -1;
    } else {
        sgn = 1;
    }
    
    if (lambda > 1) {
        eta = sqrt(-2 * log1pmx(sigma));
    } else if (lambda < 1) {
        eta = -sqrt(-2 * log1pmx(sigma));
    } else {
        eta = 0;
    }
    res = 0.5 * erfc(sgn * eta * sqrt(a / 2));
    
    for (k = 0; k < KIC; k++) {
        ck = d[k][0];
        for (n = 1; n < NIC; n++) {
            if (n > maxpow) {
                etapow[n] = eta * etapow[n-1];
                maxpow += 1;
            }
            ckterm = d[k][n]*etapow[n];
            ck += ckterm;
            if (fabs(ckterm) < MACHEP * fabs(ck)) {
                break;
            }
        }
        term = ck * afac;
        absterm = fabs(term);
        if (absterm > absoldterm) {
            break;
        }
        sum += term;
        if (absterm < MACHEP * fabs(sum)) {
            break;
        }
        absoldterm = absterm;
        afac /= a;
    }
    res += sgn * exp(-0.5 * a * eta * eta) * sum / sqrt(2 * M_PI * a);
    
    return res;
}







double ratevl(double x, const double num[], int M,
              const double denom[], int N)
{
    int i, dir;
    double y, num_ans, denom_ans;
    double absx = fabs(x);
    const double *p;
    
    if (absx > 1) {
        /* Evaluate as a polynomial in 1/x. */
        dir = -1;
        p = num + M;
        y = 1 / x;
    } else {
        dir = 1;
        p = num;
        y = x;
    }
    
    /* Evaluate the numerator */
    num_ans = *p;
    p += dir;
    for (i = 1; i <= M; i++) {
        num_ans = num_ans * y + *p;
        p += dir;
    }
    
    /* Evaluate the denominator */
    if (absx > 1) {
        p = denom + N;
    } else {
        p = denom;
    }
    
    denom_ans = *p;
    p += dir;
    for (i = 1; i <= N; i++) {
        denom_ans = denom_ans * y + *p;
        p += dir;
    }
    
    if (absx > 1) {
        i = N - M;
        return(R_pow(x, i) * num_ans / denom_ans);
    } else {
        return(num_ans / denom_ans);
    }
}




double lanczos_sum_expg_scaled(double x)
{
    return ratevl(x, lanczos_sum_expg_scaled_num,
                  sizeof(lanczos_sum_expg_scaled_num) / sizeof(lanczos_sum_expg_scaled_num[0]) - 1,
                  lanczos_sum_expg_scaled_denom,
                  sizeof(lanczos_sum_expg_scaled_denom) / sizeof(lanczos_sum_expg_scaled_denom[0]) - 1);
}

double log1p(double x)
{
    double z;
    
    z = 1.0 + x;
    if ((z < NPY_SQRT1_2) || (z > NPY_SQRT2))
        return (log(z));
    z = x * x;
    z = -0.5 * z + x * (z * polevl(x, LP, 6) / p1evl(x, LQ, 6));
    return (x + z);
}


/* log(1 + x) - x */
double log1pmx(double x)
{
    if (fabs(x) < 0.5) {
        int n;
        double xfac = x;
        double term;
        double res = 0;
        
        for(n = 2; n < MAXITER; n++) {
            xfac *= -x;
            term = xfac / n;
            res += term;
            if (fabs(term) < MACHEP * fabs(res)) {
                break;
            }
        }
        return res;
    }
    else {
        return log1p(x) - x;
    }
}



double expm1(double x)
{
    double r, xx;
    
    if (!isinf(x)) {
        if (isnan(x)) {
            return x;
        }
        else if (x > 0) {
            return x;
        }
        else {
            return -1.0;
        }
        
    }
    if ((x < -0.5) || (x > 0.5))
        return (exp(x) - 1.0);
    xx = x * x;
    r = x * polevl(xx, EP, 2);
    r = r / (polevl(xx, EQ, 3) - r);
    return (r + r);
}



double cosm1(double x)
{
    double xx;
    
    if ((x < -NPY_PI_4) || (x > NPY_PI_4))
        return (cos(x) - 1.0);
    xx = x * x;
    xx = -0.5 * xx + xx * xx * polevl(xx, coscof, 6);
    return xx;
}


/* Compute lgam(x + 1) around x = 0 using its Taylor series. */
double lgam1p_taylor(double x)
{
    int n;
    double xfac, coeff, res;
    
    if (x == 0) {
        return 0;
    }
    res = -NPY_EULER * x;
    xfac = -x;
    for (n = 2; n < 42; n++) {
        xfac *= -x;
        coeff = zeta(n, 1) * xfac / n;
        res += coeff;
        if (fabs(coeff) < MACHEP * fabs(res)) {
            break;
        }
    }
    
    return res;
}


/* Compute lgam(x + 1). */
double lgam1p(double x)
{
    if (fabs(x) <= 0.5) {
        return lgam1p_taylor(x);
    } else if (fabs(x - 1) < 0.5) {
        return log(x) + lgam1p_taylor(x - 1);
    } else {
        return lgam(x + 1);
    }
}




double zeta(x, q)
double x, q;
{
    int i;
    double a, b, k, s, t, w;
    
    if (x == 1.0)
        goto retinf;
    
    if (x < 1.0) {
    domerr:
        //sf_error("zeta", SF_ERROR_DOMAIN, NULL);
        //printf("zeta  SF_ERROR_DOMAIN\n");
        return (NAN);
    }
    
    if (q <= 0.0) {
        if (q == floor(q)) {
            //sf_error("zeta", SF_ERROR_SINGULAR, NULL);
            //printf("zeta  SF_ERROR_SINGULAR\n");
        retinf:
            return (NPY_INFINITY);
        }
        if (x != floor(x))
            goto domerr;    /* because q^-x not defined */
    }
    
    /* Asymptotic expansion
     * https://dlmf.nist.gov/25.11#E43
     */
    if (q > 1e8) {
        return (1/(x - 1) + 1/(2*q)) * R_pow(q, 1 - x);
    }
    
    /* Euler-Maclaurin summation formula */
    
    /* Permit negative q but continue sum until n+q > +9 .
     * This case should be handled by a reflection formula.
     * If q<0 and x is an integer, there is a relation to
     * the polyGamma function.
     */
    s = R_pow(q, -x);
    a = q;
    i = 0;
    b = 0.0;
    while ((i < 9) || (a <= 9.0)) {
        i += 1;
        a += 1.0;
        b = R_pow(a, -x);
        s += b;
        if (fabs(b / s) < MACHEP)
            goto done;
    }
    
    w = a;
    s += b * w / (x - 1.0);
    s -= 0.5 * b;
    a = 1.0;
    k = 0.0;
    for (i = 0; i < 12; i++) {
        a *= x + k;
        b /= w;
        t = a * b / AA[i];
        s = s + t;
        t = fabs(t / s);
        if (t < MACHEP)
            goto done;
        k += 1.0;
        a *= x + k;
        b /= w;
        k += 1.0;
    }
done:
    return (s);
}





//************************************* END igam.c*****************************************







/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/* Wendland covariance */

/* integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp)
{
    double res=0.0,y;
    y=lag/supp;
    res=R_pow(1-x,mu-1)*R_pow(x*x-y*y,alpha)/beta(2*alpha+1,mu);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}
void integr_gen(double *x, int n, void *ex){
    int i;double mu,alpha,beta,y;
    mu =    ((double*)ex)[0];  //mu
    alpha = ((double*)ex)[1];  //alpha
    beta =     ((double*)ex)[2];  //csupp
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_gen(x[i],mu,alpha,y,beta);}
    return;
}
// function computing generalized wendland
double wendintegral(double x, double *param) {
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;		     /* as instructed in WRE */
    iwork =   (int *) Calloc(subdiv, int);  /* idem */
    work = (double *) Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=x/param[2];
    upper=1;
    // Compute the integral
    if(x<=param[2]) {
    Rdqags(integr_gen, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);

    }else   {result=0;}
    Free(iwork);Free(work);
    return(result);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
// Integrand function (derivatives of the Student-t cdf):
double int_pt(double x, double df)
  {
    double res=0.0, x2=0.0, y=0.0;
    x2=R_pow(x,2);
    y=1+x2/df;
    res=0.5*dt(x,df,0)*((df+1)*x2/R_pow(df,2)/y-log(y));
    return res;
  }
// Vectorised integrand function:
void integr_pt(double *x, int n, void *ex)
{
  int i=0;
  double d=0.0;
  d=*((double *)ex);
  for(i=0;i<n;i++)
    x[i]=int_pt(x[i],d);
  return;
}




// compute  bivariate normal standard pdf:
double d2norm(double x, double y, double rho)
{
  double res=0.0, omr=1-R_pow(rho,2);
  res=(1/(2*M_PI))*exp(-0.5*(R_pow(x,2)-2*rho*x*y+R_pow(y,2))/omr)/sqrt(omr);
  return(res);
}


double d2lognorm(double x, double y, double sill,double nugget, double mux,double muy,double rho)
{
  rho=(1-nugget)*rho;
  double KK=exp(sill/2);
  x=x*KK; y=y*KK;
  double res=0.0, q=0.0, omr=R_pow(sill,2)-R_pow(rho*sill,2);
  q=(sill*R_pow((log(x)-mux),2) + sill*R_pow((log(y)-muy),2)
    -2*rho*sill*(log(x)-mux)*(log(y)-muy))/omr;
  res=exp(-q/2)/(2*x*y*M_PI*sqrt(omr));
  return(res*R_pow(KK,2));
}

double biv_sinh(double corr,double zi,double zj,double mi,double mj,double skew,double tail,double vari)
{
double b1=0.0,b2=0.0,A=0.0,B=0.0,k=0.0,res=0.0,Z1,Z2;
double xi=(zi-mi)/(sqrt(vari));
double xj=(zj-mj)/(sqrt(vari));
  b1=tail * asinh(xi)-skew;
  b2=tail * asinh(xj)-skew;
  k=1-R_pow(corr,2);
  A=R_pow(2 * M_PI * R_pow(k,0.5) * (vari),-1) * cosh(b1) * cosh(b2) * R_pow(tail,2)/sqrt((R_pow(xi,2)+1) * (R_pow(xj,2)+1));
  Z1=sinh(b1);Z2=sinh(b2);
  B=exp(- (Z1*Z1 + Z2*Z2 - 2*corr*Z1*Z2)/(2*k)  );
  res=A*B;                     
  return(res);
 } 


// pochammer factorial
double ff(double a,int k) 
{
return(gammafn(k+a)/gammafn(a));
}




double biv_chisqu2(double corr,double zi,double zj, double shape)
{
double KK1,KK2,KK3,rr;
   rr=1-corr*corr;  
   KK1=R_pow(gammafn(shape/2),2)*R_pow(rr,shape/2);
   KK2=R_pow(2,-shape)*R_pow(zi*zj,shape/2-1)*exp(-0.5*(zi+zj)/rr)/KK1;
   KK3= KK2* gammafn(shape/2)*R_pow(0.5*fabs(corr)*sqrt(zi*zj)/rr,1-shape/2)*bessel_i(fabs(corr)*sqrt(zi*zj)/rr ,shape/2-1,1);
   return(KK3);
}


/*********** bivariate bimodal********************/ 
double biv_two_piece_bimodal(double rho,double zi,double zj,double sill,double nuu,double eta,
             double p11,double mui,double muj)
{
double res;  
double nu=1/nuu;
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);


if(zi>=mui&&zj>=muj)
{res=          (zistd*zjstd*p11/R_pow(etamos,4))*biv_chisqu2(rho,R_pow(zistd/etamos,2),R_pow(zjstd/etamos,2),nu);  }
if(zi>=mui&&zj<muj)
{res=-zistd*zjstd*(1-eta-2*p11)/(R_pow(etamos,2)*R_pow(etamas,2))*biv_chisqu2(rho,R_pow(zistd/etamos,2),R_pow(zjstd/etamas,2),nu);}
if(zi<mui&&zj>=muj)
{res=-zistd*zjstd*(1-eta-2*p11)/(R_pow(etamos,2)*R_pow(etamas,2))*biv_chisqu2(rho,R_pow(zistd/etamas,2),R_pow(zjstd/etamos,2),nu);}
if(zi<mui&&zj<muj)
{res=     (zistd*zjstd*(p11+eta)/R_pow(etamas,4))*biv_chisqu2(rho,R_pow(zistd/etamas,2),R_pow(zjstd/etamas,2),nu);}
return(4*res/sill);
}







double  biv_Weibull2(double rho12,double zi,double zj,double mi,double mj, double shape1,double shape2)

{
double a1=0.0,a2=0.0,a=0.0,b=0.0,c=0.0,dens=.0,k=.0,mui=.0,muj=.0;
a1=gammafn(1+1/shape1);
a2=gammafn(1+1/shape2);
mui=exp(mi);muj=exp(mj);
k=1-rho12*rho12;
a=(shape1*shape2*R_pow(a1,shape1)*R_pow(a2,shape2)*R_pow(zi,shape1-1)*R_pow(zj,shape2-1))/(R_pow(mui,shape1)*R_pow(muj,shape2)*k);
b=exp(-(R_pow(a1*zi/mui,shape1) + R_pow(a2*zj/muj,shape2))/k);
c=bessel_i(2*fabs(rho12)*R_pow(a1*zi/mui,shape1/2)*R_pow(a2*zj/muj,shape2/2)/k,0,1);
dens=a*b*c;
return(dens);
}

double asy_log_besselI(double z,double nu)
{
     double val;
     double K=4*pow(nu,2);
     val =  (z-0.5*log(2*M_PI*z))+
           log((1-(K-1)/(8*z)*(1-(K-9)/(2*8*z)*(1-(K-25)/(3*8*z)))));
     return(val);
}


double biv_Weibull(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double res1=0.0,ui=0.0, uj=0.0,z=0.0,a=0.0,A=0.0,k=0.0,res=0.0,B=0.0;double ci=exp(mui);double cj=exp(muj);;
    k=pow(gammafn(1+1/shape),-1);
    ui=zi/ci;uj=zj/cj;
        a=1-R_pow(corr,2);
        z=2*fabs(corr)*pow(ui*uj,shape/2)*pow(k,-shape)/a;
        A=pow(shape,2)*pow(k,-2*shape)*pow(ui*uj,shape-1)/a;
        B= exp(-pow(k,-shape)*(pow(ui,shape)+pow(uj,shape))/a);
        if(z<700) 
               res=A*B*bessel_i(z,0,1)/(ci*cj);
        else{
               B=-pow(k,-shape)*(pow(ui,shape)+pow(uj,shape))/a;
               res1=(log(A)+B-(log(ci)+log(cj))) + asy_log_besselI(z,0);
               res=exp(res1);
             }
    return(res);
}



/*******************************************/
/*
double igam(double x, double a){
  double sum=0;
  double term=1.0/a;
  int n=1;
  while (term != 0){
    sum = sum + term;
    term = term*(x/(a+n));
    n++;
  }
  return (R_pow(x,a)*exp(-1*x)*sum);
}

double corr_pois1(double rho,double lambda1,double lambda2){
  double corr=0.0;
  double rho2=rho*rho;
        double aux1=lambda1/(1-rho2);
        double aux2=lambda2/(1-rho2);
        double aux3=0;
        int r=0;
        while (r <= 2000){
            aux3 = aux3 + exp(log(igam(aux1,r+1))+log(igam(aux2,r+1))-2*(lgammafn(r+1)));
            r++;
        }
        corr = rho2*(1-rho2)*aux3/sqrt(lambda1*lambda2);
    return(corr);
}
*/

/*
double e_n(int k,double x)
{ double sum=0.0;int i;
  for(i=0;i<=k;i++) sum=sum+R_pow(x,i)/gamma(i+1);
    return(sum);
}
double corr_pois(double rho,double mi,double mj)
{
int r=0; double res0=0.0,sum=0.0,pi,pj;
double rho2=rho*rho;
double ki=mi/(1-rho2);
double kj=mj/(1-rho2);
double K=rho2*(1-rho2)/sqrt(mi*mj);
while(r<4000){
  pi=exp( -ki+log(e_n(r,ki)));
  pj=exp( -kj+log(e_n(r,kj)));
  sum=sum+(1-pi)*(1-pj);
if((fabs(sum-res0)<1e-10)  ) {break;}
else {res0=sum;}
        r++;
    }
return(sum*K);
}*/

double corr_pois(double rho,double mi,double mj)
{
int r=0; double res0=0.0,sum=0.0;
double rho2=rho*rho;
double ki=mi/(1-rho2);
double kj=mj/(1-rho2);
double K=rho2*(1-rho2)/sqrt(mi*mj);
while(r<4000){
  sum=sum+ exp( log(igam(r+1,ki))+log(igam(r+1,kj)));
if((fabs(sum-res0)<1e-10)  ) {break;}
else {res0=sum;}
        r++;
    }
return(sum*K);
}


double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double res1=0.0,a=0.0,A=0.0,D=0.0,z=0.0,res=0.0,B=0.0,C=0.0;
    double ci=zi/exp(mui);double cj=zj/exp(muj);
    double gam = gammafn(shape/2);
        a=1-R_pow(corr,2);  

        z=shape*fabs(corr)*sqrt((ci)*(cj))/a;

        A=R_pow((ci)*(cj),shape/2-1) * R_pow(z/2,1-shape/2) ; ///ok
        C= exp(-shape*((ci)+(cj))/(2*a));//ok
        B=gam*R_pow(a,shape/2)*R_pow(2,shape)*R_pow(shape,-shape);  
        D=bessel_i(z,shape/2-1,1); //ok
        if(z<700) 
             res=(A*C*D)/(exp(muj)*exp(muj)*B);
        else{
            C=-shape*((ci)+(cj))/(2*a);
            res1=log(A)+C-(mui-muj-log(B))+asy_log_besselI(z,shape/2-1);
            res=exp(res1);
        }
        return(res);
}
/*******************************
double biv_gamma(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double a=0.0,A=0.0,D=0.0,z=0.0,res=0.0,B=0.0,C=0.0;
    double ci=zi/exp(mui);double cj=zj/exp(muj);
    double gam = gammafn(shape/2);
    //if(corr)   {
        a=1-R_pow(corr,2);  
        z=shape*fabs(corr)*sqrt((ci)*(cj))/a;
        A=R_pow((ci)*(cj),shape/2-1) * exp(-shape*((ci)+(cj))/(2*a)); ///ok
        C=R_pow(z/2,1-shape/2); //ok
        B=gam*R_pow(a,shape/2)*R_pow(2,shape)*R_pow(shape,-shape);  
        D=bessel_i(z,shape/2-1,1); //ok
        res=(A*C*D)/(exp(muj)*exp(muj)*B);
        return(res);
}*/


void biv_gamma_call(double *corr,double *zi,double *zj,double *mui, double *muj, double *shape, double *res)
{
    *res = biv_gamma(*corr,*zi,*zj,*mui,*muj,*shape);
}

/***********/
double Hyp_conf_laplace_approx(double a, double b, double z)
{
double y=0.0,r=0.0,A=0.0;
  y=2*a/(b-z+sqrt(R_pow(z-b,2)+ 4*a*z ));

  r=(R_pow(y,2)/a)+(R_pow(1-y,2)/(b-a));

  A=R_pow(b,b-0.5)*R_pow(r,-0.5)*R_pow((y/a),a)*R_pow((1-y)/(b-a),b-a)*exp(z*y);
  return(A);
}
/*****************************************************/
double asym_aprox_1F1(double a, double b, double z)
{
double k1,k2,k3,res;
double alpha=(b-a-1)/z;
if(b<z+a+1)
{
  k1=(gammafn(b)*exp(z)*R_pow(1-alpha,a-1))/(gammafn(a)*gammafn(b-a));
  k2=sqrt(2*M_PI*alpha/z)*R_pow(alpha/exp(1),alpha*z);
  k3=(2-a*alpha)*(1-a)/(2*R_pow(1-alpha,2));
  res=k1*k2*(1+ (k3+ 1/(12*alpha))/z);  
}
if(b>z+a+1)
{
k1=gammafn(b)/(gammafn(b-a)*R_pow(b-a-z-1,a));
k2=a*(a+1)*(a+1-b)/(2*R_pow(b-a-z-1,2));
res=k1*(1-k2);
}
if(b==z+a+1)
{k1=gammafn(b)/(2*gammafn(a)*gammafn(b-a));
 k2=R_pow(2/(b-a-1),a/2);
 k3=2*gammafn(0.5*(a+3))*sqrt(2/(b-a-1))/3;
res=k1*k2*(gammafn(a/2)-k3);}
return(res);
}



double biv_gamma2(double rho,double x,double y,double shape1, double shape2, double rate)
{
double  ss,yr,xr,xi;  
double  ar,ai=0.0;
double  br,bi=0.0;
double  r1,r2,ri=0.0;
int len=1;int lnchf=0;int ip=0;

double res,beta,b,c,f,rho2,a;
int k;

rho2=R_pow(rho,2);

beta=shape1+shape2;

if(rho){
c=rate/(2*(1-rho2));
   
ar=shape2/2;
a=0.0;
for(k=0;k<=20;k++)
{
br=beta/2+k;
xr=c*rho2*x;
yr=c*rho2*y;
//r1=Hyp_conf_laplace_approx(ar,br,c*rho2*x);
//r2=Hyp_conf_laplace_approx(ar,br,c*rho2*y);
F77_NAME(chfm)(&xr,&xi,&ar,&ai,&br,&bi,&r1,&ri,&len,&lnchf,&ip);
F77_NAME(chfm)(&yr,&xi,&ar,&ai,&br,&bi,&r2,&ri,&len,&lnchf,&ip);
if(!R_FINITE(r1*r2)||ISNAN(r1*r2)){
 r1=asym_aprox_1F1(ar,br,c*rho2*x);
 r2=asym_aprox_1F1(ar,br,c*rho2*y);}

ss=ff(shape1/2,k)*r1*r2*R_pow((R_pow(c,2)*rho2*x*y),k)/(gammafn(k+1)*R_pow(ff(beta/2,k),2));
a=a+ss;
}
b=R_pow(rate/2,beta)*R_pow(x*y,beta/2-1)*exp(-c*(x+y));
f=R_pow(gammafn(beta/2),2)*R_pow(1-rho2,shape1/2);

res=b*a/f;}
else
{
  b=(R_pow(rate,beta/2)*R_pow(x,beta/2-1)*exp(-(rate*x/2))) /(R_pow(2,beta/2)*gammafn(beta/2));
  f=(R_pow(rate,beta/2)*R_pow(y,beta/2-1)*exp(-(rate*y/2)) )/(R_pow(2,beta/2)*gammafn(beta/2));
  res=b*f;
}
  return(res);
}


//bivariate skew gaussian distribution
double biv_skew(double corr,double z1,double z2,double mi,double mj,double vari,double skew)
{
   double aux1=0.0,aux11=0.0, aux2=0.0,aux21 = 0.0, aux22=0.0, pdf1=0,pdf2=0,cdf1=0,cdf2=0,quadr,zi,zj;
    zi=z1-mi;
    zj=z2-mj;
    double det,dens,det1,det2,lim1,lim2,a11,a22;
    double nu2  = R_pow(skew,2);
    double om2  = R_pow(vari,1);
    double cor = corr;
                                      // pdf 1
                                       aux1  =  om2 + nu2 ;
                                       aux2  =  cor * aux1;
                                       det1  =  R_pow(aux1,2) - R_pow(aux2,2) ; 
                                       det2  =  aux1 - aux1 * R_pow(cor,2)  ;
                                       quadr =  (1.0/det2)*( (R_pow(zj,2.0)+R_pow(zi,2.0)) - 2.0*cor*zi*zj  ) ;  
                                       pdf1  =  (0.5/M_PI) * (1/sqrt(det1)) * exp(- 0.5 * quadr  ) ;
                                       lim1  = (zi * skew)/(aux1);
                                       lim2  = (zj * skew)/(aux1);
                                       a11   = (om2)/(aux1);
                                       a22   = (cor * om2)/(aux1);
                                       cdf1  =  cdf_norm(lim1,lim2,a11,a22) ;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                       // pdf 2
                                       aux1  =  om2 + nu2 ;
                                       aux2  =  om2 * cor - nu2 * cor ;
                                       det   =  R_pow(aux1,2) - R_pow(aux2,2) ;   
                                       quadr =  (1.0/det)*( (aux1)*(R_pow(zj,2.0)+R_pow(zi,2.0)) - 2.0*aux2*zi*zj  ) ;  
                                       pdf2  =  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
                                       aux1  =  om2 + nu2 ;
                                       aux2  =  om2 * cor - nu2 * cor ;
                                       aux11 = (skew)/( R_pow(aux1,2.0) - R_pow(aux2,2));
                                       aux21 = nu2 * (1 - R_pow(cor,2.0));
                                       aux22 = om2 * (1 + R_pow(cor,2.0));
                                      
                                       lim1  = aux11 * (zi*aux21 + zi*aux22 - 2*zj*om2*cor);
                                       lim2  = aux11 * (zj*aux21 + zj*aux22 - 2*zi*om2*cor); // bien hasta aquí
                                       a11   =  1 - skew*aux11 * (nu2 + om2 - nu2 * R_pow(cor,2.0) + 3*om2 * R_pow(cor,2.0));
                                       a22   =  - cor + (skew) * aux11 * ( nu2  *(cor - R_pow(cor,3.0)) + om2*(cor + R_pow(cor,3.0)) + 2 * cor * om2  ); 
                                       cdf2  =  cdf_norm(lim1,lim2,a11,a22) ; 
  
dens = 2*(pdf1 * cdf1 + pdf2 * cdf2);
return(dens);
}


double biv_wrapped (double alfa,double u, double v,double mi,double mj,double nugget,double sill,double corr)
{
    double x,y,s1=0.0,s12=0.0,quadr=0.0,det=0.0,wrap_gauss=0.0; // 5???
    //2*atan(mean[i])-M_PI
    x=u-2*atan(mi)-M_PI; y=v-2*atan(mj)-M_PI;
    s1=nugget+sill;
    s12=sill*corr; //sill * corr
    det=R_pow(s1,2)-R_pow(s12,2);
    double k1=-alfa,k2=-alfa;
     while(k1<=alfa){
     while(k2<=alfa){
     quadr = -0.5*(1.0/det)*(s1*R_pow((y+2*k1*M_PI),2.0)+s1*R_pow((x+2*k2*M_PI),2.0)
                                                    -2.0*s12*(x+2*k2*M_PI)*(y+2*k1*M_PI));
     wrap_gauss +=  (1/2.0*M_PI)*(1/sqrt(det)*exp(quadr)) ;
     k2 = k2+1;}
     k1 = k1+1;k2 = -alfa;}
return(wrap_gauss);
}



//bivariate skew gaussian bivariate  distribution
double biv_skew2(double corr,double zi,double zj,double vari1,double vari2,double rho ,double skew1,double skew2)
{
   
    double aux1,aux2,aux3,pdf1,pdf2,quadr;
    double dens,det,fact1,fact2, fact3,fact4,factor,c1,lim1,lim2,cdf1,cdf2,a11,a12,a22;
       double taul=R_pow(vari1,0.5); double taur=R_pow(vari2,0.5);
                                       c1= 1.0 /  (1-R_pow(corr,2)) ;
                                       aux1 = R_pow(taul,2) + R_pow(skew1,2) ;
                                       aux2 =  taul*taur*corr + skew1*skew2*corr ;
                                       aux3 = R_pow(taur,2) + R_pow(skew2,2) ;
                                       det =   aux1*aux3 - R_pow(aux2,2) ;   

                                       quadr= (1.0/det)*( aux1*R_pow(zj,2.0)+aux3*R_pow(zi,2.0) - 2.0*aux2*zi*zj  ) ;  
                                       pdf1=  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
                                       //////
                                       aux1 = (R_pow(skew1/taul,2)+1)*c1;
                                       aux2 = ((skew1*skew2)/(taul*taur)+1)*c1*corr;
                                       aux3 = (R_pow(skew2/taur,2)+1)*c1;

                                       det =   aux1*aux3 - R_pow(aux2,2) ;  
                                       a11  =  (1/det)* aux3   ;
                                       a12  =  (1/det)* aux2   ;
                                       a22  =  (1/det)* aux1   ;
                                       factor = c1 /  det ;
                                       fact1 = factor * (aux3*(skew1/R_pow(taul,2)) - aux2*(skew2/(taul*taur))*corr) ;
                                       fact2 = factor * (-aux3*(skew1/(taul*taur))*corr + aux2*(skew2/R_pow(taur,2)) ) ;
                                       fact3 = factor * (aux2*(skew1/R_pow(taul,2)) - aux1*(skew2/(taul*taur))*corr) ;
                                       fact4 = factor * (-aux2*(skew1/(taul*taur))*corr + aux1*(skew2/R_pow(taur,2)) ) ;

                                       lim1 =   fact1*zi +  fact2*zj   ;
                                       lim2 =   fact3*zi +  fact4*zj   ;
                                       cdf1 = cdf_norm2(lim1,lim2,a11,a12,a22) ;
                                       //////////////////////////////////////////////////////////////////////////////////
                                       aux1 = R_pow(taul,2) + R_pow(skew1,2) ;
                                       aux2 =  taul*taur*corr - skew1*skew2*corr ;
                                       aux3 = R_pow(taur,2) + R_pow(skew2,2) ;
                                       det =   aux1*aux3 - R_pow(aux2,2) ;   
                                       quadr= (1.0/det)*( aux1*R_pow(zj,2.0)+aux3*R_pow(zi,2.0) - 2.0*aux2*zi*zj  ) ;  
                                       pdf2=  (0.5/M_PI) * (1/sqrt(det)) * exp(- 0.5 * quadr) ;
                                       //////
                                       aux1 = R_pow(skew1/taul,2)*c1 + c1 ;
                                       aux2 = (skew1*skew2)/(taul*taur)*c1*corr - c1*corr;
                                       aux3 = R_pow(skew2/taur,2)*c1 + c1 ;
                                       det =   aux1*aux3 - R_pow(aux2,2) ;  
                                       a11  =  (1/det)* aux3   ;
                                       a12  =  (1/det)* aux2   ;
                                       a22  =  (1/det)* aux1   ;
                                       factor = c1 /  det ;
                                       fact1 = factor * (aux3*(skew1/R_pow(taul,2)) - aux2*(skew2/(taul*taur))*corr) ;
                                       fact2 = factor * (-aux3*(skew1/(taul*taur))*corr + aux2*(skew2/R_pow(taur,2)) ) ;
                                       fact3 = factor * (aux2*(skew1/R_pow(taul,2)) - aux1*(skew2/(taul*taur))*corr) ;
                                       fact4 = factor * (-aux2*(skew1/(taul*taur))*corr + aux1*(skew2/R_pow(taur,2) )) ;
                                       lim1 =   fact1*zi +  fact2*zj   ;
                                       lim2 =   fact3*zi +  fact4*zj   ;
                                       cdf2 = cdf_norm2(lim1,lim2,a11,a12,a22) ; 
                            dens=2*(pdf1*cdf1+pdf2*cdf2);

return(dens);
}


// trivariate skew gaussian distribution/
double triv_skew(double x,double c_0i,double c_0j, double rho,double data_i,double data_j,double *nuis)
{

    int N=3,m,n,i, maxpts=3*2000,fail=100;
    double esterror=10.0,abseps=1.0e-6, releps=0.0;
    double *K;K=(double *) Calloc(N,double);
    double *lower;lower=(double *) Calloc(N,double);
    double *upper;upper=(double *) Calloc(N,double);
    double *corr;corr=(double *) Calloc(N,double);
    int *infin;infin=(int *) Calloc(N,int);
    lower[0]=0;lower[1]=0;lower[2]=0;
    infin[0]=0;infin[1]=0;infin[2]=0;

    double **Om;
    Om= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){Om[i]=(double *) Calloc(N,double);}

    double **C;
    C= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){C[i]=(double *) Calloc(N,double);}

    double **M;
    M= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){M[i]=(double *) Calloc(N,double);}

    double **cov;
    cov= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){cov[i]=(double *) Calloc(N,double);}
    
    double **inverse;
    inverse= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){inverse[i]=(double *) Calloc(N,double);}  

    double **inverse1;
    inverse1= (double **) Calloc(N,double *);
    for(i=0;i<N;i++){inverse1[i]=(double *) Calloc(N,double);} 

    int *indx; indx=(int *) Calloc(N,int);
    double *col; col=(double *) Calloc(N,double);
    double *data;  data=(double *) Calloc(N,double);

    double det,pdf1,pdf2,pdf3,pdf4,dd;
    double dens,cdf1,cdf2,cdf3,cdf4;
    double nu2=R_pow(nuis[3],2.0),tau2 =nuis[2];
    double kk=tau2+nu2,co=tau2/nu2;
    double co1=co+1,com1=co-1;
    double ratio=nuis[3]/tau2;
  

  /**********    omega and its inverse*/
    Om[0][0]=1;   Om[0][1]=rho ;Om[0][2]=c_0i;
    Om[1][0]=Om[0][1];   Om[1][1]=1;     Om[1][2]=c_0j;
    Om[2][0]=Om[0][2];   Om[2][1]=Om[1][2];  Om[2][2]=1;
     ludcmp(Om,N,indx,&dd);    //Lu decomposition
     //for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
     for(n=0;n<N;n++) {
     for(m=0;m<N;m++) col[m]=0.0;
     col[n]=1.0;
     lubksb(Om,N,indx,col);
    for(m=0;m<N;m++) inverse1[m][n]=col[m];  // inverse using LU decomposition
      }
     /**********/

    data[0]=data_i-nuis[0] ;data[1]=data_j-nuis[0] ;data[2]=x-nuis[0] ;
    det=co*(R_pow(c_0j,2)-2*c_0i*rho*c_0j+R_pow(rho,2)+R_pow(c_0i,2)-1);
/********************************************************************************************/
/********************************************************************************************/
                    M[0][0]=kk     ;   M[0][1]=kk*rho ;M[0][2]=kk*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=kk*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf1=  dNnorm(N,M,data);
                    cov[0][0]= (co1)*(c_0j-1)*(c_0j+1)/det;   
                    cov[0][1]=-(co1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=-(co1)*(rho*c_0j-c_0i)/det;
                    cov[1][1]= (co1)*(c_0i-1)*(c_0i+1)/det;
                    cov[1][2]=(co1)*(c_0j-c_0i*rho)/det;
                    cov[2][2]=(co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  

                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                   // for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }

                    Matrix_prod(inverse,inverse1,C,N);
                    Matrix_prod_vec(C,data,K,N);

                    upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                    upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                    upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                    corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];
                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf1, &fail );
/********************************************************************************************/
/********************************************************************************************/
                    M[0][0]=kk     ;   M[0][1]=(tau2-nu2)*rho ;M[0][2]=(tau2-nu2)*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=kk*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf2= dNnorm(N,M,data);
                    cov[0][0]=(co1)*(c_0j-1)*(c_0j+1)/det;
                    cov[0][1]=(com1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=(com1)*(rho*c_0j-c_0i)/det;
                    cov[1][2]=(co1)*(c_0j-c_0i*rho)/det;
                    cov[1][1]=(co1)*(c_0i-1)*(c_0i+1)/det; 
                    cov[2][2]= (co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  
                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                    //for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }
                     Matrix_prod(inverse,inverse1,C,N);
                     Matrix_prod_vec(C,data,K,N);
                   
                   upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                   upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                   upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                   corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];

                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf2, &fail );
/********************************************************************************************/
/********************************************************************************************/
                    M[0][0]=kk     ;   M[0][1]=(tau2-nu2)*rho ;M[0][2]=kk*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=(tau2-nu2)*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf3= dNnorm(N,M,data);

                    cov[0][0]=(co1)*(c_0j-1)*(c_0j+1)/det;   
                    cov[0][1]=(com1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=-(co1)*(rho*c_0j-c_0i)/det;
                    cov[1][1]=(co1)*(c_0i-1)*(c_0i+1)/det;
                    cov[1][2]=-(com1)*(c_0j-c_0i*rho)/det;
                    cov[2][2]=(co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  

                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                   // for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }
                     Matrix_prod(inverse,inverse1,C,N);
                     Matrix_prod_vec(C,data,K,N);
                   upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                   upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                   upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                   corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];

                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf3, &fail );
/********************************************************************************************/
/********************************************************************************************/

                    M[0][0]=kk;   M[0][1]=kk*rho; M[0][2]=(tau2-nu2)*c_0i;
                    M[1][0]=M[0][1];   M[1][1]=kk;     M[1][2]=(tau2-nu2)*c_0j;
                    M[2][0]=M[0][2];   M[2][1]=M[1][2];  M[2][2]=kk;
                    pdf4= dNnorm(N,M,data);

                    cov[0][0]=(co1)*(c_0j-1)*(c_0j+1)/det;   
                    cov[0][1]=-(co1)*(c_0i*c_0j-rho)/det;
                    cov[0][2]=(com1)*(rho*c_0j-c_0i)/det ;
                    cov[1][1]=(co1)*(c_0i-1)*(c_0i+1)/det; 
                    cov[1][2]=-(com1)*(c_0j-c_0i*rho)/det;
                    cov[2][2]=(co1)*(rho-1)*(rho+1)/det;
                    cov[1][0]=cov[0][1];cov[2][0]=cov[0][2];cov[2][1]=cov[1][2];  
                    ludcmp(cov,N,indx,&dd);    //Lu decomposition
                    //for(m=0;m<N;m++){ dd *= cov[m][m];}  //determinant using LU
                    for(n=0;n<N;n++) {
                    for(m=0;m<N;m++) col[m]=0.0;
                    col[n]=1.0;
                    lubksb(cov,N,indx,col);
                    for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
                    }
                    
                     Matrix_prod(inverse,inverse1,C,N);
                     Matrix_prod_vec(C,data,K,N);
                 upper[0]=K[0]*ratio/sqrt(inverse[0][0]);
                 upper[1]=K[1]*ratio/sqrt(inverse[0][0]);
                 upper[2]=K[2]*ratio/sqrt(inverse[0][0]);
                 corr[0]=inverse[0][1]/inverse[0][0]; corr[1]=inverse[0][2]/inverse[0][0]; corr[2]=inverse[1][2]/inverse[0][0];
                    mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &cdf4, &fail );
/********************************************************************************************/
/********************************************************************************************/
dens=2*(pdf1*cdf1 + pdf2*cdf2 + pdf3*cdf3 + pdf4*cdf4);
 Free(lower);Free(upper);Free(corr);Free(infin);Free(K);
 Free(indx);Free(col);
 for(i=0;i<N;i++)  {Free (inverse[i]);Free(inverse1[i]);Free(Om[i]);Free (M[i]);Free (cov[i]);Free (C[i]);}
Free(inverse);Free(M);Free(cov);Free(C);Free(Om);Free(inverse1);
 return(dens);
}


double log_biv_Norm(double corr,double zi,double zj,double mi,double mj,double vari, double nugget)
{
    double u,v,u2,v2,det,s1,s12,dens;
    u=zi-mi;
    v=zj-mj;
    u2=R_pow(u,2);v2=R_pow(v,2);
    s1=vari+nugget;
    s12=vari*corr;
    det=R_pow(s1,2)-R_pow(s12,2);
    dens=(-0.5*(2*log(2*M_PI)+log(det)+(s1*(u2+v2)-2*s12*u*v)/det));
return(dens);
}

/*multivariate gaussian PDF*/
double dNnorm(int N,double **M, double *dat)
{
 double dd,pdf; int m=0,n=0,i=0;
 int *indx; indx=(int *) Calloc(N,int);
 double *col; col=(double *) Calloc(N,double);
 double **inverse;
 inverse= (double **) Calloc(N,double *);
 for(i=0;i<N;i++){inverse[i]=(double *) Calloc(N,double);} 

 ludcmp(M,N,indx,&dd);    //Lu decomposition
 for(m=0;m<N;m++){ dd *= M[m][m];}  //determinant using LU
 for(n=0;n<N;n++) {
 for(m=0;m<N;m++) col[m]=0.0;
 col[n]=1.0;
 lubksb(M,N,indx,col);
 for(m=0;m<N;m++) inverse[m][n]=col[m];  // inverse using LU decomposition
 }
 pdf=R_pow(2*M_PI,-N/2)*R_pow(dd,-0.5)*exp(-0.5*QFORM(inverse,dat,dat,N)) ;   
 Free(indx);Free(col);
 for(i=0;i<N;i++)  {Free (inverse[i]);}Free(inverse);
return(pdf);
}

void mult_pmnorm( int *nvar , double *lower , double *upper , int *infin , double *corr , int *maxpts , double *abseps , double *releps , double *esterror , double *result , int *fail )
{ 
  F77_CALL(sadmvn)( nvar, lower, upper, infin, corr, maxpts, abseps, releps, esterror, result, fail ) ;
} 

// compute the trivariate normal cdf  for bernoulli  RF :
double ptnorm(int which,int *cormod, double h0,double h1,double h2, double u0, double u1,double u2, 
                 double *nuis, double *par, double thr)
{  
  int N=3; int maxpts=3*2000;
  double esterror=10.0,abseps=1.0e-6, releps=0.0,res2=0.0;int fail=100;
  double a=(nuis[0]-thr)/(sqrt(nuis[2]+nuis[1]));
  double *lower;double *upper; int *infin;
  lower=(double *) Calloc(3,double);
  upper=(double *) Calloc(3,double);
  infin=(int *) Calloc(3,int);
  switch(which) //
              {
         case 1: /// 111
         lower[0]=0;lower[1]=0; lower[2]=0; 
         upper[0]=a;upper[1]=a; upper[2]=a;
         infin[0]=0;infin[1]=0;infin[2]=0;
            break;
            case 2: /// 110
         lower[0]=0;lower[1]=0; lower[2]=a; 
         upper[0]=a;upper[1]=a; upper[2]=0 ;
         infin[0]=0;infin[1]=0;infin[2]=1;     
            break;
            case 3: /// 101
         lower[0]=0;lower[1]=a; lower[2]=0; 
         upper[0]=a;upper[1]=0 ; upper[2]=a;
         infin[0]=0;infin[1]=1;infin[2]=0;      
            break;
            case 4: /// 011
         lower[0]=a;lower[1]=0; lower[2]=0; 
         upper[0]=0;upper[1]=a; upper[2]=a;
         infin[0]=1;infin[1]=0;infin[2]=0;        
            break;
 }
  double corr[3]={nuis[2]*CorFct(cormod,h0,u0,par,0,0),
                  nuis[2]*CorFct(cormod,h1,u1,par,0,0),
                  nuis[2]*CorFct(cormod,h2,u2,par,0,0)};


  mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &res2, &fail );
  //Free(lim_inf);Free(lim_sup);Free(infin);
  Free(lower);Free(upper);Free(infin);
  return(res2);
}


// trivariate CDF gaussian 
double p3norm(int *cormod, double h0,double h1,double h2, double u0, double u1,double u2, 
                 double *nuis, double *par)
{  
  int N=3; int maxpts=3*2000;
  double esterror=10.0,abseps=1.0e-6, releps=0.0,res2=0.0;int fail=100;
  double a=1;
  double *lower;double *upper; int *infin;
  lower=(double *) Calloc(3,double);
  upper=(double *) Calloc(3,double);
  infin=(int *) Calloc(3,int);

 lower[0]=0;lower[1]=0; lower[2]=0; 
 upper[0]=a;upper[1]=a; upper[2]=a;
 infin[0]=0;infin[1]=0;infin[2]=0;

 double corr[3]={nuis[2]*CorFct(cormod,h0,u0,par,0,0),
                  nuis[2]*CorFct(cormod,h1,u1,par,0,0),
                  nuis[2]*CorFct(cormod,h2,u2,par,0,0)};

  mult_pmnorm( &N, lower, upper, infin, corr, &maxpts, &abseps, &releps, &esterror, &res2, &fail );
  Free(lower);Free(upper);Free(infin);
  return(res2);
}



// cdf of  a bivariate Gausssian distribution
double cdf_norm(double lim1,double lim2,double a11,double a12)
{
    double res=0;
    double  lowe[2]={0,0}, uppe[2]={lim1/sqrt(a11),lim2/sqrt(a11)}, corre[1] ={a12/a11};
    int infin[2]={0,0};
    double auxil=1-R_pow(corre[0],2);
    double det=R_pow(a11,2)-R_pow(a12,2) ;
    res= a11* sqrt(auxil/det) *  F77_CALL(bvnmvn)(lowe,uppe,infin,corre);
    return(res);
}
// cdf of  a bivariate Gausssian distribution
double cdf_norm2(double lim1,double lim2,double a11,double a12, double a22)
{
  double res=0;

  double  lowe[2]={0,0}, uppe[2]={lim1/sqrt(a11),lim2/sqrt(a22)}, corre[1] ={a12/sqrt(a11*a22)};
  int infin[2]={0,0};
  double auxil=1-R_pow(corre[0],2);
  double det=a11*a22-R_pow(a12,2) ;
  res= sqrt(a11*a22)* sqrt(auxil/det) *  F77_CALL(bvnmvn)(lowe,uppe,infin,corre);
  return(res);
}


// compute the bivariate normal cdf for the bernoulli RF:
// compute the bivariate normal cdf for the bernoulli RF:
double pbnorm(int *cormod, double h, double u, double mean1, double mean2, 
  double nugget, double var,double *par, double thr)
{
  double res=0;
  double lim_inf[2]={0,0};//lower bound for the integration
  double lim_sup[2]={mean1,mean2};
  int infin[2]={0,0};//set the bounds for the integration
  double corr[1]={(1-nugget)*CorFct(cormod,h,u,par,0,0)};
    res=F77_CALL(bvnmvn)(lim_inf,lim_sup,infin,corr);
    //res = Phi2(lim_sup[0],lim_sup[1],corr[0]);
  return(res);
}



// compute a sequence (depending on space and time) of bivariate normal cdf:
void vpbnorm(int *cormod, double *h, double *u, int *nlags, int *nlagt,
	     double *nuis, double *par, double *rho, double *thr)
{
  int i,j,t=0;
  for(j=0;j<*nlagt;j++)
    for(i=0;i<*nlags;i++){
      rho[t]=pbnorm(cormod,h[i],u[j],nuis[0],nuis[0],nuis[1],nuis[2],par,thr[0]);
      t++;}
  return;
}

/*

double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11)
{
  int a=0,i=0;
  double kk1=0.0,kk2=0.0,dens1=0.0,dens2=0.0;

    for(a=fmax_int(0,u-v+NN-1);a<=NN-2;a++){
   for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=fac(NN-1+u,1)/(fac(i,1)*fac(NN-1-i,1)*fac(a-i,1)*fac(u-a+i,1));
    kk2=fac(v-u-1,1)/(fac(v+a-NN-u+1,1)*fac(NN-a-2,1));
    dens1+=kk1*kk2*R_pow(p11,i+1)*R_pow(1+p11-(x+y),u-a+i)*
             R_pow(x-p11,NN-i-1)*R_pow(y-p11,a-i)*R_pow(1-y,v-u-NN+a+1)*R_pow(y,NN-a-1);
  }}

    for(a=fmax_int(0,u-v+NN);a<=NN-1;a++){
    for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=fac(NN-1+u,1)/(fac(i,1)*fac(NN-1-i,1)*fac(a-i,1)*fac(u-a+i,1));
    kk2=fac(v-u-1,1)/(fac(v+a-NN-u,1)*fac(NN-a-1,1));
    dens2+=kk1*kk2*R_pow(p11,i)*R_pow(1+p11-(x+y),u-a+i)*
                 R_pow(x-p11,NN-i)*R_pow(y-p11,a-i)*R_pow(1-y,v-u-NN+a)*R_pow(y,NN-a);
  }}
  return(dens1+dens2);
}

// bivariate negative  binomial
double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
double kk1=0.0,dens=0.0;int i=0;
if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);

if(u==v)            {
    for(i=fmax_int(0,NN-u-1);i<=NN-1;i++){
      kk1=fac(NN-1+u,1)/(fac(i,1)*R_pow(fac(NN-1-i,1),2)*fac(u-NN+1+i,1));
      dens+=kk1*R_pow(p11,i+1)*R_pow(1+p11-(p01+p10),u-NN+1+i)*R_pow(p01-p11,NN-1-i)*R_pow(p10-p11,NN-1-i); }
}

if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
return(dens);
}
*/


double aux_biv_binomneg_simple(int NN, int u, double p01,double p10,double p11)
{
          int i=0;
  double kk1=0.0,dens=0.0;

    for(i=fmax_int(0,NN-u-1);i<=NN-1;i++){
      kk1=exp(lgammafn(NN-1+u+1)-(lgammafn(i+1)+lgammafn(NN-i)+lgammafn(NN-i)+lgammafn(u-NN+2+i)));
      dens+=kk1*pow(p11,i+1)*pow(1+p11-(p01+p10),u-NN+1+i)*pow(p01-p11,NN-1-i)*pow(p10-p11,NN-1-i);
       }
         return(dens);
}



double aux_biv_binomneg (int NN, int u, int v, double x,double y,double p11)
{
  int a=0,i=0;
  double kk1=0.0,kk2=0.0,dens1=0.0,dens2=0.0;

    for(a=fmax_int(0,u-v+NN-1);a<=NN-2;a++){
   for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=exp(lgammafn(NN+u)-(lgammafn(i+1)+lgammafn(NN-i)+lgammafn(a-i+1)+lgammafn(u-a+i+1)));         
    kk2=exp(lgammafn(v-u)-(lgammafn(v+a-NN-u+2)+lgammafn(NN-a-1)));                        
    dens1+=kk1*kk2*pow(p11,i+1)*pow(1+p11-(x+y),u-a+i)*
             pow(x-p11,NN-i-1)*pow(y-p11,a-i)*pow(1-y,v-u-NN+a+1)*pow(y,NN-a-1);
  }}

    for(a=fmax_int(0,u-v+NN);a<=NN-1;a++){
    for(i=fmax_int(0,a-u);i<=fmin_int(a,NN-1);i++){
    kk1=exp(lgammafn(NN+u)-(lgammafn(i+1)+lgammafn(NN-1-i+1)+lgammafn(a-i+1)+lgammafn(u-a+i+1)));
    kk2=exp(lgammafn(v-u)-(lgammafn(v+a-NN-u+1)+lgammafn(NN-a)));                      
    dens2+=kk1*kk2*pow(p11,i)*pow(1+p11-(x+y),u-a+i)*
                 pow(x-p11,NN-i)*pow(y-p11,a-i)*pow(1-y,v-u-NN+a)*pow(y,NN-a);
  }}
  return(dens1+dens2);
}

/*
double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
double kk1=0.0,dens=0.0;int i=0;
if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);

if(u==v)            {
    for(i=fmax_int(0,NN-u-1);i<=NN-1;i++){
      kk1=exp(lgammafn(NN-1+u+1)-(lgammafn(i+1)+lgammafn(NN-i)+lgammafn(NN-i)+lgammafn(u-NN+2+i)));
      dens+=kk1*pow(p11,i+1)*pow(1+p11-(p01+p10),u-NN+1+i)*pow(p01-p11,NN-1-i)*pow(p10-p11,NN-1-i); }
}

if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
return(dens);
}*/


double biv_binomneg (int NN, int u, int v, double p01,double p10,double p11)
{
double dens=0.0;
if(u<v)    dens=aux_biv_binomneg(NN,u,v,p01,p10,p11);

if(u==v)      dens=aux_biv_binomneg_simple(NN,v,p10,p01,p11);

if(u>v)    dens=aux_biv_binomneg(NN,v,u,p10,p01,p11);
return(dens);
}

double bin_aux(int a,int NN,int u,int v,double p1, double p2,double p11)
{
  double kk,dens;
kk=exp(lgammafn(NN+1)-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)+lgammafn(NN-u-v+a+1)));
dens=kk*(R_pow(p11,a)*R_pow(p1-p11,u-a)*R_pow(p2-p11,v-a)*R_pow(1+p11-(p1+p2),NN-u-v+a));
return(dens);
}

double biv_binom(int NN, int u, int v, double p01,double p10,double p11)
{
    
int a;
double kk=0.0,dens=0.0;
for(a=fmax_int(0,u+v-NN);a<=fmin_int(u,v);a++){
//kk=exp(logfac(NN)-(logfac(a)+logfac(u-a)+logfac(v-a)+logfac(NN-u-v+a)));
kk=exp(lgammafn(NN+1)-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)+lgammafn(NN-u-v+a+1)));
dens+=kk*(R_pow(p11,a)*R_pow(p01-p11,u-a)*R_pow(p10-p11,v-a)*R_pow(1+p11-(p01+p10),NN-u-v+a));
 }
    return(dens);
}


/*
double biv_binom (int NN, int u, int v, double p01,double p10,double p11)
{
    
int a;
double kk=0,dens=0.0;
for(a=fmax_int(0,u+v-NN);a<=fmin_int(u,v);a++){
kk=fac(NN,1)/(fac(a,1)*fac(u-a,1)*fac(v-a,1)*fac(NN-u-v+a,1));
dens+=kk*(R_pow(p11,a)*R_pow(p01-p11,u-a)*R_pow(p10-p11,v-a)*R_pow(1+p11-(p01+p10),NN-u-v+a));
 }
    return(dens);
}*/



/// biv binomial type II
/*
double  biv_binom2 (int NN_i,int NN_j, int k, int u, int v, double p01,double p10,double p11)
{

int a,i,j;
double const1=0.0,const2=0.0,const3=0.0,dens=0.0,dens1;
double P10=p01-p11;
double P01=p10-p11;
double P00=1+p11-(p01+p10);
double P11=p11;

for(i=0;i<=fmin_int(NN_i-k,u);i++){
   for(j=0;j<=fmin_int(NN_j-k,v);j++){
     
for(a=fmax_int(0,u+v-k-i-j);a<=fmin_int(u-i,v-j);a++){         
       const1=fac(k,1)/(fac(a,1)*fac(u-i-a,1)*fac(v-j-a,1)*fac(k-u-v+i+j+a,1));
        dens1=const1*R_pow(P11,a)*R_pow(P00,k-u-v+i+j+a)*
             R_pow(P10,u-i-a)*R_pow(P01,v-j-a);
   
       const2=fac(NN_i-k,1)/(fac(NN_i-k-i,1)*fac(i,1));
       const3=fac(NN_j-k,1)/(fac(NN_j-k-j,1)*fac(j,1));
       //Rprintf("%f %f %f\n",const2,const3,dens1);
       dens+=dens1*const2*const3*
             R_pow(P11+P10,i) * R_pow(P11+P01,j) *
             R_pow(P00+P01,NN_i-k-i) * R_pow(P00+P10,NN_j-k-j);
  }}}
    return(dens);
}
*/


double  biv_binom2 (int NN_i,int NN_j, int k, int u, int v, double p01,double p10,double p11)
{

int a,i,j;
double const1=0.0,const2=0.0,const3=0.0,dens=0.0,dens1;
double P10=p01-p11;
double P01=p10-p11;
double P00=1+p11-(p01+p10);
double P11=p11;

for(i=0;i<=fmin_int(NN_i-k,u);i++){
   for(j=0;j<=fmin_int(NN_j-k,v);j++){
     
for(a=fmax_int(0,u+v-k-i-j);a<=fmin_int(u-i,v-j);a++){         
       const1=exp(lgammafn(k+1)-(lgammafn(a+1)+lgammafn(u-i-a+1)+lgammafn(v-j-a+1)+lgammafn(k-u-v+i+j+a+1)));
        dens1=const1*R_pow(P11,a)*R_pow(P00,k-u-v+i+j+a)*
             R_pow(P10,u-i-a)*R_pow(P01,v-j-a);
   
       const2=exp(lgammafn(NN_i-k+1)-(lgammafn(NN_i-k-i+1)+lgammafn(i+1)));
       const3=exp(lgammafn(NN_j-k+1)-(lgammafn(NN_j-k-j+1)+lgammafn(j+1)));

       dens+=dens1*const2*const3*
             R_pow(P11+P10,i) * R_pow(P11+P01,j) *
             R_pow(P00+P01,NN_i-k-i) * R_pow(P00+P10,NN_j-k-j);
  }}}
    return(dens);
}

/*
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0;
    a_u=fmin(u,v);
    dens=0;
    for(a=0;a<=a_u;a++){
        dens+=(R_pow(NN*p11,a)*R_pow(NN*(p01-p11),u-a)*R_pow(NN*(p10-p11),v-a))/(fac(a,1)*fac(u-a,1)*fac(v-a,1));
    }
    return(exp(-NN*(p01+p10-p11))*dens);
}



double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0;
    a_u=fmin_int(u,v);
    dens=0;
    pp=p11*p01*p10;
    bb=(1-p11)/p11;
    for(a=0;a<=a_u;a++){
       dens+=(R_pow(NN*((p11*(p01+p10)-p01*p10*(p11+1))/pp),a)*R_pow(NN*((p10-p11)/(p10*p11)),u-a)*R_pow(NN*bb,v-a))/(fac(a,1)*fac(u-a,1)*fac(v-a,1));
    }
    return(exp(-NN*bb)*dens);
}
*/

// bivariate pois-binomial 
double biv_poisbin (int NN, int u, int  v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0, kk=0.0;
    a_u=fmin(u,v);
    dens=0;
    for(a=0;a<=a_u;a++){

    kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
        dens+=kk*(R_pow(NN*p11,a)*R_pow(NN*(p01-p11),u-a)*R_pow(NN*(p10-p11),v-a));
    }
    return(exp(-NN*(p01+p10-p11))*dens);
}



// bivariate pois-bineg 
double biv_poisbinneg (int NN, int u, int v, double p01,double p10,double p11)
{
    int a=0,a_u=0;
    double dens=0.0,pp=0.0,bb=0.0, kk=0.0;
    a_u=fmin_int(u,v);
    dens=0;
    pp=p11*p01*p10;
    bb=(1-p11)/p11;
    for(a=0;a<=a_u;a++){

       kk=exp(-(lgammafn(a+1)+lgammafn(u-a+1)+lgammafn(v-a+1)));
       dens+=kk*(R_pow(NN*((p11*(p01+p10)-p01*p10*(p11+1))/pp),a)*R_pow(NN*((p10-p11)/(p10*p11)),u-a)*R_pow(NN*bb,v-a));
    }
    return(exp(-NN*bb)*dens);
}



// marginal binomial
double marg_binom(int n,double x,double p)
{
 double pr=0.0;double kk;
 kk=fac(n,n-x+1)/fac(x,1);
 pr=kk*R_pow(p,x)*R_pow(1-p,n-x);
return(pr);
}
//marginal geom
double marg_geom(int x,double p)
{
 double pr=0.0;
 pr=p*R_pow(1-p,x);
    return(pr);
}
//marginal binom neg
double marg_binomneg(int n,int x,double p)
{
 double pr=0.0;
 pr=(fac(x+n-1,1)/(fac(n-1,1)*fac(x,1)))*R_pow(p,n)*R_pow(1-p,x);
    return(pr);
}
//marginal poisson
double marg_pois(int n,double x,double p)
{
 double pr=0.0,lambda=0.0;
 lambda=n*p;
 pr=R_pow(lambda,x)*exp(-lambda)/fac(x,1);
    return(pr);
}
                    

double marg_p(double categ_0,double psm,int *model,int n)
{
    double res=0.0;
    if(model[0]==2||model[0]==11) res=marg_binom (n-1,categ_0,psm);
    if(model[0]==14)           res=marg_geom (categ_0,psm);
    if(model[0]==15)           res=marg_pois (n-1,categ_0,psm);
    return(res);
}

/************ pochammer symbols*******************/
double Poch(double q,double n)
{
    return(gammafn(q+n)/gammafn(q));
}
/*********** ratio of gamma function type 1********************/
double aprox(int k,double a, double b) 
{
return(R_pow(k,a-b)*(1+ (a-b)*(a+b-1)/(2*k)));
}
/************ ratio of gamma function type 2 *******************/
double aprox1(double a) {
  double res;
    res=R_pow(R_pow(a,5)+(5*R_pow(a,4)/4)+(25*R_pow(a,3)/32)+(35*R_pow(a,2)/128)+(75*a/2048)+3/8192,1/10);
    return(res);
 }





/**********************************************************/

/* Logarithm of Gamma function */
double lgam(double x)
{
    int sign;
    return lgam_sgn(x, &sign);
}

double lgam_sgn(double x, int *sign)
{
    double p, q, u, w, z;
    int i;
    
    *sign = 1;
    
    if (!isfinite(x))
        return x;
    
    if (x < -34.0) {
        q = -x;
        w = lgam_sgn(q, sign);
        p = floor(q);
        if (p == q) {
        lgsing:
            //mtherr("lgam", SING);
           // printf("SIGN\n");
            return (NPY_INFINITY);
        }
        i = p;
        if ((i & 1) == 0)
            *sign = -1;
        else
            *sign = 1;
        z = q - p;
        if (z > 0.5) {
            p += 1.0;
            z = p - q;
        }
        z = q * sin(NPY_PI * z);
        if (z == 0.0)
            goto lgsing;
        /*     z = log(NPY_PI) - log( z ) - w; */
        z = LOGPI - log(z) - w;
        return (z);
    }
    
    if (x < 13.0) {
        z = 1.0;
        p = 0.0;
        u = x;
        while (u >= 3.0) {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while (u < 2.0) {
            if (u == 0.0)
                goto lgsing;
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if (z < 0.0) {
            *sign = -1;
            z = -z;
        }
        else
            *sign = 1;
        if (u == 2.0)
            return (log(z));
        p -= 2.0;
        x = x + p;
        p = x * polevl(x, B, 5) / p1evl(x, C, 6);
        return (log(z) + p);
    }
    
    if (x > MAXLGM) {
        return (*sign * NPY_INFINITY);
    }
    
    q = (x - 0.5) * log(x) - x + LS2PI;
    if (x > 1.0e8)
        return (q);
    
    p = 1.0 / (x * x);
    if (x >= 1000.0)
        q += ((7.9365079365079365079365e-4 * p
               - 2.7777777777777777777778e-3) * p
              + 0.0833333333333333333333) / x;
    else
        q += polevl(p, A, 4) / x;
    return (q);
}







double hyp2f1( double a,double b,double c,double x)
{
    double d, d1, d2, e;
    double p, q, r, s, y, ax;
    double ia, ib, ic, id, err;
    double t1;
    int i, aid;
    int neg_int_a = 0, neg_int_b = 0;
    int neg_int_ca_or_cb = 0;
    
    err = 0.0;
    ax = fabs(x);
    s = 1.0 - x;
    ia = round(a);        /* nearest integer to a */
    ib = round(b);
    
    if (x == 0.0) {
        return 1.0;
    }
    
    d = c - a - b;
    id = round(d);
    
    if ((a == 0 || b == 0) && c != 0) {
        return 1.0;
    }
    
    if (a <= 0 && fabs(a - ia) < EPS) {    /* a is a negative integer */
        neg_int_a = 1;
    }
    
    if (b <= 0 && fabs(b - ib) < EPS) {    /* b is a negative integer */
        neg_int_b = 1;
    }
    
    if (d <= -1 && !(fabs(d - id) > EPS && s < 0)
        && !(neg_int_a || neg_int_b)) {
        return pow(s, d) * hyp2f1(c - a, c - b, c, x);
    }
    if (d <= 0 && x == 1 && !(neg_int_a || neg_int_b))
        goto hypdiv;
    
    if (ax < 1.0 || x == -1.0) {
        /* 2F1(a,b;b;x) = (1-x)**(-a) */
        if (fabs(b - c) < EPS) {    /* b = c */
            if (neg_int_b) {
                y = hyp2f1_neg_c_equal_bc(a, b, x);
            } else {
                y = pow(s, -a);    /* s to the -a power */
            }
            goto hypdon;
        }
        if (fabs(a - c) < EPS) {    /* a = c */
            y = pow(s, -b);    /* s to the -b power */
            goto hypdon;
        }
    }
    
    
    
    if (c <= 0.0) {
        ic = round(c);        /* nearest integer to c */
        if (fabs(c - ic) < EPS) {    /* c is a negative integer */
            /* check if termination before explosion */
            if (neg_int_a && (ia > ic))
                goto hypok;
            if (neg_int_b && (ib > ic))
                goto hypok;
            goto hypdiv;
        }
    }
    
    if (neg_int_a || neg_int_b)    /* function is a polynomial */
        goto hypok;
    
    t1 = fabs(b - a);
    if (x < -2.0 && fabs(t1 - round(t1)) > EPS) {
        /* This transform has a pole for b-a integer, and
         * may produce large cancellation errors for |1/x| close 1
         */
        p = hyp2f1(a, 1 - c + a, 1 - b + a, 1.0 / x);
        q = hyp2f1(b, 1 - c + b, 1 - a + b, 1.0 / x);
        p *= pow(-x, -a);
        q *= pow(-x, -b);
        t1 = gamma(c);
        s = t1 * gamma(b - a) / (gamma(b) * gamma(c - a));
        y = t1 * gamma(a - b) / (gamma(a) * gamma(c - b));
        return s * p + y * q;
    }
    else if (x < -1.0) {
        if (fabs(a) < fabs(b)) {
            return pow(s, -a) * hyp2f1(a, c - b, c, x / (x - 1));
        }
        else {
            return pow(s, -b) * hyp2f1(b, c - a, c, x / (x - 1));
        }
    }
    
    if (ax > 1.0)        /* series diverges  */
        goto hypdiv;
    
    p = c - a;
    ia = round(p);        /* nearest integer to c-a */
    if ((ia <= 0.0) && (fabs(p - ia) < EPS))    /* negative int c - a */
        neg_int_ca_or_cb = 1;
        
        r = c - b;
        ib = round(r);        /* nearest integer to c-b */
        if ((ib <= 0.0) && (fabs(r - ib) < EPS))    /* negative int c - b */
            neg_int_ca_or_cb = 1;
            
            id = round(d);        /* nearest integer to d */
            q = fabs(d - id);
            
        /* Thanks to Christian Burger <BURGER@DMRHRZ11.HRZ.Uni-Marburg.DE>
         * for reporting a bug here.  */
            if (fabs(ax - 1.0) < EPS) {    /* |x| == 1.0   */
                if (x > 0.0) {
                    if (neg_int_ca_or_cb) {
                        if (d >= 0.0)
                            goto hypf;
                        else
                            goto hypdiv;
                    }
                    if (d <= 0.0)
                        goto hypdiv;
                    y = gamma(c) * gamma(d) / (gamma(p) * gamma(r));
                    goto hypdon;
                }
                if (d <= -1.0)
                    goto hypdiv;
            }
    
    /* Conditionally make d > 0 by recurrence on c
     * AMS55 #15.2.27
     */
    if (d < 0.0) {
        /* Try the power series first */
        y = hyt2f1(a, b, c, x, &err);
        if (err < ETHRESH)
            goto hypdon;
        /* Apply the recurrence if power series fails */
        err = 0.0;
        aid = 2 - id;
        e = c + aid;
        d2 = hyp2f1(a, b, e, x);
        d1 = hyp2f1(a, b, e + 1.0, x);
        q = a + b + 1.0;
        for (i = 0; i < aid; i++) {
            r = e - 1.0;
            y = (e * (r - (2.0 * e - q) * x) * d2 +
                 (e - a) * (e - b) * x * d1) / (e * r * s);
            e = r;
            d1 = d2;
            d2 = y;
        }
        goto hypdon;
    }
    
    
    if (neg_int_ca_or_cb)
        goto hypf;        /* negative integer c-a or c-b */
    
hypok:
    y = hyt2f1(a, b, c, x, &err);
    
    
hypdon:
    if (err > ETHRESH) {
        //mtherr("hyp2f1", PLOSS);
            //  printf( "Estimated err = %.2e\n", err );
    }
    return (y);
    
    /* The transformation for c-a or c-b negative integer
     * AMS55 #15.3.3
     */
hypf:
    y = pow(s, d) * hys2f1(c - a, c - b, c, x, &err);
    goto hypdon;
    
    /* The alarm exit */
hypdiv:
    //mtherr("hyp2f1", OVERFLOW);
    //printf( "Estimated err = %.2e\n", err );
    return NPY_INFINITY;
}






/* Apply transformations for |x| near 1
 * then call the power series
 */
double hyt2f1( double a, double b, double c, double x, double *loss )
{
    double p, q, r, s, t, y, w, d, err, err1;
    double ax, id, d1, d2, e, y1;
    int i, aid, sign;
    
    int ia, ib, neg_int_a = 0, neg_int_b = 0;
    
    ia = round(a);
    ib = round(b);
    
    if (a <= 0 && fabs(a - ia) < EPS) {    /* a is a negative integer */
        neg_int_a = 1;
    }
    
    if (b <= 0 && fabs(b - ib) < EPS) {    /* b is a negative integer */
        neg_int_b = 1;
    }
    
    err = 0.0;
    s = 1.0 - x;
    if (x < -0.5 && !(neg_int_a || neg_int_b)) {
        if (b > a)
            y = pow(s, -a) * hys2f1(a, c - b, c, -x / s, &err);
        
        else
            y = pow(s, -b) * hys2f1(c - a, b, c, -x / s, &err);
        
        goto done;
    }
    
    d = c - a - b;
    id = round(d);        /* nearest integer to d */
    
    if (x > 0.9 && !(neg_int_a || neg_int_b)) {
        if (fabs(d - id) > EPS) {
            int sgngam;
            
            /* test for integer c-a-b */
            /* Try the power series first */
            y = hys2f1(a, b, c, x, &err);
            if (err < ETHRESH)
                goto done;
            /* If power series fails, then apply AMS55 #15.3.6 */
            q = hys2f1(a, b, 1.0 - d, s, &err);
            sign = 1;
            w = lgam_sgn(d, &sgngam);
            sign *= sgngam;
            w -= lgam_sgn(c-a, &sgngam);
            sign *= sgngam;
            w -= lgam_sgn(c-b, &sgngam);
            sign *= sgngam;
            q *= sign * exp(w);
            r = pow(s, d) * hys2f1(c - a, c - b, d + 1.0, s, &err1);
            sign = 1;
            w = lgam_sgn(-d, &sgngam);
            sign *= sgngam;
            w -= lgam_sgn(a, &sgngam);
            sign *= sgngam;
            w -= lgam_sgn(b, &sgngam);
            sign *= sgngam;
            r *= sign * exp(w);
            y = q + r;
            
            q = fabs(q);    /* estimate cancellation error */
            r = fabs(r);
            if (q > r)
                r = q;
            err += err1 + (MACHEP * r) / y;
            
            y *= gamma(c);
            goto done;
        }
        else {
            /* Psi function expansion, AMS55 #15.3.10, #15.3.11, #15.3.12
             *
             * Although AMS55 does not explicitly state it, this expansion fails
             * for negative integer a or b, since the psi and Gamma functions
             * involved have poles.
             */
            
            if (id >= 0.0) {
                e = d;
                d1 = d;
                d2 = 0.0;
                aid = id;
            }
            else {
                e = -d;
                d1 = 0.0;
                d2 = d;
                aid = -id;
            }
            
            ax = log(s);
            
            /* sum for t = 0 */
            y = digamma(1.0) + digamma(1.0 + e) - digamma(a + d1) - digamma(b + d1) - ax;
            y /= gamma(e + 1.0);
            
            p = (a + d1) * (b + d1) * s / gamma(e + 2.0);    /* Poch for t=1 */
            t = 1.0;
            do {
                r = digamma(1.0 + t) + digamma(1.0 + t + e) - digamma(a + t + d1)
                - digamma(b + t + d1) - ax;
                q = p * r;
                y += q;
                p *= s * (a + t + d1) / (t + 1.0);
                p *= (b + t + d1) / (t + 1.0 + e);
                t += 1.0;
                if (t > MAX_ITERATIONS) {    /* should never happen */
                    //mtherr("hyp2f1", TOOMANY);
                   // printf( "Estimated err = %.2e\n", err );
                    *loss = 1.0;
                    return NPY_NAN;
                }
            }
            while (y == 0 || fabs(q / y) > EPS);
            
            if (id == 0.0) {
                y *= gamma(c) / (gamma(a) * gamma(b));
                goto psidon;
            }
            
            y1 = 1.0;
            
            if (aid == 1)
                goto nosum;
            
            t = 0.0;
            p = 1.0;
            for (i = 1; i < aid; i++) {
                r = 1.0 - e + t;
                p *= s * (a + t + d2) * (b + t + d2) / r;
                t += 1.0;
                p /= t;
                y1 += p;
            }
        nosum:
            p = gamma(c);
            y1 *= gamma(e) * p / (gamma(a + d1) * gamma(b + d1));
            
            y *= p / (gamma(a + d2) * gamma(b + d2));
            if ((aid & 1) != 0)
                y = -y;
            
            q = pow(s, id);    /* s to the id power */
            if (id > 0.0)
                y *= q;
            else
                y1 *= q;
            
            y += y1;
        psidon:
            goto done;
        }
        
    }
    
    /* Use defining power series if no special cases */
    y = hys2f1(a, b, c, x, &err);
    
done:
    *loss = err;
    return (y);
}





/* Defining power series expansion of Gauss hypergeometric function */

double hys2f1( double a,double b,double c,double x,double *loss )
{
    double f, g, h, k, m, s, u, umax;
    int i;
    int ib, intflag = 0;
    
    if (fabs(b) > fabs(a)) {
        /* Ensure that |a| > |b| ... */
        f = b;
        b = a;
        a = f;
    }
    
    ib = round(b);
    
    if (fabs(b - ib) < EPS && ib <= 0 && fabs(b) < fabs(a)) {
        /* .. except when `b` is a smaller negative integer */
        f = b;
        b = a;
        a = f;
        intflag = 1;
    }
    
    if ((fabs(a) > fabs(c) + 1 || intflag) && fabs(c - a) > 2
        && fabs(a) > 2) {
        /* |a| >> |c| implies that large cancellation error is to be expected.
         *
         * We try to reduce it with the recurrence relations
         */
        return hyp2f1ra(a, b, c, x, loss);
    }
    
    i = 0;
    umax = 0.0;
    f = a;
    g = b;
    h = c;
    s = 1.0;
    u = 1.0;
    k = 0.0;
    do {
        if (fabs(h) < EPS) {
            *loss = 1.0;
            return NPY_INFINITY;
        }
        m = k + 1.0;
        u = u * ((f + k) * (g + k) * x / ((h + k) * m));
        s += u;
        k = fabs(u);        /* remember largest term summed */
        if (k > umax)
            umax = k;
        k = m;
        if (++i > MAX_ITERATIONS) {    /* should never happen */
            *loss = 1.0;
            return (s);
        }
    }
    while (s == 0 || fabs(u / s) > MACHEP);
    
    /* return estimated relative error */
    *loss = (MACHEP * umax) / fabs(s) + (MACHEP * i);
    
    return (s);
}


/*
 * Evaluate hypergeometric function by two-term recurrence in `a`.
 *
 * This avoids some of the loss of precision in the strongly alternating
 * hypergeometric series, and can be used to reduce the `a` and `b` parameters
 * to smaller values.
 *
 * AMS55 #15.2.10
 */
double hyp2f1ra(double a, double b, double c, double x,
                       double *loss)
{
    double f2, f1, f0;
    int n;
    double t, err, da;
    
    /* Don't cross c or zero */
    if ((c < 0 && a <= c) || (c >= 0 && a >= c)) {
        da = round(a - c);
    }
    else {
        da = round(a);
    }
    t = a - da;
    
    *loss = 0;
    
    assert(da != 0);
    
    if (fabs(da) > MAX_ITERATIONS) {
        /* Too expensive to compute this value, so give up */
        //mtherr("hyp2f1", TLOSS);
        //printf( "TLOSS\n");
        *loss = 1.0;
        return NPY_NAN;
    }
    
    if (da < 0) {
        /* Recurse down */
        f2 = 0;
        f1 = hys2f1(t, b, c, x, &err);
        *loss += err;
        f0 = hys2f1(t - 1, b, c, x, &err);
        *loss += err;
        t -= 1;
        for (n = 1; n < -da; ++n) {
            f2 = f1;
            f1 = f0;
            f0 = -(2 * t - c - t * x + b * x) / (c - t) * f1 - t * (x -
                                                                    1) /
            (c - t) * f2;
            t -= 1;
        }
    }
    else {
        /* Recurse up */
        f2 = 0;
        f1 = hys2f1(t, b, c, x, &err);
        *loss += err;
        f0 = hys2f1(t + 1, b, c, x, &err);
        *loss += err;
        t += 1;
        for (n = 1; n < da; ++n) {
            f2 = f1;
            f1 = f0;
            f0 = -((2 * t - c - t * x + b * x) * f1 +
                   (c - t) * f2) / (t * (x - 1));
            t += 1;
        }
    }
    
    return f0;
}


/*
 15.4.2 Abramowitz & Stegun.
 */
double hyp2f1_neg_c_equal_bc(double a, double b, double x)
{
    double k;
    double collector = 1;
    double sum = 1;
    double collector_max = 1;
    
    if (!(fabs(b) < 1e5)) {
        return NPY_NAN;
    }
    
    for (k = 1; k <= -b; k++) {
        collector *= (a + k - 1)*x/k;
        collector_max = fmax(fabs(collector), collector_max);
        sum += collector;
    }
    
    if (1e-16 * (1 + collector_max/fabs(sum)) > 1e-7) {
        return NPY_NAN;
    }
    
    return sum;
}



double hypergeo(double a,double b,double c,double x)
{
    double sol;
    sol =hyp2f1( a, b, c, x );
    return(sol);
}

void hypergeo_call(double *a,double *b,double *c,double *x, double *res)
{
    *res = hypergeo(*a,*b,*c,*x);
}





/***********************************************************/
/*********** bivariate T distribution********************/ 
/***********************************************************/
/*********** bivariate T distribution********************/ 
double biv_T(double rho,double zi,double zj,double nuu)
{
  int k=0; 
  double nu=1/nuu;
  double res0=0.0,RR=0.0,pp1=0.0,pp2=0.0;
  double bb1,bb2;
  double x=zi;double y=zj;        
  double cc=(nu+1)/2; double nu2=nu/2;
  double x1=(x*x+nu); double y1=(y*y+nu);
  double rho2=R_pow(1-rho*rho,-cc);

  double b1 = (R_pow(nu,nu))*R_pow(x1*y1,-cc)*R_pow(gammafn(cc),2);
  double c1 = M_PI*R_pow(gammafn(nu/2),2)*rho2;
  double b2 = rho*x*y*R_pow(nu,nu+2)*R_pow(x1*y1,-nu2-1);
  double c2 = 2*M_PI*rho2;

  double a1 = 0; double a2 = 0;
  double aux  = R_pow(rho*x*y,2)/(x1*y1);
  double aux1 = R_pow(rho*nu,2)/(x1*y1);
 //if(fabs(rho)<=EPS1)
  /*if(!fabs(rho))  
  {
    C = lgammafn(cc)+log(R_pow((1+x*x/nu),-cc))-log(sqrt(M_PI*nu))-lgammafn(nu/2);
    B = lgammafn(cc)+log(R_pow((1+y*y/nu),-cc))-log(sqrt(M_PI*nu))-lgammafn(nu/2);
    return(exp(B)*exp(C));
  }*/
  while( k<=6000 )
    {
   // pp1=hypergeo(cc+k,cc+k,0.5,aux);
    pp1=(0.5-2*(cc+k))*log(1-aux)+log(hypergeo(0.5-(cc+k),0.5-(cc+k),0.5,aux)); //euler
    bb1=pp1+k*log(aux1)+2*(lgammafn(cc+k)-lgammafn(cc))-lgammafn(k+1)-lgammafn(nu2+k)+lgammafn(nu2);
    a1 = a1 + exp(bb1);
   // pp2=hypergeo(nu2+1+k,nu2+1+k,1.5,aux);
    pp2=(1.5-2*(nu2+1+k))*log(1-aux)+log(hypergeo(1.5-(nu2+1+k),1.5-(nu2+1+k),1.5,aux));//euler
    bb2=pp2+k*log(aux1)+2*log((1+k/nu2))+lgammafn(nu2+k)-lgammafn(k+1)-lgammafn(nu2);
    a2 = a2 + exp(bb2);
    RR=(b1/c1)*a1+(b2/c2)*a2;
   if(!R_FINITE(RR)) return(res0);
    if((fabs(RR-res0)<1e-10)  ) {break;}
    else {res0=RR;}
        k++;
    }
return(RR);
}
/*********** Appell F4 function ********/
double appellF4(double a,double b,double c,double d,double x,double y)
{
double RR=0.0,bb=0.0,res0=0.0;int k=0;
  while( k<=6000 )
    {
    bb=k*log(y)+(lgammafn(a+k)+lgammafn(b+k)+lgammafn(d))
               -(lgammafn(a)+lgammafn(b)+lgammafn(d+k)+lgammafn(k+1))
               +(c-(a+k)-(b+k))*log(1-x)+log(hypergeo(c-a-k,c-b-k,c,x)); //euler
              // +log(hypergeo(a+k,b+k,c,x));
    RR=RR+exp(bb);
 if(!R_FINITE(RR)) return(res0);
    if((fabs(RR-res0)<1e-10)  ) {break;}
    else {res0=RR;}
        k++;
    }
return(RR);
}




/****************************************/
double appellF4_mod(double nu,double rho2,double x,double y)
{
  double xx,yy,x2,y2,arg,arg1,pp1,pp2,app;
xx=x*x;yy=y*y;
x2=xx+nu;
y2=yy+nu;
arg=(nu+1)/2; 
arg1=nu/2;
pp1=R_pow(nu,nu)*R_pow(x2*y2,-arg)*R_pow(gammafn(arg),2);
pp2=M_PI*R_pow(gammafn(arg1),2)*R_pow(1-rho2,-arg);
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
{res=          (p11/R_pow(etamos,2))*appellF4_mod(nu,rho2,zistd/etamos,zjstd/etamos);}
if(zi>=mui&&zj<muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho2,zistd/etamos,zjstd/etamas);}
if(zi<mui&&zj>=muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*appellF4_mod(nu,rho2,zistd/etamas,zjstd/etamos);}
if(zi<mui&&zj<muj)
{res=    ((p11+eta)/R_pow(etamas,2))*appellF4_mod(nu,rho2,zistd/etamas,zjstd/etamas);}
return(res/sill);
}
/***** bivariate half gaussian ****/     
double biv_half_Gauss(double rho,double zi,double zj)
{
double kk=0, dens=0,a=0,b=0,rho2=1-rho*rho;
  kk=(M_PI)*sqrt(rho2);
  a=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)-2*rho*zi*zj));
  b=exp(- (1/(2*(rho2)))*(R_pow(zi,2)+R_pow(zj,2)+2*rho*zi*zj));
  dens=(a+b)/kk;
  return(dens);
}
/***** bivariate two piece gaussian ****/ 
double biv_two_pieceGaussian(double rho,double zi,double zj,double sill,double eta,
             double p11,double mui,double muj)
{
double res;  
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);
if(zi>=mui&&zj>=muj)
{res=          (p11/R_pow(etamos,2))*biv_half_Gauss(rho,zistd/etamos,zjstd/etamos);}
if(zi>=mui&&zj<muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Gauss(rho,zistd/etamos,zjstd/etamas);}
if(zi<mui&&zj>=muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Gauss(rho,zistd/etamas,zjstd/etamos);}
if(zi<mui&&zj<muj)
{res=    ((p11+eta)/R_pow(etamas,2))*biv_half_Gauss(rho,zistd/etamas,zjstd/etamas);}
return(res/sill);
}



double biv_Kumara(double rho,double zi,double zj,double ai,double aj,double shape1,double shape2)
{
  double xx=0.0,yy=0.0,ki=0.0,kj=0.0,p1=0.0,p2=0.0,rho2=0.0,res=0.0;
 ki=1-pow(zi,shape2); kj=1-pow(zj,shape2);
//if(rho) {
  rho2=rho*rho;
  xx=rho2*pow(ki*kj,shape1);
  yy=rho2*(1-pow(ki,shape1))*(1-pow(kj,shape1));
  p1=pow(shape1*shape2,2)*pow(zi*zj,shape2-1)*pow(ki*kj,shape1-1)*pow(1-rho2,2);
    p2= appellF4(2,2,1,1,xx,yy);
  res=p1*p2;
/*} else  {p1=shape1*shape2*pow(zi,shape2-1)* pow(ki,shape1-1);
         p2=shape1*shape2*pow(zj,shape2-1)* pow(kj,shape1-1);
         res=p1*p2;}*/
return(res);
}

/***********************************************************************************/
/************ functions for binomial or negative binomial  two piece *********/
/***********************************************************************************/
/***********************************************************************************/
/************ functions for binomial or negative binomial  two piece *********/
/***********************************************************************************/

double pbnorm22(double lim1,double lim2,double corr,double nugget)
{
    double  lowe[2]  = {0,0}, uppe[2] = {lim1,lim2}, corre[1] = {(1-nugget)*corr};
    double  value;
    int     infin[2] = {0,0};
    value            = F77_CALL(bvnmvn)(lowe,uppe,infin,corre); 
    return(value);
}

// cdf bivariate half-normal distribution
double pbhalf_gauss(double zi,double zj,double rho,double nugget)
{
  double dens = 0;
  dens = pbnorm22(zi,zj,rho,nugget) + pbnorm22(-zi,-zj,rho,nugget) -
              pbnorm22(-zi,zj,rho,nugget) - pbnorm22(zi,-zj,rho,nugget);
  return(dens); 
}
/****** cdf univariate two-piece gaussian distribution *****/
double pnorm_two_piece(double x, double eta)
{
    double cdf = 0;
    if (x <=  0) cdf = (1 + eta)*pnorm(x/(1 + eta),0,1,1,0);
    if (x > 0)   cdf = eta + (1 - eta)*pnorm(x/(1 - eta),0,1,1,0);
  return(cdf);
}
//***********************************************//
// cdf univariate half-gaussian distribution
double phalf_gauss (double z){
  double dens = 0;
  dens = 2 * pnorm(z,0,1,1,0) - 1;
  return(dens);
}

// cdf bivariate two_piece gaussian distribution
double pbnorm_two_piece(int *cormod, double h, double u, 
    double xi, double xj, double nugget, double var,double eta,double *par)
{
    double p11,g_eta=0.0,q_g_eta=0.0,dens=0.0,corr=0.0;
    double etamos,etamas;
    g_eta   = (1 - eta)/2;
    q_g_eta = qnorm(g_eta,0,1,1,0);
    corr    = CorFct(cormod,h,u,par,0,0);
    p11     = pbnorm22(q_g_eta,q_g_eta,corr,nugget);
    etamas = 1+eta;
    etamos = 1-eta;
    if(xi  <= 0 & xj <= 0) dens = (1 + pbhalf_gauss(-xi/etamas,-xj/etamas,corr,nugget) - phalf_gauss(-xi/etamas) - phalf_gauss(-xj/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) );
    if(xi  > 0  & xj <= 0) dens = (1 + pbhalf_gauss(0,-xj/etamas,corr,nugget) - phalf_gauss(-xj/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + 
      (phalf_gauss(xi/etamos) - pbhalf_gauss(xi/etamos,-xj/etamas,corr,nugget)) * (pnorm(q_g_eta,0,1,1,0) - p11);
    if(xi  <= 0 & xj >  0) dens =  (1 + pbhalf_gauss(-xi/etamas,0,corr,nugget) - phalf_gauss(-xi/etamas)) * ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + 
      (phalf_gauss(xj/etamos) - pbhalf_gauss(-xi/etamas,xj/etamos,corr,nugget)) * (pnorm(q_g_eta,0,1,1,0) - p11);
    if(xi  > 0  & xj >  0) dens =  ( 1 + p11 - 2*pnorm(q_g_eta,0,1,1,0) ) + (phalf_gauss(xi/etamos) - pbhalf_gauss(xi/etamos,0,corr,nugget))*(pnorm(q_g_eta,0,1,1,0)-p11) + 
     (phalf_gauss(xj/etamos) - pbhalf_gauss(0,xj/etamos,corr,nugget))* (pnorm(q_g_eta,0,1,1,0)-p11) + pbhalf_gauss(xi/etamos,xj/etamos,corr,nugget)* p11;
    return(dens);
}
/***********************************************************************************/
/***********************************************************************************/
/********** bivariate log-logistic **********/
double biv_LogLogistic(double corr,double zi,double zj,double mui, double muj, double shape)
{
    double c=gammafn(1+1/shape)*gammafn(1-1/shape);
    double A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=exp(mui);double cj=exp(muj);
    double ki=R_pow(c*zi/ci,shape)+1;
    double kj=R_pow(c*zj/cj,shape)+1;
    double rho2=R_pow(corr,2);
    double kij=ki*kj;
   // if(corr)   {
        A=(R_pow(c*shape,2)/(ci*cj))*R_pow((c*c*zi*zj)/(ci*cj),shape-1)*R_pow(1-rho2,2);
        B=R_pow(kij,-2);
        C=appellF4(2,2,1,1,
        (rho2*R_pow(c*c*zi*zj,shape))/(R_pow(ci*cj,shape)*kij),
        rho2/(kij));
        res=A*B*C;
    //}else{    B=(c*shape/ci)*R_pow((c*zi/ci),shape-1)*R_pow(ki,-2);  C=(c*shape/cj)*R_pow((c*zj/cj),shape-1)*R_pow(kj,-2);
    //    res=B*C;}
    return(res);
}

/********** bivariate logistic **********/
double biv_Logistic(double corr,double zi,double zj,double mui, double muj, double sill)
{
    double a=0.0,A=0.0,res=0.0,B=0.0,C=0.0;
    double ci=mui;double cj=muj;
    double ki=exp((zi-ci)/sqrt(sill));
    double kj=exp((zj-cj)/sqrt(sill));
    double rho2=R_pow(corr,2);
    //if(corr)   {
        a=1-rho2;
        A=(ki*kj)/(R_pow(a,-2)*sill);
        B=R_pow((ki+1)*(kj+1),-2);
        C=appellF4(2,2,1,1,
        (rho2*ki*kj)/((ki+1)*(kj+1)),
        rho2/((ki+1)*(kj+1)));
        res=A*B*C;
    //} else{ B=ki*R_pow((ki+1),-2)/sqrt(sill);C=kj*R_pow((kj+1),-2)/sqrt(sill);res=B*C;}
    return(res);
}




//********** functions for bivariate tukey h *******************//

// compute lambert w function
double LambertW(double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
  //if (dbgW) fprintf(stderr,"LambertW: z=%g\n",z);
  if (z<-em1 || isinf(z) || isnan(z)) { 
    //fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); exit(1);
      return 0.0;
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return 
     -1.0
     +2.331643981597124203363536062168*r
     -1.812187885639363490240191647568*q
     +1.936631114492359755363277457668*r*q
     -2.353551201881614516821543561516*q2
     +3.066858901050631912893148922704*r*q2
     -4.175335600258177138854984177460*q3
     +5.858023729874774148815053846119*r*q3
     -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
    p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else 
    w=log(z); /* asymptotic */
  if (z>3.0) w-=log(w); /* useful? */
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  //fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z);
    return 0.0;
  //exit(1);
}

// pdf bivariate gaussian distribution
double dbnorm(double x_i,double x_j,double mean_i,double mean_j,double sill,double corr)
{
    double  fraq = 1.0,dens = 0.0,aux1 = 1.0,z1 = 1.0,z2 = 1.0,z3 = 1.0,z = 1.0;
    fraq = 2*M_PI*sill*sqrt(1-corr*corr);
    z1   = (x_i - mean_i)*(x_i - mean_i);
    z2   = (x_j - mean_j)*(x_j - mean_j);
    z3   = 2*corr*(x_i - mean_i)*(x_j - mean_j);
    z    = (z1 + z2 - z3)/sill;
    aux1 = 2*(1-corr*corr); 
    dens = (1/fraq)*exp(-z/aux1);
    return(dens);
}


// compute the inverse lambert w transformation
double inverse_lamb(double x,double tail)
{
  double sign,value;
  value = sqrt(LambertW(tail*x*x)/tail);
  if (x > 0) sign= 1;
  if (x < 0) sign= -1;
   return(sign*value);
}


// pdf bivariate tukey h random field 
double biv_tukey_h(double corr,double data_i, double data_j, double mean_i, double mean_j, double tail, double sill)
{
  double dens = 0.0,x_i = 0.0,x_j = 0.0,est_mean_i = 0.0,est_mean_j = 0.0;
  double est_mean_ij = 1.0,extra = 1.0;

  est_mean_i = (data_i - mean_i)/sqrt(sill);
  est_mean_j = (data_j - mean_j)/sqrt(sill);
  
  x_i = inverse_lamb(est_mean_i,tail);
  x_j = inverse_lamb(est_mean_j,tail);

  est_mean_ij = 1/(est_mean_i*est_mean_j);
    extra       = 1/( (1 + LambertW(tail*est_mean_i*est_mean_i))*(1 + LambertW(tail*est_mean_j*est_mean_j)));
  dens = dbnorm(x_i,x_j,0,0,1,corr)*
              x_i*x_j * est_mean_ij * extra/sill;
  if((x_i==0.0)&&(x_j!=0.0))  dens = dbnorm(x_i,x_j,0,0,1,corr)*x_j/(est_mean_j*(1 + LambertW(tail*est_mean_j*est_mean_j)));
  if((x_j==0.0)&&(x_i!=0.0))  dens = dbnorm(x_i,x_j,0,0,1,corr)*x_i/(est_mean_i*(1 + LambertW(tail*est_mean_i*est_mean_i)));
  if((x_j==0.0)&&(x_i==0.0))  dbnorm(x_i,x_j,0,0,1,corr);
  return(dens);
}








/***** bivariate half tukey h ****/     
double biv_half_Tukeyh(double rho,double ti,double tj,double tail)
{
  double dens = 0.0;
  dens = biv_tukey_h(rho,ti,tj,0,0,tail,1) + biv_tukey_h(rho,-ti,-tj,0,0,tail,1) + biv_tukey_h(rho,-ti,tj,0,0,tail,1) + biv_tukey_h(rho,ti,-tj,0,0,tail,1);
  return(dens);
}
 

/***** bivariate two piece tukey h ****/ 
double biv_two_pieceTukeyh(double rho,double zi,double zj,double sill,double eta,double tail,
             double p11,double mui,double muj)
{
double res;  
double etamas=1+eta;
double etamos=1-eta;
double zistd=(zi-mui)/sqrt(sill);
double zjstd=(zj-muj)/sqrt(sill);

if(zi>=mui&&zj>=muj)
{res=          (p11/R_pow(etamos,2))*biv_half_Tukeyh(rho,zistd/etamos,zjstd/etamos,tail);}
if(zi>=mui&&zj<muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,zistd/etamos,zjstd/etamas,tail);}
if(zi<mui&&zj>=muj)
{res=((1-eta-2*p11)/(2*(1-eta*eta)))*biv_half_Tukeyh(rho,zistd/etamas,zjstd/etamos,tail);}
if(zi<mui&&zj<muj)
{res=    ((p11+eta)/R_pow(etamas,2))*biv_half_Tukeyh(rho,zistd/etamas,zjstd/etamas,tail);}

return(res/sill);
}


/*******************************************************************************/



double  binomialCoeff(double n, double k) 
{ 
    double res=0.0;
    res=lgammafn(n+1)-(lgammafn(k+1)+lgammafn(n-k+1));
    return(exp(res)); 
} 

/*************************************************************/
double Prt(double corr,double r, double t, double mean_i, double mean_j){
    double rho2= pow(corr,2);
    double prt;
    double term =0, term1 =0 ;
    double value = 0, value1 = 0;
    double auxi= mean_i/(1-rho2);
    double auxj= mean_j/(1-rho2);
    double n= r-t;
    double aux=0, aux1=0;
    int k=0,j=0,l1=0,l2=0,m=0;
    for(k=0; k<50; ++k){
        for(int j=0; j<=n; ++j){
            aux= binomialCoeff(n, j)*pow((1-rho2)/rho2,j+t);
            aux1= exp(lgammafn(t+j+k)+(n-j)*log(mean_i)-lgammafn(k+1))*pow(-1,j);
            term1= aux*aux1*igam(t+j+k, rho2*auxi)*igam(t+k, auxj);
            //if(term1<1e-10) break;
            value1 =value1+ term1;
        }
    }
    int iter=45;
    double aux2=0, aux3=0, aux4=0;
    for(k=0; k<iter; ++k){
        for(j=0; j<=n-1; ++j){
            for(m=0; m<iter; ++m){
                for(l1=0; l1<iter; ++l1){
                    for(l2=0; l2<iter; ++l2){
                        aux2= binomialCoeff(n-1, j)*R_pow(rho2,k+m+l1)*R_pow(1-rho2,-2*k-2*m-l1-l2-1)*R_pow(-1,j+l1+l2);
                        aux3= exp(lgammafn(t+m)+2*lgammafn(k+1)+lgammafn(m+1)+lgammafn(l1+1)+lgammafn(l2+1)+log((j+k+l1+1)*(k+l2+1))-lgammafn(t));
                        aux4= R_pow(mean_i,m+k+n+l1)*R_pow(mean_j,m+k+l2+1);
                        term= (aux2/aux3)*aux4*exp(lbeta(n+k+l1+1,t+m)+lbeta(k+l2+2,t+m))*hyperg(t+m, n+k+l1+1+t+m,-rho2*auxi)*hyperg(t+m, k+l2+2+t+m,-auxj);
                        value =value+ term;
                    }
                }
            }
        }
    }

     prt= exp(-mean_i+log(value1)-(lgammafn(n+1)+lgammafn(t)))- 
        exp(-mean_i+t*log(mean_i*mean_j/(1-rho2))+log(value)-(2*lgammafn(t)+lgammafn(n)));
    return(prt);
}

//*****************************************************************************/

double Prr(double corr,double r, double t, double mean_i, double mean_j){
    double rho2= pow(corr,2);
    double prr;
    
    double term, term1=0.0, term2=0.0, term3=0.0;
    double value = 0, value1 = 0, value2 = 0, value3 = 0 ;
    int k = 0, m=0, l1=0, l2=0 ;
    double auxi= mean_i/(1-rho2);
    double auxj= mean_j/(1-rho2);
    
    while ( term1>1e-10 && term2> 1e-10 && term3>1e-10 )
    {
      term1 = pow(rho2,k)* gammafn(r+k) * igam(r+k, auxi) * igam(r+k, auxj)/gammafn(k+1) ;
      term2 =  exp(lgammafn(r+k) + log(igam(r+k, rho2*auxi))+log(igam(r+k, auxj))-lgammafn(k+1)) ;
      term3 =  exp(lgammafn(r+k) + log(igam(r+k, auxi)) + log(igam(r+k, rho2*auxj))-lgammafn(k+1) ) ;
      value1 += term1;
      value2 += term2;
      value3 += term3;
      k++;
     } 
    
    double aux=0, aux1=0, aux2=0; 
     int iter=50;
      for(k=0; k<iter; ++k){
          for(m=0; m<iter; ++m){
              for(l1=0; l1<iter; ++l1){
                  for( l2=0; l2<iter; ++l2){
                      aux=pow(rho2,k+m)*pow(1-rho2,-2*k-2*m-l1-l2)*pow(-1,l1+l2);
                      aux1= exp(lgammafn(r+m)+2*lgammafn(k+1)+lgammafn(m+1)+lgammafn(l1+1)+lgammafn(l2+1)+log((k+l1+1)*(k+l2+1))-lgammafn(r));
                      aux2=pow(mean_i,m+k+l1)*pow(mean_j,m+k+l2);
                      term=(aux/aux1)*aux2*exp(lbeta(k+l1+2,r+m)+lbeta(k+l2+2,r+m))*hyperg(r+m, k+l1+2+r+m,-1*auxi)*hyperg(r+m, k+l2+2+r+m,-1*auxj);
                      value += term;
                  }
              }
          }
      }
    
    prr =  (pow(1-rho2,r)/gammafn(r))*(-1*value1+exp(-1*mean_i)*pow(rho2,-1*r)*value2+exp(-1*mean_j)*pow(rho2,-1*r)*value3)+ pow(mean_i*mean_j/(1-rho2),r+1)*value/pow(gammafn(r),2);
    return prr;
    
}

/*******+++++++++++++++++++++++++*********************************************/
double Pr0(double corr,double r, double t, double mean_i, double mean_j){
    double rho2= pow(corr,2);
    double pr0;
    
    double term =0 ;
    double value = 0;
    
    double auxi= mean_i/(1-rho2);
    double auxj= mean_j/(1-rho2);
   
    double aux=0, aux1=0;

    int iter=100;
    for(int k=0; k<iter; ++k){
        for(int j=0; j<=r-1; ++j){
            aux= binomialCoeff(r-1, j)*pow((1-rho2)/rho2,j+1);
            aux1= exp(lgammafn(j+k+1)+ (r-j-1)*log(mean_i)-lgammafn(k+1))*pow(-1,j);
            term= aux*aux1*igam(j+k+1, rho2*auxi)*igam(k+1, auxj);
            value =value+ term;
        }
    }
    pr0= exp(r*log(mean_i)-mean_i-lgammafn(r+1))-
    exp(-mean_i + log(value)-lgammafn(r));
    return(pr0);
}

/*******************************************************************************/
double P00(double corr,double r, double t, double mean_i, double mean_j){
    double rho2= R_pow(corr,2);
    double p00;int k = 0;
    double sum = 0.0,res0=0.0;
    double auxi= mean_i/(1-rho2);
    double auxj= mean_j/(1-rho2);
    
    while(k<1000){
             sum =sum+ exp( k*log(rho2) + log(igam(k+1, auxi)) + log(igam(k+1, auxj) )) ;
    if((fabs(sum-res0)<1e-10)  ) {break;}
else {res0=sum;}
        k++;}

    p00 = -1+ exp(-mean_i)+ exp(-mean_j)+(1-rho2)*sum;
    return(p00);
}




double biv_Poisson(double corr,double r, double t, double mean_i, double mean_j)
{
double dens;
if(r==t)
{if(r==0) dens=P00(corr,r,r,mean_i,mean_j);
     if(r>0)  dens=Prr(corr,r,r,mean_i,mean_j);
}

if(r==0&&t>0) dens=Pr0(corr,t,r,mean_j,mean_i);
if(r>0&&t==0) dens=Pr0(corr,r,t,mean_i,mean_j);

if(r>0&&t>0)
{  
if(r>t) dens=Prt(corr,r,t,mean_i,mean_j);
if(t>r) dens=Prt(corr,t,r,mean_j,mean_i);
}
return(dens);

}



