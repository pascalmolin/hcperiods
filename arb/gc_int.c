/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
#include "complex_extras.h"

double
distance_ellipse_d(double x, double y, double a)
{
    double b, t, t1, st, ct;
    b = sqrt(a*a-1);
    t = atan(a*y/(b*x));

    do
    {
        double ft, fpt;
        st = sin(t);
        ct = cos(t);
        ft = ct*st-x*a*st+y*b*ct;
        fpt = ct*(ct-x*a) - st*(st + y*b);
        t1 = t;
        t -= ft / fpt;
    }
    while (t != t1);

    x -= a * ct;
    y -= b * st;

    return sqrt(x*x + y*y);
}

double
constant_m_d(cdouble * w, slong len, double r, slong d)
{
    /* 2*Pi*r^d/sqrt(prod d(w_i,e_r)) */
    slong k;
    double p;
    const double eps = 1.e-5;
    p = 1;
    for (k = 0; k < len; k++)
    {
        double d = distance_ellipse_d(creal(w[k]), cimag(w[k]), r);
        p *= (d < eps) ? r : d;
    }
    return 2 * PI * pow(r, d) / sqrt(p);
}

slong
gc_int_params(double r, sec_t c, slong prec)
{
    slong k, n;
    double M, A, B, rho;
    double mult = .25;
    cdouble * w;

    /* small approx of roots */
    w = malloc(c.d * sizeof(cdouble));
    for (k = 0; k < c.d; k++)
        w[k] = acb_get_cdouble(c.roots + k);

    /* first estimate */
    M = constant_m_d(w, c.d, r, c.d);
    A = prec*log(10) + log(2*PI*M) + 1;
    rho = r + sqrt(r*r-1);
    B = 2*log(rho);
    printf("n1 = %lf\n", A/B);
    n = ceil(A / B);

    /* second = should be exact */
    r *= 1-mult/n;
    M = constant_m_d(w, c.d, r, c.d);
    A = prec*log(10) + log(2*PI*M) + 1;
    rho = r + sqrt(r*r-1);
    B = 2*log(rho);
    printf("n2 = %lf\n", A/B);

    free(w);

    return ceil(A/B);
}

void
gc_distance_ellipse(mag_t d, const acb_t p, mag_t r)
{
    return;
}

void
gc_constant_M(mag_t M, acb_srcptr w, slong len, slong d, mag_t r)
{
    slong i;
    mag_t di, p;
    /* 2*Pi*r^d/sqrt(prod d(w_i,e_r)) */
    mag_const_pi(M);
    mag_mul_2exp_si(M, M, 1);
    mag_pow_ui(p, r, d);
    mag_mul(M, M, p);
    mag_one(p);
    for (i = 0; i < len; i++)
    {
        gc_distance_ellipse(di, w + i, r);
        mag_mul(p, p, di);
    }
    mag_rsqrt(p, p);
    mag_mul(M, M, p);
}

/* error for computation with n points */
void
gc_error(mag_t e, slong n, acb_srcptr w, slong d, mag_t r)
{
    return;
}
