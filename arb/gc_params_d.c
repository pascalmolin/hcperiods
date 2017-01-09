/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

#define PARAMS 0

double
distance_ellipse_d(double x, double y, double a)
{
    const double eps = 1.e-10;
    double b, t, ft, st, ct;
    b = sqrt(a*a-1);
    t = atan(a*y/(b*x));

    do
    {
        double fpt;
        st = sin(t);
        ct = cos(t);
        ft = ct*st - x*a*st + y*b*ct;
        fpt = ct*(ct-x*a) - st*(st + y*b);
        t -= ft / fpt;
    }
    while (fabs(ft) > eps);

    x -= a * ct;
    y -= b * st;

    return sqrt(x*x + y*y);
}

double
constant_m_d(const cdouble * w, slong len, double r, slong d)
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
gc_params_d(const cdouble * w, slong len, double r, slong i, slong prec)
{
    slong n;
    double M, A, B, rho;
    double mult = .25;

    if (r <= 0)
    {
        slong k;
        for (k = 0; k < len; k++)
        {
            double t = (cabs(w[k] - 1) + cabs(w[k] + 1)) / 2;
            if (!k || t < r)
                r = t;
        }
    }

    if (r <= 1)
    {
        slong k;
        flint_printf("gc int: r must be > 1 (r = %lf)\n", r);
        for (k = 0; k < len; k++)
            printf("u[%ld] = %f + %f*I\n", k, creal(w[k]), cimag(w[k]));
        abort();
    }
#if PARAMS
    {
        /* zero estimate */
        double r1 = (1 + r) / 2;
        M = constant_m_d(w, len, r1, i);
        A = prec*LOG2 + log(2*PI*M) + 1;
        rho = r1 + sqrt(r1*r1 - 1);
        B = 2*log(rho);
        printf("r = %lf -> n0 = %lf [M = %lf, A = %lf, B = %lf]\n", r1, A/B, M, A, B);
    }
#endif

    /* first estimate */
    M = constant_m_d(w, len, r, i);
    A = prec*LOG2 + log(2*PI*M) + 1;
    rho = r + sqrt(r*r-1);
    B = 2*log(rho);
#if PARAMS
    printf("r = %lf -> n1 = %lf [M = %lf, A = %lf, B = %lf]\n", r, A/B, M, A, B);
#endif
    n = ceil(A / B);

    /* second = should be exact */
    r *= 1-mult/n;
    M = constant_m_d(w, len, r, i);
    A = prec*LOG2 + log(2*PI*M) + 1;
    rho = r + sqrt(r*r-1);
    B = 2*log(rho);
#if PARAMS
    printf("r = %lf -> n2 = %lf [M = %lf, A = %lf, B = %lf]\n", r, A/B, M, A, B);
#endif
    n = ceil(A / B);

    return n;
}
