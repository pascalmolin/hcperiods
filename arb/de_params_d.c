/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* rough estimates using doubles */
double
de_distance_d(double x, double y, double r)
{
    /* FIXME: this is not the real distance... */
    cdouble xItau;
    xItau = casinh(catanh(x+I*y)/LAMBDA);
    return
    (fabs(cimag(xItau))-r) * LAMBDA * cabs(ccosh(x+I*y)) / pow(cabs(ccosh(LAMBDA * csinh(x+I*y))), 2);
}

double
de_constant_m_d(const cdouble * w, slong len, double r, slong d, slong m)
{
    /* 2*Pi*r^d/sqrt(prod d(w_i,e_r)) */
    slong k;
    double p;
    const double eps = 1.e-5;
    p = 1;
    for (k = 0; k < len; k++)
    {
        double d = de_distance_d(creal(w[k]), cimag(w[k]), r);
        p *= (d < eps) ? r : d;
    }
    return pow(r, d) / pow(p, 1/m);
}

static double
de_constant_b(double r, double alpha)
{
    double xr, cr, y0, a2, b;
    a2 = 2*alpha;
    y0 = LAMBDA * sin(r);
    cr = cos(r);
    xr = cr*sqrt(PI2/y0-1);
    b = xr / 2 * (pow(cos(y0), -a2) + pow(xr, -a2))
        + pow(sinh(xr), -a2) / a2;
    return 2/cr*b;
}

slong
de_params_d(double * h, double * lambda, double * r, const cdouble * w, slong len, slong i, slong m, slong prec)
{
    slong n;
    double D, M1, M2, B, alpha;

    *lambda = LAMBDA;

    D = prec * LOG2;
    alpha = (double)(m-1) / (double)m;
    M1 = 1;
    M2 = 1;
    if (*r <= 0)
    {
        slong k;
        *r = 1.4;
        for (k = 0; k < len; k++)
        {
            double rk = fabs(cimag(casinh(catanh(w[k])/LAMBDA)));
            if (rk < *r)
                *r = rk;
        }
        *r *= .9;
        //flint_printf("### r was not set, choose r = %lf\n",r);
    }
    B = de_constant_b(*r, alpha);

    *h = 2*PI*(*r) / (D + log(2*M2*B+1));
    n = (slong)(ceil(asinh((D+log(32*M1/alpha))/(2*alpha*LAMBDA)) / *h));

    return n;
}
