/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
#include <complex.h>
typedef complex double cdouble;

double
distance_ellipse(cdouble * w, slong i, slong j, slong len)
{
    return 0.;
}

slong
gc_int_params(double r, sec_t c, slong prec)
{
    return 100;
}

void
gc_distance_ellipse(mag_t d, acb_t p, mag_t r)
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
