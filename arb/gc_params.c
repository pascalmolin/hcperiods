/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

slong
gc_params(acb_srcptr u, slong len, double r, slong i, slong prec)
{
    slong n, k;
    cdouble * w;
    w = malloc(len * sizeof(cdouble));
    for (k = 0; k < len; k++)
        w[k] = acb_get_cdouble(u + k);

    n = gc_params_d(w, len, r, i, prec);

    free(w);
    return n;
}

slong
gc_params_tree(const tree_t tree, sec_t c, slong prec)
{
    return gc_params(tree->data[tree->min].u, c.n, tree->r, c.n - 1, prec);
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
