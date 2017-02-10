/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* returns x0, c.x0 + x1, c^2.x0 + 2.c.x1 + x2, ... */
void
acb_vec_polynomial_shift(acb_ptr x, slong len, const acb_t c, slong prec)
{
    slong k, l;
    for (k = 1; k < len; k++)
        for (l = len - 1; l >= k; l--)
            acb_addmul(x + l, x + l - 1, c, prec);
}

/* returns c0.x0, c0.c.x1, c0.c^2.x2, ... */
void
acb_vec_mul_geom(acb_ptr x, slong len, acb_t c0, const acb_t c, slong prec)
{
    slong k;
    acb_mul(x + 0, x + 0, c0, prec);
    for (k = 1; k < len; k++)
    {
        acb_mul(c0, c0, c, prec);
        acb_mul(x + k, x + k, c0, prec);
    }
}
/* returns c0+x0, c0+c.x1, c0+c^2.x2, ... */
void
acb_vec_add_geom_arb(acb_ptr x, slong len, acb_t c0, const arb_t c, slong prec)
{
    slong k;
    acb_add(x + 0, x + 0, c0, prec);
    for (k = 1; k < len; k++)
    {
        acb_mul_arb(c0, c0, c, prec);
        acb_add(x + k, x + k, c0, prec);
    }
}
/* returns x0+c0, x1-c0.c, x2+c0.c^2,... */
void
acb_vec_sub_geom_arb(acb_ptr x, slong len, acb_t c0, const arb_t c, slong prec)
{
    slong k;
    acb_add(x + 0, x + 0, c0, prec);
    for (k = 1; k < len; k++)
    {
        acb_mul_arb(c0, c0, c, prec);
        if (k % 2)
                acb_sub(x + k, x + k, c0, prec);
        else
                acb_add(x + k, x + k, c0, prec);
    }
}

void
integrals_edge_factors_gc(acb_ptr res, const acb_t cab, const acb_t ba2, sec_t c, slong prec)
{
    slong i;
    acb_t cj, ci;

    acb_init(cj);
    acb_init(ci);

    /* polynomial shift */
    acb_vec_polynomial_shift(res, c.g, cab, prec);

    /* constants cj, j = 1 */
    /* c_1 = (1-zeta^-1) ba2^(-d/2) (-I)^i
     *     = 2 / ba2^(d/2) */

    acb_pow_ui(cj, ba2, c.n / 2, prec);
    if (c.n % 2)
    {
        acb_t t;
        acb_init(t);
        acb_sqrt(t, ba2, prec);
        acb_mul(cj, cj, t, prec);
        acb_clear(t);
    }
    acb_inv(cj, cj, prec);
    acb_mul_2exp_si(cj, cj, 1);

    _acb_vec_scalar_mul(res, res, c.g, cj, prec);

    /* constant ci = -I * ba2*/
    acb_one(ci);
    for (i = 1; i < c.g; i++)
    {
        acb_mul(ci, ci, ba2, prec);
        acb_div_onei(ci, ci);
        acb_mul(res + i, res + i, ci, prec);
    }

    acb_clear(ci);
    acb_clear(cj);
}

void
integrals_edge_factors(acb_ptr res, const acb_t cab, const acb_t ba2, sec_t c, const cohom_t dz, slong prec)
{
    slong j;
    acb_t z, zj, cj;

    acb_init(z);
    acb_init(zj);
    acb_init(cj);

    acb_unit_root(z, c.m, prec);

    for (j = jmin(c.m, c.n, c.delta); j < c.m; j++)
    {
        slong ni = imax(j, c.m, c.n, c.delta);

        /* polynomial shift */
        acb_vec_polynomial_shift(res, ni, cab, prec);

        /* constant cj */
        acb_pow_ui(zj, z, j, prec);
        acb_sub_si(zj, zj, -1, prec);
        acb_neg(zj, zj);
        acb_pow_ui(cj, cab, j, prec);
        acb_div(cj, zj, cj, prec);

        _acb_vec_scalar_mul(res, res, ni, cj, prec);

        /* constant ci */
        acb_set(cj, ba2);
        acb_vec_mul_geom(res, ni, cj, ba2, prec);
        res += ni;
    }
    acb_clear(z);
    acb_clear(zj);
    acb_clear(cj);
}
