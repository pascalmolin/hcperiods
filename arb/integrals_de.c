/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
de_integrals_precomp(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c,
        const cohom_t dz, const de_int_t de, slong prec)
{
    slong k, l;
    arb_t x;
    acb_t y, wy, wyx;

    arb_init(x);
    acb_init(y);
    acb_init(wy);
    acb_init(wyx);

    /* compute integral */
    _acb_vec_zero(res, c.g);
    for (l = 0; l < de->n; l++)
    {
        slong ix, iy;

        /* compute 1/y(x) */
        mth_root_pol_def(y, u, d1, d, de->x + l, c.m, prec);
        acb_set_arb(wy, de->ch2m + l);
        acb_div(y, wy, y, prec);

        /* all differentials for x */
        iy = 1; ix = 0;
        acb_mul_arb(wy, y, de->dx + l, prec);
        acb_set(wyx, wy);

        for (k = 0; k < c.g; k++)
        {
            for (; iy < dz[k].y; ix = 0, iy++)
            {
                acb_mul(wy, wy, y, prec);
                acb_set(wyx, wy);
            }
            for (; ix < dz[k].x; ix++)
                acb_mul_arb(wyx, wyx, de->x + l, prec);
            acb_add(res + k, res + k, wyx, prec);

        }

        if (l == 0)
            continue;

        /* now on -x */
        arb_neg(x, de->x + l);
        mth_root_pol_def(y, u, d1, d, x, c.m, prec);
        acb_set_arb(wy, de->ch2m + l);
        acb_div(y, wy, y, prec);

        iy = 1; ix = 0;
        acb_mul_arb(wy, y, de->dx + l, prec);
        acb_set(wyx, wy);
        for (k = 0; k < c.g; k++)
        {
            for (; iy < dz[k].y; ix = 0, iy++)
            {
                acb_mul(wy, wy, y, prec);
                acb_set(wyx, wy);
            }
            for (; ix < dz[k].x; ix++)
                acb_mul_arb(wyx, wyx, de->x + l, prec);
            if (ix % 2)
                acb_sub(res + k, res + k, wyx, prec);
            else
                acb_add(res + k, res + k, wyx, prec);
        }
    }

    _acb_vec_scalar_mul_arb(res, res, c.g, de->factor, prec);

    arb_clear(x);
    acb_clear(y);
    acb_clear(wy);
    acb_clear(wyx);
}

void
de_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c, slong prec)
{
    slong n;
    double h;
    cohom_t dz;
    de_int_t de;
    dz = malloc(c.g * sizeof(dform_t));
    holomorphic_differentials(dz, c.n, c.m);
    n = de_params(&h, u, d, 0, c.n - 1, c.m, prec);
    de_int_init(de, h, n, c.m, prec);
    de_integrals_precomp(res, u, d1, d, c, dz, de, prec);
    de_int_clear(de);
    free(dz);
}

void
integrals_edge_de(acb_ptr res, sec_t c, edge_t e, const cohom_t dz, const de_int_t de, slong prec)
{
    slong n1;
    acb_ptr u, cab, ab2;

    /* reduce roots */
    u = _acb_vec_init(c.n);
    n1 = ab_points(u, c.roots, e, c.n, c.m, prec);
    ab2 = u + c.n - 2;
    cab = u + c.n - 1;

    de_integrals_precomp(res, u, n1, c.n - 2, c, dz, de, prec);
    integrals_edge_factors(res, cab, ab2, c, dz, prec);

    _acb_vec_clear(u, c.n);
}

void
integrals_tree_de(acb_mat_t integrals, const data_t data, sec_t c, const tree_t tree, const cohom_t dz, slong prec)
{
    slong k;
    ulong n;
    double h;
    de_int_t de;

    n = de_params_tree(&h, tree, c, prec);
    de_int_init(de, h, n, c.m, prec);

    for (k = 0; k < c.n - 1; k++)
    {
        slong n1;
        acb_ptr res, u, cab, ab2;

        res = integrals->rows[k];
        u = data->upoints->rows[k];
        n1 = data->n1[k];
        ab2 = u + c.n - 2;
        cab = u + c.n - 1;

        de_integrals_precomp(res, u, n1, c.n - 2, c, dz, de, prec);
        integrals_edge_factors(res, cab, ab2, c, dz, prec);
    }

    de_int_clear(de);
}
