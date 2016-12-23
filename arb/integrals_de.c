/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
de_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c,
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
        acb_inv(y, y, prec);

        /* all differentials for x */
        iy = 1; ix = 0;
        acb_set_arb(wy, de->dx + l);
        acb_mul_arb(wy, y, de->dx + l, prec);
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
        mth_root_pol_def(y, u, d1, d, de->x + l, c.m, prec);
        acb_inv(y, y, prec);
        iy = 1; ix = 0;
        acb_set_arb(wy, de->dx + l);
        acb_mul_arb(wy, y, de->dx + l, prec);
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

    arb_clear(x);
    acb_clear(y);
    acb_clear(wy);
    acb_clear(wyx);
}

void
integrals_edge_de(acb_ptr res, sec_t c, edge_t e, const cohom_t dz, const de_int_t de, slong prec)
{
    slong d, d1;
    acb_ptr u, cab, ab2;

    /* reduce roots */
    d = c.d;
    u = _acb_vec_init(d);
    d1 = ab_points(u, c.roots, e, d, prec);
    ab2 = u + d - 2;
    cab = u + d - 1;

    de_integrals(res, u, d1, d - 2, c, dz, de, prec);
    integrals_edge_factors(res, cab, ab2, c, dz, prec);

    _acb_vec_clear(u, d);
}

void
integrals_tree_de(acb_mat_t integrals, sec_t c, const tree_t tree, const cohom_t dz, slong prec)
{
    slong k;
    ulong n;
    arf_t h;
    de_int_t de;

    de_int_params(h, &n, tree, c, prec);
    de_int_init(de, h, n, prec);

    for (k = 0; k < c.d - 1; k++)
        integrals_edge_de(integrals->rows[k], c, tree->e[k], dz, de, prec);

    de_int_clear(de);
}
