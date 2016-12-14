/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
ab_points(acb_ptr u, acb_srcptr x, edge_t e, slong d, slong prec)
{
    slong k, l;
    acb_t mab, ba;
    acb_init(mab);
    acb_init(ba);
    acb_set(ba, x + e.b);
    acb_sub(ba, ba, x + e.a, prec);
    acb_set(mab, x + e.a);
    acb_add(mab, mab, x + e.b, prec);
    acb_neg(mab, mab);
    for (k = 0, l = 0; k < d; k++)
    {
        if (k == e.a || k == e.b)
            continue;
        acb_mul_2exp_si(u + l, x + k, 1);
        acb_add(u + l, u + l, mab, prec);
        acb_div(u + l, u + l, ba, prec);
        l++;
    }
    acb_clear(mab);
    acb_clear(ba);
}

static void
integrals_edge(acb_ptr res, sec_t c, edge_t e, const cohom_t dz,
        const de_int_t de, slong prec)
{
    slong k, l, d = c.d;
    arb_t x;
    acb_t y, wy, wyx;
    acb_ptr u;

    u = _acb_vec_init(d - 2);
    arb_init(x);
    acb_init(y);
    acb_init(wy);
    acb_init(wyx);

    /* reduce roots */
    ab_points(u, c.roots, e, d, prec);

    /* compute integral */
    _acb_vec_zero(res, c.g);
    for (l = 0; l < de->n; l++)
    {
        slong ix, iy;
        /* compute 1/y(x) */
        nth_root_pol_def(y, u, de->x + l, d - 2, c.m, prec);
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
        nth_root_pol_def(y, u, x, d - 2, c.m, prec);
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
    /* periods a->b of all differentials */
    _acb_vec_clear(u, d - 2);
    arb_clear(x);
    acb_clear(y);
    acb_clear(wy);
    acb_clear(wyx);
}

void
integrals_tree(acb_mat_t integrals, sec_t c, const tree_t tree, const cohom_t dz, slong prec)
{
    slong k;
    ulong n;
    arf_t h;
    de_int_t de;

    de_int_params(h, &n, tree->tau, 1., 1., prec);
    de_int_init(de, h, n, prec);

    for (k = 0; k < c.d - 1; k++)
        integrals_edge(integrals->rows[k], c, tree->e[k], dz, de, prec);

    de_int_clear(de);
}
