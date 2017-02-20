/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
de_integrals_precomp(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c,
        const cohom_t dz, const de_int_t de, slong prec)
{
    slong l;
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
        slong j;
        acb_ptr r;

        /* compute 1/y(x) */
        mth_root_pol_def(y, u, d1, d, de->x + l, c.m, prec);
        acb_set_arb(wy, de->ch2m + l);
        acb_div(y, wy, y, prec);

        /* all differentials for x */
        j = jmin(c.m, c.n, c.delta);
        if (j > 1)
        {
            acb_pow_ui(wy, y, j, prec);
            acb_mul_arb(wy, wy, de->dx + l, prec);
        }
        else
            acb_mul_arb(wy, y, de->dx + l, prec);

        for (r = res; j < c.m; j++)
        {
            slong ni = imax(j, c.m, c.n, c.delta);
            acb_set(wyx, wy);
            acb_vec_add_geom_arb(r, ni, wyx, de->x + l, prec);
            r += ni;
        }

        if (l == 0)
            continue;

        /* now on -x */
        arb_neg(x, de->x + l);
        mth_root_pol_def(y, u, d1, d, x, c.m, prec);
        acb_set_arb(wy, de->ch2m + l);
        acb_div(y, wy, y, prec);

        j = jmin(c.m, c.n, c.delta);
        if (j > 1)
        {
            acb_pow_ui(wy, y, j, prec);
            acb_mul_arb(wy, wy, de->dx + l, prec);
        }
        else
            acb_mul_arb(wy, y, de->dx + l, prec);

        for (r = res; j < c.m; j++)
        {
            slong ni = imax(j, c.m, c.n, c.delta);
            acb_set(wyx, wy);
            acb_vec_sub_geom_arb(r, ni, wyx, de->x + l, prec);
            r += ni;
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
integrals_edge_de(acb_ptr res, ydata_t ye, sec_t c, const cohom_t dz, const de_int_t de, slong prec)
{
    de_integrals_precomp(res, ye->u, ye->n1, c.n - 2, c, dz, de, prec);
    integrals_edge_factors(res, ye->c, ye->ba2, c, dz, prec);
}

void
integrals_tree_de(acb_mat_t integrals, sec_t c, const tree_t tree, const cohom_t dz, slong prec)
{
    slong k;
    ulong n;
    double h;
    de_int_t de;

    n = de_params_tree(&h, tree, c, prec);
    de_int_init(de, h, n, c.m, prec);

    for (k = 0; k < c.n - 1; k++)
        integrals_edge_de(integrals->rows[k], tree->data + k, c, dz, de, prec);

    de_int_clear(de);
}
