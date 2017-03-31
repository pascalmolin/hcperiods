/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
de_integrals_precomp(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c,
        const de_int_t de, int flag, slong prec)
{
    slong l;
    arb_t x;
    acb_t y, wy, wyx;
    acb_ptr z;
    void (*mth_root_pol) (acb_t y, acb_srcptr u, slong d1, slong d,
            const arb_t x, acb_srcptr z, slong m, slong prec);

    arb_init(x);
    acb_init(y);
    acb_init(wy);
    acb_init(wyx);

    z = _acb_vec_init(c.m);
    _acb_vec_unit_roots(z, c.m, prec);

#if DEBUG
    flint_printf("\nde integral, d1=%ld, d=%ld, prec=%ld\n", d1, d, prec);
#endif

    mth_root_pol = &mth_root_pol_def;
    if (prec > 1000 && d > 5)
        mth_root_pol = &mth_root_pol_turn;

    if (flag & AJ_ROOT_DEF)
        mth_root_pol = &mth_root_pol_def;
    else if (flag & AJ_ROOT_PROD)
        mth_root_pol = &mth_root_pol_prod;
    else if (flag & AJ_ROOT_TURN)
        mth_root_pol = &mth_root_pol_turn;

    /* compute integral */
    _acb_vec_zero(res, c.g);
    for (l = 0; l < de->n; l++)
    {
        slong j;
        acb_ptr r;

        /* compute 1/y(x) */

        mth_root_pol(y, u, d1, d, de->x + l, z, c.m, prec);
        acb_set_arb(wy, de->ch2m + l);
        acb_div(y, wy, y, prec);

        /* all differentials for x */
        if (c.j1 > 1)
        {
            acb_pow_ui(wy, y, c.j1, prec);
            acb_mul_arb(wy, wy, de->dx + l, prec);
        }
        else
            acb_mul_arb(wy, y, de->dx + l, prec);

        for (r = res, j = 0; j < c.nj; j++)
        {
            if (j)
                acb_mul(wy, wy, y, prec);
            acb_set(wyx, wy);
            acb_vec_add_geom_arb(r, c.ni[j], wyx, de->x + l, prec);
            r += c.ni[j];
        }
#if DEBUG > 2
        flint_printf("\nl = %ld, res = ", l);
        _acb_vec_printd(res, c.g, 30, "\n");
#endif

        if (l == 0)
            continue;

        /* now on -x */
        arb_neg(x, de->x + l);
        mth_root_pol(y, u, d1, d, x, z, c.m, prec);
        acb_set_arb(wy, de->ch2m + l);
        acb_div(y, wy, y, prec);

        if (c.j1 > 1)
        {
            acb_pow_ui(wy, y, c.j1, prec);
            acb_mul_arb(wy, wy, de->dx + l, prec);
        }
        else
            acb_mul_arb(wy, y, de->dx + l, prec);

        for (r = res, j = 0; j < c.nj; j++)
        {
            if (j)
                acb_mul(wy, wy, y, prec);
            acb_set(wyx, wy);
            acb_vec_sub_geom_arb(r, c.ni[j], wyx, de->x + l, prec);
            r += c.ni[j];
        }
    }

    _acb_vec_scalar_mul_arb(res, res, c.g, de->factor, prec);

#if DEBUG
        flint_printf("\nend integration ");
        _acb_vec_printd(res, c.g, 30, "\n");
#endif

    arb_clear(x);
    acb_clear(y);
    acb_clear(wy);
    acb_clear(wyx);
    _acb_vec_clear(z, c.m);
}

void
de_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c, int flag, slong prec)
{
    slong n;
    arf_t h, l;
    mag_t e;
    de_int_t de;

    arf_init(h);
    arf_init(l);
    mag_init(e);

    n = de_params(e, h, l, u, d, 0, c.n - 1, c.j1, c.m, prec);
    de_int_init(de, h, l, n, e, c.m, prec);
    de_integrals_precomp(res, u, d1, d, c, de, flag, prec);

    de_int_clear(de);
    arf_clear(h);
    arf_clear(l);
    mag_clear(e);
}

void
integrals_edge_de(acb_ptr res, ydata_t ye, sec_t c, const de_int_t de, int flag, slong prec)
{
    de_integrals_precomp(res, ye->u, ye->n1, c.n - 2, c, de, flag, prec);
    if (0 && !mag_is_zero(de->e))
        _acb_vec_add_error_mag(res, c.g, de->e);
    integrals_edge_factors(res, ye->ba2, ye->ab, ye->c, c, prec);
}

void
integrals_tree_de(acb_mat_t integrals, sec_t c, const tree_t tree, int flag, slong prec)
{
    slong k, paramprec = 64;
    ulong n;
    arf_t h, l;
    mag_t e;
    de_int_t de;

    arf_init(h);
    arf_init(l);
    mag_init(e);

    n = de_params_tree(e, h, l, tree, c, paramprec);
#if DEBUG
    flint_printf("\nprecomputed DE, n = %ld, h = %lf, prec=%ld\n", n, h, prec);
#endif
    de_int_init(de, h, l, n, e, c.m, prec);

    for (k = 0; k < tree->n; k++)
        integrals_edge_de(integrals->rows[k], tree->data + k, c, de, flag, prec);

#if DEBUG > 1
    flint_printf("\ntree integrals\n");
    acb_mat_printd(integrals, 30);
#endif

    de_int_clear(de);
    arf_clear(h);
    arf_clear(l);
    mag_clear(e);
}

void
integral_d2(acb_ptr res, ydata_t ye, sec_t c, slong prec)
{
    slong j;
    fmpq_t h, jm;
    arb_t g;
    acb_ptr p = res;

    arb_init(g);
    fmpq_init(h);
    fmpq_set_si(h, 1, 2);
    fmpq_init(jm);
    _acb_vec_zero(p, c.g);
    for (j = c.g ; j > 0; j--, p++)
    {
        arb_ptr r = acb_realref(p);
        arb_const_sqrt_pi(r, prec);
        fmpq_set_si(jm, j, c.m);
        arb_gamma_fmpq(g, jm, prec);
        arb_mul(r, r, g, prec);
        fmpq_add(jm, jm, h);
        arb_gamma_fmpq(g, jm, prec);
        arb_div(r, r, g, prec);
    }
    fmpq_clear(jm);
    fmpq_clear(h);
    arb_clear(g);
    integrals_edge_factors(res, ye->ba2, ye->ab, ye->c, c, prec);
}
