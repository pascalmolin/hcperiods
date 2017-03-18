/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
gc_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, slong g, slong n, slong prec)
{
    slong i, l;
    fmpq_t ln;
    arb_t w, x;
    acb_t y, yxi;

    fmpq_init(ln);
    arb_init(w);
    arb_init(x);
    acb_init(y);
    acb_init(yxi);
#if DEBUG
    flint_printf("\ngc integral : d1 = %ld, d = %ld, g = %ld, n = %ld, prec = %ld",
            d1, d, g, n, prec);
#endif

    /* compute integral */
    _acb_vec_zero(res, g);
    for (l = 0; l < n; l++)
    {
        fmpq_t ln;

        /* compute x */
        fmpq_init(ln);
        fmpq_set_si(ln, 2 * l + 1, 2 * n);
        arb_cos_pi_fmpq(x, ln, prec);
        fmpq_clear(ln);

        /* compute 1/y(x) */
        mth_root_pol_def(y, u, d1, d, x, 2, prec);
        acb_inv(y, y, prec);

        /* differentials : j = 1 && i < g */
        acb_set(yxi, y);
        acb_add(res + 0, res + 0, yxi, prec);
#if 0 && DEBUG
        flint_printf("\nl = %ld, res[%ld] = ", l, 0);
        acb_printd(res + 0, 20);
#endif

        for (i = 1; i < g; i++)
        {
            acb_mul_arb(yxi, yxi, x, prec);
            acb_add(res + i, res + i, yxi, prec);
#if 0 && DEBUG
            flint_printf("\nl = %ld, res[%ld] = ", l, i);
            acb_printd(res + i, 20);
#endif
        }

        continue;
        /* TODO: could reuse -x, but not a big deal */

        if (l == 0)
            continue;

        /* now on -x */
        arb_neg(x, x);

        mth_root_pol_def(y, u, d1, d, x, 2, prec);
        acb_inv(y, y, prec);

        /* differentials : j = 1 && i < g */
        acb_set(yxi, y);
        acb_add(res + 0, res + 0, yxi, prec);

        for (i = 1; i < g; i++)
        {
            acb_mul_arb(yxi, yxi, x, prec);
            if (i % 2)
                acb_add(res + i, res + i, yxi, prec);
            else
                acb_sub(res + i, res + i, yxi, prec);
        }
    }

    /* multiply by weight = Pi / n */

    arb_const_pi(w, prec);
    arb_div_ui(w, w, n, prec);

    for (i = 0; i < g; i++)
        acb_mul_arb(res + i, res + i, w, prec);

    fmpq_clear(ln);
    arb_clear(x);
    arb_clear(w);
    acb_clear(y);
    acb_clear(yxi);
}

void
integrals_edge_gc(acb_ptr res, ydata_t ye, sec_t c, slong n, mag_t e, slong prec)
{
    gc_integrals(res, ye->u, ye->n1, c.n - 2, c.g, n, prec);
    if (0 && e)
        _acb_vec_add_error_mag(res, c.g, e);
    /*integrals_edge_factors_gc(res, ye->ba2, ye->ab, ye->c, c, prec);*/
    integrals_edge_factors(res, ye->ba2, ye->ab, ye->c, c, prec);
}

void
integrals_tree_gc(acb_mat_t integrals, sec_t c, const tree_t tree, slong prec)
{
    slong k, n;
    mag_t e;

    mag_init(e);
    n = gc_params_tree(e, tree, c, prec);

    for (k = 0; k < tree->n; k++)
        integrals_edge_gc(integrals->rows[k], tree->data + k, c, n, e, prec);

    mag_clear(e);
}
