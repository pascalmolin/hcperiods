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

        for (i = 1; i < g; i++)
        {
            acb_mul_arb(yxi, yxi, x, prec);
            acb_add(res + i, res + i, yxi, prec);
        }

        continue;
        /* could reuse -x, but not a big deal */

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
    acb_clear(y);
    acb_clear(yxi);
}

void
integrals_edge_gc(acb_ptr res, sec_t c, edge_t e, slong n, slong prec)
{
    slong n1;
    acb_ptr u, cab, ab2;

    /* reduce roots */
    u = _acb_vec_init(c.n);
    n1 = ab_points(u, c.roots, e, c.n, c.m, prec);
    ab2 = u + c.n - 2;
    cab = u + c.n - 1;

    gc_integrals(res, u, n1, c.n - 2, c.g, n, prec);
    integrals_edge_factors_gc(res, cab, ab2, c, prec);

    _acb_vec_clear(u, c.n);
}

void
integrals_tree_gc(acb_mat_t integrals, const data_t data, sec_t c, const tree_t tree, slong prec)
{
    slong k, n;

    n = gc_params_tree(tree, c, prec);

    for (k = 0; k < c.n - 1; k++)
    {
        slong n1;
        acb_ptr res, u, cab, ab2;

        res = integrals->rows[k];
        u = data->upoints->rows[k];
        n1 = data->n1[k];
        ab2 = u + c.n - 2;
        cab = u + c.n - 1;

        gc_integrals(res, u, n1, c.n - 2, c.g, n, prec);
        integrals_edge_factors_gc(res, cab, ab2, c, prec);
    }


}
