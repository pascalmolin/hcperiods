/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

static void
integrals_edge_gc(acb_ptr res, sec_t c, edge_t e, const cohom_t dz,
        slong n, slong prec)
{
    slong i, l, d = c.d;
    fmpq_t ln;
    arb_t w, x;
    acb_t y, yxi;
    acb_ptr u;

    u = _acb_vec_init(d - 2);
    fmpq_init(ln);
    arb_init(w);
    arb_init(x);
    acb_init(y);
    acb_init(yxi);

    /* reduce roots */
    ab_points(u, c.roots, e, d, prec);

    /* compute integral */
    _acb_vec_zero(res, c.g);
    for (l = 0; l < n; l++)
    {
        fmpq_t ln;

        /* compute x */
        fmpq_init(ln);
        fmpq_set_si(ln, 2 * l - 1, 2 * n);
        arb_cos_pi_fmpq(x, ln, prec);
        fmpq_clear(ln);

        /* compute 1/y(x) */
        nth_root_pol_def(y, u, x, d - 2, 2, prec);
        acb_inv(y, y, prec);

        /* differentials : j = 1 && i < g */
        acb_set(yxi, y);
        acb_add(res + 0, res + 0, yxi, prec);

        for (i = 1; i < c.g; i++)
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

        nth_root_pol_def(y, u, x, d - 2, 2, prec);
        acb_inv(y, y, prec);

        /* differentials : j = 1 && i < g */
        acb_set(yxi, y);
        acb_add(res + 0, res + 0, yxi, prec);

        for (i = 1; i < c.g; i++)
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

    for (i = 0; i < c.g; i++)
        acb_mul_arb(res + i, res + i, w, prec);

    _acb_vec_clear(u, d - 2);
    fmpq_clear(ln);
    arb_clear(x);
    acb_clear(y);
    acb_clear(yxi);
}

void
integrals_tree_gc(acb_mat_t integrals, sec_t c, const tree_t tree, const cohom_t dz, slong prec)
{
    slong k;
    ulong n;

    n = gc_int_params(tree, c, prec);

    for (k = 0; k < c.d - 1; k++)
        integrals_edge_gc(integrals->rows[k], c, tree->e[k], dz, n, prec);

}
