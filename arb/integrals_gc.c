/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

static void
integrals_edge_gc(acb_ptr res, sec_t c, edge_t e, const cohom_t dz,
        slong n, slong prec)
{
    slong k, l, d = c.d;
    fmpq_t ln;
    arb_t w, x;
    acb_t y, yj, yjxi;
    acb_ptr u;

    u = _acb_vec_init(d - 2);
    fmpq_init(ln);
    arb_init(w);
    arb_init(x);
    acb_init(y);
    acb_init(yj);
    acb_init(yjxi);

    /* reduce roots */
    ab_points(u, c.roots, e, d, prec);

    /* compute integral */
    _acb_vec_zero(res, c.g);
    for (l = 0; l < n; l++)
    {
        fmpq_t ln;
        slong i, j;

        /* compute x */
        fmpq_set_si(ln, 2 * l - 1, 2 * n);
        arb_cos_pi_fmpq(x, ln, prec);

        /* compute 1/y(x) */
        nth_root_pol_def(y, u, x, d - 2, 2, prec);
        acb_inv(y, y, prec);

        /* differentials, jd + mi >= delta */
        acb_one(yj);
        for (k = 0, j = 1; j < c.m; j++)
        {
            slong lim = j * c.d - c.delta;

            acb_mul(yj, yj, y, prec);

            if (c.m > lim)
                continue;
            
            acb_set(yjxi, yj);
            acb_add(res + k, res + k, yjxi, prec);
            k++;

            for (i = 2; i < c.d && c.m * i < j * c.d - c.delta; i++)
            {
                acb_mul_arb(yjxi, yjxi, x, prec);
                acb_add(res + k, res + k, yjxi, prec);
                k++;
            }
        }

        if (l == 0)
            continue;

        /* now on -x */
        arb_neg(x, x);

        nth_root_pol_def(y, u, x, d - 2, 2, prec);
        acb_inv(y, y, prec);

        /* differentials, jd + mi >= delta */
        acb_one(yj);
        for (k = 0, j = 1; j < c.m; j++)
        {
            slong lim = j * d - c.delta;

            acb_mul(yj, yj, y, prec);

            if (c.m > lim)
                continue;
            
            acb_set(yjxi, yj);
            acb_add(res + k, res + k, yjxi, prec);
            k++;

            for (i = 2; i < c.d && c.m * i < j * c.d - c.delta; i++)
            {
                acb_mul_arb(yjxi, yjxi, x, prec);
                if (i % 2)
                    acb_add(res + k, res + k, yjxi, prec);
                else
                    acb_sub(res + k, res + k, yjxi, prec);
                k++;
            }
        }
    }

    /* multiply by weight = Pi / n */

    arb_const_pi(w, prec);
    arb_div_ui(w, w, n, prec);

    for (k = 0; k < c.g; k++)
        acb_mul_arb(res + k, res + k, w, prec);

    _acb_vec_clear(u, d - 2);
    fmpq_clear(ln);
    arb_clear(x);
    acb_clear(y);
    acb_clear(yj);
    acb_clear(yjxi);
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
