/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
gc_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, slong g, slong n, int flag, slong prec)
{
    slong l;
    fmpq_t ln;
    arb_t w, x;
    acb_t y, yxi;
    void (*sqrt_pol) (acb_t y, acb_srcptr u, slong d1, slong d,
            const arb_t x, slong prec);

    fmpq_init(ln);
    arb_init(w);
    arb_init(x);
    acb_init(y);
    acb_init(yxi);
#if DEBUG
    flint_printf("\ngc integral : d1 = %ld, d = %ld, g = %ld, n = %ld, prec = %ld",
            d1, d, g, n, prec);
#endif

    sqrt_pol = &sqrt_pol_def;
    if (prec > 1000 && d > 5)
        sqrt_pol = &sqrt_pol_turn;

    if (flag & AJ_ROOT_DEF)
        sqrt_pol = &sqrt_pol_def;
    else if (flag & AJ_ROOT_TURN)
        sqrt_pol = &sqrt_pol_turn;

    /* compute integral */
    _acb_vec_zero(res, g);

    for (l = 0; 2 * l + 1 < n; l++)
    {

        /* compute x */
        fmpq_set_si(ln, 2 * l + 1, 2 * n);
        arb_cos_pi_fmpq(x, ln, prec);

        /* compute 1/y(x) */
        sqrt_pol(y, u, d1, d, x, prec);
        acb_inv(y, y, prec);

        /* differentials */
        acb_vec_add_geom_arb(res, g, y, x, prec);

        /* now on -x */
        arb_neg(x, x);

        sqrt_pol(y, u, d1, d, x, prec);
        acb_inv(y, y, prec);
        acb_vec_add_geom_arb(res, g, y, x, prec);
    }

    if (n % 2)
    {
        arb_zero(x);
        sqrt_pol(y, u, d1, d, x, prec);
        acb_inv(y, y, prec);
        acb_add(res + 0, res + 0, y, prec);
    }

    /* multiply by weight = Pi / n */
    arb_const_pi(w, prec);
    arb_div_ui(w, w, n, prec);
    _acb_vec_scalar_mul_arb(res, res, g, w, prec);

    fmpq_clear(ln);
    arb_clear(x);
    arb_clear(w);
    acb_clear(y);
    acb_clear(yxi);
}

void
integrals_edge_gc(acb_ptr res, ydata_t ye, sec_t c, slong n, mag_t e, int flag, slong prec)
{
    gc_integrals(res, ye->u, ye->n1, c.n - 2, c.g, n, flag, prec);
    if (0 && e)
        _acb_vec_add_error_mag(res, c.g, e);
    /*integrals_edge_factors_gc(res, ye->ba2, ye->ab, ye->c, c, prec);*/
    integrals_edge_factors(res, ye->ba2, ye->ab, ye->c, c, prec);
}

void
integrals_tree_gc(acb_mat_t integrals, sec_t c, const tree_t tree, int flag, slong prec)
{
    slong k, n;
    mag_t e;

    mag_init(e);
    n = gc_params_tree(e, tree, c, prec);

    for (k = 0; k < tree->n; k++)
        integrals_edge_gc(integrals->rows[k], tree->data + k, c, n, e, flag, prec);

    mag_clear(e);
}
