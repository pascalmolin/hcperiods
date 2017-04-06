/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
gc_integrals_precomp(acb_ptr res, acb_srcptr u, slong d1, slong d, slong g, const gc_int_t gc, int flag, slong prec)
{
    slong l;
    arb_t w, x;
    acb_t y, yxi;
    void (*sqrt_pol) (acb_t y, acb_srcptr u, slong d1, slong d,
            const arb_t x, slong prec);

    arb_init(w);
    arb_init(x);
    acb_init(y);
    acb_init(yxi);
#if DEBUG
    flint_printf("\ngc integral : d1 = %ld, d = %ld, g = %ld, n = %ld, prec = %ld",
            d1, d, g, gc->n, prec);
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

    for (l = 0; l < gc->len; l++)
    {

        /* compute 1/y(x) */
        sqrt_pol(y, u, d1, d, gc->x + l, prec);
        acb_inv(y, y, prec);

        /* differentials */
        acb_vec_add_geom_arb(res, g, y, gc->x + l, prec);

        /* now on -x */
        arb_neg(x, gc->x + l);

        sqrt_pol(y, u, d1, d, x, prec);
        acb_inv(y, y, prec);
        acb_vec_add_geom_arb(res, g, y, x, prec);
    }

    if (gc->n % 2)
    {
        arb_zero(x);
        sqrt_pol(y, u, d1, d, x, prec);
        acb_inv(y, y, prec);
        acb_add(res + 0, res + 0, y, prec);
    }

    /* multiply by weight = Pi / n */
    arb_const_pi(w, prec);
    arb_div_ui(w, w, gc->n, prec);
    _acb_vec_scalar_mul_arb(res, res, g, w, prec);

    arb_clear(x);
    arb_clear(w);
    acb_clear(y);
    acb_clear(yxi);
}

void
gc_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, slong g, slong n, int flag, slong prec)
{
    gc_int_t gc;
    gc_int_init(gc, n, prec);
    gc_integrals_precomp(res, u, d1, d, g, gc, flag, prec);
    gc_int_clear(gc);
}

void
integrals_edge_gc(acb_ptr res, ydata_t ye, sec_t c, int flag, slong prec)
{
    slong n;
    mag_t e;
    mag_init(e);
    n = gc_params(e, ye->u, c.n - 2, c.g, prec);
    gc_integrals(res, ye->u, ye->n1, c.n - 2, c.g, n, flag, prec);
    if (0 && e)
        _acb_vec_add_error_mag(res, c.g, e);
    /*integrals_edge_factors_gc(res, ye->ba2, ye->ab, ye->c, c, prec);*/
    integrals_edge_factors(res, ye->ba2, ye->ab, ye->c, c, prec);
    mag_clear(e);
}

void
integrals_tree_gc(acb_mat_t integrals, sec_t c, const tree_t tree, int flag, slong prec)
{
    slong k;

    for (k = 0; k < tree->n; k++)
        integrals_edge_gc(integrals->rows[k], tree->data + k, c, flag, prec);
}
