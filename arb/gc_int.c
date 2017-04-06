/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* cos((2k+1)/(2n)),0 <= k < len, len = n / 2 */
static void
_arb_vec_cos_pi(arb_ptr c, slong len, slong n, slong prec)
{
    slong k, s, wp;
    acb_ptr x;
    acb_t t;

    x = _acb_vec_init(len + 1);
    wp = prec + 6 + 2 * FLINT_BIT_COUNT(len + 1);

    acb_init(t);
    acb_unit_root(t, 4 * n, wp);
    _acb_vec_set_powers(x, t, len + 1, wp);
    acb_clear(t);

    s = n % 2;

    for (k = 0; k < len / 2; k++)
    {
        arb_set_round(c + k, acb_realref(x + 2 * k + 1), prec);
        arb_set_round(c + len - k - 1, acb_imagref(x + 2 * k + 1 + s), prec);
    }
    if (n % 4 == 2)
    {
        arb_sqrt_ui(c + (n / 4), 2, prec);
        arb_mul_2exp_si(c + (n / 4), c + (n / 4), -1);
    }
    else if (n % 4 == 3)
    {
        arb_set_round(c + (n / 4), acb_realref(x + (n / 2)), prec);
    }

    _acb_vec_clear(x, len);
}

void
gc_int_init(gc_int_t gc, slong n, slong prec)
{
    gc->n = n;
    gc->len = n / 2;
    gc->x = _arb_vec_init(gc->len);
    _arb_vec_cos_pi(gc->x, gc->len, n, prec);
}

void
gc_int_clear(gc_int_t gc)
{
    _arb_vec_clear(gc->x, gc->len);
}
