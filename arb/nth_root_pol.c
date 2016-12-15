/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
nth_root_pol_def(acb_t y, acb_srcptr u, const arb_t x, slong d, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_one(y);
    for (k = 0; k < d; k++)
    {
        acb_set_arb(t, x);
        acb_sub(t, t, u + k, prec);
        acb_root_ui(t, t, m, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
nth_root_pol_prod(acb_t y, acb_srcptr u, const arb_t x, slong d, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_zero(y);
    for (k = 0; k < d; k++)
    {
        acb_set_arb(t, x);
        acb_sub(t, t, u + k, prec);
        acb_log(t, t, prec + 4);
        acb_add(y, y, t, prec + 4);
    }
    acb_div_ui(y, y, m, prec + 4);
    acb_exp(y, y, prec);

    acb_clear(t);
}

void
nth_root_pol_turn(acb_t y, acb_srcptr u, const arb_t x, acb_srcptr z, slong d, slong m, slong prec)
{
    slong q; /* integer mod 2m */

    /* make and product and update q */

    /* multiply by (z^1/2)^q */
    if (q)
        acb_mul(y, y, z + q, prec);
}
