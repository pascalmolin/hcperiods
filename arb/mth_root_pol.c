/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
sqrt_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_one(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k < d1)
            acb_neg(t, t);
        acb_sqrt(t, t, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
mth_root_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_one(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k < d1)
            acb_neg(t, t);
        acb_root_ui(t, t, m, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
mth_root_pol_prod(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_zero(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k < d1)
            acb_neg(t, t);
        acb_log(t, t, prec + 4);
        acb_add(y, y, t, prec + 4);
    }
    acb_div_ui(y, y, m, prec + 4);
    acb_exp(y, y, prec);

    acb_clear(t);
}

void
mth_root_pol_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec)
{
    slong q; /* integer mod 2m */

    /* make and product and update q */

    /* multiply by (z^1/2)^q */
    if (q)
        acb_mul(y, y, z + q, prec);
}
