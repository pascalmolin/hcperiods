/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

#define FAIL 2147483647

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
        if (k >= d1)
            acb_neg(t, t);
        acb_sqrt(t, t, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
mth_root_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_one(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k >= d1)
            acb_neg(t, t);
        acb_root_ui(t, t, m, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
mth_root_pol_prod(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_zero(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k >= d1)
            acb_neg(t, t);
        acb_log(t, t, prec + 4);
        acb_add(y, y, t, prec + 4);
    }
    acb_div_ui(y, y, m, prec + 4);
    acb_exp(y, y, prec);

    acb_clear(t);
}

static int
acb_sgn_imag(const acb_t x)
{
    if (arb_contains_zero(acb_imagref(x)))
        return arb_is_nonnegative(acb_realref(x));
    else
        return arf_sgn(arb_midref(acb_imagref(x)));
}

static slong
acb_prod_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong prec)
{
    int isgn_y,  isgn_t;
    slong q = 0, k;
    acb_t t;

    acb_init(t);
    acb_sub_arb(y, u + 0, x, prec);
    if (!d1)
        acb_neg(y, y);

    isgn_y = acb_sgn_imag(y);
    for (k = 1; k < d; k++)
    {
        if (!isgn_y)
        {
#if DEBUG
            flint_printf("[y=%ld]",k);
#endif
            return FAIL;
        }

        acb_sub_arb(t, u + k, x, prec);
        if (k >= d1)
            acb_neg(t, t);

        acb_mul(y, y, t, prec);

        isgn_t = acb_sgn_imag(t);
        if (!isgn_t)
        {
#if DEBUG
            flint_printf("[t=%ld]",k);
#endif
            return FAIL;
        }

        if (isgn_y == isgn_t)
        {
            isgn_y = acb_sgn_imag(y);
            if (isgn_y != isgn_t)
            {
                if (isgn_t > 0)
                    q++;
                else
                    q--;
            }
        }
        else
            isgn_y = acb_sgn_imag(y);
    }
    acb_clear(t);
    return q;
}

void
sqrt_pol_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong prec)
{
    slong q;

    q = acb_prod_turn(y, u, d1, d, x, prec);
    if (q == FAIL)
        sqrt_pol_def(y, u, d1, d, x, prec);
    else
    {
        acb_sqrt(y, y, prec);
        if (q % 2)
            acb_neg(y, y);
    }
}

void
mth_root_pol_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec)
{
    slong q;

    q = acb_prod_turn(y, u, d1, d, x, prec);
    if (q == FAIL)
        mth_root_pol_def(y, u, d1, d, x, z, m, prec);
    else
    {
        acb_root_ui(y, y, m, prec);
        q = (q % m + m) % m;
        if (q)
            acb_mul(y, y, z + q, prec);
    }
}
