/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* u[0..l1[ contains roots re(ui)<=0
   u[l1..d-2[ roots with re(ui) > 0
   the last two components are set to
   (b-a)/2 and (a+b)/(b-a)
   returns l1
*/
slong
ab_points(acb_ptr u, acb_srcptr x, edge_t e, slong d, slong prec)
{
    slong k, l;
    acb_t ab, ba; /* a + b and b - a */

    acb_init(ab);
    acb_init(ba);

    acb_set(ba, x + e.b);
    acb_sub(ba, ba, x + e.a, prec);
    acb_set(ab, x + e.a);
    acb_add(ab, ba, x + e.b, prec);

    for (k = 0, l = 0; k < d; k++)
    {
        if (k == e.a || k == e.b)
            continue;
        acb_mul_2exp_si(u + l, x + k, 1);
        acb_sub(u + l, u + l, ab, prec);
        acb_div(u + l, u + l, ba, prec);
        l++;
    }

    /* now l = d - 2, set last two */

    acb_mul_2exp_si(u + l, ba, -1);
    acb_div(u + l + 1, ab, ba, prec);

    /* reorder */
    for (k = 0; k < l; k++)
        if (arb_is_positive(acb_realref(u + k)))
            acb_swap(u + k--, u + l--);
    acb_clear(ab);
    acb_clear(ba);

    return l;
}

/* returns x0, c.x0 + x1, c^2.x0 + c.x1 + x2, ... */
void
acb_vec_polynomial_shift(acb_ptr x, const acb_t c, slong len, slong prec)
{
    slong k;
    for (k = 1; k < len; k++)
        acb_addmul(x + k, x + k - 1, c, prec);
}

void
integrals_edge_factors_gc(acb_ptr res, const acb_t cab, const acb_t ba2, sec_t c, slong prec)
{
    slong i;
    acb_t cj, ci;

    acb_init(cj);
    acb_init(ci);

    /* polynomial shift */
    acb_vec_polynomial_shift(res, cab, c.g, prec);

    /* constants cj, j = 1 */
    /* c_1 = (1-zeta^-1) ba2^(-d/2) (-I)^i
     *     = 2 / ba2^(d/2) */

    acb_pow_ui(cj, ba2, c.d / 2, prec);
    if (c.d % 2)
    {
        acb_t t;
        acb_init(t);
        acb_sqrt(t, ba2, prec);
        acb_mul(cj, cj, t, prec);
        acb_clear(t);
    }
    acb_inv(cj, cj, prec);
    acb_mul_2exp_si(cj, cj, 1);

    _acb_vec_scalar_mul(res, res, c.g, cj, prec);

    /* constant ci = -I * ba2*/
    acb_one(ci);
    for (i = 1; i < c.g; i++)
    {
        acb_mul(ci, ci, ba2, prec);
        acb_div_onei(ci, ci);
        acb_mul(res + i, res + i, ci, prec);
    }

    acb_clear(ci);
    acb_clear(cj);
}

void
integrals_edge_factors(acb_ptr res, const acb_t cab, const acb_t ba2, sec_t c, const cohom_t dz, slong prec)
{
    return;
}
