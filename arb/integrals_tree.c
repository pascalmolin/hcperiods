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
