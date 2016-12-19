/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
ab_points(acb_ptr u, acb_srcptr x, edge_t e, slong d, slong prec)
{
    slong k, l;
    acb_t mab, ba;
    acb_init(mab);
    acb_init(ba);
    acb_set(ba, x + e.b);
    acb_sub(ba, ba, x + e.a, prec);
    acb_set(mab, x + e.a);
    acb_add(mab, mab, x + e.b, prec);
    acb_neg(mab, mab);
    for (k = 0, l = 0; k < d; k++)
    {
        if (k == e.a || k == e.b)
            continue;
        acb_mul_2exp_si(u + l, x + k, 1);
        acb_add(u + l, u + l, mab, prec);
        acb_div(u + l, u + l, ba, prec);
        l++;
    }
    acb_clear(mab);
    acb_clear(ba);
}

void
integrals_tree(acb_mat_t integrals, sec_t c, const tree_t tree, const cohom_t dz, slong prec)
{
    if (c.m == 2)
        integrals_tree_gc(integrals, c, tree, dz, prec);
    else
        integrals_tree_de(integrals, c, tree, dz, prec);
}
