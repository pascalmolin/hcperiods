/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

static void
powers_x(arb_ptr vx, const arb_t x, const arb_t dx, slong xmax, slong prec)
{
    return;
}

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
mth_root_pol_affine_product(acb_t y, acb_srcptr u, const arb_t x, slong d, slong m, slong prec)
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

static void
powers_inv_y(acb_ptr vy, const arb_t x, acb_srcptr u, const slong * yvec, slong n, slong m, slong prec)
{
    acb_t yp, ym;

    acb_init(yp);
    acb_init(ym);
    acb_clear(yp);
    acb_clear(ym);
}

static void
integrals_edge(acb_ptr res, sec_t c, edge_t e, const cohom_t dz,
        const slong * yvec, slong ny, const de_int_t de, slong prec)
{
    slong k, l, d = c.d;
    arb_ptr x;
    acb_ptr u, y;
    x = _arb_vec_init(d - 1);
    y = _acb_vec_init(ny);
    u = _acb_vec_init(d - 2);

    /* reduce roots */
    ab_points(u, c.roots, e, d, prec);

    /* compute integral */
    _acb_vec_zero(res, c.g);
    for (l = 0; l < de->n; l++)
    {
        /* powers dx^k/k */
        powers_x(x, de->x + l, de->dx + l, d - 1, prec);
        /* powers 1/y(x)^k+(-1)^k/y(-x)^k */
        powers_inv_y(y, de->x + l, u, yvec, ny, c.m, prec);
        /* all differentials */
        for (k = 0; k < c.g; k++)
            acb_addmul_arb(res + l, y + dz[k].y, x + dz[k].x, prec);
    }
    /* periods a->b of all differentials */
    _arb_vec_clear(x, d - 1);
    _acb_vec_clear(y, ny);
    _acb_vec_clear(u, d - 2);
    return;
}

static slong *
select_ypowers(const cohom_t dz, slong g, slong * len)
{
    slong * yvec, ny, k;

    for (ny = 1, k = 1; k < g; k++)
        if (dz[k].y > dz[k-1].y)
            ny++;

    * len = ny;
    yvec = flint_malloc(ny * sizeof(int));

    for (ny = 0, k = 1; k < g; k++)
        if (dz[k].y > dz[k-1].y)
            yvec[ny++] = dz[k].y;

    return yvec;
}

void
integrals_tree(acb_mat_t periods, sec_t c, const tree_t tree, const cohom_t dz, slong prec)
{
    slong k;
    slong * yvec, ny;
    de_int_t de;

    /*de_int_init(de, tree->tau, prec);*/

    yvec = select_ypowers(dz, c.g, &ny);

    for (k = 0; k < c.d - 1; k++)
        integrals_edge(periods->rows[k], c, tree->e[k], dz, yvec, ny, de, prec);

    flint_free(yvec);

    de_int_clear(de);
}
