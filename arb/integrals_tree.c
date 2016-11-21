/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"

static void
powers_x(arb_ptr vx, const arb_t x, const arb_t dx, slong xmax, slong prec)
{
    return;
}

static void
nth_root_pol_affine(acb_t y, acb_srcptr u, const arb_t x, slong d, slong prec)
{
    return;
}

static void
powers_inv_y(acb_ptr vy, const arb_t x, acb_srcptr u, const slong * yvec, slong n, slong prec)
{
    return;
}

static void
integrals_edge(acb_ptr res, edge_t e, const cohom_t dz, slong g,
        acb_srcptr x, slong d, const slong * yvec, slong ny, const de_int_t de, slong prec)
{
    slong k, l;
    arb_ptr xk;
    acb_ptr u, yk;
    xk = _arb_vec_init(d - 1);
    yk = _acb_vec_init(ny);
    u = _acb_vec_init(d - 2);

    /* reduce roots */
    affine_transform(u, x, e, prec);

    /* compute integral */
    _acb_vec_zero(res, g);
    for (l = 0; l < de->n; l++)
    {
        /* powers dx^k/k */
        powers_x(xk, de->x + l, de->dx + l, d - 1, prec);
        /* powers 1/y(x)^k+(-1)^k/y(-x)^k */
        powers_inv_y(yk, de->x + l, u, yvec, ny, prec);
        /* all differentials */ 
        for (k = 0; k < g; k++)
            acb_addmul_arb(res + l, yk + dz[k].y, xk + dz[k].x, prec);
    }
    /* periods a->b of all differentials */
    _arb_vec_clear(xk, d - 1);
    _acb_vec_clear(yk, ny);
    _acb_vec_clear(u, d - 2);
    return;
}

void
integrals_tree(acb_mat_t periods, const tree_t tree, const cohom_t dz, slong g, acb_srcptr x, slong d, slong prec)
{
    slong k;
    slong * yvec, ny;
    de_int_t de;

    /*de_int_init(de, tree->tau, prec);*/

    for (k = 0; k < d - 1; k++)
        integrals_edge(periods->rows[k], tree->e[k], dz, g, x, d, yvec, ny, de, prec);

    de_int_clear(de);
}
