/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"

static void
periods_edge(acb_ptr res, edge_t e, const cohom_t dz, slong g,
        acb_srcptr x, slong d, de_int_t de, slong prec)
{
    /* periods a->b of all differentials */
    return; 
}

void
periods_tree(acb_mat_t periods, const tree_t tree, const cohom_t dz, slong g, acb_srcptr x, slong d, slong prec)
{
    slong k;
    de_int_t de;

    /*de_int_init(de, tree->tau, prec);*/

    for (k = 0; k < d - 1; k++)
        periods_edge(periods->rows[k], tree->e[k], dz, g, x, d, de, prec);

    de_int_clear(de);
}
