/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"


void compute_periods_tree(acb_ptr res, const tree_edge[], acb_srcptr x, slong len)
void compute_spanning_tree(tree_edge[], acb_srcptr x, slong d)

void
se_curve_init(se_curve_t c, slong n, acb_ptr x, slong d)
{
    slong k, g;

    g = (n-1)*(d-1) / 2;
    c->n = n;
    c->d = d;
    c->g = g;
    c->roots = acb_vec_init_copy(x, d);

    c->dx = malloc(g * sizeof(slong));
    c->dy = malloc(g * sizeof(slong));

    c->tree = malloc((d-1) * sizeof(slong[2]));

    c->inter = malloc( (d-1) * sizeof (* int_tree) );
    for (k = 0; k < d-1; k++)
        c->inter[k] = malloc( (d-1) * sizeof(int_tree));

    c->ABtoC = malloc( 2 * g * sizeof(loop));
    for (k = 0; k < 2 * g; k++)
        c->ABtoC[k] = malloc( (d-1) * sizeof(int_tree));

    acb_mat_init(c->bigperiods, g, 2 * g);
    acb_mat_init(c->tau, g, g);

}
void
se_curve_clear(se_curve_t c)
{
    _acb_vec_clear(c->roots);
}

