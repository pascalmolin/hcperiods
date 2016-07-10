/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"


void
se_curve_init(se_curve_t c, slong n, acb_ptr x, slong d)
{
    slong k, g;

    g = (n-1)*(d-1) / 2;
    c->n = n;
    c->d = d;
    c->g = g;
    c->roots = x;
    c->roots = _acb_vec_init(d);

    tree_init(c->tree, d);
    c->dz = malloc(g * sizeof(dform_t));

    c->inter = malloc( (d-1) * sizeof (inter_t *) );
    for (k = 0; k < d-1; k++)
        c->inter[k] = malloc( (d-1) * sizeof(inter_t));

    c->loop_a = malloc(g * sizeof(loop_t));
    c->loop_b = malloc(g * sizeof(loop_t));
    /*
    for (k = 0; k < g; k++)
    {
        c->loop_a[k] = malloc( (d-1) * sizeof(int_tree));
        c->loop_b[k] = malloc( (d-1) * sizeof(int_tree));
    }
    */

    acb_mat_init(c->omega0, g, g);
    acb_mat_init(c->omega1, g, g);
    acb_mat_init(c->tau, g, g);

    arb_mat_init(c->proj, 2*g, 2*g);
}

void
se_curve_clear(se_curve_t c)
{
    _acb_vec_clear(c->roots, c->d);
    acb_mat_clear(c->omega0);
    acb_mat_clear(c->omega1);
    acb_mat_clear(c->tau);
    arb_mat_clear(c->proj);
    tree_clear(c->tree);
    free(c->loop_a);
    free(c->loop_b);
    free(c->inter);
    free(c->dz);
}

void
se_curve_compute(se_curve_t c, slong prec)
{
    acb_mat_t periods;

    /* homology */
    spanning_tree(c->tree, c->roots, c->d);
    intersection_tree(c->inter, c->tree, c->roots, c->d);
    symplectic_basis(c->loop_a, c->loop_b, c->inter, c->g, c->d);

    /* cohomology */
    differentials(c->dz, c->d, c->n);

    /* integration */
    acb_mat_init(periods, c->d-1, c->g);
    periods_tree(periods, c->tree, c->dz, c->g, c->roots, c->d, prec);

    /* period matrices */
    period_matrix(c->omega0, c->loop_a, periods, c->g, c->d, prec);
    period_matrix(c->omega1, c->loop_b, periods, c->g, c->d, prec);

    acb_mat_clear(periods);

    tau_matrix(c->tau, c->omega0, c->omega1, prec);

    /* abel-jacobi map */
    c->p0 = 0;
}
