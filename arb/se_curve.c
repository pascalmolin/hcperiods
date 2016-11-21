/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"

void
se_curve_init_roots(se_curve_t c, slong n, acb_srcptr x, slong d)
{
    slong k, g;

    g = ((n-1)*(d-1) - n_gcd(n,d) + 1)/ 2;
    c->n = n;
    c->d = d;
    c->g = g;
    c->roots = _acb_vec_init(d);
    for (k = 0; k < d; k++)
        acb_set(c->roots + k, x + k);

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
se_curve_init_poly(se_curve_t c, slong n, acb_srcptr f, slong len, slong prec)
{
    slong d = len - 1;
    acb_ptr x;
    x = _acb_vec_init(d);
    /* isolate roots */
    if (_acb_poly_find_roots(x, f, NULL, len, 0, prec) < d)
    {
        flint_printf("missing roots, abort.\n");
        abort();
    }
    se_curve_init_roots(c, n, x, d);
    _acb_vec_clear(x, d);
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
    acb_mat_t integrals;

    /* homology */
    spanning_tree(c->tree, c->roots, c->d);
    intersection_tree(c->inter, c->tree, c->roots, c->d);
    symplectic_basis(c->loop_a, c->loop_b, c->inter, c->g, c->d, c->n);

    /* cohomology */
    differentials(c->dz, c->d, c->n);

    /* integration */
    acb_mat_init(integrals, c->d-1, c->g);
    integrals_tree(integrals, c->tree, c->dz, c->g, c->roots, c->d, prec);

    /* period matrices */
    period_matrix(c->omega0, c->loop_a, integrals, c->g, c->d, prec);
    period_matrix(c->omega1, c->loop_b, integrals, c->g, c->d, prec);

    acb_mat_clear(integrals);

    tau_matrix(c->tau, c->omega0, c->omega1, prec);

    /* abel-jacobi map */
    c->p0 = 0;
}
