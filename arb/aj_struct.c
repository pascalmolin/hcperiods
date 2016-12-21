/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
abel_jacobi_init_roots(abel_jacobi_t aj, slong m, acb_srcptr x, slong d)
{
    slong g;

    sec_init(&aj->c, m, x, d);
    g = aj->c.g;

    tree_init(aj->tree, d - 1);
    aj->dz = malloc(g * sizeof(dform_t));

    homol_init(&aj->loop_a, g);
    homol_init(&aj->loop_b, g);

    acb_mat_init(aj->omega0, g, g);
    acb_mat_init(aj->omega1, g, g);
    acb_mat_init(aj->tau, g, g);

    arb_mat_init(aj->proj, 2*g, 2*g);
}

void
abel_jacobi_init_poly(abel_jacobi_t aj, slong n, acb_srcptr f, slong len, slong prec)
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
    abel_jacobi_init_roots(aj, n, x, d);
    _acb_vec_clear(x, d);
}

void
abel_jacobi_clear(abel_jacobi_t aj)
{
    sec_clear(aj->c);
    acb_mat_clear(aj->omega0);
    acb_mat_clear(aj->omega1);
    acb_mat_clear(aj->tau);
    arb_mat_clear(aj->proj);
    tree_clear(aj->tree);
    homol_clear(aj->loop_a, aj->c.g);
    homol_clear(aj->loop_b, aj->c.g);
    free(aj->dz);
}

void
abel_jacobi_compute(abel_jacobi_t aj, slong prec)
{
    sec_t c = aj->c;
    acb_mat_t integrals;

    /* homology */
    spanning_tree(aj->tree, c.roots, c.d, (c.m == 2) ? INT_GC : INT_DE);
    flint_printf("spanning tree\n");
    symplectic_basis(aj->loop_a, aj->loop_b, aj->tree, c);
    flint_printf("symplectic basis\n");

    /* cohomology */
    holomorphic_differentials(aj->dz, c.d, c.m);
    flint_printf("differentials\n");

    /* integration */
    acb_mat_init(integrals, c.d-1, c.g);
    integrals_tree(integrals, c, aj->tree, aj->dz, prec);
    flint_printf("integrals\n");

    /* period matrices */
    period_matrix(aj->omega0, aj->loop_a, integrals, c.g, c.d, prec);
    period_matrix(aj->omega1, aj->loop_b, integrals, c.g, c.d, prec);
    flint_printf("periods\n");

    acb_mat_clear(integrals);

    tau_matrix(aj->tau, aj->omega0, aj->omega1, prec);

    /* abel-jacobi map */
    aj->p0 = 0;
}
