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

    aj->type = (m == 2) ? INT_GC : INT_DE;

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
    acb_mat_t u;
    acb_mat_t integrals;
    slong * d1;

    /* homology */
    spanning_tree(aj->tree, c.roots, c.d, aj->type);

    acb_mat_init(u, c.d-1, c.d);
    d1 = flint_malloc((c.d-1) * sizeof(slong));
    ab_points_tree(u, d1, aj->tree, c, prec);

    symplectic_basis(aj->loop_a, aj->loop_b, aj->tree, c);

    /* cohomology */
    holomorphic_differentials(aj->dz, c.d, c.m);

    /* integration */
    acb_mat_init(integrals, c.d-1, c.g);
    if (aj->type == INT_GC)
        integrals_tree_gc(integrals, c, aj->tree, prec);
    else
        integrals_tree_de(integrals, c, aj->tree, aj->dz, prec);

    /* period matrices */
    period_matrix(aj->omega0, aj->loop_a, aj->dz, integrals, c, prec);
    period_matrix(aj->omega1, aj->loop_b, aj->dz, integrals, c, prec);

    acb_mat_clear(integrals);
    acb_mat_clear(u);
    flint_free(d1);

    tau_matrix(aj->tau, aj->omega0, aj->omega1, prec);

    /* abel-jacobi map */
    aj->p0 = 0;
}
