/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
void
acb_branch_points(acb_ptr x, slong n, const acb_poly_t pol, slong prec)
{
    /* isolate roots */
    if (acb_poly_find_roots(x, pol, NULL, 0, prec) < n)
    {
        acb_poly_printd(pol, 20);
        flint_printf("missing roots, abort.\n");
        abort();
    }
    /* set real roots to exact real */
    if (acb_poly_validate_real_roots(x, pol, prec))
    {
        slong k;
        for (k = 0; k < n; k++)
            if (arb_contains_zero(acb_imagref(x + k)))
                arb_zero(acb_imagref(x + k));
    }
    /* and order them */
    _acb_vec_sort_lex(x, n);
}

void
abel_jacobi_init_poly(abel_jacobi_t aj, slong m, const acb_poly_t pol)
{
    slong n, g;

    sec_init_poly(&aj->c, m, pol);
    n = aj->c.n;
    g = aj->c.g;
    if (g < 1)
    {
        flint_printf("\nno periods if genus < 1 (g = %ld)\n", g);
        abort();
    }

    aj->type = 0;

    aj->roots = _acb_vec_init(n);

    tree_init(aj->tree, n - 1 - (aj->c.delta == m));

    homol_init(&aj->loop_a, g);
    homol_init(&aj->loop_b, g);

    acb_mat_init(aj->integrals, aj->tree->n, g);
    acb_mat_init(aj->omega0, g, g);
    acb_mat_init(aj->omega1, g, g);
    acb_mat_init(aj->tau, g, g);

    arb_mat_init(aj->proj, 2*g, 2*g);

    /* abel-jacobi map */
    aj->p0 = 0;
}

void
abel_jacobi_clear(abel_jacobi_t aj)
{
    _acb_vec_clear(aj->roots, aj->c.n);
    acb_mat_clear(aj->integrals);
    acb_mat_clear(aj->omega0);
    acb_mat_clear(aj->omega1);
    acb_mat_clear(aj->tau);
    arb_mat_clear(aj->proj);
    tree_clear(aj->tree);
    homol_clear(aj->loop_a, aj->c.g);
    homol_clear(aj->loop_b, aj->c.g);
    sec_clear(aj->c);
}

void
abel_jacobi_compute(abel_jacobi_t aj, int flag, slong prec)
{
    sec_t c = aj->c;

    aj->type = (flag & AJ_USE_DE || (c.m > 2 && c.n > 2)) ? INT_DE : (c.m == 2) ? INT_GC : INT_D2;
#if DEBUG
    acb_poly_printd(c.pol, 20);
#endif

    /* branch points */
    progress("## branch points\n");
    acb_branch_points(aj->roots, c.n, c.pol, prec);
#if DEBUG
    flint_printf("\nuse points X = ");
    _acb_vec_printd(aj->roots, c.n, 30, "\n");
#endif

    /* homology */
    progress("## spanning tree\n");
    spanning_tree(aj->tree, aj->roots, c.n, aj->type);
#if DEBUG
    tree_print(aj->tree);
#endif

    tree_ydata_init(aj->tree, aj->roots, c.n, c.m, prec);

    progress("## symplectic basis\n");
    symplectic_basis(aj->loop_a, aj->loop_b, aj->tree, c);

    if (flag & AJ_NO_INT)
        return;

    /* integration */
    progress("## integrals %s\n", (aj->type == INT_GC) ? "gc" : "de");
    acb_mat_init(aj->integrals, aj->tree->n, c.g);
    if (aj->type == INT_D2)
        integral_d2(aj->integrals->rows[0], aj->tree->data + 0, c, prec);
    else if (aj->type == INT_GC)
        integrals_tree_gc(aj->integrals, c, aj->tree, flag, prec);
    else
        integrals_tree_de(aj->integrals, c, aj->tree, flag, prec);
    tree_ydata_clear(aj->tree);

#if DEBUG > 2
    flint_printf("\n\ntree integrals\n");
    acb_mat_printd(aj->integrals, 30);
    flint_printf("\n\n");
#endif
    if (flag & AJ_NO_AB)
        return;

    /* period matrices */
    progress("## periods\n");
    period_matrix(aj->omega0, aj->loop_a, aj->integrals, c, prec);
    period_matrix(aj->omega1, aj->loop_b, aj->integrals, c, prec);
#if DEBUG > 2
    progress("\n\nperiods A\n");
    acb_mat_printd(aj->omega0, 20);
    progress("\n\nperiods B\n");
    acb_mat_printd(aj->omega1, 20);
#endif

    if (flag & AJ_TRIM)
    {
        acb_mat_trim(aj->omega0);
        acb_mat_trim(aj->omega1);
    }


    if (flag & AJ_NO_TAU)
        return;

    progress("## tau\n");
    tau_matrix(aj->tau, aj->omega0, aj->omega1, prec);

    if (flag & AJ_TRIM)
        acb_mat_trim(aj->tau);
}
