/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
abel_jacobi_init_roots(abel_jacobi_t aj, slong m, acb_srcptr x, slong n, int flag)
{
    slong g;

    sec_init(&aj->c, m, x, n);
    g = aj->c.g;
    if (g < 1)
    {
        flint_printf("\nno periods if genus < 1 (g = %ld)\n", g);
        abort();
    }

    aj->type =  (flag & AJ_USE_DE || (m > 2 && n > 2)) ? INT_DE : (m == 2) ? INT_GC : INT_D2;

    tree_init(aj->tree, n - 1 - (aj->c.delta == m));

    homol_init(&aj->loop_a, g);
    homol_init(&aj->loop_b, g);

    acb_mat_init(aj->omega0, g, g);
    acb_mat_init(aj->omega1, g, g);
    acb_mat_init(aj->tau, g, g);

    arb_mat_init(aj->proj, 2*g, 2*g);

    /* abel-jacobi map */
    aj->p0 = 0;
}

void
abel_jacobi_init_poly(abel_jacobi_t aj, slong m, acb_poly_t f, int flag, slong prec)
{
    slong n;
    acb_ptr x;
    n = acb_poly_degree(f);
    x = _acb_vec_init(n);
    /* isolate roots */
    if (acb_poly_find_roots(x, f, NULL, 0, prec) < n)
    {
        flint_printf("missing roots, abort.\n");
        abort();
    }
    /* set real roots to exact real */
    if (acb_poly_validate_real_roots(x, f, prec))
    {
        slong k;
        for (k = 0; k < n; k++)
            if (arb_contains_zero(acb_imagref(x + k)))
                arb_zero(acb_imagref(x + k));
    }
    /* and order them */
    _acb_vec_sort_lex(x, n);
    abel_jacobi_init_roots(aj, m, x, n, flag);
    _acb_vec_clear(x, n);
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
}

void
abel_jacobi_compute(abel_jacobi_t aj, int flag, slong prec)
{
    sec_t c = aj->c;
    acb_mat_t integrals;
#if DEBUG
    flint_printf("\nuse points X = ");
    _acb_vec_printd(aj->c.roots, c.n, 30, "\n");
#endif

    /* homology */
    progress("## spanning tree\n");
    spanning_tree(aj->tree, c.roots, c.n, aj->type);
#if DEBUG
    tree_print(aj->tree);
#endif

    tree_ydata_init(aj->tree, c.roots, c.n, c.m, prec);

    progress("## symplectic basis\n");
    symplectic_basis(aj->loop_a, aj->loop_b, aj->tree, c);

    if (flag & AJ_NO_AB)
        return;

    /* integration */
    progress("## integrals %s\n", (aj->type == INT_GC) ? "gc" : "de");
    acb_mat_init(integrals, aj->tree->n, c.g);
    if (aj->type == INT_D2)
        integral_d2(integrals->rows[0], aj->tree->data + 0, c, prec);
    else if (aj->type == INT_GC)
        integrals_tree_gc(integrals, c, aj->tree, prec);
    else
        integrals_tree_de(integrals, c, aj->tree, prec);
    tree_ydata_clear(aj->tree);
#if DEBUG
    flint_printf("\n\ntree integrals\n");
    acb_mat_printd(integrals, 30);
    flint_printf("\n\n");
#endif

    /* period matrices */
    progress("## periods\n");
    period_matrix(aj->omega0, aj->loop_a, integrals, c, prec);
    period_matrix(aj->omega1, aj->loop_b, integrals, c, prec);
#if DEBUG > 1
    progress("\n\nperiods A\n");
    acb_mat_printd(aj->omega0, 20);
    progress("\n\nperiods B\n");
    acb_mat_printd(aj->omega1, 20);
#endif

    acb_mat_clear(integrals);

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
