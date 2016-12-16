/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

static slong
count_coeffs(fmpz * row, slong len)
{
    slong j, c = 0;
    for (j = 0; j < len; j++)
        if (!fmpz_is_zero(row + j))
            c++;
    return c;
}

static void
set_loops(gamma_t * l, slong n, fmpz * row, slong m1, slong len)
{
    slong i, j;
    /* coeff k * (m-1) + s -> gamma_k^s */
    for (i = 0, j = 0; i < n; i++)
    {
        for (; fmpz_is_zero(row + j); j++);
        if (!fmpz_fits_si(row + j))
        {
            flint_printf("error: coefficient too large in basis\n");
            abort();
        }
        l[i].coeff = fmpz_get_si(row + j);
        l[i].shift = j % m1;
        l[i].index = j / m1;
    }
}

void
symplectic_basis(homol_t alpha, homol_t beta, const tree_t tree, sec_t c)
{
    slong i, len = (c.d-1)*(c.m-1);
    si_mat_t m, p;
    si_mat_init(m, len, len);
    si_mat_init(p, len, len);
    /* compute big intersection matrix, size len = (d-1)*(m-1) */
    intersection_tree(m, tree, c.d, c.m);
    symplectic_reduction(p, m, c.g, len);
    si_mat_clear(m, len, len);

    for (i = 0; i < c.g; i++)
    {
        slong n;
        fmpz * row;
        /* loops alpha = even coordinates */
        row = p->rows[2 * i];
        n = count_coeffs(row, len);
        alpha[i].n = n;
        alpha[i].l = flint_malloc(n * sizeof(gamma_t));
        set_loops(alpha[i].l, n, row, c.m - 1, len);
        /* loops beta = odd coordinates */
        row = p->rows[2 * i + 1];
        n = count_coeffs(row, len);
        beta[i].n = n;
        beta[i].l = flint_malloc(n * sizeof(gamma_t));
        set_loops(beta[i].l, n, row, c.m - 1, len);
    }
    si_mat_clear(p, len, len);
}
