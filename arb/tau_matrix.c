/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include <perm.h>
#include <acb_poly.h>
#include "abel_jacobi.h"

static void
acb_mat_real(arb_mat_t r, const acb_mat_t x)
{
    slong i, j;
    for (i = 0; i < acb_mat_nrows(x); i++)
        for (j = 0; j < acb_mat_ncols(x); j++)
            arb_set(arb_mat_entry(r, i, j), acb_realref(acb_mat_entry(x, i, j)));
}

static void
acb_mat_imag(arb_mat_t r, const acb_mat_t x)
{
    slong i, j;
    for (i = 0; i < acb_mat_nrows(x); i++)
        for (j = 0; j < acb_mat_ncols(x); j++)
            arb_set(arb_mat_entry(r, i, j), acb_imagref(acb_mat_entry(x, i, j)));
}

static void
acb_mat_set_real(acb_mat_t x, const arb_mat_t r)
{
    slong i, j;
    for (i = 0; i < acb_mat_nrows(x); i++)
        for (j = 0; j < acb_mat_ncols(x); j++)
            arb_set(acb_realref(acb_mat_entry(x, i, j)), arb_mat_entry(r, i, j));
}

static void
acb_mat_set_imag(acb_mat_t x, const arb_mat_t r)
{
    slong i, j;
    for (i = 0; i < acb_mat_nrows(x); i++)
        for (j = 0; j < acb_mat_ncols(x); j++)
            arb_set(acb_imagref(acb_mat_entry(x, i, j)), arb_mat_entry(r, i, j));
}

int
acb_mat_solve_real(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    int result;
    arb_mat_t a, b, t, c, d, x;
    slong n, m, *perm;
    arb_mat_t LU;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(a, n, n);
    arb_mat_init(b, n, n);
    arb_mat_init(t, n, n);
    arb_mat_init(c, n, m);
    arb_mat_init(d, n, m);
    arb_mat_init(x, n, m);

    acb_mat_real(a, A);
    acb_mat_imag(b, A);
    acb_mat_real(c, B);
    acb_mat_imag(d, B);

    perm = _perm_init(n);
    arb_mat_init(LU, n, n);

    if ((result = arb_mat_lu(perm, LU, a, prec)))
    {
        /* d <- a^(-1)d */
        arb_mat_solve_lu_precomp(d, perm, LU, d, prec);
        /* c <- c + ba^(-1)d */
        arb_mat_mul(x, b, d, prec);
        arb_mat_add(c, c, x, prec);
        /* t <- a^(-1)b */
        arb_mat_solve_lu_precomp(t, perm, LU, b, prec);
        /* b <- ba^(-1)b */
        arb_mat_mul(b, b, t, prec);
        /* a <- a + ba^(-1)b */
        arb_mat_add(a, a, b, prec);
        if ((result = arb_mat_lu(perm, LU, a, prec)))
        {
            /* (a + ba^(-1)b) x = c + ba^(-1)d */
            arb_mat_solve_lu_precomp(x, perm, LU, c, prec);
            acb_mat_set_real(X, x);
            /* y = a^(-1)d - a^(-1)b x */
            arb_mat_mul(c, t, x, prec);
            arb_mat_sub(x, d, c, prec);
            acb_mat_set_imag(X, x);
        }
    }
    arb_mat_clear(a);
    arb_mat_clear(b);
    arb_mat_clear(t);
    arb_mat_clear(c);
    arb_mat_clear(d);
    arb_mat_init(x, n, m);
    arb_mat_clear(LU);
    _perm_clear(perm);

    return result;
}

int
acb_mat_solve_charpoly(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    slong k;
    acb_poly_t chi;
    acb_mat_t AB;
    acb_poly_init(chi);
    acb_mat_init(AB, acb_mat_nrows(B), acb_mat_ncols(B));
    acb_mat_charpoly(chi, A, prec);
    acb_poly_scalar_div(chi, chi, acb_poly_get_coeff_ptr(chi, 0), prec);
    /* evaluate on B */
    acb_mat_scalar_mul_acb(X, B, acb_poly_get_coeff_ptr(chi, 1), prec);
    acb_mat_set(AB, B);
    for (k = 2; k <= acb_poly_degree(chi); k++)
    {
        acb_mat_mul(AB, A, AB, prec);
        acb_mat_scalar_addmul_acb(X, AB, acb_poly_get_coeff_ptr(chi, k), prec);
    }
    acb_mat_clear(AB);
    acb_poly_clear(chi);
    return 1;
}

/* try several strategies */
void
acb_mat_solve_multi(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, slong prec)
{
    if (!acb_mat_solve(X, A, B, prec))
    {
        acb_t det;
        acb_init(det);
        acb_mat_det(det, A, prec);
        if (acb_contains_zero(det))
        {
            flint_printf("ERROR: matrix is not invertible\n");
            flint_printf("det(A) = ");
            acb_printd(det, 30);
            flint_printf("\n");
            acb_mat_printd(A,20);
            flint_printf("\n");
            acb_mat_printd(B,20);
            flint_printf("\n");
            abort();
        }
        flint_printf("acb matrix inversion failed\n");
        if (!acb_mat_solve_real(X, A, B, prec))
        {
            acb_t i;
            acb_mat_t iA, iB;
            flint_printf("arb matrix inversion A failed\n");

            acb_init(i);
            acb_onei(i);
            acb_mat_init(iA, acb_mat_nrows(A), acb_mat_ncols(A));
            acb_mat_init(iB, acb_mat_nrows(B), acb_mat_ncols(B));
            acb_mat_scalar_mul_acb(iA, A, i, prec);
            acb_mat_scalar_mul_acb(iB, B, i, prec);
            if (!acb_mat_solve_real(X, iA, iB, prec))
            {
                flint_printf("arb matrix inversion B failed\n");
                acb_mat_solve_charpoly(X, A, B, prec);
            }
            acb_clear(i);
            acb_mat_clear(iA);
            acb_mat_clear(iB);
        }
    }
}

void tau_matrix(acb_mat_t tau, const acb_mat_t omega0, const acb_mat_t omega1, slong prec)
{
    acb_mat_solve_multi(tau, omega0, omega1, prec);
}
