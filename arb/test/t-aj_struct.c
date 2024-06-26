#include "abel_jacobi.h"

/* check positivity */
void
acb_mat_imag(arb_mat_t im, const acb_mat_t m)
{
    slong i, n;
    n = acb_mat_nrows(m) * acb_mat_ncols(m);
    for (i = 0; i < n; i++)
        arb_set(im->entries + i, acb_imagref(m->entries + i));
}

int
check_pos_imag(const acb_mat_t tau, flint_rand_t state, slong prec)
{
    slong k, g, r = 1;
    arb_mat_t im;
    arb_mat_t x, y, z;

    g = acb_mat_nrows(tau);
    arb_mat_init(im, g, g);
    arb_mat_init(x, 1, g);
    arb_mat_init(y, g, 1);
    arb_mat_init(z, 1, 1);

    acb_mat_imag(im, tau);
    for (k = 0; k < 30; k++)
    {
        arb_vec_set_random_nonzero(x->entries, g, state, prec, 4);

        arb_mat_transpose(y, x);
        arb_mat_mul(x, x, im, prec);
        arb_mat_mul(z, x, y, prec);

        if (!arb_is_positive(acb_mat_entry(z, 0, 0)))
        {
            r = 0;
            break;
        }
    }

    arb_mat_clear(im);
    arb_mat_clear(x);
    arb_mat_clear(y);
    arb_mat_clear(z);
    return r;
}

int main() {

    slong n, i, m;
    int flag = 0;
    slong prec = 200;
    flint_rand_t state;
    
    flint_printf("abel_jacobi struct...");
    fflush(stdout);
    flint_rand_init(state);

    for (n = 3; n < 10; n++)
    {
        fmpz_poly_t pol;
        fmpz_poly_init(pol);

        for (i = 0; i < 5; i++)
        {

            do {
                fmpz_poly_randtest(pol, state, n + 1, 10);
            } while (fmpz_poly_degree(pol) < n || !fmpz_poly_is_squarefree(pol));
 
            for (m = 2; m < 7; m++)
            {
                abel_jacobi_t aj;

                progress("[n=%ld, m=%ld]",n, m);
                abel_jacobi_init_poly(aj, m, pol);
                abel_jacobi_compute(aj, flag, prec);

                if (!check_pos_imag(aj->tau, state, prec))
                {
                    flint_printf("FAIL: im(tau) non positive definite\n");
                    flint_printf("curve y^%ld = ", m);
                    fmpz_poly_print_pretty(pol, "x");
                    flint_printf("\n");
                    acb_mat_printd(aj->tau, 5);
                    abort();
                }

                abel_jacobi_clear(aj);
            }

        }

        fmpz_poly_clear(pol);
    }

    flint_rand_clear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
