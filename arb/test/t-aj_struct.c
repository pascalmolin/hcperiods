#include "abel_jacobi.h"

void
acb_poly_random(acb_poly_t pol, flint_rand_t state, slong n, slong prec)
{
    slong i;
    do
        acb_poly_randtest(pol, state, n + 1, prec, 4);
    while (acb_poly_degree(pol) < n);

    for (i = 0; i <= n; i++)
    {
        mag_zero(arb_radref(acb_realref(pol->coeffs + i)));
        mag_zero(arb_radref(acb_imagref(pol->coeffs + i)));
    }
}


int main() {

    slong n, i, m;
    int flag = 0;
    slong prec = 200;
    flint_rand_t state;
    
    flint_printf("abel_jacobi struct...");
    fflush(stdout);
    flint_randinit(state);

    for (n = 3; n < 10; n++)
    {
        acb_poly_t pol;
        acb_poly_init(pol);

        for (i = 0; i < 5; i++)
        {

            acb_poly_random(pol, state, n, prec);
 
            for (m = 2; m < 7; m++)
            {
                abel_jacobi_t aj;

                progress("[n=%ld, m=%ld]",n,m);
                abel_jacobi_init_poly(aj, m, pol);
                abel_jacobi_compute(aj, flag, prec);

                abel_jacobi_clear(aj);
            }

        }

        acb_poly_clear(pol);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
