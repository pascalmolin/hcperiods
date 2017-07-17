#include "abel_jacobi.h"

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

                abel_jacobi_clear(aj);
            }

        }

        fmpz_poly_clear(pol);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
