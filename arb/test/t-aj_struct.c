#include "abel_jacobi.h"

int main() {

    slong n, i, m;
    slong prec = 200;
    flint_rand_t state;
    
    flint_printf("abel_jacobi struct...");
    fflush(stdout);
    flint_randinit(state);

    for (n = 3; n < 10; n++)
    {
        for (i = 0; i < 5; i++)
        {
            acb_ptr x;
            slong i;
 
            x = _acb_vec_init(n);
            for (i = 0; i < n; i++)
                acb_randtest_precise(x + i, state, prec, 4);

            for (m = 2; m < 7; m++)
            {
                abel_jacobi_t aj;

                progress("[n=%ld, m=%ld]",n,m);
                abel_jacobi_init_roots(aj, m, x, n);
                abel_jacobi_compute(aj, prec);

                abel_jacobi_clear(aj);
            }

            _acb_vec_clear(x, n);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
