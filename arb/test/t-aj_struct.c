#include "abel_jacobi.h"

int main() {

    slong n, i, m;
    int flag = 0;
    slong prec = 200;
    flint_rand_t state;
    
    flint_printf("abel_jacobi struct...");
    fflush(stdout);
    flint_randinit(state);

    for (n = 3; n < 5; n++)
    {
        acb_ptr x;
        x = _acb_vec_init(n);

        for (i = 0; i < 3; i++)
        {
            slong i;
 
            for (i = 0; i < n; i++)
                acb_randtest_precise(x + i, state, prec, 4);

            for (m = 2; m < 7; m++)
            {
                abel_jacobi_t aj;

                progress("[n=%ld, m=%ld]",n,m);
                abel_jacobi_init_roots(aj, m, x, n, flag);
                abel_jacobi_compute(aj, flag, prec);

                abel_jacobi_clear(aj);
            }

        }

        _acb_vec_clear(x, n);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
