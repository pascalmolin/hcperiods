#include "abel_jacobi.h"

int main() {

    slong d, i, m;
    slong prec = 200;
    flint_rand_t state;
    
    flint_printf("abel_jacobi struct...");
    fflush(stdout);
    flint_randinit(state);

    for (d = 3; d < 10; d++)
    {
        for (i = 0; i < 5; i++)
        {
            acb_ptr x;
            slong i;
 
            x = _acb_vec_init(d);
            for (i = 0; i < d; i++)
            {
                acb_randtest_precise(x + i, state, prec, 4);
                /*
                flint_printf("\nx[%ld] = ", i);
                acb_printd(x + i, 10);
                */
            }

            for (m = 2; m < 7; m++)
            {
                abel_jacobi_t aj;
                abel_jacobi_init_roots(aj, m, x, d);

                abel_jacobi_compute(aj, prec);

                abel_jacobi_clear(aj);
            }

            _acb_vec_clear(x, d);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
