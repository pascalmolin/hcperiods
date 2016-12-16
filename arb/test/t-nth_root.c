#include "abel_jacobi.h"

int main() {

    slong i, prec;
    acb_ptr u;
    acb_t y1, y2;
    arb_t x;

    flint_rand_t state;
    
    flint_printf("nth_root...");
    fflush(stdout);
    flint_randinit(state);

    for (prec = 60; prec < 1000; prec *= 2)
    {
        for (i = 0; i < 10; i++)
        {
            slong m, d, j;

            d = 3 + n_randint(state, 20);
            m = n_randint(state, 12);
            u = _acb_vec_init(d);
            arb_init(x);
            acb_init(y1);
            acb_init(y2);

            for (j = 0; j < d; j++)
                acb_randtest_precise(u + j, state, prec, 10);

            for (j = 0; j < 10; j++)
            {
                arb_randtest_precise(x, state, prec, 1);
                nth_root_pol_def(y1, u, x, d, m, prec);
                nth_root_pol_prod(y2, u, x, d, m, prec);
                if (!acb_overlaps(y1, y2))
                {
                    flint_printf("FAIL:\n\n");
                    flint_printf("d = %ld, m = %ld, prec = %ld\n", d, m, prec);
                    flint_printf("nth_root_def = ");
                    acb_printd(y1, 20);
                    flint_printf("\nnth_root_pol = ");
                    acb_printd(y2, 20);
                    flint_printf("\n\n");
                    abort();
                }
            }

            arb_clear(x);
            acb_clear(y1);
            acb_clear(y2);
            _acb_vec_clear(u, d);

        }
    }
    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
