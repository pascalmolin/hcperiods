#include "abel_jacobi.h"

int main() {

    slong i, prec;
    acb_ptr u, z;
    acb_t y1, y2, y3, y4;
    arb_t x;

    flint_rand_t state;

    flint_printf("mth_root...");
    fflush(stdout);
    flint_randinit(state);

    for (prec = 60; prec < 500; prec *= 2)
    {
        for (i = 0; i < 10; i++)
        {
            slong m, d, d1, j;

            d = 3 + n_randint(state, 20);
            m = 2 + n_randint(state, 12);
            u = _acb_vec_init(d);
            z = _acb_vec_init(m);
            _acb_vec_unit_roots(z, m, prec);
            arb_init(x);
            acb_init(y1);
            acb_init(y2);
            acb_init(y3);
            acb_init(y4);

            acb_vec_set_random_u(u, d, state, prec, 4, .01);

            for (j = 0; j < 10; j++)
            {
                arb_randtest_precise(x, state, prec, 3);
                arb_tanh(x, x, prec);

                for (d1 = 0; d1 <= d; d1++)
                {
                    mth_root_pol_def(y1, u, d1, d, x, NULL, m, prec);
                    mth_root_pol_prod(y2, u, d1, d, x, NULL, m, prec);
                    mth_root_pol_turn(y3, u, d1, d, x, z, m, prec);
                    if (!acb_overlaps(y1, y2) || !acb_overlaps(y2,y3) || !acb_overlaps(y1,y3))
                    {
                        flint_printf("FAIL:\n\n");
                        flint_printf("d = %ld, m = %ld, prec = %ld\n", d, m, prec);
                        flint_printf("mth_root_def = ");
                        acb_printd(y1, 20);
                        flint_printf("\nmth_root_pol = ");
                        acb_printd(y2, 20); 
                        flint_printf("\nmth_root_turn = ");
                        acb_printd(y3, 20);
                        flint_printf("\n\n");
                        abort();
                    }
                    if (m == 2)
                    {
                        sqrt_pol_def(y2, u, d1, d, x, prec);
                        sqrt_pol_turn(y3, u, d1, d, x, prec);
                        if (!acb_overlaps(y1, y2) || !acb_overlaps(y2,y3) || !acb_overlaps(y1,y3))
                        {
                            flint_printf("FAIL:\n\n");
                            flint_printf("d = %ld, m = %ld, prec = %ld\n", d, m, prec);
                            flint_printf("mth_root_def = ");
                            acb_printd(y1, 20);
                            flint_printf("\nsqrt_root_def ");
                            acb_printd(y2, 20);
                            flint_printf("\nsqrt_root_turn ");
                            acb_printd(y3, 20);
                            flint_printf("\n\n");
                            abort();
                        }
                    }
                }
            }

            arb_clear(x);
            acb_clear(y1);
            acb_clear(y2);
            acb_clear(y3);
            _acb_vec_clear(u, d);
            _acb_vec_clear(z, m);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
