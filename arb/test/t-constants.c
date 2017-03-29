#include "abel_jacobi.h"

int main() {

    slong n, i, m;
    slong prec = 200;
    flint_rand_t state;
    
    flint_printf("constants...");
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
                acb_randtest_precise(x + i, state, prec, 5);

            for (m = 2; m < 7; m++)
            {
                slong k;
                sec_t c;
                tree_t tree;

                sec_init(&c, m, n);
                tree_init(tree, n - 1);
                spanning_tree(tree, x, n, INT_DE);
                tree_ydata_init(tree, x, n, m, prec);

                /* check some roots */
                for (k = 0; k < n - 1; k++)
                {
                    /* f((a+b)/2) */
                    slong i;
                    slong nab = tree->data[k].n1;
                    slong a = tree->e[k].a, b = tree->e[k].b;
                    acb_ptr uab = tree->data[k].u;
                    acb_t z, fx, y;
                    acb_init(z);
                    acb_init(y);
                    acb_init(fx);
                    acb_one(fx);
                    for (i = 0; i < n; i++)
                    {
                        acb_add(z, x + a, x + b, prec);
                        acb_mul_2exp_si(z, z, -1);
                        acb_sub(z, z, x + i, prec);
                        acb_mul(fx, fx, z, prec);
                    }
                    acb_zero(z);
                    mth_root_pol_def(y, uab, nab, n - 2, acb_realref(z), m, prec);
                    acb_mul(y, y, tree->data[k].c, prec);
                    acb_pow_ui(y, y, m, prec);

                    if (!acb_overlaps(y, fx))
                    {
                        flint_printf("FAIL\n");
                        flint_printf("m = %ld, n = %ld, k = %ld (a,b = %ld,%ld)\n",m, n, k, a, b);

                        for (i = 0; i < n; i++)
                        {
                            flint_printf("\nx[%ld] = ", i);
                            acb_printd(x + i, 20);
                        }
                        for (i = 0; i < n - 2; i++)
                        {
                            flint_printf("\nu[%ld] = ", i);
                            acb_printd(uab + i, 20);
                        }
                        flint_printf("\n(b-a)/2 = ");
                        acb_printd(uab + n - 2, 20);
                        flint_printf("\n(a+b)/(b-a) = ");
                        acb_printd(uab + n - 1, 20);
                        flint_printf("\nCab = ");
                        acb_printd(uab + n, 20);

                        flint_printf("\nf(x) = ");
                        acb_printd(fx, 20);
                        flint_printf("\ny^m = ");
                        acb_printd(y, 20);
                        abort();
                    }

                    acb_clear(fx);
                    acb_clear(y);
                    acb_clear(z);
                }

                tree_ydata_clear(tree);
                tree_clear(tree);
                sec_clear(c);

            }

            _acb_vec_clear(x, n);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
