#include "abel_jacobi.h"

int main() {

    slong n, i, m, prec = 100;
    flint_rand_t state;
    
    flint_printf("intersection matrix...");
    fflush(stdout);
    flint_randinit(state);

    for (n = 3; n < 10; n++)
    {
        for (i = 0; i < 5; i++)
        {
            tree_t tree;
            acb_ptr x;
            slong i;
 
            x = _acb_vec_init(n);
            for (i = 0; i < n; i++)
            {
                acb_randtest_precise(x + i, state, prec, 4);
                /*
                flint_printf("\nx[%ld] = ", i);
                acb_printd(x + i, 10);
                */
            }

            tree_init(tree, n - 1);

            spanning_tree(tree, x, n, INT_DE);

            for (m = 2; m < 7; m++)
            {
                sec_t c;
                fmpz_mat_t ca;
                slong len = (n-1)*(m-1);

                sec_init(&c, m, x, n);
                tree_ydata_init(tree, x, n, m, prec);

                fmpz_mat_init(ca, len, len);

                intersection_tree(ca, tree, n, m);

#if 0
                flint_printf("\n\n ======= d = %ld, m = %ld ===== \n", d, m);
                fmpz_mat_print_pretty(c);
#endif

                tree_ydata_clear(tree);
                fmpz_mat_clear(ca);
            }

            tree_clear(tree);
            _acb_vec_clear(x, n);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
