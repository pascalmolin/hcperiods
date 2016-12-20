#include "abel_jacobi.h"

int main() {

    slong d, i, m;
    flint_rand_t state;
    
    flint_printf("intersection matrix...");
    fflush(stdout);
    flint_randinit(state);

    for (d = 3; d < 10; d++)
    {
        for (i = 0; i < 5; i++)
        {
            tree_t tree;
            acb_ptr x;
            slong i;
 
            x = _acb_vec_init(d);
            for (i = 0; i < d; i++)
            {
                acb_randtest_precise(x + i, state, 100, 4);
                /*
                flint_printf("\nx[%ld] = ", i);
                acb_printd(x + i, 10);
                */
            }

            tree_init(tree, d - 1);

            spanning_tree(tree, x, d, INT_DE);

            for (m = 2; m < 7; m++)
            {
                fmpz_mat_t c;
                slong len = (d-1)*(m-1);
                fmpz_mat_init(c, len, len);

                intersection_tree(c, tree, d, m);

#if 0
                flint_printf("\n\n ======= d = %ld, m = %ld ===== \n", d, m);
                fmpz_mat_print_pretty(c);
#endif

                fmpz_mat_clear(c);
            }

            tree_clear(tree);
            _acb_vec_clear(x, d);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
