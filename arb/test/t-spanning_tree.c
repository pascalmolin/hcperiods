#include "abel_jacobi.h"

int main() {

    slong d, i;
    flint_rand_t state;
    
    flint_printf("spanning_tree...");
    fflush(stdout);
    flint_randinit(state);

    for (d = 3; d < 30; d++)
    {
        for (i = 0; i < 10; i++)
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
            spanning_tree(tree, x, d, INT_GC);

            tree_clear(tree);
            _acb_vec_clear(x, d);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
