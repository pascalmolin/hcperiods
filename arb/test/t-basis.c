#include "abel_jacobi.h"

int main() {

    slong d, i, m;
    flint_rand_t state;
    
    flint_printf("symplectic basis...");
    fflush(stdout);
    flint_randinit(state);

    for (d = 3; d < 10; d++)
    {
        for (i = 0; i < 5; i++)
        {
            acb_ptr x;
            tree_t tree;
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
                sec_t c;
                homol_t alpha, beta;

                sec_init(&c, m, x, d);

                alpha = malloc(c.g * sizeof(loop_t));
                beta = malloc(c.g * sizeof(loop_t));

                symplectic_basis(alpha, beta, tree, c);


                homol_clear(alpha, c.g);
                homol_clear(beta, c.g);

                sec_clear(c);
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
