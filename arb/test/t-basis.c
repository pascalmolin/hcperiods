#include "abel_jacobi.h"

int main() {

    slong n, i, m, prec = 60;
    flint_rand_t state;
    
    flint_printf("symplectic basis...");
    fflush(stdout);
    flint_randinit(state);

    for (n = 3; n < 10; n++)
    {
        for (i = 0; i < 5; i++)
        {
            acb_ptr x;
            tree_t tree;
 
            x = _acb_vec_init(n);
            acb_vec_set_random(x, n, state, prec, 4);

            tree_init(tree, n - 1);
            spanning_tree(tree, x, n, INT_DE);

            for (m = 2; m < 7; m++)
            {
                sec_t c;
                homol_t alpha, beta;

                sec_init(&c, m, n);
                tree_ydata_init(tree, x, n, m, prec);

                alpha = flint_malloc(c.g * sizeof(loop_t));
                beta = flint_malloc(c.g * sizeof(loop_t));

                symplectic_basis(alpha, beta, tree, c);

                homol_clear(alpha, c.g);
                homol_clear(beta, c.g);

                tree_ydata_clear(tree);
                sec_clear(c);
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
