#include "abel_jacobi.h"

int main() {

    slong g;
    flint_rand_t state;
    
    flint_printf("symplectic_reduction...");
    fflush(stdout);
    flint_randinit(state);

    for (g = 2; g < 20; g++)
    {       
        si_mat_t m1, m2, p;
        slong i, j, g2, l;

        g2 = 2 * g;
        m1 = si_mat_init(g2, g2);
        m2 = si_mat_init(g2, g2);
        p = si_mat_init(g2, g2);
        flint_printf("\n========= init completed =====\n");

        for (l = 0; l < 10; l++)
        { 

            /* random matrix */
            for (i = 0; i < g2; i++)
            {
                m1[i][i] = m2[i][i] = 0;
                for (j = i + 1; j < g2; j++)
                {
                    int t = n_randint(state, 3) - 1;
                    m1[i][j] = m2[i][j] = t;
                    m1[j][i] = m2[j][i] = -t;
                }
            }

            symplectic_reduction(p, m2, g, g2);

            if (!is_symplectic_j(m2, g, g2))
            {
                flint_printf("\n== initial m ==\n");
                si_mat_print_gp(m1, g2, g2);

                flint_printf("\n== reduced m ==\n");
                si_mat_print(m2, g2, g2);

                flint_printf("\n== p ==\n");
                si_mat_print(p, g2, g2);

                abort();
            }
        }

        si_mat_clear(m1, g2, g2);
        si_mat_clear(m2, g2, g2);
        si_mat_clear(p, g2, g2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
