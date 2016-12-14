#include "abel_jacobi.h"

#define m1(i,j) fmpz_mat_entry(m1,i,j)
#define m2(i,j) fmpz_mat_entry(m2,i,j)

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
        si_mat_init(m1, g2, g2);
        si_mat_init(m2, g2, g2);
        si_mat_init(p, g2, g2);
        flint_printf("\n========= init completed =====\n");

        for (l = 0; l < 10; l++)
        { 

            /* random matrix */
            for (i = 0; i < g2; i++)
            {
                fmpz_zero(m1(i, i));
                fmpz_zero(m2(i, i));
                for (j = i + 1; j < g2; j++)
                {
                    slong t = n_randint(state, 3) - 1;
                    fmpz_set_si(m1(i, j), t);
                    fmpz_set_si(m2(i, j), t);
                    fmpz_set_si(m1(j, i), -t);
                    fmpz_set_si(m2(j, i), -t);
                }
            }

        flint_printf("\n== initial m ==\n");
        si_mat_print(m2, g2, g2);

            symplectic_reduction(p, m2, g, g2);

            if (!is_symplectic_j(m2, g, g2))
            {
                flint_printf("\n== initial m ==\n");
                si_mat_print(m1, g2, g2);

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
