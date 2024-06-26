#include "abel_jacobi.h"

#define m1(i,j) fmpz_mat_entry(m1,i,j)
#define m2(i,j) fmpz_mat_entry(m2,i,j)

int main() {

    slong g;
    flint_rand_t state;
    
    flint_printf("symplectic_reduction...");
    fflush(stdout);
    flint_rand_init(state);

    for (g = 2; g < 30; g++)
    {       
        fmpz_mat_t m1, m2, p;
        slong i, j, g2, l;

        g2 = 2 * g;
        fmpz_mat_init(m1, g2, g2);
        fmpz_mat_init(m2, g2, g2);
        fmpz_mat_init(p, g2, g2);

        for (l = 0; l < 10; l++)
        { 

            /* random matrix with 0,1,-1,
             * density decreases with l */
            for (i = 0; i < g2; i++)
            {
                fmpz_zero(m1(i, i));
                fmpz_zero(m2(i, i));
                for (j = i + 1; j < g2; j++)
                {
#if 0
                    slong t;
                    if (n_randint(state, 2) > 0)
                        t = 0;
                    else
                        t = n_randint(state, 3) - 1;
#else
                        slong t = n_randint(state, 3) - 1;
#endif
                    fmpz_set_si(m1(i, j), t);
                    fmpz_set_si(m2(i, j), t);
                    fmpz_set_si(m1(j, i), -t);
                    fmpz_set_si(m2(j, i), -t);
                }
            }

#if 0
                flint_printf("\n== initial m ==\n");
                fmpz_mat_print(m1);
#endif

            symplectic_reduction(p, m2, g);

            if (!is_symplectic_j(m2, g, g2))
            {
                flint_printf("\n== initial m ==\n");
                fmpz_mat_print(m1);

                flint_printf("\n== reduced m ==\n");
                fmpz_mat_print(m2);

                flint_printf("\n== p ==\n");
                fmpz_mat_print(p);

                abort();
            }
        }

        fmpz_mat_clear(m1);
        fmpz_mat_clear(m2);
        fmpz_mat_clear(p);
    }

    flint_rand_clear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
