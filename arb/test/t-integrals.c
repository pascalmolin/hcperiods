#include "abel_jacobi.h"

/* check with gp */
/*
 x = [5, 8 , -9, 8, 1 , 7, -4 , 3 , -5, 0, 0 , -9, 4 , 4, 4, 7 , 10, 4, -3, -9];
 y = [9, -5, 1 , 0, -4, 8, -10, -8, 0 , 6, -8, -9, -6, 0, 5, -7, 4 , 7, 4 , 5];
 u = vector(#x,k,(x[k]+I*y[k])/sqrt(2));
 { for(d=1,3, print( vector((d+1)\2,i,
     intnum(t=[-1,-1/2],[1,-1/2],t^(i-1)/(sqrt(1-t^2)*prod(k=1,d,sqrt(t-u[k]))))
     ))) }
 */
#define dmax 20
int x[dmax] = {5, 8 , -9, 8, 1 , 7, -4 , 3 , -5, 0, 0 , -9, 4 , 4, 4, 7 , 10, 4, -3, -9};
int y[dmax] = {9, -5, 1 , 0, -4, 8, -10, -8, 0 , 6, -8, -9, -6, 0, 5, -7, 4 , 7, 4 , 5};

int main() {

    slong prec;

    flint_rand_t state;

    flint_printf("integrals...");
    fflush(stdout);
    flint_randinit(state);

    for (prec = 60; prec < 150; prec *= 2)
    {
        /* start with m = 2 only */
        slong j, d; /*, m = 2;*/
        acb_ptr u;
        arb_t s2;
        u = _acb_vec_init(dmax);

        arb_init(s2);
        arb_set_si(s2, 2);
        arb_sqrt(s2, s2, prec);

        for (j = 0; j < dmax; j++)
        {
            arb_set_si(acb_realref(u + j), x[j]);
            arb_set_si(acb_imagref(u + j), y[j]);
            acb_div_arb(u + j, u + j, s2, prec);
        }
        arb_clear(s2);


        for (d = 1; d < 4; d ++)
        {
            slong g = (d + 1) / 2;
            slong n, j;
            acb_ptr res;

            res = _acb_vec_init(g);

            n = gc_params(u, d, -1, g - 1, prec);

            /* take d1 = d ie sqrt(x-u[i]) for all */
            gc_integrals(res, u, d, d, g, n, prec);

#if 0
            for (j = 0; j < g; j++)
            {
                flint_printf("\nI[%ld,%ld] = ", d, j);
                acb_printd(res + j, 20);
            }
#endif

#if 0
            for (j = 0; j < 0; j++)
            {
                if (!acb_overlaps(res y1, y2))
                {
                    flint_printf("FAIL:\n\n");
                    flint_printf("d = %ld, m = %ld, prec = %ld\n", d, m, prec);
                    flint_printf("mth_root_def = ");
                    acb_printd(y1, 20);
                    flint_printf("\nmth_root_pol = ");
                    acb_printd(y2, 20);
                    flint_printf("\n\n");
                    abort();
                }
            }
#endif

            _acb_vec_clear(res, g);
        }
        _acb_vec_clear(u, dmax);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
