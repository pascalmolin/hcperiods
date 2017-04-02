#include "abel_jacobi.h"

#define dmax 30
#define imax 40
#define nf 2

typedef slong (*param_f) (acb_srcptr x, slong d, slong prec);

int main() {

    slong d, i;
    slong prec = 60;
    flint_rand_t state;

    flint_printf("integration parameters...");
    fflush(stdout);
    flint_randinit(state);

    for (d = 3; d < 30; d++)
    {
        slong f;
        acb_mat_t pols;
        slong nmin[nf], nmax[nf];
        double nmed[nf];
        mag_t e_de, e_gc;
        arf_t h, l;
        arf_init(h);
        arf_init(l);
        mag_init(e_gc);
        mag_init(e_de);

        acb_mat_init(pols, imax, d);

        /* create imax random examples */
        for (i = 0; i < imax; i++)
            acb_vec_set_random_u(pols->rows[i], d, state, prec, 4, .01);

        /* compute integration parameters */
        for (f = 0; f < nf; f++)
        {
            nmed[f] = 0.; nmin[f] = LONG_MAX; nmax[f] = 0;
            for (i = 0; i < imax; i++)
            {
                slong n;
#if 0
                flint_printf("\nd = %ld, i = %ld, %6s\n", d, i, f ? "de" : "gauss");
                for (n = 0; n < d; n++)
                    flint_printf("\nu_%ld = ", n),
                    acb_printd(acb_mat_entry(pols, i, n), 10);
#endif
                if (f == 0)
                    n = gc_params(e_gc, pols->rows[i], d, 0, prec);
                else
                    n = de_params(e_de, h, l, pols->rows[i], d, 0., 0, 1, 2, prec);
                nmed[f] += n;
                if (n < nmin[f]) nmin[f] = n;
                if (n > nmax[f]) nmax[f] = n;
            }
            nmed[f] /= imax;
            flint_printf("\n%6s: d = %3ld, min, max, med = %3ld, %8ld, %8.3lf",
                    f ? "de" : "gauss", d, nmin[f], nmax[f], nmed[f]);
        }

        acb_mat_clear(pols);
        arf_clear(h);
        arf_clear(l);
        mag_clear(e_gc);
        mag_clear(e_de);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}
