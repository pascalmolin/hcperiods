#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "abel_jacobi.h"

/* precisions */
#define np 1
slong vp[np] = { 256 };

/* families of polynomials */
#define nf 2

typedef void (*pol_f) (fmpz_poly_t poly, ulong n);

/* x^n+1 */
void
pol_xn1(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_fit_length(poly, n + 1);
    fmpz_poly_set_coeff_si(poly, 0, 1);
    fmpz_poly_set_coeff_si(poly, n, 1);
}

pol_f vf[nf] = { pol_xn1, fmpz_poly_fibonacci };

/* degrees */
#define nn 5
slong vn[nn] = { 5, 6, 9, 11, 12 };

/* sheets */
#define nm 4
slong vm[nm] = { 2, 3, 4, 5 };

/* checked values */
#define nc 4

int
acb_get_file(acb_t z, FILE * fp, slong prec)
{
    int i;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    arb_ptr x = acb_realref(z);

    for (i = 0; i < 2; i++)
    {
        if ((read = getline(&line, &len, fp)) < 0 || arb_set_str(x, line, prec))
        {
            flint_printf("error while setting value ref <- %s\n", line);
            abort();
        }
        x = acb_imagref(z);
    }
    free(line);
    return 0;
}


void
generate()
{
    int ip, f, in, im;
    int flag = 0;

    flint_printf("[");
    for (ip = 0; ip < np; ip++)
    {
        slong prec = vp[ip] + 50;
        slong digits = ((double)prec * .30103) + 5;

        flint_printf("< %ld, [", digits);
        for (f = 0; f < nf; f++)
        {
            for (in = 0; in < nn; in++)
            {
                fmpz_poly_t pol;
                abel_jacobi_t aj;

                fmpz_poly_init(pol);
                vf[f](pol, vn[in]);

                abel_jacobi_init_poly(aj, 2, pol);
                /* should not be necessary */
                abel_jacobi_compute(aj, flag, prec);

                /* print roots */
                flint_printf("< [");
                _acb_vec_arf_printd(aj->roots, vn[in], digits, ", ");
                flint_printf("]");

                /* print exponents */
                flint_printf(",[%ld",vm[0]);
                for (im = 1; im < nm; im++)
                    flint_printf(",%ld",vm[im]);
                flint_printf("]>,\n");

                fmpz_poly_clear(pol);
                abel_jacobi_clear(aj);
            }
        }
        flint_printf("] >\n");
    }
    flint_printf("]\n");
}

void
check_periods()
{
    FILE * fp;
    int ip, f, in, im;
    int flag = 0;

    flint_printf("periods...");
    fflush(stdout);

    fp = fopen("ref_periods.txt", "r");
    if (fp == NULL)
        abort();

    for (ip = 0; ip < nm; ip++)
    {
        slong prec = vp[ip];
        acb_t ref;
        acb_init(ref);

        for (f = 0; f < nf; f++)
        {
            for (in = 0; in < nn; in++)
            {
                fmpz_poly_t pol;

                fmpz_poly_init(pol);

                vf[f](pol, vn[in]);

                for (im = 0; im < nm; im++)
                {
                    slong c;
                    abel_jacobi_t aj;

                    abel_jacobi_init_poly(aj, vm[im], pol);
                    abel_jacobi_compute(aj, flag, prec);

                    /* check 4 entries of each matrix */
                    for (c = 0; c < nc; c++)
                    {
                        acb_ptr coeff;

                        acb_get_file(ref, fp, prec);

                        coeff = acb_mat_entry(aj->omega0, c % 2, c / 2);
                        if (!acb_overlaps(ref, coeff))
                        {
                            flint_printf("FAIL:\n\n");
                            flint_printf("pol = ");
                            fmpz_poly_print(pol);
                            flint_printf("\nn = %ld, m = %ld, prec = %ld\n", vn[in], vm[im], prec);
                            flint_printf("\nref = ");
                            acb_printd(ref, 20);
                            flint_printf("\nperiod = ");
                            acb_printd(coeff, 20);
                            flint_printf("\n\n");
                            abort();
                        }
                    }
                    abel_jacobi_clear(aj);
                }

                fmpz_poly_clear(pol);
            }
        }
        acb_clear(ref);
    }

    fclose(fp);
    flint_cleanup();
    printf("PASS\n");
}

int main(int argc, char * argv[])
{
    if (argc > 1 && strcmp(argv[1], "gen") == 0)
        generate();
    else
        return 0; /*check_periods();*/
   return 0;
}
