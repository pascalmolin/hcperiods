#include "abel_jacobi.h"

#define OUTMAGMA 0

/* families of polynomials */
#define nf 2

typedef void (*pol_f) (acb_poly_t poly, ulong n, slong prec);

/* x^n+1 */
void
pol_xn1(acb_poly_t poly, ulong n, slong prec)
{
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, n, 1);
}

/* sum x^k/k! */
void
pol_exp(acb_poly_t poly, ulong n, slong prec)
{
    acb_poly_one(poly);
    acb_poly_shift_left(poly, poly, 1);
    acb_poly_exp_series(poly, poly, n, prec);
}

pol_f vf[nf] = { pol_xn1, pol_exp };

/* degrees */
#define nn 5
slong vn[nn] = { 3, 6, 9, 11, 12 };

/* sheets */
#define nm 4
slong vm[nm] = { 2, 3, 4, 5 };

/* checked values */
#define nc 4

/* here need 160 values */
#define nref nf * nn * nm * nc

const char * ref_r[nref] = {
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"
};

const char * ref_i[nref] = {
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0",
    "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"
};

int main() {

    slong prec;

    flint_printf("integrals...");
    fflush(stdout);

    for (prec = 128; prec < 100; prec *= 2)
    {
        slong f, k;
#if !OUTMAGMA
        const char ** real = ref_r, ** imag = ref_i;
#endif
        acb_t ref;
        acb_init(ref);

        k = 0;
        for (f = 0; f < nf; f++)
        {
            slong in;

            for (in = 0; in < nn; in++)
            {
                slong im;
                acb_poly_t pol;

                acb_poly_init(pol);

                vf[f](pol, vn[in], prec);

                for (im = 0; im < nm; im++)
                {
#if !OUTMAGMA
                    slong c;
#endif
                    abel_jacobi_t aj;
                    abel_jacobi_init_poly(aj, vm[im], pol, prec);
                    abel_jacobi_compute(aj, prec);
#if OUTMAGMA
                    /* output parameters for magma */
                    flint_printf("\nn := %ld; m := %ld; prec := %ld;\n", vn[in], vm[im], prec);
                    flint_printf("roots := [");
                    _acb_vec_arf_printd(aj->c.roots, vn[in], 30, ", ");
                    flint_printf("];\n\n");
#else
                    /* check first entry */
                    for (c = 0; c < nc; c++)
                    {
                        acb_ptr coeff;

                        if (arb_set_str(acb_realref(ref), *(real++), prec) ||
                                arb_set_str(acb_imagref(ref), *(imag++), prec) )
                        {
                            flint_printf("error while setting ref <- %s+I*%s\n", real, imag);
                            abort();
                        }

                        coeff = acb_mat_entry(aj->omega0, c % 2, c / 2);
                        if (!acb_overlaps(ref, coeff))
                        {
                            flint_printf("FAIL:\n\n");
                            flint_printf("pol = ");
                            acb_poly_printd(pol, 10);
                            flint_printf("\nn = %ld, m = %ld, prec = %ld\n", vn[in], vm[im], prec);
                            flint_printf("\nref = ");
                            acb_printd(ref, 20);
                            flint_printf("\nperiod = ");
                            acb_printd(coeff, 20);
                            flint_printf("\n\n");
                            abort();
                        }
                    }
#endif
                    abel_jacobi_clear(aj);
                }

                acb_poly_clear(pol);
            }
        }
        acb_clear(ref);
    }

    flint_cleanup();
    printf("PASS\n");
    return 0;
}
