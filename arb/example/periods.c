#include <string.h>

#include "flint/arith.h"
#include "acb_poly.h"
#include "abel_jacobi.h"

/* x^n-1 */
void
pol_xn1(acb_poly_t poly, slong n, slong prec)
{
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, n, -1);
}
/* sum x^k/k! */
void
pol_exp(acb_poly_t poly, slong n, slong prec)
{
    acb_poly_one(poly);
    acb_poly_shift_left(poly, poly, 1);
    acb_poly_exp_series(poly, poly, n, prec);
}
void
pol_bern(acb_poly_t pol, slong n, slong prec)
{
    fmpq_poly_t h;
    fmpq_poly_init(h);
    arith_bernoulli_polynomial(h, n);
    acb_poly_set_fmpq_poly(pol, h, prec);
    fmpq_poly_clear(h);
}

int main(int argc, char * argv[])
{
    int i, print = 1, flag = 0, run = 1;
    slong n = 5, m = 2, prec = 128, digits = 0;
    void (*f_print) (const acb_mat_t, slong) = &acb_mat_printd;
    void (*f_pol) (acb_poly_t pol, slong n, slong prec) = &pol_xn1;
    acb_poly_t poly;
    abel_jacobi_t aj;

    acb_poly_init(poly);

    if (argc < 2)
    {
        flint_printf("periods [-m m] [--prec p] poly \n");
        flint_printf("Print period matrix of curve y^m = f_n(x).\n");
        flint_printf("Default m = 2, f_n(x) = x^5 + 1, prec = 128.\n");
        flint_printf("Polynomials f_n(x):\n");
        flint_printf("  --xn1 n: x^n + 1\n");
        flint_printf("  --bern n: Bernoulli polynomial Bn(x)\n");
        flint_printf("  --cheb n: Chebychev polynomial Tn(x)\n");
        flint_printf("  --pol n cn ... c1 c0 : cn x^n + ... + c1 x + c0\n");
        flint_printf("Output options:\n");
        flint_printf("  --quiet: no printing\n");
        flint_printf("  --trim: reduce to obtained precision\n");
        flint_printf("  --big: big period matrix\n");
        flint_printf("  --gp: output for pari/gp (discard error balls)\n");
        flint_printf("  --de: force use of DE integration\n");
        return 1;
    }

    for (i = 1; i < argc;)
    {
        /* parameters */
        if (!strcmp(argv[i], "-m"))
        {
            i++;
            m = atol(argv[i++]);
        }
        else if (!strcmp(argv[i], "--prec"))
        {
            i++;
            prec = atol(argv[i++]);
        }
        else if (!strcmp(argv[i], "--de"))
        {
            i++;
            flag |= AJ_USE_DE;
        }
        /* mth root */
        else if (!strcmp(argv[i], "--rootdef"))
        {
            i++;
            flag |= AJ_ROOT_DEF;
        }
        else if (!strcmp(argv[i], "--rootprod"))
        {
            i++;
            flag |= AJ_ROOT_PROD;
        }
        else if (!strcmp(argv[i], "--rootturn"))
        {
            i++;
            flag |= AJ_ROOT_TURN;
        }
        /* families */
        else if (!strcmp(argv[i], "--xn1"))
        {
            i++;
            n = atol(argv[i++]);
            f_pol = &pol_xn1;
        }
        else if (!strcmp(argv[i], "--exp"))
        {
            i++;
            n = atol(argv[i++]);
            f_pol = &pol_exp;
        }
        else if (!strcmp(argv[i], "--bern"))
        {
            i++;
            n = atol(argv[i++]);
            f_pol = &pol_bern;
        }
        else if (!strcmp(argv[i], "--pol"))
        {
            int j;
            i++;
            n = atol(argv[i++]);
            for (j = 0; j <= n && i < argc; j++)
                acb_poly_set_coeff_si(poly, n - j, atol(argv[i++]));
            f_pol = NULL;
        }
        /* restrict computations / output */
        else if (!strcmp(argv[i], "--big"))
        {
            i++;
            flag |= AJ_NO_TAU;
        }
        else if (!strcmp(argv[i], "--int"))
        {
            i++;
            flag |= AJ_NO_AB;
        }
        else if (!strcmp(argv[i], "--tree"))
        {
            i++;
            flag |= AJ_NO_INT;
        }
        else if (!strcmp(argv[i], "--gp"))
        {
            i++;
            f_print = &acb_mat_print_gp;
        }
        else if (!strcmp(argv[i], "--trim"))
        {
            i++;
            flag |= AJ_TRIM;
        }
        else if (!strcmp(argv[i], "--digits"))
        {
            i++;
            digits = atol(argv[i++]);
        }
        /* benchmark mode */
        else if (!strcmp(argv[i], "--bench"))
        {
            i++;
            run = atol(argv[i++]);
            print = 0;
        }
         else
            i++;
    }
    if (!digits)
        digits = (slong)(prec * .301);

    /* compute pol to actual accuracy */
    if (f_pol)
        f_pol(poly, n, prec);

    abel_jacobi_init_poly(aj, m, poly);

    for (i = 0; i < run; i++)
        abel_jacobi_compute(aj, flag, prec);

    if (print)
    {
        if (flag & AJ_NO_INT)
        {
            tree_print(aj->tree);
        }
        else if (flag & AJ_NO_AB)
        {
            f_print(aj->integrals, digits);
        }
        else if (flag & AJ_NO_TAU)
        {
            f_print(aj->omega0, digits);
            flint_printf("\n");
            f_print(aj->omega1, digits);
        }
        else
            f_print(aj->tau, digits);
        flint_printf("\n");
    }

    abel_jacobi_clear(aj);
    acb_poly_clear(poly);
    flint_cleanup();
}
