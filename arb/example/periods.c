#include <string.h>

#include "flint/arith.h"
#include "fmpz_poly.h"
#include "abel_jacobi.h"
#include "parse.h"

void
fmpz_poly_numer(fmpz_poly_t num, fmpq_poly_t pol)
{
    int k, n = fmpq_poly_degree(pol);
    fmpz * coeffs = fmpq_poly_numref(pol);
    for (k = 0; k <= n; k++)
        fmpz_poly_set_coeff_fmpz(num, k, coeffs + k);
}

/* x^n-1 */
void
pol_xn1(fmpz_poly_t poly, slong n, slong prec)
{
    fmpz_poly_set_coeff_si(poly, 0, -1);
    fmpz_poly_set_coeff_si(poly, n, 1);
}
void
pol_bern(fmpz_poly_t pol, slong n, slong prec)
{
    fmpq_poly_t h;
    fmpq_poly_init(h);
    arith_bernoulli_polynomial(h, n);
    fmpz_poly_numer(pol, h);
    fmpq_poly_clear(h);
}

int
usage()
{
    flint_printf("periods [-v] [-m m] [--prec p] poly \n");
    flint_printf("Print period matrix of curve y^m = f_n(x).\n");
    flint_printf("Default m = 2, f_n(x) = x^5 + 1, prec = 128.\n");
    flint_printf("Polynomials f_n(x):\n");
    flint_printf("  --xn1 n: x^n - 1\n");
    flint_printf("  --bern n: Bernoulli polynomial Bn(x)\n");
    flint_printf("  --exp n: exponential polynomial sum x^k/k!\n");
    flint_printf("  --coeffs n cn ... c1 c0 : cn x^n + ... + c1 x + c0\n");
    flint_printf("  --bernrev n: reverse Bernoulli x^nBn(1/x)\n");
    flint_printf("  --exprev n: reverse exponential sum x^k/(n-k)!\n");
    flint_printf("  --stoll: 82342800*x^6 - 470135160*x^5 + 52485681*x^4 + 2396040466*x^3 + 567207969*x^2 - 985905640*x + 247747600\n");
    flint_printf("  --pol '<string>': polynomial\n");
    flint_printf("Output options:\n");
    flint_printf("  --quiet: no printing\n");
    flint_printf("  --trim: reduce to obtained precision\n");
    flint_printf("  --big: big period matrix\n");
    flint_printf("  --gp: output for pari/gp (discard error balls)\n");
    flint_printf("  --error: output worst error\n");
    flint_printf("  --de: force use of DE integration\n");
    flint_printf("  --desame: keep integration points\n");
    return 1;
}

int main(int argc, char * argv[])
{
    int i, print = 1, flag = 0, run = 1, rev = 0;
    slong n = 5, m = 2, prec = 128, digits = 0;
    void (*f_print) (const acb_mat_t, slong) = &acb_mat_printd;
    void (*f_pol) (fmpz_poly_t pol, slong n, slong prec) = &pol_xn1;
    fmpz_poly_t poly;
    abel_jacobi_t aj;

    fmpz_poly_init(poly);

    if (argc < 2)
        return usage();

    for (i = 1; i < argc;)
    {
        /* parameters */
        if (!strcmp(argv[i], "-v"))
            i++, flag += AJ_VERBOSE;
        else if (!strcmp(argv[i], "-m"))
            i++, m = atol(argv[i++]);
        else if (!strcmp(argv[i], "--prec"))
            i++, prec = atol(argv[i++]);
        else if (!strcmp(argv[i], "--de"))
            i++, flag |= AJ_USE_DE;
        else if (!strcmp(argv[i], "--desame"))
            i++, flag |= AJ_DE_SAME;
        /* mth root */
        else if (!strcmp(argv[i], "--rootdef"))
            i++, flag |= AJ_ROOT_DEF;
        else if (!strcmp(argv[i], "--rootprod"))
            i++, flag |= AJ_ROOT_PROD;
        else if (!strcmp(argv[i], "--rootturn"))
            i++, flag |= AJ_ROOT_TURN;
        /* families */
        else if (!strcmp(argv[i], "--xn1"))
        {
            i++;
            n = atol(argv[i++]);
            f_pol = &pol_xn1;
        }
        else if (!strcmp(argv[i], "--bern"))
        {
            i++, n = atol(argv[i++]), f_pol = &pol_bern;
        }
        else if (!strcmp(argv[i], "--bernrev"))
        {
            i++, n = atol(argv[i++]), f_pol = &pol_bern, rev = 1;
        }
        else if (!strcmp(argv[i], "--pol"))
        {
            i++;
            if (!fmpz_poly_parse(poly, argv[i++]))
                abort();
            f_pol = NULL;
        }
        else if (!strcmp(argv[i], "--coeffs"))
        {
            int j;
            i++;
            n = atol(argv[i++]);
            for (j = 0; j <= n && i < argc; j++)
                fmpz_poly_set_coeff_si(poly, n - j, atol(argv[i++]));
            f_pol = NULL;
        }
        else if (!strcmp(argv[i], "--stoll"))
        {
            i++;
            fmpz_poly_set_coeff_si(poly, 6, 82342800);
            fmpz_poly_set_coeff_si(poly, 5, - 470135160);
            fmpz_poly_set_coeff_si(poly, 4, + 52485681);
            fmpz_poly_set_coeff_si(poly, 3, + 2396040466);
            fmpz_poly_set_coeff_si(poly, 2, + 567207969);
            fmpz_poly_set_coeff_si(poly, 1, - 985905640);
            fmpz_poly_set_coeff_si(poly, 0, 247747600);
            f_pol = NULL;
        }
        /* restrict computations / output */
        else if (!strcmp(argv[i], "--big"))
            i++, flag |= AJ_NO_TAU;
        else if (!strcmp(argv[i], "--int"))
            i++, flag |= AJ_NO_AB;
        else if (!strcmp(argv[i], "--tree"))
            i++, flag |= AJ_NO_INT;
        else if (!strcmp(argv[i], "--gp"))
            i++, f_print = &acb_mat_print_gp;
        else if (!strcmp(argv[i], "--error"))
            i++, f_print = &acb_mat_print_error;
        else if (!strcmp(argv[i], "--trim"))
            i++, flag |= AJ_TRIM;
        else if (!strcmp(argv[i], "--digits"))
            i++, digits = atol(argv[i++]);
        /* benchmark mode */
        else if (!strcmp(argv[i], "--bench"))
        {
            i++;
            run = atol(argv[i++]);
            print = 0;
        }
        else
        {
            flint_printf("unrecognized argument %s\n", argv[i]);
            return usage();
        }
    }
    if (!digits)
        digits = (slong)(prec * .301);

    /* compute pol to actual accuracy */
    if (f_pol)
        f_pol(poly, n, prec + n + 40);
    if (rev)
        fmpz_poly_reverse(poly, poly, n);

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
    fmpz_poly_clear(poly);
    flint_cleanup();
}
