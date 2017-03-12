#include <string.h>

#include "flint/arith.h"
#include "acb_poly.h"
#include "abel_jacobi.h"

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
void
pol_bern(acb_poly_t pol, slong n, slong prec)
{
    fmpq_poly_t h;
    fmpq_poly_init(h);
    arith_bernoulli_polynomial(h, n);
    acb_poly_set_fmpq_poly(pol, h, prec);
    fmpq_poly_clear(h);
}
enum { XN1, EXP, BERN, COEFF };

int main(int argc, char * argv[])
{
    int i, print = 1, type = XN1;
    slong n = 9, m = 4, prec = 250, digits = 0;
    acb_poly_t poly;
    abel_jacobi_t aj;

    acb_poly_init(poly);

    if (argc < 2)
    {
        flint_printf("periods [-m m] [-n n] [--prec p]\n");
        flint_printf("Print period matrix of curve y^m = f_n(x).\n");
        flint_printf("Default m = 4, n = 9, f_n(x) = x^n + 1, prec = 250.\n");
        flint_printf("Other polynomials f_n(x):\n");
        flint_printf("  --bern n: Bernoulli polynomial Bn(x)");
        flint_printf("  --cheb n: Chebychev polynomial Tn(x)");
        flint_printf("  --coeff n cn ... c1 c0 : cn x^n + ... + c1 x + c0");
        flint_printf("Output options:\n");
        flint_printf("  --quiet: no printing\n");
        flint_printf("  --mid: do not print error balls\n");
        flint_printf("  --big: big period matrix\n");
        return 1;
    }

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--quiet"))
        {
            print = 0;
        }
        else if (!strcmp(argv[i], "-m"))
        {
            m = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-n"))
        {
            n = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "--prec"))
        {
            prec = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "--digits"))
        {
            digits = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "--exp"))
        {
            type = EXP;
            n = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "--bern"))
        {
            type = BERN;
            n = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "--coeff"))
        {
            int j;
            type = COEFF;
            i++;
            n = atol(argv[i++]);
            for (j = 0; j <= n && i < argc; j++, i++)
                acb_poly_set_coeff_si(poly, j, atol(argv[i]));
        }
    }
    if (!digits)
        digits = (slong)(prec * .301);

    if (type == XN1)
        pol_xn1(poly, n, prec);
    else if (type == EXP)
        pol_exp(poly, n, prec);
    else if (type == BERN)
        pol_bern(poly, n, prec);
    else if (type != COEFF)
        abort();

    abel_jacobi_init_poly(aj, m, poly, prec);
    abel_jacobi_compute(aj, prec);

    if (print)
        acb_mat_printd(aj->tau, digits);

    abel_jacobi_clear(aj);
    acb_poly_clear(poly);
    flint_cleanup();
}
