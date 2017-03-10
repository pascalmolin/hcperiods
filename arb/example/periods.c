#include <string.h>
#include "abel_jacobi.h"

int main(int argc, char * argv[])
{
    int i, print = 1;
    slong n = 9, m = 4, prec = 250, digits = 0;
    acb_poly_t poly;
    abel_jacobi_t aj;

    if (0 && argc < 2)
    {
        flint_printf("periods [--quiet] [-m m] [-n n] [--prec p] [-c c0 c1 .. cn] [--big]\n");
        flint_printf("Print period matrix of curve y^m = f_n(x).\n");
        flint_printf("Default m = 4, n = 9, p = 250, f_n = x^n + 1");
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
        }
        else if (!strcmp(argv[i], "-n"))
        {
            n = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "--prec"))
        {
            prec = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "--digits"))
        {
            digits = atol(argv[i+1]);
        }
    }
    if (!digits)
        digits = (slong)(prec * .301);

    acb_poly_init(poly);
    acb_poly_fit_length(poly, n + 1);
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, n, 1);

    abel_jacobi_init_poly(aj, m, poly, prec);
    abel_jacobi_compute(aj, prec);

    if (print)
        acb_mat_printd(aj->tau, digits);

    abel_jacobi_clear(aj);
    acb_poly_clear(poly);
    flint_cleanup();
}
