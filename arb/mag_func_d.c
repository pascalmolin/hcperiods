#include <math.h>

typedef double (* max_func_d)(void * p, double x);

double
golden_search(max_func_d f, void * p, double a, double b)
{
    const double eps = 1e-5;
    const double phi = (1+sqrt(5))/2;
    double c, d;
    c = b - (b - a) / phi;
    d = a + (b - a) / phi;
    while (fabs(c - d) > eps)
    {
        if (f(p, c) < f(p, d))
            b = d;
        else
            a = c;

        c = b - (b - a) / phi;
        d = a + (b - a) / phi;
    }
    return (b + a) / 2;
}

double
max_func_d_gs(max_func_d f, void * p, double a, double b, int steps)
{
    int k;
    double h, m, r = - 1;

    h = (b-a)/steps;

    b = a + h;
    m = golden_search(f, p, a, (b = a+h));
    for (k = 0; k < steps; k++)
    {
        b = a + h;
        m = golden_search(f, p, a, b);
        if (!k || m > r)
            r = m;
    }
    return r;
}
