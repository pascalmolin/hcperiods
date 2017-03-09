#include "abel_jacobi.h"

void
acb_randtest_exclude(acb_t z, acb_t b, flint_rand_t state, slong prec, slong mag_bits)
{
    do
    {
        acb_randtest_precise(z, state, prec, mag_bits);
    }
    while(acb_overlaps(z, b));
}

void
acb_vec_set_random_u(acb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits, double eps)
{
    slong k;

    acb_t b;
    acb_init(b);
    mag_one(arb_radref(acb_realref(b)));
    mag_set_d(arb_radref(acb_imagref(b)), eps);
    for (k = 0; k < len; k++)
        acb_randtest_exclude(u + k, b, state, prec, mag_bits);
}

void
acb_vec_set_random(acb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits)
{
    slong k;
    for (k = 0; k < len; k++)
        acb_randtest_precise(u + k, state, prec, mag_bits);
}

void
_acb_vec_printd(acb_srcptr u, slong len, slong d, const char * sep)
{
    slong k;

    for (k = 0; k < len; k++)
    {
        if (k)
            flint_printf(sep);
        acb_printd(u + k, d);
    }
}

void
_acb_vec_arf_printd(acb_srcptr u, slong len, slong d, const char * sep)
{
    slong k;

    for (k = 0; k < len; k++)
    {
#define re arb_midref(acb_realref(u + k))
#define im arb_midref(acb_imagref(u + k))
        if (k)
            flint_printf(sep);
        arf_printd(re, d);
        if (arf_sgn(im) >= 0)
            flint_printf(" + ");
        else
            flint_printf(" ");
        arf_printd(im, d);
        flint_printf(" * I");
    }
}
