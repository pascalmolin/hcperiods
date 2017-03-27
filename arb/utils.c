#include "abel_jacobi.h"

void
_acb_vec_add_error_mag(acb_ptr res, slong len, const mag_t e)
{
    slong k;
    for (k = 0; k < len; k++)
        acb_add_error_mag(res + k, e);
}

void
acb_randtest_exclude(acb_t z, acb_t b, flint_rand_t state, slong prec, slong mag_bits)
{
    do
        acb_randtest_precise(z, state, prec, mag_bits);
    while(acb_overlaps(z, b));
}

void
acb_vec_set_random_u(acb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits, double eps)
{
    slong k;

    acb_t b;
    acb_init(b);
    mag_set_d(arb_radref(acb_realref(b)), 1 + eps);
    mag_set_d(arb_radref(acb_imagref(b)), eps);
    for (k = 0; k < len; k++)
        acb_randtest_exclude(u + k, b, state, prec, mag_bits);
    acb_clear(b);
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
        int zr, zi;
#define re arb_midref(acb_realref(u + k))
#define im arb_midref(acb_imagref(u + k))
        if (k)
            flint_printf(sep);
        zr = arf_sgn(re);
        zi = arf_sgn(im);
        if (!zr && !zi)
            flint_printf("0");
        if (zr)
        {
            arf_printd(re, d);
            if (zi > 0)
                flint_printf(" + ");
            else if (zi < 0)
                flint_printf(" ");
        }
        if (zi)
        {
            arf_printd(im, d);
            flint_printf(" * I");
        }
    }
}

void
acb_mat_print_gp(const acb_mat_t m, slong digits)
{
    long i, nr, nc;
    nr = acb_mat_nrows(m);
    nc = acb_mat_ncols(m);
    flint_printf("[");
    for (i = 0; i < nr; i++)
    {
        if (i)
            printf(";");
        _acb_vec_arf_printd(m->rows[i], nc, digits, ", ");
    }
    flint_printf("]");
}

/* adapted from acb_vec_sort_pretty */
typedef int(*__compar_fn_t) (const void *, const void *);
int acb_cmp_lex(const acb_t a, const acb_t b)
{
    arb_t t;
    int res;
    arb_init(t);
    res = 0;
	arb_sub(t, acb_realref(a), acb_realref(b), MAG_BITS);
    if (arb_contains_zero(t))
		arb_sub(t, acb_imagref(a), acb_imagref(b), MAG_BITS);
	res = arb_is_positive(t) ? 1 : -1;
    arb_clear(t);
    return res;
}
void _acb_vec_sort_lex(acb_ptr vec, slong len)
{
    qsort(vec, len, sizeof(acb_struct), (__compar_fn_t) acb_cmp_lex);
}
