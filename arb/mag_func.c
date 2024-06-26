/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "mag_func.h"

/* b subset union( t + k * step, k = 0..n-1 ) */
void
arb_subdivide(arb_t t, arf_t step, const arb_t b, slong n, slong prec)
{
    double r;

#define mb arb_midref(b)
#define mt arb_midref(t)
#define rb arb_radref(b)
#define rt arb_radref(t)

#if 0
    arf_set_mag(mt, rb);
    arf_div_ui(mt, mt, n, prec, ARF_RND_NEAR);
    arf_mul_i(mt, mt, n - 1, prec, ARF_RND_NEAR);
    arf_sub(mt, mb, mt, prec, ARF_RND_NEAR);
#else
    r = mag_get_d(rb) / n;
    arf_set_d(mt, r * (-n+1));
    arf_add(mt, mt, mb, prec, ARF_RND_NEAR);
    arf_set_d(step, 2*r);
#endif

    /*
     choose radius s.t.
     * mt - r < mb - R
     * mt + (n-1) * step + r > mb + R
     * 2*r > step
     */
     mag_set_d(rt, r);
}

void
arb_bound_func_arb(arb_t m, arb_func_t f, void * params, const arb_t b, slong n, slong depth, slong prec)
{
    const int tolerance = -2;
    slong k;
    arb_t t, abs;
    arf_t step;

    arb_init(t);
    arb_init(abs);
    arf_init(step);
    if (mag_is_special(arb_radref(b)))
    {
        flint_printf("\nERROR: illegal ball in mag_func ");
        arb_printd(b, 20); flint_printf(" [depth = %d]\n",depth);
        flint_abort();
    }
    if (!depth)
    {
        flint_printf("\nERROR: limit depth exceeded in mag_func.\n Ball ");
        arb_printd(b, 20); flint_printf("\n");
        flint_abort();
    }

    arb_subdivide(t, step, b, n, prec);

    for (k = 0; k < n; k++)
    {
        if (f(abs, t, params, prec) == 0
                || !arb_is_finite(abs)
                || ( arb_rel_error_bits(abs) > tolerance
                     && mag_cmp_2exp_si(arb_radref(abs), 2) > 0)
           )
        {
#if VERBOSE > 1
            flint_printf("\n  ... "); arb_printd(t, 20);
            flint_printf(" failed "); arb_printd(abs, 20);
#endif
            arb_bound_func_arb(abs, f, params, t, 3, depth - 1, prec);
        }
        if (k == 0)
            arb_set(m, abs);
        else
            arb_union(m, m, abs, prec);
        arb_add_arf(t, t, step, prec);
    }

#if VERBOSE > 0
    flint_printf("\n  [%d] mag on ",depth); arb_printd(b, 10);
    flint_printf(" -> "); arb_printd(m, 10);
#endif

    arf_clear(step);
    arb_clear(t);
    arb_clear(abs);

}

void
arb_bound_func_arf(arb_t m, arb_func_t f, void * params, const arf_t tmin, const arf_t tmax, slong n, slong depth, slong prec)
{
    arb_t b;
    arb_init(b);
    arb_set_interval_arf(b, tmin, tmax, prec);
    arb_bound_func_arb(m, f, params, b, n, depth, prec);
#if VERBOSE > 0
    flint_printf("\n## arb bound "); arb_printd(b, 10);
    flint_printf(" -> "); arb_printd(m, 10);
#endif
    arb_clear(b);
}

slong
mag_func_arb(mag_t m, arb_func_t f, void * params, const arb_t b, slong n, slong prec)
{
    slong k, count;
    arb_t t, abs;
    arf_t step;
    mag_t m2;

    mag_init(m2);
    arb_init(t);
    arb_init(abs);
    arf_init(step);

    if (mag_is_special(arb_radref(b)) || mag_get_d(arb_radref(b)) < .00001)
        flint_abort();

    arb_subdivide(t, step, b, n, prec);

    count = n;
    for (k = 0; k < n; k++)
    {
#if VERBOSE
        int r;
        if ((r = f(abs, t, params, prec)) == 0 || !arb_is_finite(abs))
        {
            flint_printf("\n## refine ");
            arb_printd(t, 5);
            if (r == 0)
            {
                flint_printf(" -- function returned 0");
            }
            else
            {
                flint_printf(" -- eval =  ");
                arb_printd(abs, 5);
            }
            count += mag_func_arb(m2, f, params, t, 5, prec);
        }
#else
        if (f(abs, t, params, prec) == 0 || !arb_is_finite(abs))
            count += mag_func_arb(m2, f, params, t, 5, prec);
#endif
        else
            arb_get_mag(m2, abs);
        if (k == 0)
            mag_set(m, m2);
        else
            mag_max(m, m, m2);
        arb_add_arf(t, t, step, prec);
    }

    arf_clear(step);
    arb_clear(t);
    arb_clear(abs);
    mag_clear(m2);

    return count;
}
    
slong
mag_func_arf(mag_t m, arb_func_t f, void * params, const arf_t tmin, const arf_t tmax, slong n, slong prec)
{
    slong count;
    arb_t b;
    arb_init(b);
    arb_set_interval_arf(b, tmin, tmax, prec);
    count = mag_func_arb(m, f, params, b, n, prec);
    arb_clear(b);
    return count;
}
