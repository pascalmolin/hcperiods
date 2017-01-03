/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include <arb.h>

void
arb_max_func_arf(arb_t m, int (* f)(arb_t abs, const arb_t, void * params, slong prec), void * params, arf_t tmin, arf_t tmax, slong n, slong prec)
{
    arb_t t, abs;
    arf_t step, tmp;
    slong k;

    arb_init(t);

    arf_init(step);
    arf_sub(step, tmax, tmin, prec, ARF_RND_NEAR);
    arf_div_ui(step, step, n, prec, ARF_RND_NEAR);

    arf_init(tmp);
    arf_sub(tmp, tmax, tmin, prec, ARF_RND_UP);
    arf_div_ui(tmp, tmp, n, prec, ARF_RND_UP);

    /* t0 = tmin + t/2 + B(eps)*/
#if 1
    arf_mul_2exp_si(arb_midref(t), step, -1);
    arf_add(arb_midref(t), tmin, arb_midref(t), prec, ARF_RND_NEAR);
    arf_get_mag(arb_radref(t), tmp);
#else
    arf_init(t0);
    arf_mul_2exp_si(t0, step, -1);
    arf_add(t0, t0, tmin, prec, ARF_RND_NEAR);
    mag_init(eps);
    arb_get_mag(eps, tmp);
#endif

    for (k = 0; k < n; k++)
    {
        if (!f(abs, t, params, prec) || !arb_is_finite(abs))
        {
            arf_t tmin, tmax;
            arb_get_lbound_arf(tmin, t, prec);
            arb_get_ubound_arf(tmax, t, prec);
            arb_max_func_arf(abs, f, params, tmin, tmax, 3, prec);
        }
        if (k == 0)
            arb_set(m, abs);
        else
            arb_union(m, m, abs, prec);
        arb_add_arf(t, t, step, prec);
    }
}

void
arb_max_func_arb(arb_t m, int (* f)(arb_t abs, const arb_t, void * params, slong prec), void * params, arb_t tmin, arb_t tmax, slong n, slong prec)
{
    arf_t t0, t1;
    arf_init(t0);
    arf_init(t1);
    arb_get_lbound_arf(t0, tmin, prec);
    arb_get_ubound_arf(t1, tmax, prec);
    arb_max_func_arf(m, f, params, t0, t1, n, prec);
}
