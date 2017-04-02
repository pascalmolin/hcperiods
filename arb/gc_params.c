/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
arb_func_r_gc(arb_t r, const acb_t u, arb_srcptr l, slong prec)
{
    acb_t z;
    arb_t t;
    acb_init(z);
    arb_init(t);
    acb_add_si(z, u, 1, prec);
    acb_abs(t, z, prec);
    acb_add_si(z, u, -1, prec);
    acb_abs(r, z, prec);
    arb_add(r, r, t, prec);
    arb_mul_2exp_si(r, r, -1);
    /* check */
    arb_sub_ui(t, r, 1, prec);
    if (!arb_is_positive(t))
    {
        flint_printf("\n\nERROR: not enough precision in parameters, r <= 1 for ");
        acb_printd(u, 10);
        abort();
    }
    arb_clear(t);
    acb_clear(z);
}

void
arb_vec_r(arb_t r0, arb_ptr r, void (*f)(arb_t r, const acb_t u, arb_srcptr l, slong prec), acb_srcptr u, slong len, slong prec)
{
    slong k;
    arf_t t;
    arf_init(t);
    for (k = 0; k < len; k++)
    {
        f(r + k, u + k, NULL, prec);
        arb_get_abs_lbound_arf(t, r + k, prec);
        if (!k || arf_cmp(t, arb_midref(r0)) < 0)
            arf_set(arb_midref(r0), t);
    }
    arf_clear(t);
}

slong
maj_yr(arb_t m, const arb_t m0, const arb_t r, arb_srcptr vr, slong len, slong prec)
{
    slong k, K;

    arb_t d;
    arb_init(d);

    if (!arb_is_finite(r))
        abort();

    arb_one(m);
    for (K = 0, k = 0; k < len; k++)
    {
        if (arb_overlaps(r, vr + k))
            K++;
        else
        {
            arb_sub_arf(d, vr + k, arb_midref(r), prec);
            arb_mul(m, m, d, prec);
        }
    }
    if (!arb_is_finite(m))
        abort();
    /* should mult m by sh(r0)^K = derivative */
    arb_log(m, m, prec);
    arb_mul_2exp_si(m, m, -1);
    arb_sub(m, m0, m, prec);

    arb_clear(d);

    return K;
}

void
choose_r(arb_t r, const arb_t r0, const arb_t m, slong K, slong prec)
{
#if 0
    double R, A, M, eps;
    R = acosh( arf_get_d(arb_midref(r0), ARF_RND_DOWN) );
    flint_printf("\n## R0 = %lf from r0 = ", R); arb_printd(r0, 10);
    M = arf_get_d(arb_midref(m), ARF_RND_UP);
    A = 1 + 2 * M / K;
    eps = R / (A + log( A / R ));
    flint_printf("\n## choose eps = %lf", eps);
    R = cosh(R * (1. - eps));
    flint_printf("\n## new R = %lf", R);
    if (R < 1)
        abort();
    arb_zero(r);
    arf_set_d(arb_midref(r), R);
#else
    arb_t A;
    arb_init(A);
    /* r = r0 (1-eps), eps*(A-log(eps)) = r0, A = 1+2*m/K */
    /* A = 1+2m/K */
    arb_get_ubound_arf(arb_midref(A), m, prec);
    arb_div_si(A, A, K, prec);
    arb_mul_2exp_si(A, A, 1);
    arb_add_si(A, A, 1, prec);
#if 0
    flint_printf("\n## A = 1+2m/K = "); arb_printd(A, 10);
    flint_printf("\n## r0 = "); arb_printd(r0, 10);
    /* eps = r0 / W(r0*exp(-A)) */
    arb_acosh(r, r0, prec);
    arb_exp(A, A, prec);
    arb_div(A, r, A, prec);
    flint_printf("\n## take W of r/e^A = "); arb_printd(A, 10);
    arb_lambertw(A, A, 0, prec);
    flint_printf("\n## W = "); arb_printd(A, 10);
    arb_div(A, r, A, prec); /* epsilon */
    arb_submul(r, r, A, prec);
    flint_printf("\n## R = "); arb_printd(r, 10);
    arb_cosh(r, r, prec);
    mag_zero(arb_radref(r));
    arb_clear(A);
#else
    /* A = r0/(A+log(A/r0)) */
    arb_acosh(r, r0, prec);
    arb_div(r, A, r, prec);
    arb_log(r, r, prec);
    arb_add(A, r, A, prec);
    arb_div(A, r0, A, prec);
    arb_sub(r, r0, A, prec);
    arb_cosh(r, r, prec);
    mag_zero(arb_radref(r));
    arb_clear(A);
#endif
#endif
}

slong
gc_params(mag_t e, acb_srcptr u, slong len, slong d, slong prec)
{
    slong n, K;
    slong pp = 32;
    acb_t z;
    arf_t t;
    arb_t r0, r, m0, m;
    arb_ptr vr;

    acb_init(z);
    arf_init(t);
    arb_init(r0);
    arb_init(r);
    arb_init(m0);
    arb_init(m);

    /* compute rk and r0 */
    vr = _arb_vec_init(len);
    arb_vec_r(r0, vr, arb_func_r_gc, u, len, pp);

    /* m0 = D + log(2pi) + d*log(r0) */
    arb_const_log_sqrt2pi(m0, pp);
    arb_mul_2exp_si(m0, m0, 1);
    arb_log(m, r0, pp);
    arb_addmul_si(m0, m, d, pp);
    arb_const_log2(m, pp);
    arb_addmul_si(m0, m, prec, pp);

    /* fixme: let the user provide some r0? */
    /* choose r */
#if 1
    arb_add_si(r, r0, 1, prec);
    arb_mul_2exp_si(r, r, -1);
    mag_zero(arb_radref(r));
#else
    arb_set(r, r0);
    mag_set_d(arb_radref(r), 2./prec);
#endif

    while ((K = maj_yr(m, m0, r, vr, len, pp)))
    {
#if DEBUG > 1
        flint_printf("\nr = "); arb_printd(r, 10);
        flint_printf(" -> K = %ld && m = ", K); arb_printd(m, 10);
#endif
        choose_r(r, r0, m, K, pp);
    }
#if DEBUG > 1
    flint_printf("\nchoose r = "); arb_printd(r, 10);
    flint_printf(" -> m = "); arb_printd(m, 10);
#endif

    /* final result */
    arb_acosh(r, r, pp); /* r' = log(r+sqrt(r^2-1)) */
    arb_div(m, m, r, pp);
    arb_mul_2exp_si(m, m, -1);
    arb_get_ubound_arf(t, m, pp);
    n = arf_get_si(t, ARF_RND_CEIL);

    /* error */
    if (e)
    {
        /* exp(m)/(exp(2nr')-1) */
        arb_exp(m, m, pp);
        arb_mul_si(r, r, 2 * n, pp);
        arb_exp(r, r, pp);
        arb_sub_ui(r, r, 1, pp);
        arb_div(m, m, r, pp);
        arb_get_mag(e, m);
    }

    _arb_vec_clear(vr, len);
    acb_clear(z);
    arf_clear(t);
    arb_clear(r0);
    arb_clear(r);
    arb_clear(m0);
    arb_clear(m);

    return n;
}

slong
gc_params_tree(mag_t e, const tree_t tree, sec_t c, slong prec)
{
    return gc_params(e, tree->data[tree->min].u, c.n - 2, c.n - 1, prec);
}
