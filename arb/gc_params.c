/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

slong
maj_yr(arb_t m, const arb_t m0, const arb_t r, arb_srcptr vr, slong len, slong prec)
{
    slong k, K;

    arb_t chr0, chr;
    arb_init(chr0);
    arb_init(chr);

    arb_one(m);
    arb_set_arf(chr0, arb_midref(r));
    arb_cosh(chr0, chr0, prec);
    for (K = 0, k = 0; k < len; k++)
    {
        if (arb_overlaps(r, vr + k))
            K++;
        else
        {
            arb_cosh(chr, vr + k, prec);
            arb_sub(chr, chr, chr0, prec);
            arb_mul(m, m, chr, prec);
        }
    }
    if (K)
    {
        arb_set_arf(chr, arb_midref(r));
        arb_sinh(chr, chr, prec);
        arb_pow_ui(chr, chr, K, prec);
        arb_mul(m, m, chr, prec);
    }
    arb_log(m, m, prec);
    arb_mul_2exp_si(m, m, -1);
    arb_sub(m, m0, m, prec);

    arb_clear(chr0);
    arb_clear(chr);

    return K;
}

void
choose_r(arb_t r, const arb_t r0, const arb_t m, slong K, slong prec)
{
#if 1
    arb_zero(r);
    arf_set_d(arb_midref(r), .9 * arf_get_d(arb_midref(r0), ARF_RND_DOWN));
#else
    arb_t A;
    arb_init(A);
    /* r = r0 (1-eps), eps*(A-log(eps)) = r0, A = 1+2*m/K */
    arb_div_si(A, m, K, prec);
    arb_mul_2exp_si(A, A, 1);
    arb_add_si(A, A, 1, prec);
    arb_div(r, A, m, prec);
    arb_add(r, r, A, prec);
    arb_div(A, r0, r, prec);
    arb_sub(r, r0, A, prec);
    /* mag_zero(arb_radref(r)); */
    arb_clear(A);
#endif
}

slong
gc_params(mag_t e, acb_srcptr u, slong len, slong d, slong prec)
{
    slong n, k, K;
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
    for (k = 0; k < len; k++)
    {
#if 0
        acb_asin(z, u + k, pp);
        arb_abs(vr + k, acb_imagref(z));
#else
        acb_add_si(z, u + k, 1, pp);
        acb_abs(vr + k, z, pp);
        acb_add_si(z, u + k, -1, pp);
        acb_abs(r, z, pp);
        arb_add(vr + k, vr + k, r, pp);
        arb_mul_2exp_si(vr + k, vr + k, -1);
        arb_acosh(vr + k, vr + k, pp);
#endif
        if (!arb_is_finite(vr + k))
        {
            flint_printf("\nERROR: could not compute r for u = "); acb_printd(u + k, 20);
            abort();
        }
        arb_get_abs_lbound_arf(t, vr + k, pp);
        if (!k || arf_cmp(t, arb_midref(r0)) < 0)
            arf_set(arb_midref(r0), t);
    }
#if DEBUG > 1
    flint_printf("\nr0 = "); arb_printd(r0, 10);
#endif

    /* m0 = D + log(2pi) + d*r0 */
    arb_const_log_sqrt2pi(m0, pp);
    arb_mul_2exp_si(m0, m0, 1);
    arb_addmul_si(m0, r0, d, pp);
    arb_const_log2(m, pp);
    arb_addmul_si(m0, m, prec, pp);

    /* fixme: let the user provide some r0? */
    /* choose r */
    arb_set(r, r0);
    mag_set_d(arb_radref(r), 2./prec);

    while ((K = maj_yr(m, m0, r, vr, len, pp)))
        choose_r(r, r0, m, K, pp);

    /* final result */
    arb_div(m, m, r, pp);
    arb_mul_2exp_si(m, m, -1);
    arb_get_ubound_arf(t, m, pp);
    n = arf_get_si(t, ARF_RND_CEIL);

    /* error */
    if (e)
    {
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
    return gc_params(e, tree->data[tree->min].u, c.n, c.n - 1, prec);
}
