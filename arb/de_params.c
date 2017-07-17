/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
acb_func_r_de(acb_t z, const acb_t u, arb_srcptr l, slong prec)
{
    acb_atanh(z, u, prec);
    acb_div_arb(z, z, l, prec);
    acb_asinh(z, z, prec);
}

void
acb_vec_r(arb_t r0, acb_ptr r, void (*f)(acb_t r, const acb_t u, arb_srcptr l, slong prec), acb_srcptr u, slong len, slong prec)
{
    slong k;
    arf_t t;
    arf_init(t);
    for (k = 0; k < len; k++)
    {
        f(r + k, u + k, NULL, prec);
        arb_get_abs_lbound_arf(t, acb_imagref(r + k), prec);
        if (!k || arf_cmp(t, arb_midref(r0)) < 0)
            arf_set(arb_midref(r0), t);
    }
    arf_clear(t);
}

/* precise bounds for actual integral */
typedef struct
{
    arb_t r;
    arb_t lambda;
    acb_srcptr u;
    slong len;
    slong i;
    slong m;
} params_t;

static int
max_f(arb_t m, const arb_t t, params_t * p, slong prec)
{
    slong k;
    acb_t z, zu;
    arb_t abs;

    arb_init(abs);
    acb_init(z);
    acb_init(zu);

    arb_set(acb_realref(z), t);
    arb_set(acb_imagref(z), p->r);
    acb_sinh(z, z, prec);
    acb_mul_arb(z, z, p->lambda, prec);
    acb_tanh(z, z, prec);

    arb_one(m);
    for (k = 0; k < p->len; k++)
    {
        acb_sub(zu, z, p->u + k, prec);
        if (acb_contains_zero(zu))
        {
            arb_clear(abs);
            acb_clear(zu);
            acb_clear(z);
            return 0;
        }
        acb_abs(abs, zu, prec);
        arb_mul(m, m, abs, prec);
    }
#if 1
    arb_root_ui(m, m, p->m, prec);
    acb_abs(abs, z, prec);
    arb_pow_ui(abs, abs, p->i, prec);
    arb_div(m, abs, m, prec);
#else
    acb_abs(m, z, prec);
#endif

    arb_clear(abs);
    acb_clear(zu);
    acb_clear(z);
    return 1;
}

void
de_constant_m2(arb_t m2, const arf_t l, acb_srcptr u, slong len, double r, slong i, slong m, slong prec)
{
    arb_t tmp;
    arf_t tmin;
    arf_t tmax;
    params_t p;
    p.m = m;
    p.i = i;
    p.u = u;
    p.len = len;
    arb_init(p.r);
    arb_set_d(p.r, r);
    arb_init(p.lambda);
    arb_set_arf(p.lambda, l);

    arf_init(tmin);
    arf_init(tmax);
    arf_zero(tmin);
    /* tmax = acosh(Pi/(2*l*sin(r))) */
    arb_init(tmp);
    arb_set_d(tmp, r);
    arb_sin(tmp, tmp, prec);
    arb_mul_arf(tmp, tmp, l, prec);
    arb_mul_2exp_si(tmp, tmp, 1);
    arb_const_pi(m2, prec);
    arb_div(m2, m2, tmp, prec);
    arb_acosh(m2, m2, prec);
    arb_get_ubound_arf(tmax, m2, prec);

#define MAXDEPTH 20
    arb_bound_func_arf(m2, (arb_func_t)&max_f, &p, tmin, tmax, 30, MAXDEPTH, prec);

    arb_clear(tmp);
    arb_clear(p.r);
    arb_clear(p.lambda);
}

void
de_constant_m1(arb_t m1, acb_srcptr u, slong len, slong m, slong prec)
{
    slong k;
    arb_t z1;
    arb_t tmp;

    arb_init(z1);
    arb_one(m1);
    arb_init(tmp);

    /* [0,1] */
    arb_union(z1, z1, m1, prec);

    for (k = 0; k < len; k++)
    {
        if (arb_overlaps(acb_realref(u + k), z1))
        {
            arb_mul(m1, m1, acb_imagref(u + k), prec);
        }
        else
        {
            acb_abs(tmp, u + k, prec);
            arb_mul(m1, m1, tmp, prec);
        }
    }
    arb_root_ui(m1, m1, m, prec);
    arb_inv(m1, m1, prec);

    arb_clear(tmp);
    arb_clear(z1);
}

void
de_constant_b(arb_t b, const arf_t l, double r, slong j, slong m, slong prec)
{
    fmpq_t a2;
    arb_t pi2, y0, cr, xr, t;

    fmpq_init(a2);
    arb_init(pi2);
    arb_init(y0);
    arb_init(cr);
    arb_init(xr);
    arb_init(t);

    /* a2 = 2(1-j/m) */
    fmpq_init(a2);
    fmpq_set_si(a2, 2*(m-j), m);

    /* y0 = lambda * sin (r) */
    arb_const_pi(pi2, prec);
    arb_mul_2exp_si(pi2, pi2, -1);

    arb_set_d(y0, r);
    arb_cos(cr, y0, prec);
    arb_sin(y0, y0, prec);
    arb_mul_arf(y0, y0, l, prec);

    /* xr = cos(r)*sqrt(Pi/(2y0)-1) */
    arb_div(xr, pi2, y0, prec);
    arb_sub_ui(xr, xr, 1, prec);
    arb_sqrt(xr, xr, prec);

    arb_mul(xr, cr, xr, prec);

    /* b = ( xr(1/cos(y0)^2a+1/xr^2a)+1/a*sinh(xr)^2a ) / cos(r) */

    arb_cos(t, y0, prec);
    arb_pow_fmpq(t, t, a2, prec);
    arb_inv(t, t, prec);
    arb_pow_fmpq(b, xr, a2, prec);
    arb_inv(b, b, prec);
    arb_add(b, b, t, prec);
    arb_mul(b, b, xr, prec);

    arb_sinh(t, xr, prec);
    arb_pow_fmpq(t, t, a2, prec);
    arb_mul_ui(t, t, m - j, prec);
    arb_div_ui(t, t, m, prec);
    arb_inv(t, t, prec);
    arb_mul_2exp_si(t, t, 1);
    arb_add(b, b, t, prec);

    arb_div(b, b, cr, prec);

    arb_clear(pi2);
    arb_clear(y0);
    arb_clear(cr);
    arb_clear(xr);
    arb_clear(t);
    fmpq_clear(a2);
}

void
de_error(mag_t e, const arf_t h, const arf_t l, slong n, double r, acb_srcptr u, slong len, slong d, slong j, slong m, slong prec)
{
    arb_t m1, m2, b, t;

    arb_init(m1);
    arb_init(m2);
    arb_init(b);
    arb_init(t);

    /* m2 */
    de_constant_m2(m2, l, u, len, r, d, m, prec);
    de_constant_b(b, l, r, j, m, prec);
    arb_mul(m2, m2, b, prec);

    arb_set_d(b, r);
    arb_const_pi(t, prec);
    arb_mul(t, t, b, prec);
    arb_mul_2exp_si(t, t, 1);
    arb_div_arf(t, t, h, prec);
    arb_exp(t, t, prec);
    arb_sub_ui(t, t, 1, prec);
    arb_div(m2, m2, t, prec);

    /* m1 */
    arb_set_arf(t, h);
    arb_mul_ui(t, t, n, prec);
    arb_sinh(t, t, prec);
    arb_mul_arf(t, t, l, prec);
    arb_const_log2(m1, prec);
    arb_sub(t, t, m1, prec);

    arb_mul_ui(t, t, 2*(m-j), prec);
    arb_div_ui(t, t, m, prec);
    arb_exp(t, t, prec);
    arb_mul_ui(t, t, (m-j), prec);
    arb_div_ui(t, t, m, prec);
    arb_mul_arf(t, t, l, prec);

    de_constant_m1(m1, u, len, m, prec);
    arb_div(m1, m1, t, prec);

    arb_add(m1, m1, m2, prec);
    arb_get_mag(e, m1);

    arb_clear(m1);
    arb_clear(m2);
    arb_clear(b);
    arb_clear(t);
}

slong
de_params(mag_t e, arf_t h, arf_t l, acb_srcptr u, slong len, double r, slong d, slong j, slong m, slong prec)
{
    slong n, k;
    cdouble * w;
    double hh, ll = 1.57, rr = .9 * r;
    const slong paramprec = 64;
    w = malloc(len * sizeof(cdouble));
    for (k = 0; k < len; k++)
        w[k] = acb_get_cdouble(u + k);

    n = de_params_d(&hh, &ll, &rr, w, len, d, m, prec);
    arf_set_d(h, hh);
    arf_set_d(l, ll);

    de_error(e, h, l, n, rr, u, len, d, j, m, paramprec);

    free(w);
    return n;
}


slong
de_params_tree(mag_t e, arf_t h, arf_t l, const tree_t tree, sec_t c, slong prec)
{
    return de_params(e, h, l, tree->data[tree->min].u, c.n - 2, tree->r, c.n - 1, c.j1, c.m, prec);
}
