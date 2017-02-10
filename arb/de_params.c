/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

slong
de_params(double * h, acb_srcptr u, slong len, double r, slong i, slong m, slong prec)
{
    slong n, k;
    cdouble * w;
    w = malloc(len * sizeof(cdouble));
    for (k = 0; k < len; k++)
        w[k] = acb_get_cdouble(u + k);

    n = de_params_d(h, w, len, r, i, m, prec);

    free(w);
    return n;
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
    arb_root_ui(m, m, p->m, prec);
    acb_abs(abs, z, prec);
    arb_pow_ui(abs, abs, p->i, prec);
    arb_div(abs, abs, m, prec);

    arb_clear(abs);
    acb_clear(zu);
    acb_clear(z);
    return 1;
}

void
de_constant_m2(mag_t b, acb_srcptr u, slong len, double r, slong i, slong m, slong prec)
{
    arb_t abs, tmp;
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
    arb_const_pi(p.lambda, prec);
    arb_mul_2exp_si(p.lambda, p.lambda, -1);
    arb_init(abs);

    arf_init(tmin);
    arf_init(tmax);
    arf_zero(tmin);
    /* tmax = acosh(Pi/(2*l*sin(r))) */
    arb_init(tmp);
    arb_set_d(tmp, r);
    arb_sin(tmp, tmp, prec);
    arb_mul(tmp, tmp, p.lambda, prec);
    arb_mul_2exp_si(tmp, tmp, 1);
    arb_const_pi(abs, prec);
    arb_div(abs, abs, tmp, prec);
    arb_acosh(abs, abs, prec);
    arb_get_ubound_arf(tmax, abs, prec);

    arb_bound_func_arf(abs, (arb_func_t)&max_f, &p, tmin, tmax, 50, prec);
    arb_get_mag(b, abs);

    arb_clear(abs);
    arb_clear(p.r);
    arb_clear(p.lambda);
}

void
de_constant_m1(mag_t b, acb_srcptr u, slong len, double r, slong i, slong m, slong prec)
{
    slong k;
    arb_t abs;
    arb_t z1;
    arb_t tmp;

    arb_init(abs);
    arb_init(z1);
    arb_one(abs);
    arb_init(tmp);

    /* [0,1] */
    arb_union(z1, z1, abs, prec);

    for (k = 0; k < len; k++)
    {
        if (arb_overlaps(acb_realref(u + k), z1))
        {
            arb_mul(abs, abs, acb_imagref(u + k), prec);
        }
        else
        {
            acb_abs(tmp, u + k, prec);
            arb_mul(abs, abs, tmp, prec);
        }
    }
    arb_root_ui(abs, abs, m, prec);
    arb_inv(abs, abs, prec);
    arb_get_mag(b, abs);
    arb_clear(tmp);
    arb_clear(abs);
    arb_clear(z1);
}

#if 0
void
de_int_cst_b(mag_t b, mag_t tau, slong j, slong m)
{
    mag_t pi2, y0, xt, a2;
    mag_init(pi2);
    mag_init(y0);
    mag_init(xt);
    mag_init(a2);
    /* y0 = lambda * sin (tau) */
    mag_pi(pi2);
    mag_div_ui(pi2,pi2,2);
    mag_sin(y0, tau);
    mag_mul(y0, y0, pi2);
    /* xt = cos(tau)sqrt(Pi/(2y0)-1) */
    mag_div(y0,pi2,y0);
    mag_sub_ui(y0,y0,1);
    mag_sqrt(y0,y0);
    mag_cos(xt,tau);
    mag_mul(xt,y0);
    /* b = ( xt(1/cos(y0)^2a+1/xt^2a)+1/a*sinh(xt)^2a ) / cos(tau) */

    mag_clear(pi2);
    mag_clear(y0);
    mag_clear(xt);
    mag_clear(a2);
}
#endif

slong
de_params_tree(double * h, const tree_t tree, sec_t c, slong prec)
{
    return de_params(h, tree->data[tree->min].u, c.n, tree->r, c.n - 1, c.m, prec);
}
