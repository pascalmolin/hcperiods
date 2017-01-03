/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
#include "complex_extras.h"


/* rough estimates using doubles */
double
de_distance_d(double x, double y, double r)
{
    /* FIXME: this is not the real distance... */
    cdouble xItau;
    xItau = casinh(catanh(x+I*y)/LAMBDA);
    return fabs(cimag(xItau))-r;
}

double
de_constant_m_d(const cdouble * w, slong len, double r, slong d, slong m)
{
    /* 2*Pi*r^d/sqrt(prod d(w_i,e_r)) */
    slong k;
    double p;
    const double eps = 1.e-5;
    p = 1;
    for (k = 0; k < len; k++)
    {
        double d = de_distance_d(creal(w[k]), cimag(w[k]), r);
        p *= (d < eps) ? r : d;
    }
    return pow(r, d) / pow(p, 1/m);
}

static double
de_constant_b(double r, double alpha)
{
    double xr, cr, y0, a2, b;
    a2 = 2*alpha;
    y0 = LAMBDA * sin(r);
    cr = cos(r);
    xr = cr*sqrt(PI2/y0-1);
    b = xr / 2 * (pow(cos(y0), -a2) + pow(xr, -a2))
        + pow(sinh(xr), -a2) / a2;
    return 2/cr*b;
}

slong
de_params_d(double *ph, const cdouble * w, slong len, double r, slong i, slong m, slong prec)
{
    slong n;
    double h, D, M1, M2, B, alpha;

    D = prec * LOG2;
    alpha = (m-1) / m;
    M1 = 1;
    M2 = 1;
    /* TODO: heuristic choice of r */
    B = de_constant_b(r, alpha);

    h = 2*PI*r / (D + log(2*M2*B+1));
    n = asinh((D+log(32*M1/alpha))/(2*alpha*LAMBDA));

    * ph = h;
    return n;
}

slong
de_params(double *ph, acb_srcptr u, slong len, double r, slong i, slong m, slong prec)
{
    slong n, k;
    cdouble * w;
    w = malloc(len * sizeof(cdouble));
    for (k = 0; k < len; k++)
        w[k] = acb_get_cdouble(u + k);

    n = de_params_d(ph, w, len, r, i, m, prec);

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
            return 0;
        acb_abs(abs, zu, prec);
        arb_mul(m, m, abs, prec);
    }
    arb_root_ui(m, m, p->m, prec);
    acb_abs(abs, z, prec);
    arb_pow_ui(abs, abs, p->i, prec);
    arb_div(abs, abs, m, prec);
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

    arb_max_func_arf(abs, &max_f, &p, tmin, tmax, 50, prec);
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
de_int_params(arf_t h, acb_srcptr u, slong len, double r, sec_t c, slong prec)
{
    arf_set_d(h, .1);
    return 100;
}

slong
de_int_params_tree(arf_t h, const tree_t tree, sec_t c, slong prec)
{
    arf_set_d(h, .1);
    return 100;
}

void
de_int_init(de_int_t de, arf_t h, ulong n, slong prec)
{
    slong k;
    arb_t ah, lambda;
    arb_t kh, sh, ch, shsh, chsh;

    de->n = n;
    de->prec = prec;

    arf_init(de->h);
    arf_set(de->h, h);
    arb_init(lambda);
    arb_const_pi(lambda, prec);
    arb_mul_2exp_si(lambda, lambda, -1);

    arb_init(de->factor);
    arb_mul_arf(de->factor, lambda, h, prec);

    de->x = _arb_vec_init(n);
    de->dx = _arb_vec_init(n);


    arb_init(sh);
    arb_init(ch);
    arb_init(shsh);
    arb_init(chsh);
    arb_init(kh);

    arb_zero(de->x + 0);
    arb_set(de->dx + 0, de->factor);

    for (k = 1; k < n; k++)
    {
        arb_mul_si(kh, ah, k, prec);
        arb_sinh_cosh(sh, ch, kh, prec);
        arb_mul(sh, sh, lambda, prec);
        arb_sinh_cosh(shsh, chsh, sh, prec);
        arb_div(de->x + k, shsh, chsh, prec);
        arb_div(de->dx, ch, shsh, prec);
    }

    arb_clear(sh);
    arb_clear(ch);
    arb_clear(shsh);
    arb_clear(chsh);
    arb_clear(kh);
    arb_clear(lambda);
    arb_clear(ah);
}

void
de_int_clear(de_int_t de)
{
    _arb_vec_clear(de->x, de->n);
    _arb_vec_clear(de->dx, de->n);
    arb_clear(de->factor);
}
