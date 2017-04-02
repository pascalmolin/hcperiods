/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

static void
arb_tanh_cosh2(arb_t t, arb_t c, const arb_t x, slong prec)
{
    arb_t e2x, p, m;
    arb_init(e2x);
    arb_init(m);
    arb_init(p);
    arb_mul_2exp_si(e2x, x, 1);
    arb_exp(e2x, e2x, prec);
    arb_add_si(p, e2x, 1, prec);

    arb_inv(e2x, e2x, prec);
    arb_one(m);
    arb_sub(m, m, e2x, prec);

    arb_add(c, p, m, prec);
    arb_div(t, m, p, prec);

    arb_clear(e2x);
    arb_clear(m);
    arb_clear(p);
}

void
de_int_init(de_int_t de, const arf_t h, const arf_t l, ulong n, const mag_t e, ulong m, slong prec)
{
    slong k;
    arb_t ah, lambda;
    arb_t kh, sh, ch, shsh, chsh;

    de->n = n;
    de->prec = prec;

    arf_init(de->l);
    arf_init(de->h);
    mag_init(de->e);

    arf_set(de->h, h);
    arf_set(de->l, l);
    mag_set(de->e, e);

    arb_init(lambda);
    arb_set_arf(lambda, l);

    arb_init(de->factor);
    arb_mul_arf(de->factor, lambda, de->h, prec);

    de->x = _arb_vec_init(n);
    de->dx = _arb_vec_init(n);
    de->ch2m = _arb_vec_init(n);

    arb_init(ah);
    arb_init(sh);
    arb_init(ch);
    arb_init(shsh);
    arb_init(chsh);
    arb_init(kh);

    arb_set_arf(ah, de->h);
    arb_zero(de->x + 0);
    arb_one(de->dx + 0);
    arb_one(de->ch2m + 0);

    for (k = 1; k < n; k++)
    {
        arb_mul_si(kh, ah, k, prec);
        arb_sinh_cosh(sh, ch, kh, prec);
        arb_mul(sh, sh, lambda, prec);
#if 1
        arb_sinh_cosh(shsh, chsh, sh, prec);
        arb_div(de->x + k, shsh, chsh, prec);
        arb_mul(chsh, chsh, chsh, prec);
        arb_div(de->dx + k, ch, chsh, prec);
        arb_root_ui(de->ch2m + k, chsh, m, prec);
#else
        arb_tanh_cosh2(de->x + k, de->ch2m + k, sh, prec);
        arb_div(de->dx + k, ch, de->ch2m + k, prec);
        arb_root_ui(de->ch2m + k, de->ch2m + k, m, prec);
#endif
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
    arf_clear(de->h);
    arf_clear(de->l);
    mag_clear(de->e);
    arb_clear(de->factor);
    _arb_vec_clear(de->x, de->n);
    _arb_vec_clear(de->dx, de->n);
    _arb_vec_clear(de->ch2m, de->n);
}
