/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/*
void

de_params(double * h, ulong *n,
*/

void
de_int_init(de_int_t de, double h, ulong n, slong prec)
{
    slong k;
    arb_t ah, lambda;
    arb_t kh, sh, ch, shsh, chsh;

    de->n = n;
    de->prec = prec;

    arb_init(ah);
    arb_set_d(ah, h);
    arb_init(lambda);
    arb_const_pi(lambda, prec);
    arb_mul_2exp_si(lambda, lambda, -1);

    arb_init(de->factor);
    arb_mul(de->factor, lambda, ah, prec);

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
de_init_clear(de_int_t de)
{
    _arb_vec_clear(de->x, de->n);
    _arb_vec_clear(de->dx, de->n);
    arb_clear(de->factor);
}
