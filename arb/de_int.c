/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

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

void
de_int_params_se(arf_t h, slong * n, double tau, mag_t M1, mag_t M2, slong prec)
{
    mag_t D;

}
#endif

void
de_int_params(arf_t h, ulong *n, const tree_t tree, sec_t c, slong prec)
{
    arf_set_d(h, .1);
    *n = 100;
}

void
de_int_init(de_int_t de, arf_t h, ulong n, slong prec)
{
    slong k;
    arb_t ah, lambda;
    arb_t kh, sh, ch, shsh, chsh;

    de->n = n;
    de->prec = prec;

    arb_init(de->h);
    arb_set_arf(de->h, h);
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
