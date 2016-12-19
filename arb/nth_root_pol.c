/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

void
nth_root_pol_def(acb_t y, acb_srcptr u, const arb_t x, slong d, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_one(y);
    for (k = 0; k < d; k++)
    {
        acb_set_arb(t, x);
        acb_sub(t, t, u + k, prec);
        acb_root_ui(t, t, m, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
nth_root_pol_prod(acb_t y, acb_srcptr u, const arb_t x, slong d, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_zero(y);
    for (k = 0; k < d; k++)
    {
        acb_set_arb(t, x);
        acb_sub(t, t, u + k, prec);
        acb_log(t, t, prec + 4);
        acb_add(y, y, t, prec + 4);
    }
    acb_div_ui(y, y, m, prec + 4);
    acb_exp(y, y, prec);

    acb_clear(t);
}

void
nth_root_pol_turn(acb_t y, acb_srcptr u, const arb_t x, acb_srcptr z, slong d, slong m, slong prec)
{
  
    slong q, k isgn_s, isgn_t, sgn_res;  
    acb_t s, t;
    acb_init(s);
    acb_init(t);
    q = 0;
    sgn_res = 0;
    acb_one(y);
    acb_set_arb(s, x);
    acb_sub(s, s, u + 0, prec);

    if ( acb_is_real(s) && arb_is_negative(acb_realref(s)) )
    {
        sgn_res++;
        q++;
    }
    arb_sgn(isgn_s,acb_imagref(s));

    for (k = 1; k < d; k++)
    { 
      acb_set_arb(t, x);
      acb_sub(t, t, u + k, prec);

      if ( acb_is_real(t) && arb_is_negative(acb_realref(t)) )
      {
        sgn_res++;
        q++;
      }
      arb_sgn(isgn_t,acb_imagref(t));
      acb_mul(s, s, t, prec);
      if ( arb_eq(isgn_s,isgn_t) )
      {
        arb_sgn(isgn_s,acb_imagref(s));
        if ( arb_ne(isgn_s,isgn_t) )
        {
          if ( arb_contains_nonnegative(isgn_t) )
          {
            q = q + 2;
          }
          else
          {
            q = q - 2;
          }
        }
      }
      else
      {
         arb_sgn(isgn_s,acb_imagref(s));
      }

      if ( acb_is_real(s) && arb_is_negative(acb_realref(s)) )
      {
        sgn_res++;
        q++;
      }
    } 
     
    acb_root_ui(s, s, m, prec);
    acb_set(t, z);
    acb_root_ui(t, t, 2, prec);
    acb_pow_ui(t, t, q, prec);
    acb_mul(y, y, t, prec);
    acb_mul(y, y, s, prec);
    acb_mul(y, y, sgn_res % 2);

    acb_clear(s);
    acb_clear(t);
}
