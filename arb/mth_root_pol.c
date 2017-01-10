/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

void
sqrt_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_one(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k < d1)
            acb_neg(t, t);
        acb_sqrt(t, t, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
mth_root_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_one(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k < d1)
            acb_neg(t, t);
        acb_root_ui(t, t, m, prec);
        acb_mul(y, y, t, prec);
    }

    acb_clear(t);
}

void
mth_root_pol_prod(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong m, slong prec)
{
    slong k;
    acb_t t;
    acb_init(t);

    acb_zero(y);
    for (k = 0; k < d; k++)
    {
        acb_sub_arb(t, u + k, x, prec);
        if (k < d1)
            acb_neg(t, t);
        acb_log(t, t, prec + 4);
        acb_add(y, y, t, prec + 4);
    }
    acb_div_ui(y, y, m, prec + 4);
    acb_exp(y, y, prec);

    acb_clear(t);
}


int
im_sgn(const acb_t x)
{
        if ( arb_is_nonnegative(acb_imagref(x)) )
		return 1;
	else if ( arb_is_negative(acb_imagref(x)) )
		return -1;
	else
		flint_printf("sign cannot be determined\n"); abort();
}

void
mth_root_pol_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec)
{
    slong q, k, isgn_s,  isgn_t;
    acb_t s, t;
    q = 0;
    acb_init(s);
    acb_init(t);
    
    acb_one(y);
    acb_sub_arb(s, u + 0, x, prec); 
    if ( d1 != 0 )
          acb_neg(s, s);
    isgn_s = im_sgn(s);
    for (k = 1; k < d; k++)
    { 
      acb_sub_arb(t, u + k, x, prec);
      if (k < d1)
            acb_neg(t, t);
      isgn_t = im_sgn(t);
      acb_mul(s, s, t, prec);
      if ( isgn_s == isgn_t )
      {
        isgn_s = im_sgn(s);
        if ( isgn_s != isgn_t )
        {
          if ( isgn_t > 0 )
            q++;
          else
            q--;
        }
      }
      else
         isgn_s = im_sgn(s);
    } 

    q = (q % m + m) % m;;
    acb_root_ui(s, s, m, prec);
    acb_mul(y, s, z + q, prec);

    acb_clear(s);
    acb_clear(t);
}
