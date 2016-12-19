/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
sec_init(sec_t * c, slong m, acb_srcptr x, slong d)
{
    slong k, g;

    g = ((m-1)*(d-1) - n_gcd(m,d) + 1)/ 2;

    c->m = m;
    c->d = d;
    c->g = g;
    c->roots = _acb_vec_init(d);
    for (k = 0; k < d; k++)
        acb_set(c->roots + k, x + k);
}

void
sec_clear(sec_t c)
{
    _acb_vec_clear(c.roots, c.d);
}
