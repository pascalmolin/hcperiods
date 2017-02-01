/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
sec_init(sec_t * c, slong m, acb_srcptr x, slong n)
{
    slong k;

    c->m = m;
    c->n = n;
    c->delta = n_gcd(m, n);
    c->g = ((m-1)*(n-1) - c->delta + 1)/ 2;
    c->roots = _acb_vec_init(n);
    for (k = 0; k < n; k++)
        acb_set(c->roots + k, x + k);
}

void
sec_clear(sec_t c)
{
    _acb_vec_clear(c.roots, c.n);
}
