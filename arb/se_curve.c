/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
#define TMP 0

void
sec_init(sec_t * c, slong m, acb_srcptr x, slong n)
{
#if TMP
    flint_printf("\n\ninit curve y^%ld = f_%ld(x) =",m,n);
#endif
    if (m < 2)
    {
        flint_printf("\nneed exponent m > 1 in se curve, got %ld\n", m);
        abort();
    }
    c->m = m;
    c->n = n;
    c->delta = n_gcd(m, n);
    c->g = ((m-1)*(n-1) - c->delta + 1)/ 2;
    c->roots = _acb_vec_init(n);
    _acb_vec_set(c->roots, x, n);
#if TMP
    _acb_vec_printd(c->roots, n, 30, ", ");
#endif
}

void
sec_clear(sec_t c)
{
    _acb_vec_clear(c.roots, c.n);
}
