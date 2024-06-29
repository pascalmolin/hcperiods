/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
#include <gr_poly.h>
#include <ulong_extras.h>
#define TMP 0

void
sec_init_poly(sec_t * c, slong m, const gr_poly_t pol, gr_ctx_t ctx)
{
    sec_init(c, m, gr_poly_length(pol, ctx) - 1);
    GR_MUST_SUCCEED(gr_poly_set(c->pol, pol, ctx)); /* keep a copy */
    c->ctx = &ctx;
}

void
sec_init(sec_t *c, slong m, slong n)
{
    slong j, delta;
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
    c->delta = delta = n_gcd(m, n);
    c->g = ((m-1)*(n-1) - delta + 1)/ 2;

    gr_poly_init(c->pol, c->ctx);

    /* differentials */
    c->j1 = jmin(m, n, delta);
    c->nj = m - c->j1;
    c->ni = flint_malloc(c->nj * sizeof(slong));
    for (j = 0; j < c->nj; j++)
        c->ni[j] = imax(c->j1 + j, m, n, delta);

#if DEBUG > 1
    flint_printf("\ndifferentials: j1 = %ld, ni = %ld", c->j1, c->ni[0]);
    for (j = 1; j < c->nj; j++)
        flint_printf(", %ld", c->ni[j]);
    flint_printf("\n");
#endif
}

void
sec_clear(sec_t c)
{
    gr_poly_clear(c.pol, c.ctx);
    flint_free(c.ni);
}
