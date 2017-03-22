/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
loop_init(loop_t * l, slong len)
{
    l->n = len;
    l->l = flint_malloc(len * sizeof(gamma_t));
}

void
loop_clear(loop_t l)
{
    if (l.n)
        flint_free(l.l);
}

void
loop_print(const loop_t l)
{
    slong k;
    for (k = 0; k < l.n; k++)
        flint_printf("%s%ld*g_%ld^(%ld)",k?"+ ":"",l.l[k].coeff,l.l[k].index,l.l[k].shift);
}

void
homol_init(homol_t * cyc, slong len)
{
    slong k;
    * cyc = flint_malloc(len * sizeof(loop_t));
    for (k = 0; k < len; k++)
        (*cyc)[k].n = 0;
}

void
homol_clear(homol_t cyc, slong len)
{
    slong k;

    for (k = 0; k < len; k++)
        loop_clear(cyc[k]);

    flint_free(cyc);
}
