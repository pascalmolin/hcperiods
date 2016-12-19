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
    flint_free(l.l);
}

void
homol_clear(homol_t l, slong len)
{
    slong k;

    for (k = 0; k < len; k++)
        loop_clear(l[k]);

    flint_free(l);
}
