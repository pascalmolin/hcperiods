/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"

void
ab_points_worst(cdouble * w, const tree_t tree, sec_t c)
{
    slong k, l, i, j;
    cdouble a, b;

    /* small approx of points */
    for (k = 0; k < c.d; k++)
        w[k] = acb_get_cdouble(c.roots + k);

    /* get worst edge */
    k = tree->min;

    i = tree->e[k].a;
    j = tree->e[k].b;

    a = w[i];
    b = w[j];

    for (k = 0, l = 0; k < c.d; k++)
    {
        if (k == i || k == j)
            continue;
        else
            w[l++] = (2*w[k]-a-b)/(b-a);
    }
}

