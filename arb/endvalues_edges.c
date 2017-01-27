/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
endvalues_edge(double * va, double * vb, double * dir,
        const cdouble * w, slong ia, slong ib, slong d)
{
    slong k;
    double a, b;
    cdouble ba = w[ib] - w[ia];

    *dir = carg(ba);
    a = b = (d - 1) * (*dir);
    
    for (k = 0; k < d; k++)
    {
        if (k == ia || k == ib)
            continue;
        a += carg((w[ia] - w[k]) / ba);
        b += carg((w[ib] - w[k]) / ba);
    }
    *va = a;
    *vb = b;
}

void
shift_info_tree(tree_t tree, cdouble * w, slong d)
{
    slong k;

    for (k = 0; k < tree->n; k++)
         /* compute endvalues for shifting numbers */
        endvalues_edge(&tree->e[k].va, &tree->e[k].vb, &tree->e[k].dir, w, tree->e[k].a, tree->e[k].b, d);

}
