/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
symplectic_basis(homol_t loop_a, homol_t loop_b, const tree_t tree, sec_t c)
{
    slong a, len = (c.d-1)*(c.m-1);
    si_mat_t m, p;
    si_mat_init(m, len, len);
    si_mat_init(p, len, len);
    /* compute big intersection matrix, size len = (d-1)*(m-1) */
    intersection_tree(m, tree, c.d, c.m);
    symplectic_reduction(p, m, c.g, len);
    for (a = 0; a < c.g; a++)
    {
        slong len, i, s; /* index and shift */
        /* let loop_a[a] be the appropriate sum */
    }
    si_mat_clear(m, len, len);
    si_mat_clear(p, len, len);
}
