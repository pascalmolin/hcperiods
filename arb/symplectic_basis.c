/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* compute big intersection matrix, size len = (d-1)*(m-1) */
void
big_intersection_matrix(si_mat_t c, inter_mat inter, slong d, slong m)
{
    slong i, ni, s, ns; /* index and shift */
    ni = d - 1;
    ns = m - 1;
    /* entry i * ni + s corresponds to loop alpha_i^(s) */
    for (i = 0; i < ni; i++)
    {
        /* intersection with self shifts */
        /* intersection with others */
    }
}

void
symplectic_basis(homol_t loop_a, homol_t loop_b, const inter_mat inter, sec_t c)
{
    slong a, len = (c.d-1)*(c.m-1);
    si_mat_t m, p;
    si_mat_init(m, len, len);
    si_mat_init(p, len, len);
    big_intersection_matrix(m, inter, c.d, c.m);
    symplectic_reduction(p, m, c.g, len);
    for (a = 0; a < c.g; a++)
    {
        slong len, i, s; /* index and shift */
        /* let loop_a[a] be the appropriate sum */
    }
    si_mat_clear(m, len, len);
    si_mat_clear(p, len, len);
}
