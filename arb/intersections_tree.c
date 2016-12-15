/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* set
   c[i+k][j+l] = 1 if k-l = sp mod m
   c[i+k][j+l] = -1 if k-l = sm mod m
 */

#define c(i,j) fmpz_mat_entry(c,i,j)

static void
fill_block(si_mat_t c, slong i, slong j, slong sp, slong sm, slong m)
{
    slong k, l;
    for (l = 0; l < m - 1; l++)
    {
        k = l + sp % m;
        if (k < m - 1)
        {
            *c(i + k, j + l) = 1;
            *c(j + l, i + k) = -1;
        }
        k = l + sm % m;
        if (k < m - 1)
        {
            *c(i + k, j + l) = -1;
            *c(j + l, i + k) = 1;
        }
    }
}

void
intersection_tree(si_mat_t c, const tree_t tree, slong d, slong m)
{
    slong k, l, size = m - 1;

    fmpz_mat_zero(c);

    /* the entry c[ k * (m-1) + s ] corresponds to loop alpha_i^(s) */
    for (k = 0; k < d - 1; k++)
    {
        edge_t ek = tree->e[k];

        /* intersection with self shifts */

        fill_block(c, k * size, k * size, 1, -1, m);

        /* intersection with others */
        for (l = k + 1; l < d - 1; l++)
        {
            edge_t el = tree->e[l];

            if (el.a != ek.a && el.a != ek.b)
                /* no intersection */
                continue;

            /* compute angle */

            /* compute shift */

            if(el.a == ek.a)
            {
                /* case ab.ad */
                double v = el.va - ek.va;
                fill_block(c, k * size, l * size, 1, -1, m);
            }
            else if (el.a == ek.b)
            {
                /* case ab.bd */
                double v = el.va - ek.vb;
                fill_block(c, k * size, l * size, 1, -1, m);
            }
        }
    }
}
