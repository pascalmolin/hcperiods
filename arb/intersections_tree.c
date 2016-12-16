/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

#define TWOPI (2 * acos(-1.))
#define c(i,j) fmpz_mat_entry(c,i,j)

/* set
   c[i+k][j+l] = 1 if k-l = sp mod m
   c[i+k][j+l] = -1 if k-l = sm mod m
 */

static void
fill_block(fmpz_mat_t c, slong i, slong j, slong sp, slong sm, slong m)
{
    slong k, l;
    /* important: make sp and sm positive */
    sp = (sp % m + m) % m;
    sm = (sm % m + m) % m;
    for (l = 0; l < m - 1; l++)
    {
        k = (l + sp) % m;
        if (k < m - 1)
        {
            *c(i + k, j + l) = 1;
            *c(j + l, i + k) = -1;
        }
        k = (l + sm) % m;
        if (k < m - 1)
        {
            *c(i + k, j + l) = -1;
            *c(j + l, i + k) = 1;
        }
    }
}

void
intersection_tree(fmpz_mat_t c, const tree_t tree, slong d, slong m)
{
    slong k, l, size = m - 1;

    fmpz_mat_zero(c);

    /* the entry c[ k * (m-1) + s ] corresponds to the
       loop gamma_k^(s) */
    for (k = 0; k < d - 1; k++)
    {
        slong s;
        edge_t ek = tree->e[k];

        /* intersection with self shifts */

        fill_block(c, k * size, k * size, 1, -1, m);

        /* intersection with other shifts */
        for (l = k + 1; l < d - 1; l++)
        {
            edge_t el = tree->e[l];

            if (ek.a != el.a && ek.b != el.a)
                /* no intersection */
                continue;

            if(el.a == ek.a)
            {
                /* case ab.ad */
                s = lrint((el.va - ek.va ) / TWOPI);
                if (el.dir > ek.dir)
                    fill_block(c, k * size, l * size, -s, -1-s, m);
                else
                    fill_block(c, k * size, l * size, 1-s, -s, m);
            }
            else if (el.a == ek.b)
            {
                /* case ab.bd */
                s = lrint(.5 + (el.va - ek.vb ) / TWOPI);
                fill_block(c, k * size, l * size, -s, 1-s, m);
            }
            else
            {
                flint_printf("invalid tree\n");
                abort();
            }
        }
    }
}
