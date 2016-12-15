/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
intersection_tree(inter_mat inter, tree_t tree, acb_srcptr x, slong d)
{
    slong k, l;

    for (k = 0; k < d - 1; k++)
    {
        edge_t ek = tree->e[k];
        for (l = k; l < d - 1; l++)
        {
            edge_t el = tree->e[l];
            if(el.a == ek.a)
            {
                /* case ab.ad */
            }
            else if (el.a == ek.b)
            {
                /* case ab.bd */
            }
        }
    }
    return;
}
