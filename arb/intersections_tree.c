/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"
#include <complex.h>

/* set
   c[i+k][j+l] = 1 if k-l = sp mod m
   c[i+k][j+l] = -1 if k-l = sm mod m
 */

#define c(i,j) fmpz_mat_entry(c,i,j)
#define PI 2*asin(1.0)


int 
real_sgn(double x)
{
  if (x > 0.0) return 1;
  if (x < 0.0) return -1;
  return 0;
}


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
intersection_tree(si_mat_t c, const cdouble * w, const tree_t tree, slong d, slong m)
{
    slong k, l, shft, sp, sm, size = m - 1;
    double rho; /* angle arg((b-a)/(d-c)) <-- need w for this (low prec branch pts) */ 
    fmpz_mat_zero(c);

    /* the entry c[ k * (m-1) + s ] corresponds to cycle gamma_{k+1}^{(s)} */
    for (k = 0; k < d - 1; k++)
    {
        edge_t ek = tree->e[k];

        /* intersection with self shifts */
        fill_block(c, k * size, k * size, 1, -1, m);

        /* intersection with other cycles */
        for (l = k + 1; l < d - 1; l++)
        {
            edge_t el = tree->e[l];

            if (ek.a != el.a && ek.b != el.a)
                /* no intersection */
                continue;

            /* compute angle */
            rho = carg((w[ek.b]-w[ek.a])/(w[el.b]-w[el.a]));
                                        
            if (ek.a == el.a)
            {
                /* case a=c */
                shft = round((1/(2*PI)) * ( rho + el.va - ek.va ));
                if ( real_sgn(rho) > 0 )
                {
                  sp = 1 - shft;
                  sm = -shft; 
                }
                else
                {
                  sp = -shft;
                  sm = -1 - shft;
                }
               
            }
            else if (ek.b == el.b)
            {
                /* case b=d  */
                shft = round((1/(2*PI)) * ( rho + el.vb - ek.vb ));
                if ( real_sgn(rho) > 0 )
                {
                  sp = 1 - shft;
                  sm = -shft; 
                }
                else
                {
                  sp = -shft;
                  sm = -1 - shft;
                }
            }
            else if (ek.b == el.a)
            {
                /* case b=c  */
                shft = round((1/(2*PI)) * ( PI + rho + el.va - ek.vb ));
                sp = -shft;
                sm = 1 - shft; 
            }
            else if (ek.a == el.b)
            {
                /* case a=d  */
                shft = round((1/(2*PI)) * ( PI + rho + el.vb - ek.va ));
                sp = -shft;
                sm = 1 - shft; 
            }

            /* corresponding block matrix */
            fill_block(c, k * size, l * size, sp, sm, m);
        }
    }
}





