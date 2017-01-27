/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* i,j s.t. dj >= mi + delta */
void
holomorphic_differentials(cohom_t dz, slong d, slong m)
{
    slong i, j, k, delta;
    delta = n_gcd(m, d);
    k = 0;
    for (j = 1 + (m + delta - 1) / d; j < m; j++)
    {
        for (i = 1; i < d && m * i <= j * d - delta; i++)
        {
            dz[k].x = i - 1;
            dz[k].y = j;
            k++;
        }
    }
    if (2 * k != (d-1)*(m-1) - delta + 1)
        abort();
}
