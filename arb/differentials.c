/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

/* i,j s.t. dj >= mi + delta */
void
holomorphic_differentials(cohom_t dz, slong n, slong m)
{
    slong i, j, k, delta;
    delta = n_gcd(m, n);
    k = 0;
    for (j = 1 + (m + delta - 1) / n; j < m; j++)
    {
        for (i = 1; i < n && m * i <= j * n - delta; i++)
        {
            dz[k].x = i - 1;
            dz[k].y = j;
            k++;
        }
    }
    if (2 * k != (n-1)*(m-1) - delta + 1)
        abort();
}
