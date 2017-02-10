/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

void
holomorphic_differentials(cohom_t dz, slong n, slong m)
{
    slong i, j, k, delta;
    delta = n_gcd(m, n);
    k = 0;
    for (j = jmin(m, n, delta); j < m; j++)
    {
        for (i = 1; i < n && m * i <= n * j - delta; i++)
        {
            dz[k].x = i - 1;
            dz[k].y = j;
            k++;
        }
    }
    if (2 * k != (n-1)*(m-1) - delta + 1)
        abort();
}
