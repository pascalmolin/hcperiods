/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void differentials(cohom_t dz, slong d, slong m)
{
    slong i, j, k, delta;
    delta = n_gcd(m, d);
    k = 0;
    for (i = 1; i < d; i++)
    {
        for (j = 1; j < m; j++)
        {
            if (j * d >= m * i + delta && m * i > j)
            {
                dz[k].x = i - 1;
                dz[k].y = j;
                k++;
            }
        }
    }
}
