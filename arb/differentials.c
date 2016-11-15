/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"

void differentials(cohom_t dz, slong d, slong n)
{
    slong i, j, k, m;
    m = n_gcd(n, d);
    k = 0;
    for (i = 1; i < d; i++)
    {
        for (j = 1; j < n; j++)
        {
            if (j * d >= n * i + m && n * i > j)
            {
                dz[k].x = i - 1;
                dz[k].y = j;
                k++;
            }
        }
    }
}
