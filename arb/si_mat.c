/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "flint.h"
#include "si_mat.h"

si_mat_t
si_mat_init(long nr, long nc)
{
    si_mat_t m;
    long i, j;
    m = flint_malloc(nr * sizeof(long *));
    for (i = 0; i < nr; i++)
        m[i] = flint_malloc(nc * sizeof(long));
    for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
            m[i][j] = 0;
    return m;
}

void
si_mat_clear(si_mat_t m, long nr, long nc)
{
    long i;
    for (i = 0; i < nr; i++)
        flint_free(m[i]);
    flint_free(m);
}

void
row_swap(si_mat_t m, long i1, long i2, long len)
{
    long * c = m[i1];
    m[i1] = m[i2];
    m[i2] = c;
}

void
col_swap(si_mat_t m, long j1, long j2, long len)
{
    long i;
    for (i = 0; i < len; i++)
    {
        int t = m[i][j1];
        m[i][j1] = m[i][j2];
        m[i][j2] = t;
    }
}

void
col_neg(si_mat_t m, long j, long len)
{
    long i;
    for (i = 0; i < len; i++)
        m[i][j] = -m[i][j];
}

void
row_neg(si_mat_t m, long i, long len)
{
    long j;
    for (j = 0; j < len; j++)
        m[i][j] = -m[i][j];
}

void
col_addmul(si_mat_t m, long j, long k, long v, long len)
{
    long i;
    for (i = 0; i < len; i++)
        m[i][j] += v * m[i][k];
}

void
row_addmul(si_mat_t m, long i, long k, long v, long len)
{
    long j;
    for (j = 0; j < len; j++)
        m[i][j] += v * m[k][j];
}

/* cj <- a cj + b ck
   ck <- c cj + d ck */
void
row_bezout(si_mat_t m, long i1, long i2, long a, long b, long c, long d, long len)
{
    long j;
    for (j = 0; j < len; j++)
    {
        long s = m[i1][j], t = m[i2][j];
        m[i1][j] = a * s + b * t;
        m[i2][j] = c * s + d * t;
    }
}
void
col_bezout(si_mat_t m, long j1, long j2, long a, long b, long c, long d, long len)
{
    long i;
    for (i = 0; i < len; i++)
    {
        long s = m[i][j1], t = m[i][j2];
        m[i][j1] = a * s + b * t;
        m[i][j2] = c * s + d * t;
    }
}

void
si_mat_set_id(si_mat_t p, long len)
{
    long i, j;
    for (i = 0; i < len; i++)
    {
        for (j = 0; j < len; j++)
            p[i][j] = (i==j);
    }
}

int
si_mat_eq(si_mat_t a, si_mat_t b, long nr, long nc)
{
    long i, j;
    for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
            if (a[i][j] != b[i][j])
                return 0;
    return 1;
}

void
si_mat_print(si_mat_t m, long nr, long nc)
{
    long i, j;
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
            printf("%2d", m[i][j]);
        printf("\n");
    }
}
    
void
si_mat_print_gp(si_mat_t m, long nr, long nc)
{
    long i, j;
    flint_printf("[");
    for (i = 0; i < nr; i++)
    {
        if (i)
            printf(";");
        for (j = 0; j < nc; j++)
        {
            if (j)
                printf(",");
            printf("%d", m[i][j]);
        }
    }
    flint_printf("]");
}

int
is_symplectic_j(si_mat_t m, long g, long g2)
{
    long i, j;
    for (i = 0; i < g; i++)
    {
        for (j = 0; j < g2; j++)
        {
            if (j!=2*i+1 && m[2*i][j])
                return 0;
            if (j!=2*i && m[2*i+1][j])
                return 0;
        }
        if (m[2*i][2*i+1]+m[2*i+1][2*i])
                return 0;
    }
    return 1;
}
