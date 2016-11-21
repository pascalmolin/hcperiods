/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

static void
col_swap(si_mat_t m, slong j1, slong j2, slong len)
{
    slong * c = m[j1];
    m[j1] = m[j2];
    m[j2] = c;
}
static void
row_swap(si_mat_t m, slong i1, slong i2, slong len)
{
    slong j, t;
    for (j = 0; j < len; j++)
    {
        t = m[j][i1];
        m[j][i1] = m[j][i2];
        m[j][i2] = t;
    }
}
static void
col_neg(si_mat_t m, slong j, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        m[j][i] = -m[j][i];
}
static void
row_neg(si_mat_t m, slong i, slong len)
{
    slong j;
    for (j = 0; j < len; j++)
        m[j][i] = -m[j][i];
}
static void
col_muladdmul(si_mat_t m, slong j, slong u, slong k, slong v, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        m[j][i] = u * m[j][i] + v * m[k][i];
}
static void
row_muladdmul(si_mat_t m, slong i, slong u, slong k, slong v, slong len)
{
    slong j;
    for (j = 0; j < len; j++)
        m[j][i] = u * m[j][i] + v * m[j][k];
}
/* assume m is a 2g*2g antisymmetric matrix,
 * find a basis such that m has standard form
 * J(2g). Modifies m */
void
symplectic_reduction(si_mat_t p, si_mat_t m, slong g, slong len)
{
    slong d, i, j, k;

    for (i = 0; i < len; i++)
        for (j = 0; j < len; j++)
            p[j][i] = (i==j);

    /* loop on symplectic 2-subspace */
    for (d = 0; d < g; d++)
    {
        slong ii, jj;
        i = 2*d; j = i+1;
        /* find pivot value i,j -> subspace c_i, c_j */
        /* TODO: select min absolute value or first +/- 1 */
        for (j = i + 1; !m[j][i]; j++);
        /* swap i,j <- 2d, 2d+1 */
        ii = 2*d; jj = 2*d+1;
        if (i > ii)
        {
            col_swap(p, i, ii, len);
            col_swap(m, i, ii, len);
            row_swap(m, i, ii, len);
            i = ii;
        }
        if (j > jj)
        {
            col_swap(p, j, jj, len);
            col_swap(m, j, jj, len);
            row_swap(m, j, jj, len);
            j = jj;
        }
        if (m[j][i] < 0)
        {
            /* c_j <- - c_j */
            col_neg(p, j, len);
            col_neg(m, j, len);
            row_neg(m, j, len);
        }
        /* clear row i */
        for (k = j + 1; k < len; k++)
        {
            ulong u = 1, v;
            if (!m[i][k])
                continue;
            if (m[i][k]%m[i][j] == 0)
                v = m[i][k] / m[i][j];
            /* ck <- ck - v * cj */
            else
                n_xgcd(&u, &v, m[i][j], m[i][k]);

            col_muladdmul(p, k, u, j, -v, len);
            col_muladdmul(m, k, u, j, -v, len);
            row_muladdmul(m, k, u, j, -v, len);
        }
    }
}

