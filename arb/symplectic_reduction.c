/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

#define m(i,j) fmpz_mat_entry(m,i,j)

/*
 assume m is a 2g*2g antisymmetric matrix,
 find a basis such that m has standard form
 J(2g). Modifies m in place

 algorithm similar to snf:
 1. take minimal N[i,j] non-zero value above diagonal
 2. move it to front
 3. reduce lines i and j to have symplectic
    [0,1;-1,0] and zeroes outside
 the same operations are done on P to track the base change.
 */

static void
swap_step(si_mat_t p, si_mat_t m, slong i1, slong i2, slong len)
{
    if (i1 == i2)
        return;
    row_swap(p, i1, i2, len);
    row_swap(m, i1, i2, len);
    col_swap(m, i1, i2, len);
}

static void
neg_step(si_mat_t p, si_mat_t m, slong i, slong len)
{
    row_neg(p, i, len);
    row_neg(m, i, len);
    col_neg(m, i, len);
}
static void
transvect_step(si_mat_t p, si_mat_t m, slong i1, slong i2, fmpz_t q, slong len)
{
    row_addmul(p, i1, i2, q, len);
    row_addmul(m, i1, i2, q, len);
    col_addmul(m, i1, i2, q, len);
}

static void
bezout_step(si_mat_t p, si_mat_t m, slong i, slong k, fmpz_t a, fmpz_t b, slong len)
{
    fmpz_t u, v, g, u1, v1;
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(g);
    fmpz_init(u1);
    fmpz_init(v1);

    fmpz_xgcd(g, u, v, a, b);

    fmpz_divexact(u1, b, g);
    fmpz_neg(u1, u1);
    fmpz_divexact(v1, a, g);

    row_bezout(p, i, k, u, v, u1, v1, len);
    row_bezout(m, i, k, u, v, u1, v1, len);
    col_bezout(m, i, k, u, v, u1, v1, len);

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(g);
    fmpz_clear(u1);
    fmpz_clear(v1);
}

/* return 1 if 1 is on the line or smallest entry */
slong
pivot_line(si_mat_t m, slong i, slong len)
{
    slong j, jv = len;
    fmpz_t v;
    fmpz_init(v);
    for (j = i + 1; j < len; j++)
    {
        if (*m(i, j) == 0)
            continue;
        /* or return +/-1 if possible */
        if (fmpz_is_pm1(m(i, j)))
            return j;
        if (*v == 0 || fmpz_cmpabs(v, m(i, j)) > 0)
        {
            fmpz_set(v, m(i, j));
            jv = j;
        }
    }
    fmpz_clear(v);
    return jv;
}

void
symplectic_reduction(si_mat_t p, si_mat_t m, slong g, slong len)
{
    slong d, i, j, k;

    si_mat_set_id(p, len);

    /* main loop on symplectic 2-subspace */
    for (d = 0; d < g; d++)
    {
        int cleared = 0;
        i = 2 * d;
        /* lines 0..2d-1 already cleared */
        while ((j = pivot_line(m, i, len)) == len)
        {
            /* no intersection -> move ci to end */
            swap_step(p, m, i, len - 1, len);
            len--;
            if (len == 2*d)
                return;
        }

        /* move j to i + 1 */
        if (j != i+1)
            swap_step(p, m, j, i + 1, len);
        j = i + 1;

        /* make sure m(i, j) > 0 */
        if (fmpz_sgn(m(i, j)) < 0)
            /* could also neg i or swap i and j */
            neg_step(p, m, j, len);

        while(!cleared)
        {
            fmpz_t q;
            fmpz_init(q);
            /* clear row i */
            for (k = j + 1; k < len; k++)
            {
                if (*m(i, k))
                {
                    if (fmpz_divisible(m(i, k), m(i, j)))
                    {
                        fmpz_divexact(q, m(i, k), m(i, j));
                        fmpz_neg(q, q);
                        transvect_step(p, m, k, j, q, len);
                    }
                    else
                        bezout_step(p, m, j, k, m(i, j), m(i, k), len);
                }
            }
            cleared = 1;
            /* clear row j */
            for (k = j + 1; k < len; k++)
            {
                if (*m(j, k))
                {
                    if (fmpz_divisible(m(j, k), m(i, j)))
                    {
                        fmpz_divexact(q, m(j, k), m(i, j));
                        transvect_step(p, m, k, i, q, len);
                    }
                    else
                    {
                        bezout_step(p, m, i, k, m(j, i), m(j, k), len);
                        /* warning: row i now contains some ck.cl... */
                        cleared = 0;
                    }
                }
            }
            fmpz_clear(q);
        }
    }
}
