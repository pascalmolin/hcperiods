/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "flint.h"
#include "si_mat.h"

#define m(i,j) fmpz_mat_entry(m,i,j)

void
si_mat_init(si_mat_t m, long nr, long nc)
{
    fmpz_mat_init(m, nr, nc);
}

void
si_mat_clear(si_mat_t m, long nr, long nc)
{
    fmpz_mat_clear(m);
}

void
row_swap(si_mat_t m, long i1, long i2, long len)
{
    fmpz_mat_swap_rows(m, NULL, i1, i2);
}

void
col_swap(si_mat_t m, long j1, long j2, long len)
{
    long i;
    for (i = 0; i < m->r; i++)
        fmpz_swap(m(i, j1), m(i, j2));
}

void
col_neg(si_mat_t m, long j, long len)
{
    long i;
    for (i = 0; i < m->r; i++)
        fmpz_neg(m(i, j), m(i, j));
}

void
row_neg(si_mat_t m, long i, long len)
{
    long j;
    for (j = 0; j < m->c; j++)
        fmpz_neg(m(i, j), m(i, j));
}

void
col_addmul(si_mat_t m, long j, long k, fmpz_t v, long len)
{
    long i;
    for (i = 0; i < m->r; i++)
        fmpz_addmul(m(i, j), v, m(i, k));
}

void
row_addmul(si_mat_t m, long i, long k, fmpz_t v, long len)
{
    long j;
    for (j = 0; j < m->c; j++)
        fmpz_addmul(m(i, j), v, m(k, j));
}

/* cj <- a cj + b ck
   ck <- c cj + d ck */
void
row_bezout(si_mat_t m, long i1, long i2, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, long len)
{
    long j;
    fmpz_t s, t;
    fmpz_init(s);
    fmpz_init(t);
    for (j = 0; j < m->c; j++)
    {
        fmpz_set(s, m(i1, j));
        fmpz_set(t, m(i2, j));
        fmpz_mul(m(i1, j), a, s);
        fmpz_addmul(m(i1, j), b, t);
        fmpz_mul(m(i2, j), c, s);
        fmpz_addmul(m(i2, j), d, t);
    }
    fmpz_clear(s);
    fmpz_clear(t);
}
void
col_bezout(si_mat_t m, long j1, long j2, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, long len)
{
    long i;
    fmpz_t s, t;
    fmpz_init(s);
    fmpz_init(t);
    for (i = 0; i < m->r; i++)
    {
        fmpz_set(s, m(i, j1));
        fmpz_set(t, m(i, j2));
        fmpz_mul(m(i, j1), a, s);
        fmpz_addmul(m(i, j1), b, t);
        fmpz_mul(m(i, j2), c, s);
        fmpz_addmul(m(i, j2), d, t);
    }
    fmpz_clear(s);
    fmpz_clear(t);
}

void
si_mat_set_id(si_mat_t p, long len)
{
    fmpz_mat_one(p);
}

int
si_mat_eq(si_mat_t a, si_mat_t b, long nr, long nc)
{
    fmpz_mat_equal(a, b);
}

void
si_mat_print(si_mat_t m, long nr, long nc)
{
    fmpz_mat_print(m);
}
    
#if 0
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
            printf("%d", m(i, j));
        }
    }
    flint_printf("]");
}
#endif

int
is_symplectic_j(si_mat_t m, long g, long g2)
{
    long i, j;
    for (i = 0; i < g; i++)
    {
        for (j = 0; j < g2; j++)
        {
            if (j!=2*i+1 && *m(2*i, j))
                return 0;
            if (j!=2*i && *m(2*i+1, j))
                return 0;
        }
        if (*m(2*i, 2*i+1) + *m(2*i+1, 2*i))
                return 0;
    }
    return 1;
}
