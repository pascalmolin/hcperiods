/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

#define TWOPI (2 * acos(-1.))
#define c(i,j) fmpz_mat_entry(c,i,j)

/* set
   c[i+k][j+l] = 1 if l-k = sp mod m
   c[i+k][j+l] = -1 if l-k = sm mod m
 */

static void
fill_block(fmpz_mat_t c, slong i, slong j, slong sp, slong sm, slong m)
{
    slong k, l;
    /* important: make sp and sm positive */
    sp = (sp % m + m) % m;
    sm = (sm % m + m) % m;
    for (k = 0; k < m - 1; k++)
    {
        l = (k + sp) % m;
        if (l < m - 1)
        {
            *c(i + k, j + l) = 1;
            *c(j + l, i + k) = -1;
        }
        l = (k + sm) % m;
        if (l < m - 1)
        {
            *c(i + k, j + l) = -1;
            *c(j + l, i + k) = 1;
        }
    }
}

slong
shift_abbd(const acb_ptr uab, slong nab, const acb_ptr ucd, slong ncd, slong n, slong m)
{
    slong prec = 40;
    arb_t one, m_one;
    acb_t r, yab, ycd;
    arb_t a, pi;
    arb_init(one);
    arb_init(m_one);
    acb_init(r);
    acb_init(yab);
    acb_init(ycd);
    slong s;
    fmpz_t sz;

    arb_one(one);
    arb_neg(m_one, one);

    mth_root_pol_def(ycd, ucd, ncd, n, m_one, m, prec);
    acb_mul(ycd, ycd, ucd + n, prec);
    acb_div(r, uab + n - 2, ucd + n -2, prec);
    acb_sub_arb(r, r, one, prec);
    acb_neg(r, r);
    acb_root_ui(r, r, m, prec);
    acb_mul(ycd, ycd, r, prec);

    mth_root_pol_def(yab, uab, nab, n, one, m, prec);
    acb_mul(yab, yab, uab + n, prec);
    acb_div(r, ucd + n - 2, uab + n -2, prec);
    acb_sub_arb(r, r, one, prec);
    acb_neg(r, r);
    acb_root_ui(r, r, m, prec);
    acb_mul(yab, yab, r, prec);

    acb_div(r, ycd, yab, prec);

    arb_init(a);
    arb_init(pi);
    acb_arg(a, r, prec);

    arb_const_pi(pi, prec);
    arb_mul_2exp_si(pi, pi, 1);
    arb_div(a, a, pi, prec);
    arb_mul_ui(a, a, m, prec);

    fmpz_init(sz);
    if (!arb_get_unique_fmpz(sz, a))
        abort();
    s = fmpz_get_si(sz);
    fmpz_clear(sz);

    arb_clear(a);
    arb_clear(pi);

    arb_clear(one);
    arb_clear(m_one);
    acb_clear(r);
    acb_clear(yab);
    acb_clear(ycd);

    return s;
}

void
intersection_tree(fmpz_mat_t c, const tree_t tree, slong d, slong m)
{
    slong k, l, size = m - 1;
    arb_t one, m_one;
    arb_init(one);
    arb_init(m_one);

    fmpz_mat_zero(c);

    /* the entry c[ k * (m-1) + s ] corresponds to the
       loop gamma_k^(s) */
    for (k = 0; k < d - 1; k++)
    {
        slong s;
        edge_t ek = tree->e[k];

        /* intersection with self shifts */
        fill_block(c, k * size, k * size, 1, -1, m);

        /* intersection with other shifts */
        for (l = k + 1; l < d - 1; l++)
        {
            edge_t el = tree->e[l];

            if (ek.a != el.a && ek.b != el.a)
                /* no intersection */
                continue;

            if(el.a == ek.a)
            {
                /* case ab.ad */
                s = lrint((el.va - ek.va ) / TWOPI);
                if (el.dir > ek.dir)
                    fill_block(c, k * size, l * size, -s, -1-s, m);
                else
                    fill_block(c, k * size, l * size, 1-s, -s, m);
            }
            else if (el.a == ek.b)
            {
                /* case ab.bd */
                s = lrint(.5 + (el.va - ek.vb ) / TWOPI);
                fill_block(c, k * size, l * size, -s, 1-s, m);
            }
            else
            {
                flint_printf("invalid tree\n");
                abort();
            }
        }
    }
}
