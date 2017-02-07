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

#define FORMULA 1

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

static slong
shift_ratio(const acb_t x, const acb_t y, slong m, slong prec)
{
    slong s;
    fmpz_t sz;
    acb_t r;
    arb_t a, pi;

    acb_init(r);
    arb_init(a);
    arb_init(pi);

    acb_div(r, x, y, prec);
    if (!acb_is_finite(r))
    {
        flint_printf("\nERROR: infinite value\n");
        arb_printd(x, 10); flint_printf("\n");
        arb_printd(y, 10); flint_printf("\n");
        abort();
    }

    arb_const_pi(pi, prec);

    /* if re < 0, rotate first */
    if (arb_is_nonpositive(acb_realref(r)))
    {
        acb_neg(r, r);
        acb_arg(a, r, prec);
        arb_add(a, a, pi, prec);
    }
    else
        acb_arg(a, r, prec);

    arb_mul_2exp_si(pi, pi, 1);
    arb_div(a, a, pi, prec);
    arb_mul_ui(a, a, m, prec);

    fmpz_init(sz);
    if (!arb_get_unique_fmpz(sz, a))
    {
        flint_printf("\nERROR: shift not an integer\n");
        arb_printd(a, 10);
        flint_printf("\nm = %ld, r = ", m);
        acb_printd(r, 10);
        abort();
    }
    s = fmpz_get_si(sz);
    fmpz_clear(sz);

    acb_clear(r);
    arb_clear(a);
    arb_clear(pi);
    return s;
}

/* (1 + sgn*(d-c)/(b-a)) * Cab * yab(x), x = +/-1, sgn = +/- 1 */
static void
limit_edge(acb_t z, acb_srcptr uab, slong nab, slong n, slong m, const acb_t cd2, int x, int sgn, slong prec)
{
    arb_t u;
    acb_t r;

    acb_init(r);
    arb_init(u);

    /* yab(x) * Cab */
    arb_set_si(u, x);
    mth_root_pol_def(z, uab, nab, n - 2, u, m, prec);
    acb_mul(z, z, uab + n, prec);

#if FORMULA==1
    acb_div(r, cd2, uab + n -2, prec);
    acb_add_si(r, r, sgn, prec);
    if (sgn < 0)
        acb_neg(r, r);
    acb_root_ui(r, r, m, prec);
    acb_mul(z, z, r, prec);
#endif

    acb_clear(r);
    arb_clear(u);
}

static slong
shift_abad(acb_srcptr uab, slong nab, acb_srcptr ucd, slong ncd, slong n, slong m)
{
    slong s, prec = 40;
    acb_t yab, yad;

    acb_init(yab);
    acb_init(yad);

    limit_edge(yab, uab, nab, n, m, ucd + n - 2,  -1, 1, prec);
    limit_edge(yad, ucd, ncd, n, m, uab + n - 2,  -1, 1, prec);
    s = shift_ratio(yad, yab, m, prec);

    acb_clear(yab);
    acb_clear(yad);
    return s;
}

static slong
shift_abbd(acb_srcptr uab, slong nab, acb_srcptr ucd, slong ncd, slong n, slong m)
{
    slong s, prec = 40;
    acb_t yab, ybd;

    acb_init(yab);
    acb_init(ybd);

    limit_edge(yab, uab, nab, n, m, ucd + n - 2,  1, -1, prec);
    limit_edge(ybd, ucd, ncd, n, m, uab + n - 2, -1, -1, prec);
    s = shift_ratio(ybd, yab, m, prec);

    acb_clear(yab);
    acb_clear(ybd);
    return s;
}

static int
is_neg_abad(const acb_t ab2, const acb_t cd2)
{
    int s;
    slong prec = 40;
    acb_t r;

    acb_init(r);
    acb_div(r, ab2, cd2, prec);
    s = arb_is_nonpositive(acb_imagref(r));
    acb_clear(r);
    return s;
}

static void
arb_angle(arb_t rho, const acb_t ab2, const acb_t cd2)
{
    acb_t t;
    slong prec = 40;
    acb_init(t);
    acb_div(t, ab2, cd2, prec);
    acb_arg(rho, t, prec);
    acb_clear(t);
}

void
intersection_tree(fmpz_mat_t c, const data_t data, const tree_t tree, slong n, slong m)
{
    slong k, l, size = m - 1;
    arb_t one, m_one;
    arb_init(one);
    arb_init(m_one);
    acb_ptr uab, ucd;
    slong nab, ncd;

    fmpz_mat_zero(c);

    /* the entry c[ k * (m-1) + s ] corresponds to the
       loop gamma_k^(s) */
    for (k = 0; k < n - 1; k++)
    {
        slong s;
        edge_t ek = tree->e[k];

        uab = data->upoints->rows[k];
        nab = data->n1[k];

        /* intersection with self shifts */
        fill_block(c, k * size, k * size, 1, -1, m);

        /* intersection with other shifts */
        for (l = k + 1; l < n - 1; l++)
        {
            arb_t rho;
            edge_t el = tree->e[l];

            if (ek.a != el.a && ek.b != el.a)
                /* no intersection */
                continue;

            ucd = data->upoints->rows[l];
            ncd = data->n1[l];
#if FORMULA==2
            arb_init(rho);
            arb_angle(rho, uab + nab - 2, acd + ncd - 2);
#endif

            if(el.a == ek.a)
            {
                /* case ab.ad */
#if FORMULA==3
                s = lrint((el.va - ek.va ) / TWOPI);
                if (el.dir > ek.dir)
#else
                s = shift_abad(uab, nab, ucd, ncd, n, m);
                if (is_neg_abad(uab + n -2, ucd + n -2))
#endif
                    fill_block(c, k * size, l * size, -s, -1-s, m);
                else
                    fill_block(c, k * size, l * size, 1-s, -s, m);
            }
            else if (el.a == ek.b)
            {
                /* case ab.bd */
#if FORMULA==3
                s = lrint(.5 + (el.va - ek.vb ) / TWOPI);
#else
                s = shift_abbd(uab, nab, ucd, ncd, n, m);
#endif
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
