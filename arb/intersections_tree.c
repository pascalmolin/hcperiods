/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

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

static void
arb_angle(arb_t rho, const acb_t ab, const acb_t cd, slong prec)
{
    acb_t t;
    arb_t pi;
    arb_init(pi);
    acb_init(t);
    acb_div(t, ab, cd, prec);
    if (!acb_is_finite(t))
    {
        flint_printf("\nERROR: infinite value\n");
        acb_printd(ab, 20); flint_printf("\n");
        acb_printd(cd, 20); flint_printf("\n");
        abort();
    }
    arb_const_pi(pi, prec);
    /* if re < 0, rotate first */
    if (arb_is_nonpositive(acb_realref(t)))
    {
        acb_neg(t, t);
        acb_arg(rho, t, prec);
        if (arb_is_nonnegative(rho))
            arb_sub(rho, rho, pi, prec);
        else
            arb_add(rho, rho, pi, prec);
    }
    else
        acb_arg(rho, t, prec);

    arb_div(rho, rho, pi, prec);
    arb_mul_2exp_si(rho, rho, -1);

    arb_clear(pi);
    acb_clear(t);
}

static slong
arb_get_si(const arb_t a, slong prec)
{
    slong s;
    fmpz_t sz;

    fmpz_init(sz);
    if (!arb_get_unique_fmpz(sz, a))
    {
        flint_printf("\nERROR: shift not an integer\n");
        arb_printd(a, 20);
        abort();
    }
    s = fmpz_get_si(sz);
    fmpz_clear(sz);
    return s;
}

/* Cab * yab(x), x = +/-1 */
static void
limit_edge(acb_t z, acb_srcptr uab, slong nab, slong n, slong m, int x, slong prec)
{
    arb_t u;
    acb_t r;

    acb_init(r);
    arb_init(u);

    /* yab(x) * Cab */
    arb_set_si(u, x);
    mth_root_pol_def(z, uab, nab, n - 2, u, m, prec);
    acb_mul(z, z, uab + n, prec);

    acb_clear(r);
    arb_clear(u);
}

enum { ABAD, ABBD };

static slong
shift_number(acb_srcptr uab, slong nab, acb_srcptr ucd, slong ncd, slong n, slong m, int type, slong prec)
{
    slong s;
    arb_t rho, psi;
    acb_t yab, ycd;

    acb_init(yab);
    acb_init(ycd);

    if (type == ABAD)
        limit_edge(yab, uab, nab, n, m, -1, prec);
    else
        limit_edge(yab, uab, nab, n, m, 1, prec);
    limit_edge(ycd, ucd, ncd, n, m, -1, prec);

    
    //flint_printf("\nyab");
    //acb_printd(yab, 20);
    //flint_printf("\nycd");
    //acb_printd(ycd, 20);

    arb_init(rho);
    arb_init(psi);
    arb_angle(rho, uab + n -2, ucd + n - 2, prec);
    arb_angle(psi, ycd, yab, prec);
    arb_mul_ui(psi, psi, m, prec);
    arb_add(rho, rho, psi, prec);
    if (type == ABBD)
    {
        //flint_printf("\nabbd, ba = ");
        //acb_printd(uab+n-2,20);
        //flint_printf("\n      dc = ");
        //acb_printd(ucd+n-2,20);
        //flint_printf("\n      rho = ");
        //arb_printd(rho,20);
        arb_one(psi);
        arb_mul_2exp_si(psi, psi, -1);
        arb_add(rho, rho, psi, prec);
    }
    s = arb_get_si(rho, prec);

    acb_clear(yab);
    acb_clear(ycd);
    return s;
}

static int
is_neg_abad(const acb_t ab2, const acb_t cd2, slong prec)
{
    int s;
    acb_t r;

    acb_init(r);
    acb_div(r, ab2, cd2, prec);
    s = arb_is_nonpositive(acb_imagref(r));
    acb_clear(r);
    return s;
}

void
intersection_tree(fmpz_mat_t c, const data_t data, const tree_t tree, slong n, slong m)
{
    slong prec = 64;
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

            arb_init(rho);
            arb_angle(rho, uab + n - 2, ucd + n - 2, prec);

            if(el.a == ek.a)
            {
                /* case ab.ad */
                s = shift_number(uab, nab, ucd, ncd, n, m, ABAD, prec);
                if (is_neg_abad(uab + n -2, ucd + n -2, prec))
                    fill_block(c, k * size, l * size, -s, -1-s, m);
                else
                    fill_block(c, k * size, l * size, 1-s, -s, m);
            }
            else if (el.a == ek.b)
            {
                /* case ab.bd */
                s = shift_number(uab, nab, ucd, ncd, n, m, ABBD, prec);
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
