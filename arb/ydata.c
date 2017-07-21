/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

/* cab on component n */
static void
constant_cab(acb_t c, const acb_t ba2, slong n1, slong n, slong m, slong prec)
{
    slong k;

    acb_pow_ui(c, ba2, n, prec);
    if ((k = arb_is_negative(acb_realref(c))))
        acb_neg(c, c);

    acb_root_ui(c, c, m, prec);

    if ((k + n1) % 2 == 0)
    {
        acb_t z;
        acb_init(z);
        acb_unit_root(z, 2 * m, prec);
        acb_mul(c, c, z, prec);
        acb_clear(z);
    }

    if (acb_contains_zero(c))
    {
        flint_printf("\n\nERROR: Cab contains 0\nm, n = %ld, %ld\n\n", m, n);
        abort();
    }
}

/* u[0..n1[ contains roots re(ui)>0
   u[n1..n-2[ roots with re(ui) <= 0
   the last two components are set to
   (b-a)/2 and (a+b)/(b-a)
   returns n1
*/
static slong
ab_points(acb_ptr u, acb_srcptr x, edge_t e, slong n, slong m, slong prec)
{
    slong k, l;
    acb_ptr ab, ba;
    ba = u + n - 2; /* b - a */
    ab = u + n - 1; /* a + b */

    acb_set(ba, x + e.b);
    acb_add(ab, ba, x + e.a, prec);
    acb_sub(ba, ba, x + e.a, prec);

    for (k = 0, l = 0; k < n; k++)
    {
        if (k == e.a || k == e.b)
            continue;
        acb_mul_2exp_si(u + l, x + k, 1);
        acb_sub(u + l, u + l, ab, prec);
        acb_div(u + l, u + l, ba, prec);
        l++;
    }
    if (l != n - 2)
        abort();

    /* now l = n - 2, reorder */
    for (k = 0; k < l; k++)
        if (arb_is_nonpositive(acb_realref(u + k)))
            acb_swap(u + k--, u + --l);

    /* set last two constants */
    acb_div(ab, ab, ba, prec); /* (a+b)/(b-a) */
    acb_mul_2exp_si(ba, ba, -1);   /* (b-a)/2 */

    return l;
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
    mth_root_pol_def(z, uab, nab, n - 2, u, NULL, m, prec);
    acb_mul(z, z, uab + n, prec);

    acb_clear(r);
    arb_clear(u);
}

void
ydata_init_edge(ydata_t yab, acb_srcptr x, edge_t e, slong n, slong m, slong prec)
{
    slong l;
    acb_ptr u;

    u = _acb_vec_init(n - 2 + 5);
    l = ab_points(u, x, e, n, m, prec);

    yab->m = m;
    yab->n = n;
    yab->u = u;
    yab->n1 = l;

    yab->ba2 = u + n - 2;
    yab->ab = u + n - 1;
    yab->c = u + n;
    yab->ya =  u + n + 1;
    yab->yb =  u + n + 2;

    constant_cab(yab->c, yab->ba2, l, n, m, prec);

    /* limits at a and b */
    limit_edge(yab->ya, u, l, n, m, -1, prec);
    limit_edge(yab->yb, u, l, n, m, 1, prec);
}

void
ydata_clear(ydata_t yab)
{
    _acb_vec_clear(yab->u, yab->n - 2 + 5);
}

void
tree_ydata_init(tree_t tree, acb_srcptr x, slong n, slong m, slong prec)
{
    slong k;
    tree->data = flint_malloc(tree->n * sizeof(ydata_t));
    for (k = 0; k < tree->n; k++)
        ydata_init_edge(tree->data + k, x, tree->e[k], n, m, prec);
}

void
tree_ydata_clear(tree_t tree)
{
    slong k;
    for (k = 0; k < tree->n; k++)
        ydata_clear(tree->data + k);
    flint_free(tree->data);
}

/* maximum size of ((b-a)/2)^(-jn/m+i) binom(i-1,l) ((b+a)/(b-a))^(i-1-l) */
void
mag_constants_edge(mag_t m, acb_srcptr x, edge_t e, sec_t c)
{
    long prec = 64, i, j, l;
    acb_t t;
    arb_t mab, mba, mc, mi;
    mag_t ml;

    arb_init(mab);
    arb_init(mba);
    arb_init(mc);
    arb_init(mi);
    mag_init(ml);

    mag_zero(m);

    acb_init(t);
    acb_add(t, x + e.b, x + e.a, prec);
    acb_abs(mab, t, prec);
    acb_sub(t, x + e.b, x + e.a, prec);
    acb_abs(mba, t, prec);
    arb_div(mab, mab, mba, prec);  /* (b+a)/(b-a) */
    arb_mul_2exp_si(mba, mba, -1); /* (b-a)/2 */
    arb_pow_ui(mc, mba, c.n, prec);
    arb_root(mc, mc, c.m, prec);
    arb_inv(mc, mc, prec); /* (2/b-a)^(n/m) */
    acb_clear(t);

    for (j = 0; j < c.nj; j++)
    {
        arb_pow_ui(mi, mc, j + c.j1, prec);
        for (i = 1; i < 1 + c.ni[j]; i++)
        {
            arb_mul(mi, mi, mba, prec);
            for (l = 0; l <= i - 1; l++)
            {
                arb_t bin, bil;
                arb_init(bin);
                arb_init(bil);
                arb_bin_uiui(bin, i - 1, l, prec);
                arb_pow_ui(bil, mab, i - 1 - l, prec);
                arb_mul(bin, bin, bil, prec);
                arb_mul(bin, bin, mi, prec);
                arb_get_mag(ml, bin);
                mag_max(m, m, ml);
                arb_clear(bin);
                arb_clear(bil);
            }
        }
    }
    arb_clear(mab);
    arb_clear(mba);
    arb_clear(mc);
    arb_clear(mi);
    mag_clear(ml);
}

void
mag_constants_tree(mag_t mag, tree_t tree, acb_srcptr x, sec_t c)
{
    mag_t me;
    slong k;
    mag_init(me);
    mag_zero(mag);
    for (k = 0; k < tree->n; k++)
    {
        mag_constants_edge(me, x, tree->e[k], c);
        mag_max(mag, mag, me);
    }
    mag_clear(me);
}

slong
extraprec_tree(tree_t tree, acb_srcptr x, sec_t c)
{
    slong prec;
    mag_t mag;
    mag_init(mag);
    mag_constants_tree(mag, tree, x, c);
    prec = (long)mag_get_d_log2_approx(mag) + 1;
    mag_clear(mag);
    return prec;
}
