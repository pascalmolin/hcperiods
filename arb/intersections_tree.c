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

enum { ABAD, ABBD };

static slong
shift_number(ydata_t yab, ydata_t ycd, slong m, int type, slong prec)
{
    slong s;
    arb_t rho, psi;

    arb_init(rho);
    arb_init(psi);
    arb_angle(rho, yab->ba2, ycd->ba2, prec);

    if (type == ABAD)
        arb_angle(psi, ycd->ya, yab->ya, prec);
    else
        arb_angle(psi, ycd->ya, yab->yb, prec);
    
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
intersection_tree(fmpz_mat_t c, const tree_t tree, slong m)
{
    slong prec = 64;
    slong k, l, size = m - 1;

    fmpz_mat_zero(c);

    /* the entry c[ k * (m-1) + s ] corresponds to the
       loop gamma_k^(s) */
    for (k = 0; k < tree->n; k++)
    {
        slong s;
        edge_t ek = tree->e[k];

        /* intersection with self shifts */
        fill_block(c, k * size, k * size, 1, -1, m);

        /* intersection with other shifts */
        for (l = k + 1; l < tree->n; l++)
        {
            edge_t el = tree->e[l];

            if (ek.a != el.a && ek.b != el.a)
                /* no intersection */
                continue;

            if(el.a == ek.a)
            {
                /* case ab.ad */
                s = shift_number(tree->data + k, tree->data + l, m, ABAD, prec);
                /* FIXME: use (sp, sp-1) */
                if (is_neg_abad(tree->data[k].ba2, tree->data[l].ba2, prec))
                    fill_block(c, k * size, l * size, -s, -1-s, m);
                else
                    fill_block(c, k * size, l * size, 1-s, -s, m);
            }
            else if (el.a == ek.b)
            {
                /* case ab.bd */
                s = shift_number(tree->data + k, tree->data + l, m, ABBD, prec);
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
