/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
#include <complex.h>
#define PI2 asin(0.)
#define LAMBDA PI2

typedef complex double cdouble;

static cdouble
acb_get_cd(const acb_t z)
{
    return
        arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR)
    + _Complex_I * arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
}

static double
tau_3(cdouble a, cdouble b, cdouble c)
{
    cdouble z, xItau;
    z = (2*c-a-b)/(b-a);
    xItau = casinh(catanh(z)/LAMBDA);
    return fabs(cimag(xItau));
}

static double
tau_edge(const cdouble * w, slong i, slong j, slong len, slong * l)
{
    slong k;
    double tmp, tau = PI2;
    for (k = 0; k < len; k++)
    {
        if (k == i || k == j)
            continue;
        tmp = tau_3(w[i], w[j], w[k]);
        if (tmp < tau)
        {
            tau = tmp;
            if (l)
                *l = k;
        }
    }
    return tau;
}

static void
endvalues_edge(double * va, double * vb, const cdouble * w, slong ia, slong ib, slong len)
{
    slong k;
    cdouble fa, fb;

    fa = fb = w[ib] - w[ia];
    
    for (k = 0; k < len; k++)
    {
        if (k == ia || k == ib)
            continue;
        fa *= (w[ia] - w[k]);
        fb *= (w[ia] - w[k]);
    }
    *va = carg(fa);
    *vb = carg(fb);
}

static void
edges_init(edge_t * e, const cdouble * w, slong len)
{
    slong i, j, k;

    k = 0;
    for (i = 0; i < len; i++)
    {
        for (j = i + 1; j < len; j++)
        {
            e[k].a = i;
            e[k].b = j;
            e[k].tau = tau_edge(w, i, j, len, NULL);
            k++;
        }
    }
}

static edge_t
edge_flip(edge_t e)
{
    edge_t f;
    f.tau =  e.tau;
    f.a = e.b;
    f.b = e.a ;
    f.va = e.vb;
    f.vb = e.va ;
    return f;
}

int
edge_cmp(const edge_t * x, const edge_t * y)
{
        return (x->tau < y->tau) ? -1 : (x->tau > y->tau);
}

void
spanning_tree(tree_t tree, acb_srcptr x, slong len)
{
    slong k, n;
    cdouble * w;
    int * t;
    edge_t * e;

    /* small approx of roots */
    w = malloc(len * sizeof(cdouble));
    for (k = 0; k < len; k++)
        w[k] = acb_get_cd(x + k);

    n = (len * (len - 1)) / 2;
    e = malloc(n * sizeof(edge_t));
    edges_init(e, w, len);

    /* order edges */
    qsort(e, n, sizeof(edge_t), (int(*)(const void*,const void*))edge_cmp);

    t = malloc(len * sizeof(int));
    for (k = 0; k < len; k++)
        t[k] = 0;

    n--;
    for (k = 0; k < tree->n; k++)
    {

        /* discard if both left or taken */
        for (; t[e[n].a] == t[e[n].b]; n--);

        if (t[e[n].b])
            /* reorder edge */
            e[n] = edge_flip(e[n]);

        t[e[n].a] = 1;
        t[e[n].b] = 1;

        /* compute endvalues for shifting numbers */
        endvalues_edge(&e[n].va, &e[n].vb, w, e[n].a, e[n].b, len);

        tree->e[k] = e[n];
    }

    /* save complexity estimate */
    tree->tau = e[n].tau;

    free(w);
    free(e);
    free(t);
}

void
tree_init(tree_t tree, slong n)
{
    tree->n = n;
    tree->e = malloc(n * sizeof(edge_t));
}

void
tree_clear(tree_t tree)
{
    free(tree->e);
}
