/******************************************************************************

 Copyright (C) 2016 Pascal Molin, Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"
#include <complex.h>

#define PI2 acos(0.)
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
endvalues_edge(double * va, double * vb, const cdouble * w, slong ia, slong ib, slong d)
{
    slong k;
    double a, b;
    cdouble ba = w[ib] - w[ia];

    a = b = (d - 1) * carg(ba);
    
    for (k = 0; k < d; k++)
    {
        if (k == ia || k == ib)
            continue;

        a += carg((w[ia] - w[k]) / ba);
        b += carg((w[ib] - w[k]) / ba);
    }
    *va = a;
    *vb = b;
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
    slong k, i, n;
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

    for (k = 0; k < tree->n; k++)
    {

        /* start from last edge, discard if both left or taken */
        /* remark: stupid to loop many times, otherwise one can
         * take non-connected edges but one must reorder the
         * edges at the end */
        for (i = n - 1; k && t[e[i].a] == t[e[i].b]; i--);

        if (t[e[i].b])
            /* reorder edge */
            e[i] = edge_flip(e[i]);

        t[e[i].a] = 1;
        t[e[i].b] = 1;

        /* compute endvalues for shifting numbers */
        endvalues_edge(&e[i].va, &e[i].vb, w, e[i].a, e[i].b, len);

        tree->e[k] = e[i];

        /* save complexity estimate */
        if (e[i].tau < tree->tau)
            tree->tau = e[i].tau;
    }

    free(w);
    free(e);
    free(t);
}

void
tree_init(tree_t tree, slong n)
{
    tree->n = n;
    tree->e = malloc(n * sizeof(edge_t));
    tree->tau = PI2;
}

void
tree_clear(tree_t tree)
{
    free(tree->e);
}
