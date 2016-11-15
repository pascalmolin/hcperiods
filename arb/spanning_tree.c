/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "abel_jacobi.h"
#include <complex.h>
#define LAMBDA M_PI_2

typedef complex double cdouble;

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
    double tmp, tau = M_PI_2;
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

static cdouble
acb_get_cd(const acb_t z)
{
    return
        arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR)
    + _Complex_I * arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
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

int
edge_cmp(const edge_t * x, const edge_t * y)
{
        return (x->tau < y->tau) ? -1 : (x->tau > y->tau);
}

void spanning_tree(tree_t tree, acb_srcptr x, slong len)
{
    slong k, n;
    cdouble * w;
    int * t;
    edge_t * e;

    w = malloc(len * sizeof(cdouble));
    for (k = 0; k < len; k++)
        w[k] = acb_get_cd(x + k);

    n = (len * (len - 1)) / 2;
    e = malloc(n * sizeof(edge_t));
    edges_init(e, w, len);
    free(w);

    qsort(e, n, sizeof(edge_t), (int(*)(const void*,const void*))edge_cmp);
    
    t = malloc(len * sizeof(int));
    for (k = 0; k < len; k++)
        t[k] = 0;

    n--;
    for (k = 0; k < tree->n; k++)
    {
        for (; t[e[n].a] && t[e[n].b]; n--);
        tree->e[k] = e[n];
        t[e[n].a] = 1;
        t[e[n].b] = 1;
    }
    tree->tau = e[n].tau;

    free(e);
    free(t);
}

void tree_init(tree_t tree, slong n)
{
    tree->n = n;
    tree->e = malloc(n * sizeof(edge_t));
}

void tree_clear(tree_t tree)
{
    free(tree->e);
}
