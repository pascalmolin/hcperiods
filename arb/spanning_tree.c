/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"
#include "complex_extras.h"

#define edge_flip(e) do { \
    slong a = e.a; e.a = e.b; e.b = a; \
} while (0);

typedef double (* param_f)(cdouble a, cdouble b, cdouble c);

/* superelliptic case, ellipse */

static double
param_gc_r(cdouble a, cdouble b, cdouble c)
{
    return (cabs(c-a)+cabs(c-b))/cabs(b-a);
}

/* general case, tanhsinh estimate */

static double
param_de_tau(cdouble a, cdouble b, cdouble c)
{
    cdouble z, xItau;
    z = (2*c-a-b)/(b-a);
    xItau = casinh(catanh(z)/LAMBDA);
    return fabs(cimag(xItau));
}

static double
param_edge(param_f param, const cdouble * w, slong i, slong j, slong len)
{
    slong k;
    double tmp, p = 0;
    for (k = 0; k < len; k++)
    {
        if (k == i || k == j)
            continue;
        tmp = param(w[i], w[j], w[k]);
        if (!p || tmp < p)
            p = tmp;
    }
    return p;
}

static void
edges_init(edge_t * e, param_f param, const cdouble * w, slong len)
{
    slong i, j, k;

    k = 0;
    for (i = 0; i < len; i++)
    {
        for (j = i + 1; j < len; j++)
        {
            e[k].a = i;
            e[k].b = j;
            e[k].r = param_edge(param, w, i, j, len);
            k++;
        }
    }
}

int
edge_cmp(const edge_t * x, const edge_t * y)
{
        return (x->r < y->r) ? -1 : (x->r > y->r);
}

static void
endvalues_edge(double * va, double * vb, double * dir,
        const cdouble * w, slong ia, slong ib, slong d)
{
    slong k;
    double a, b;
    cdouble ba = w[ib] - w[ia];

    *dir = carg(ba);
    a = b = (d - 1) * (*dir);
    
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

void
spanning_tree(tree_t tree, acb_srcptr x, slong len, int type)
{
    slong k, i, n;
    cdouble * w;
    int * t;
    edge_t * e;

    /* small approx of roots */
    w = malloc(len * sizeof(cdouble));
    for (k = 0; k < len; k++)
        w[k] = acb_get_cdouble(x + k);

    n = (len * (len - 1)) / 2;
    e = malloc(n * sizeof(edge_t));

    if (type == INT_GC)
        edges_init(e, param_gc_r, w, len);
    else if (type == INT_DE)
        edges_init(e, param_de_tau, w, len);
    else
        printf("unknown type\n"), abort();

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

        /* ensure a already in tree */
        if (t[e[i].b])
            edge_flip(e[i]);

        t[e[i].a] = 1;
        t[e[i].b] = 1;

        /* compute endvalues for shifting numbers */
        endvalues_edge(&e[i].va, &e[i].vb, &e[i].dir, w, e[i].a, e[i].b, len);

        tree->e[k] = e[i];

        /* save complexity estimate */
        if (!tree->r || e[i].r < tree->r)
        {
            tree->r = e[i].r;
            tree->min = k;
        }
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
    tree->r = 0.;
    tree->min = 0;
}

void
tree_clear(tree_t tree)
{
    free(tree->e);
}
