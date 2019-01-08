/* Copyright (C) 2000  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. */
#include <pari/pari.h>
#include <pari/paripriv.h>

/*********************************************************************/
/**                                                                 **/
/**                PERIODS OF HYPERELLIPTIC CURVES                  **/
/**               contributed by Pascal Molin (2019)                **/
/**                                                                 **/
/*********************************************************************/

/*
 The periods are the integrals int_{c_i} w_j where
 - w_j = x^j/y dx is a basis of holomorphic differentials
 - c_i = (a_i,b_i) is a symplectic basis of loops

 the basis is expressed on tree integrals
*/


#define hc_get_pol(c)        gel(c,1)
#define hc_get_roots(c)      gel(c,2)
#define hc_get_periods(c)    gel(c,3)
#define hc_get_tree(c)       gmael(c,4,1)
#define hc_get_integrals(c)  gmael(c,4,2)
#define hc_get_hombasis(c)   gmael(c,4,3)
#define hc_get_n(c) poldegree(hc_get_pol(c))

/** everything is dynamic **/

/*********************************************************************/
/*                                                                   */
/*             Local square roots and integrals                      */
/*                                                                   */
/*********************************************************************/

/* assume -1<= x <= 1,
 * u a t_VEC of complex numbers with increasing real(u_i) */
GEN
sqrt_pol_def(GEN u, GEN x, long prec)
{
    long k, r=0;
    GEN y = gen_1;
    pari_sp av = avma;
    for (k = 1; k < lg(u); k++) 
    {
        GEN uk = gel(u,k);
        uk = signe(gel(uk,1)) > 0 ? r++, gsub(uk, x) : gsub(x,uk);
        y = gmul(y,gsqrt(uk,prec));
    }
    switch (r % 4)
    {
        case 1: y = gmul(y, gen_I()); break;
        case 2: y = gneg(y); break;
        case 3: y = gneg(gmul(y,gen_I())); break;
        default: break;
    }
    return gerepilecopy(av, y);
}

/* e = Vecsmall(i,j), a = roots[i], b = roots[j],
 * return [ [ (2x-a-b)/(b-a)...], (b-a)/2, (a+b)/(b-a), cab, fa, fb ]
 * */
GEN
ydata_init(GEN roots, GEN e, long prec)
{
    long ia = e[1], ib = e[2], k, l, d = lg(roots)-1;
    GEN a, b, ab, ba, u, cab, fa, fb;
    pari_sp av = avma;
    a = gel(roots, ia);
    b = gel(roots, ib);
    ab = gmul2n(gadd(a,b),-1);
    ba = gmul2n(gsub(b,a),-1);

    u = cgetg(lg(roots)-2, t_VEC);

    for (k = 1, l = 1; k < lg(roots); k++)
    {
        if (k == ia || k == ib)
            continue;
        gel(u, l) = gdiv(gsub(gel(roots, k), ab), ba);
        l++;
    }

    /* reorder real < 0 first */
    //for (k = 1; k < l; k++)
    //   if (signe(greal(u + k))<0)
    //       swap(gel(u, k--),gel(u, --l));
    u = lexsort(u);

    cab = gpowgs(gsqrt(ba,prec),d);
    fa = gmul(cab, sqrt_pol_def(u,gen_m1,prec));
    fb = gmul(cab, sqrt_pol_def(u,gen_1,prec));
    cab = ginv(cab);
    
    return gerepilecopy(av, mkvecn(6, u, ba, gdiv(ab,ba), cab, fa, fb));
}

GEN
ydata_tree(GEN X, GEN tree, long prec)
{
    long k;
    GEN data = cgetg(lg(tree),t_VEC);
    for (k = 1; k < lg(tree); k++)
        gel(data, k) = mkvec2(gel(tree,k),ydata_init(X, gmael(tree,k,2), prec));
    return data;
}

/* [ cos((2*k+1)Pi/(4n)), 0 <= k < n ] */
GEN
gc_init(long n, long prec)
{
    long k;
    GEN z, x;
    pari_sp av = avma;
    /* force n even? no need */
    //n = (n+1) / 2;
    z = grootsof1(8*n, prec); // need only 1/4
    x = cgetg(n + 1, t_VEC);
    //return gerepilecopy(av, vecslice(z, 1, 2 * n));
    for (k = 1; k <= n; k++)
    {
        gel(x, k) = greal(gel(z, 2*k));
        //gel(x, n+1-k) = gimag(gel(z, 2*k+1));
    }
    return gerepilecopy(av, x);
}

long
gc_cost(GEN r, long bitprec)
{
    long pr2 = ndec2prec(18);
    double dr;
    r = gacosh(r, pr2);
    //output(r);
    dr = rtodbl(r);
    return ceil((bitprec*M_LN2+log(2*M_PI/dr)+1) / (4*dr));
}

GEN
integral_edge(GEN ydata, long g, GEN gc, long prec)
{
    long l, n;
    GEN res, u = gel(ydata, 1), cab = gel(ydata,4);
    pari_sp av = avma;

    /* 2*n points (0 is not an integration point) */
    n = lg(gc) - 1;
    res = zerovec(g);

    for (l = 1; l <= n; l++)
    {
        GEN yinv, xl = gel(gc, l);
        /* on +x */
        yinv = ginv(sqrt_pol_def(u, xl, prec));
        /* differentials x^k/y */
        res = gadd(res, gpowers0(xl,g-1,yinv));
        
        /* same on -x */ 
        yinv = ginv(sqrt_pol_def(u, gneg(xl), prec));
        res = gadd(res, gpowers0(gneg(xl),g-1,yinv));
    }
    /* multiply by weight = Pi / (2n) * Cab */
    res = gmul(res,gdivgs(gmul(mppi(prec),cab),2*n));
    /* FIXME: constants x->u */
    settyp(res, t_COL);
    return gerepilecopy(av, res);
}

GEN
integrals_tree(GEN ydata, long g, long prec)
{
    GEN gc = NULL, mat, s;
    long n = 0, k, bitprec = prec2nbits(prec);
    pari_sp av = avma;

    mat = cgetg(lg(ydata), t_MAT);
    /* integrals are sorted hardest first */
    s = vecsort0(ydata,mkvecsmall(1),1);
    for (k = 1; k < lg(ydata); k++)
    {
        long nk = gc_cost(gmael3(ydata,s[k],1,1), bitprec);
        pari_printf("%ld->%ld points\n",k,nk);
        if (nk > n || (nk - n) * g > nk)
        {
            pari_warn(warner,"compute %ld integration points", nk);
            gc = gc_init(nk, prec);
            n = nk;
        }
        gel(mat, k) = integral_edge(gmael(ydata,s[k],2),g,gc,prec);
    }
    return gerepilecopy(av, mat);
}

/*********************************************************************/
/*                                                                   */
/*              Integrals along minimum spanning tree                */
/*                                                                   */
/*********************************************************************/

/* distance (ellipse excentricity) */
static GEN
r_3(GEN a, GEN b, GEN c, long prec)
{
    return gdiv(gadd(gabs(gsub(c, a), prec), gabs(gsub(c, b), prec)),
            gabs(gsub(b, a), prec));
}

static GEN
param_edge(GEN X, long i, long j, long prec)
{
    long k, l1 = lg(X) - 1;
    GEN rij = mkoo();
    for (k = 1; k <= l1; ++k)
    {
        if (k == i || k == j)
            continue;
        rij = gmin(rij, r_3(gel(X, i), gel(X, j), gel(X, k), prec));
    }
    return rij;
}

static GEN
complete_graph(GEN X, long n, long prec)
{
    long i, j, k;
    pari_sp av = avma;
    GEN G = cgetg( n*(n-1)/2+1, t_VEC);
    k = 1;
    for (i = 1; i <= n; i++)
        for (j = i + 1; j <= n; j++)
            gel(G, k++) = mkvec2(param_edge(X, i, j, prec), mkvecsmall2(i,j));
    return gerepilecopy(av, lexsort(G));
}

/* compute best spanning tree */
GEN hc_spanning_tree(GEN X, long prec) {
    long n, k, len, nedges;
    GEN G, tree, t; /* vertices already connected */
    pari_sp av = avma;
    n = lg(X) - 1;
    G = complete_graph(X, n, prec);
    //pari_printf("graph: %Ps\n",G);
    len = lg(G) - 1;
    nedges = n - 1;

    tree = cgetg(nedges + 1, t_VEC);
    t = const_vecsmall(n, 0);
    for (k = 1; k <= nedges; k++)
    {
        long a, b, i;
        /* consider next edge with exactly one vertex taken (no cycle) */
        for (i = len; k > 1 && t[gmael(G,i,2)[1]] == t[gmael(G,i,2)[2]]; i--);
        /* this is the best edge allowed */
        a = gmael(G,i,2)[1]; b = gmael(G,i,2)[2];
        /* ensure a already in tree, flip edge if needed */
        if (t[b])
        {
            gmael(G,i,2)[1] = b;
            gmael(G,i,2)[2] = a;
        }
        t[a] = t[b] = 1;
        gel(tree, k) = gel(G,i);
    }
    /* FIXME: do not order by difficulty yet */
    return gerepilecopy(av, tree);
    return gerepilecopy(av, lexsort(tree));
}

/*********************************************************************/
/*                                                                   */
/*                 Symplectic pairing and basis                      */
/*                                                                   */
/*********************************************************************/

/* intersections */
#if 0
int
shift_number(GEN yab, GEN yad, long prec)
{
   angle = gdiv(garg(gmael(yab,2, mppi(prec));
}

int
intersection(GEN ab, GEN cd, int sgn)
{
    long a=ab[1],b=ab[2],c=cd[1],d=cd[2];
    if (a==c&&b==d || a==d&&b==c)
        return 0;
    else if (a==c)
        return intersection_abad(X,ab,cd);
    else if (b==c)
        return intersection_abbd(X,ab,cd);
    else if (a==d)
        return -intersection_abbd(X,cd,ab);
    else if (b==d)
}
#endif


GEN
intersections_tree(GEN ydata)
{
    long k, l, n = lg(ydata)-1;
    GEN mat;
    pari_sp av = avma;

    mat = zeromatcopy(n,n);
    
    for (k = 1; k <= n; k++)
    {
        GEN ek,yab,fa,fb;
        ek = gmael3(ydata,k,1,2);
        yab = gmael(ydata,k,2);
        fa = gel(yab,5);
        fb = gel(yab,6);
        for (l = k+1; l <= n; l++)
        {
            GEN el,ycd;
            el = gmael3(ydata,l,1,2);
            /* FIXME */
            //if (ek[1] != el[1] && ek[2] != el[1])
            //    continue; /* no intersection */
            pari_printf("intersection %Ps . %Ps\n",ek,el);
            ycd = gmael(ydata,l,2);
            if(el[1] == ek[1])
            {
                /* case ab.ad */
                GEN fc = gel(ycd,5);
                pari_printf("[ab.ad], ratio %Ps\n", gdiv(fa,fc));
                gcoeff(mat,k,l) = stoi(signe(gimag(gdiv(fa,fc))));
                gcoeff(mat,l,k) = gneg(gcoeff(mat,k,l));
            }
            else if(el[1]==ek[2])
            {
                /* case ab.bd */
                GEN fc = gel(ycd,5);
                pari_printf("[ab.bd], ratio %Ps\n", gdiv(fb,fc));
                gcoeff(mat,k,l) = stoi(signe(gimag(gdiv(fb,fc))));
                gcoeff(mat,l,k) = gneg(gcoeff(mat,k,l));
            }
            else if (ek[1] != el[1] && ek[2] != el[1] && ek[1] != el[2] && ek[2] != el[2])
            {
                continue; /* no intersection */
            }
            else
            {
                pari_err_BUG("hyperellperiods, intersection case should not occur");
                return NULL;
            }
        }
    }
    output(mat);
    return gerepilecopy(av, mat);
}

/* compute symplectic homology basis */

/* exchange rows and columns i,j, in place */
static void
row_swap(GEN m, long i1, long i2)
{
    long k;
    for (k = 1; k < lg(m); k++)
        swap(gcoeff(m,i1,k),gcoeff(m,i2,k));
}
static void
col_swap(GEN m, long j1, long j2)
{
    swap(gel(m,j1),gel(m,j2));
}
static void
swap_step(GEN p, GEN m, long i1, long i2)
{
    if (i1==i2)
        return;
    //pari_printf("SWAP [%ld,%ld]\n",i1,i2);
    col_swap(p, i1, i2);
    row_swap(m, i1, i2);
    col_swap(m, i1, i2);
    //outmat(m);
}
/* ci1 <- u ci1 + v ci2, ci2 <- u1 ci1 + v1 ci2, in place */
static void
row_bezout(GEN m, long i1, long i2, GEN u, GEN v, GEN u1, GEN v1)
{
    long j;
    for (j = 1; j < lg(m); j++)
    {
        GEN mi1j = gcoeff(m,i1,j), mi2j = gcoeff(m,i2,j);
        gcoeff(m,i1,j) = gadd(gmul(u,mi1j),gmul(v,mi2j));
        gcoeff(m,i2,j) = gadd(gmul(u1,mi1j),gmul(v1,mi2j));
    }
}
static void
col_bezout(GEN m, long j1, long j2, GEN u, GEN v, GEN u1, GEN v1)
{
    GEN mj1 = gel(m,j1), mj2 = gel(m,j2);
    gel(m,j1) = gadd(gmul(u,mj1),gmul(v,mj2));
    gel(m,j2) = gadd(gmul(u1,mj1),gmul(v1,mj2));
}
/* (i,k) <- (u i + v * k, u1 i + v1 k)*/
static void
bezout_step(GEN p, GEN m, long i, long k, GEN a, GEN b)
{
    GEN d,u,v,u1,v1;
    d = gbezout(a,b,&u,&v);
    u1 = gneg(gdiv(b,d));
    v1 = gdiv(a,d);
    //pari_printf("BEZOUT [%ld,%ld] <- (%Ps,%Ps;%Ps,%Ps)*[%ld,%ld]\n",i,k,u,u1,v,v1,i,k);
    col_bezout(p,i,k,u,v,u1,v1);
    row_bezout(m,i,k,u,v,u1,v1);
    col_bezout(m,i,k,u,v,u1,v1);
    //outmat(m);
}
/* i <- i + q * k */
static void
transvect_step(GEN p, GEN m, long i, long k, GEN q)
{
    GEN u = gen_1, v = q, u1 = gen_0, v1 = gen_1;
    //pari_printf("TRANS [%ld] <- [%ld] - %Ps.[%ld]\n",i,i,q,k);
    col_bezout(p,i,k,u,v,u1,v1);
    row_bezout(m,i,k,u,v,u1,v1);
    col_bezout(m,i,k,u,v,u1,v1);
    //outmat(m);
}
/* choose +/-1 or smallest element */
static long
pivot_line(GEN m, long i, long len)
{
    long j, jx = 0;
    GEN x = NULL;
    for (j = 1; j <= len; j++)
    {
        GEN z = gcoeff(m,i,j);
        if (z==gen_0)
            continue;
        if (z==gen_1 || z==gen_m1)
            return j;
        if (x == NULL || abscmpii(x, z) > 0)
            x = z, jx = j;
    }
    return jx;
}
/* returns p s.t. p~ *m*p = J_g(d), d vector of diagonal coefficients */
GEN
symplectic_reduction(GEN m, long g, GEN * d)
{
    long dim, i, j, k, len;
    GEN p, diag;
    pari_sp av = avma;

    len = lg(m)-1;
    if (!gequal(m,gneg(gtrans(m))))
        pari_err_DOMAIN("symplectic_reduction","antisymmetry","fails in", m, NULL);
    if (2 * g > len)
        pari_err_DOMAIN("symplectic_reduction","dimension",">", stoi(g), stoi(len));

#define m(i,j) gcoeff(m,i,j)

    p = matid(len);
    diag = zerovec(g);
    m = shallowcopy(m);
    //pari_printf("ENTER SYMPLECTIC REDUCTION, g=%ld, m=%Ps\n",g,m);
    /* main loop on symplectic 2-subspace */
    for (dim = 0; dim < g; dim++)
    {
        int cleared = 0;
        i = 2 * dim + 1;
        //pari_printf("NOW dim=%ld,len=%ld,g=%ld\n",dim,len,g);
        //outmat(m);
        /* lines 0..2d-1 already cleared */
        while ((j = pivot_line(m, i, len)) == 0)
        {
            /* no intersection -> move ci to end */
            swap_step(p, m, i, len);
            len--;
            if (len == 2*dim)
            {
                //pari_printf("len=%ld, dim=%ld, early abort",dim,len);
                //outmat(m);
                //outmat(p);
                gerepileall(av, 2, &p, &diag);
                if (d) *d = diag;
                return p;
            }
        }
        //pari_printf("choose pivot i,j=%ld,%ld\n",i,j);

        /* move j to i + 1 */
        if (j != i+1)
            swap_step(p, m, j, i + 1);
        j = i + 1;

        /* make sure m(i, j) > 0 */
        if (signe(m(i, j)) < 0)
            swap_step(p, m, i, j);
            // instead of neg i or j...
            //neg_step(p, m, j);

        while(!cleared)
        {
            /* clear row i */
            for (k = j + 1; k <= len; k++)
            {
                if (m(i, k) != gen_0)
                {
                    if (gequal0(remii(m(i, k), m(i, j))))
                    {
                        GEN q = gneg(gdiv(m(i, k), m(i, j)));
                        transvect_step(p, m, k, j, q);
                    }
                    else
                        bezout_step(p, m, j, k, m(i, j), m(i, k));
                }
            }
            cleared = 1;
            /* clear row j */
            for (k = j + 1; k <= len; k++)
            {
                if (m(j, k) != gen_0)
                {
                    if (gequal0(remii(m(j, k), m(i, j))))
                    {
                        GEN q = gdiv(m(j, k), m(i, j));
                        transvect_step(p, m, k, i, q);
                    }
                    else
                    {
                        bezout_step(p, m, i, k, m(j, i), m(j, k));
                        /* warning: row i now contains some ck.cl... */
                        cleared = 0;
                    }
                }
            }
        }
        gel(diag, dim + 1) = m(i,j);
    }
    //pari_printf("reduction finished\n");
    //outmat(m);
    gerepileall(av, 2, &p, &diag);
    if (d) *d = diag;
    return p;
}

GEN hc_homology_basis(GEN hc) { return NULL; }

/*********************************************************************/
/*                                                                   */
/*                        Period matrices                            */
/*                                                                   */
/*********************************************************************/

/* compute big period matrix */
GEN hc_big_periods(GEN hc) { return NULL; }

/* compute small period matrix */
GEN hc_small_periods(GEN hc) { return NULL; }

/*********************************************************************/
/*                                                                   */
/*                   Hyperelliptic curve object                      */
/*                                                                   */
/*********************************************************************/

GEN
hcinit(GEN pol, long prec)
{
    GEN hc, X, tree, ydata, integrals, mat, p, jg, periods;
    pari_sp av = avma;
    long g = (poldegree(pol,-1) - 1) / 2;
    X = roots(pol, prec);
    //pari_printf("roots: %Ps\n",X);
    tree = hc_spanning_tree(X, prec);
    //pari_printf("tree: %Ps\n",tree);
    ydata = ydata_tree(X, tree, prec);
    //pari_printf("ydata: %Ps\n",ydata);
    integrals = integrals_tree(ydata, g, prec);
    mat = intersections_tree(ydata);
    pari_printf("intersections: %Ps\n",mat);
    p = symplectic_reduction(mat, g, NULL);
    pari_printf("intersections\n");
    outmat(mat);
    jg = gmul(gmul(gtrans(p),mat),p);
    outmat(jg);
    pari_printf("change of basis\n");
    outmat(p);
    periods = gmul(integrals,p);
    pari_printf("periods\n");
    outmat(periods);
    hc = mkvec4(pol, X, periods, mkvec3(tree,integrals,p));
    return gerepilecopy(av, hc);
}

GEN
hcperiods(GEN hc, long prec)
{
    if (typ(hc) == t_VEC && lg(hc) == 5 && lg(gel(hc,4)) == 3)
    {
        /* hc */
        return gel(hc,3);
    }
    else if (typ(hc) == t_POL)
    {
        pari_sp av = avma;
        hc = hcinit(hc, prec);
        return gerepilecopy(av, hcperiods(hc,prec));
    }
    pari_err_TYPE("hcperiods",hc);
    return NULL;
}
