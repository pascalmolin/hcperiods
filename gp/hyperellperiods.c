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

 periods are combination of tree integrals

 make things dynamic ?
*/

#define hc_pol(c)        gel(c,1)
#define hc_roots(c)      gel(c,2)
#define hc_periods(c)    gel(c,3)
#define hc_tree(c)       gmael(c,4,1)
#define hc_integrals(c)  gmael(c,4,2)
#define hc_hombasis(c)   gmael(c,4,3)
#define hc_n(c) poldegree(hc_pol(c))

/** everything is dynamic **/

/*********************************************************************/
/*                                                                   */
/*             Local square roots and integrals                      */
/*                                                                   */
/*********************************************************************/

/*
 On [a,b], x = m*(u+u0), u<-[-1,1], u0 = (a+b)/(b-a), m=(b-a)/2
 and f(x) = - m^n.(1-u^2).prod(u-u_i), u_i = (2x_i-b-a)/(b-a)
 so y(x) = C sqrt(1-u^2) g(u)
 with C = I * sqrt(m)^n
 and g(u)= prod sqrt(u-u_i) holomorphic near [-1,1].

 Hence the period 2*int_a^b x^{k-1}dx/y is
 2 m^k/C * int_-1^1 (u+u0)^{k-1}/g(u) du/sqrt(1-u^2)
 */

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
  return [ [ u_i ], m, u0, 2/C, fa, fb ]
  where u_i = (2x_i-a-b)/(b-a), m = (b-a)/2, u0 = (a+b)/(b-a), C = I*m^n/2, fa, fb ]
  and fa = C*g(-1), fb = C*g(1).
 */
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

    cab = gmul(gen_I(),gpowgs(gsqrt(ba,prec),d));
    fa = gmul(cab, sqrt_pol_def(u,gen_m1,prec));
    fb = gmul(cab, sqrt_pol_def(u,gen_1,prec));
    cab = gdivsg(2,cab);

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
/* y[k] = sum C_k^i c^(k-i) x[i] */
static GEN
binomial_transform(GEN x, GEN c)
{
    long k,l,n = lg(x)-1;
    pari_sp av = avma;
    GEN y;
    y = leafcopy(x);
    for (k = 2; k <= n; k++)
        for (l = n; l >= k; l--)
            gel(y, l) = gadd(gel(y,l),gmul(c,gel(y,l-1)));
    return gerepilecopy(av,y);
}
/* x[0],c.[x1],c^2x[2]... */
static GEN
geom_shift(GEN x, GEN c)
{
    long k;
    GEN d = gen_1, y = cgetg(lg(x),t_VEC);
    for (k = 1; k < lg(x); k++)
    {
        d = gmul(d,c);
        gel(y, k) = gmul(d, gel(x, k));
    }
    return y;
}
GEN
integral_edge(GEN ydata, long g, GEN gc, long prec)
{
    long l, n;
    GEN res, u = gel(ydata, 1), ba2 = gel(ydata,2);
    GEN u0 = gel(ydata,3), cab = gel(ydata,4);
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
        if (gc_needed(av, 3))
            res = gerepilecopy(av, res);
    }
    /* multiply by Pi / (2n) * Cab */
    res = gmul(res,gdivgs(gmul(mppi(prec),cab),2*n));
    //output(res);
    /* mul by ba2^k and shift by u0 */
    res = binomial_transform(res,u0);
    //output(res);
    res = geom_shift(res, ba2);
    //res = geom_shift(res, gdiv(ba2,gen_I()));
    //output(res);
    settyp(res, t_COL);
    return gerepilecopy(av, res);
}

GEN
integrals_tree(GEN ydata, long g, long prec)
{
    GEN gc = NULL, mat, s;
    long n = 0, k, bitprec = prec2nbits(prec);
    pari_sp av = avma;
    pari_timer ti;

    mat = cgetg(lg(ydata), t_MAT);
    /* integrals are sorted hardest first */
    s = vecsort0(ydata,mkvecsmall(1),1);
    if (DEBUGLEVEL)
        timer_start(&ti);
    for (k = 1; k < lg(ydata); k++)
    {
        long nk = gc_cost(gmael3(ydata,s[k],1,1), bitprec);

        if (nk > n || n > 1.3 * nk)
        {
            pari_warn(warner,"compute %ld integration points", nk);
            gc = gc_init(nk, prec);
            n = nk;
        }
        gel(mat, s[k]) = integral_edge(gmael(ydata,s[k],2),g,gc,prec);
        if (DEBUGLEVEL)
        {
            GEN e = gmael3(ydata,s[k],1,2);
            timer_printf(&ti,"integral [%ld->%ld], %ld points", e[1], e[2], nk);
        }
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
    nedges = n % 2 ? n - 1 : n - 2;

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
    /* important: do not order by complexity yet */
    return gerepilecopy(av, tree);
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
            //pari_printf("intersection %Ps . %Ps\n",ek,el);
            ycd = gmael(ydata,l,2);
            if(el[1] == ek[1])
            {
                /* case ab.ad */
                GEN fc = gel(ycd,5);
                //pari_printf("[ab.ad], ratio %Ps\n", gdiv(fa,fc));
                gcoeff(mat,k,l) = stoi(signe(gimag(gdiv(fa,fc))));
                gcoeff(mat,l,k) = gneg(gcoeff(mat,k,l));
            }
            else if(el[1]==ek[2])
            {
                /* case ab.bd */
                GEN fc = gel(ycd,5);
                //pari_printf("[ab.bd], ratio %Ps\n", gdiv(fb,fc));
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
    col_swap(p, i1, i2);
    row_swap(m, i1, i2);
    col_swap(m, i1, i2);
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
    col_bezout(p,i,k,u,v,u1,v1);
    row_bezout(m,i,k,u,v,u1,v1);
    col_bezout(m,i,k,u,v,u1,v1);
}
/* i <- i + q * k */
static void
transvect_step(GEN p, GEN m, long i, long k, GEN q)
{
    GEN u = gen_1, v = q, u1 = gen_0, v1 = gen_1;
    col_bezout(p,i,k,u,v,u1,v1);
    row_bezout(m,i,k,u,v,u1,v1);
    col_bezout(m,i,k,u,v,u1,v1);
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
    /* main loop on symplectic 2-subspace */
    for (dim = 0; dim < g; dim++)
    {
        int cleared = 0;
        i = 2 * dim + 1;
        /* lines 0..2d-1 already cleared */
        while ((j = pivot_line(m, i, len)) == 0)
        {
            /* no intersection -> move ci to end */
            swap_step(p, m, i, len);
            len--;
            if (len == 2*dim)
            {
                gerepileall(av, 2, &p, &diag);
                if (d) *d = diag;
                return p;
            }
        }

        /* move j to i + 1 */
        if (j != i+1)
            swap_step(p, m, j, i + 1);
        j = i + 1;

        /* make sure m(i, j) > 0 */
        if (signe(m(i, j)) < 0)
            swap_step(p, m, i, j);

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
    gerepileall(av, 2, &p, &diag);
    if (d) *d = diag;
    return p;
}

GEN symplectic_homology_basis(GEN mat, long g) {
    long k;
    GEN p, p1;
    pari_sp av = avma;
    p1 = symplectic_reduction(mat,g,NULL);
    p = cgetg(2*g+1,t_MAT);
    for (k=1;k<=g;k++)
    {
        gel(p,k) = gel(p1,2*k-1);
        gel(p,g+k) = gel(p1,2*k);
    }
    return gerepilecopy(av,p);
}

/*********************************************************************/
/*                                                                   */
/*                        Period matrices                            */
/*                                                                   */
/*********************************************************************/

/* big period matrix */
GEN
hc_big_periods(GEN hc) {
    return hc_periods(hc);
}
/* small period matrix */
GEN
hc_small_periods(GEN hc) {
    long g;
    GEN ab = hc_periods(hc);
    pari_sp av = avma;
    if (lg(ab) < 3) return cgetg(1,t_MAT);
    g = nbrows(ab);
    return gerepilecopy(av,gauss(vecslice(ab,g+1,2*g),vecslice(ab,1,g)));
}

/*********************************************************************/
/*                                                                   */
/*                   Hyperelliptic curve object                      */
/*                                                                   */
/*********************************************************************/

GEN
hcinit(GEN pol, long prec)
{
    GEN hc, X, tree, ydata, integrals, mat, ab, periods;
    pari_sp av = avma;
    pari_timer ti;
    long pr2 = ndec2prec(34);
    long g = (poldegree(pol,-1) - 1) / 2;
    if (g<1)
    {
        pari_warn(warner,"genus 0 curve, no period.");
        hc = mkvec4(pol,roots(pol,prec),cgetg(1,t_MAT),zerovec(3));
        return gerepilecopy(av,hc);
    }
    if (DEBUGLEVEL)
        timer_start(&ti);
    X = roots(pol, prec);
    if (DEBUGLEVEL)
        timer_printf(&ti,"roots");
    //pari_printf("roots: %Ps\n",X);
    tree = hc_spanning_tree(X, pr2);
    if (DEBUGLEVEL)
        timer_printf(&ti,"spanning tree");
    //pari_printf("tree: %Ps\n",tree);
    ydata = ydata_tree(X, tree, prec);
    if (DEBUGLEVEL)
        timer_printf(&ti,"prepare tree");
    //pari_printf("ydata: %Ps\n",ydata);
    integrals = integrals_tree(ydata, g, prec);
    integrals = gmul(integrals, gsqrt(pollead(pol,-1),prec));
    if (DEBUGLEVEL)
        timer_printf(&ti,"integrals");
    mat = intersections_tree(ydata);
    ab = symplectic_homology_basis(mat, g);
    if (DEBUGLEVEL)
        timer_printf(&ti,"symplectic basis");
    periods = gmul(integrals,ab);
    //pari_printf("periods\n");
    //outmat(periods);
    hc = mkvec4(pol, X, periods, mkvec3(tree,integrals,ab));
    return gerepilecopy(av, hc);
}

/* Richelot's algorithm
   X=[u,u',v,v',w,w']
   compute
   int_u^u' [dx,xdx]/sqrt(-(x-u)(x-u')(x-v)(x-v')(x-w)(x-w'))
*/
GEN
richelot(GEN X, long prec)
{
  GEN a, ap, b, bp, c, cp, D, ba, ca, cb, Ia,Ib,Ic;
  pari_sp av = avma;
  long bitprec;

  if (lg(X) != 7)
      pari_err_TYPE("richelot",X);
  X = sort(X);

  a = gel(X,1); ap = gel(X,2);
  b = gel(X,3); bp = gel(X,4);
  c = gel(X,5); cp = gel(X,6);
  D = gen_1;
  bitprec = prec2nbits(prec)-5;

#define neqprec(a,ap) (gexpo(gsub(a,ap)) > -bitprec)
  while (neqprec(a,ap) || neqprec(b,bp) || neqprec(c,cp))
  {
    GEN ba,ca,cb,bpa,cpa,cpb,bap,cap,cbp,bpap,cpap,cpbp;
    GEN A, B, C, dba, dca, dcb, mba, mca, mcb, aap, bbp, ccp;
    GEN num, den;

    ba = gsub(b,a); bap = gsub(b,ap); bpa = gsub(bp,a); bpap = gsub(bp,ap);
    ca = gsub(c,a); cap = gsub(c,ap); cpa = gsub(cp,a); cpap = gsub(cp,ap);
    cb = gsub(c,b); cbp = gsub(c,bp); cpb = gsub(cp,b); cpbp = gsub(cp,bp);

#define mul4(a,b,c,d) gmul(gmul(a,b),gmul(c,d))
    A = gsqrt(mul4(cb,cbp,cpb,cpbp),prec);
    B = gsqrt(mul4(ca,cap,cpa,cpap),prec);
    C = gsqrt(mul4(ba,bap,bpa,bpap),prec);

    aap = gmul(a,ap);
    bbp = gmul(b,bp);
    ccp = gmul(c,cp);

    dba = gadd(ba,bpap);
    dca = gadd(ca,cpap);
    dcb = gadd(cb,cpbp);

    mba = gsub(bbp,aap);
    mca = gsub(ccp,aap);
    mcb = gsub(ccp,bbp);

    a  = gdiv(gsub(mca,B),dca);
    ap = gdiv(gsub(mba,C),dba);

    b  = gdiv(gadd(mba,C),dba);
    bp = gdiv(gsub(mcb,A),dcb);

    c  = gdiv(gadd(mcb,A),dcb);
    cp = gdiv(gadd(mca,B),dca);

    num = gadd(gsub(gmul(aap,dcb), gmul(bbp,dca)), gmul(ccp,dba));
    den = gmul(dba, gmul(dca, dcb));
    D  = gmul2n(gdiv(gmul(D,num), den),2);

    gerepileall(av,7,&a,&ap,&b,&bp,&c,&cp,&D);
  }

  ba = gsub(b,a);
  ca = gsub(c,a);
  cb = gsub(c,b);
  D = gmul(Pi2n(1,prec),gsqrt(D,prec));

  Ia = gdiv(D,gmul(ba,ca));
  Ib = gdiv(D,gmul(ba,cb));
  Ic = gdiv(D,gmul(ca,cb));

  Ia = mkcol2(Ia,gmul(a,Ia));
  Ib = mkcol2(Ib,gmul(b,Ib));
  Ic = mkcol2(Ic,gmul(c,Ic));
  return gerepilecopy(av, mkmat3(Ia,Ib,Ic));
}

GEN
hyperellperiods(GEN hc, long flag, long prec)
{
    if (typ(hc) == t_VEC && lg(hc) == 5
            && typ(gel(hc,4)) == t_VEC && lg(gel(hc,4)) == 4)
    {
        return flag ? hc_big_periods(hc) : hc_small_periods(hc);
    }
    else if (typ(hc) == t_POL)
    {
        pari_sp av = avma;
        hc = hcinit(hc, prec);
        return gerepilecopy(av, hyperellperiods(hc,flag,prec));
    }
    pari_err_TYPE("hcperiods",hc);
    return NULL;
}

/*********************************************************************/
/*                                                                   */
/*                           SIEGEL REDUCTION                        */
/*                                                                   */
/*********************************************************************/

GEN
genus2periods(GEN C, long prec)
{
    GEN pol, m, Om, Omr, d;
    long j1, j2, pb = prec2nbits(prec) / 2;
    pari_sp av = avma;
    pol = typ(C) == t_POL ? C : gel(genus2red(C,NULL),3);
    m = greal(hyperellperiods(pol,1,prec));
    for(j1=1,j2=2;j1<=3;)
    {
        Om = mkmat2(gel(m,j1),gel(m,j2));
        if (gexpo(gabs(det(Om),prec))>-pb)
            break;
        if (++j2>4)
            j1++, j2=j1+1;
    }
    Omr = bestappr(gmul(ginv(Om),m),int2n(pb));
    d = denom_i(Omr);
    Omr = det(hnf(gmul(Omr, d)));
    Om = gabs(det(Om),prec);
    return gerepilecopy(av, gdiv(gmul(Om, Omr),gmul(d,d)));
}

#if 0
/*********************************************************************/
/*                                                                   */
/*                           SIEGEL REDUCTION                        */
/*                                                                   */
/*********************************************************************/
static GEN
test_matrix(int k)
{
    int l = 1;
    static GEN cols, L;

    if (L)
        return gel(L, k);

    cols = cgetg(10, t_VEC);
    gel(cols, k++) = mkcol2(gen_m1, gen_m1);
    gel(cols, k++) = mkcol2(gen_m1,  gen_0);
    gel(cols, k++) = mkcol2(gen_m1,  gen_1);
    gel(cols, k++) = mkcol2( gen_0, gen_m1);
    gel(cols, k++) = mkcol2( gen_0,  gen_0);
    gel(cols, k++) = mkcol2( gen_0,  gen_1);
    gel(cols, k++) = mkcol2( gen_1, gen_m1);
    gel(cols, k++) = mkcol2( gen_1,  gen_0);
    gel(cols, k++) = mkcol2( gen_1,  gen_1);

    L = cgetg(20, t_VEC);
    k = 1;
#define matcols(i,j) mkmat2(gel(cols,i),gel(cols,j))
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1,         gen_0);
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(8, 5));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(2, 5));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(5, 6));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(5, 4));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1,         gen_1);
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1,        gen_m1);
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(2, 6));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(8, 4));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(6, 8));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(4, 2));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(9, 8));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(1, 2));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(6, 9));
    gel(L, k++) = mkmat22( gen_0, gen_m1,         gen_1, matcols(4, 1));
    gel(L, k++) = mkmat22( gen_1, gen_m1, matcols(8, 5), matcols(5, 6));
    gel(L, k++) = mkmat22( gen_1, gen_m1, matcols(5, 6), matcols(8, 5));
    gel(L, k++) = mkmat22( gen_1,  gen_0, matcols(7, 3),         gen_1);
    gel(L, k++) = mkmat22(gen_m1,  gen_0, matcols(7, 3),        gen_m1);

    return gel(L,k);
}
static int
real_less12(GEN x) /* x <= 1/2 ? */
{ return (gcmp(gmul2n(gabs(greal(x)),1),gen_1) <= 0); }
static int
abs2_less(GEN x, GEN y) /* 2*|x| <= |y| ? */
{ return (gcmp(gmul2n(gabs(x),1),gabs(y)) <= 0); }
static int
abs_less(GEN x, GEN y) /* 2*|x| <= |y| ? */
{ return ((gcmp(gmul2n(gabs(x),1),gabs(y)) <= 0); }
static int
is_minkowski_reduced(GEN m)
{
    if (gcmp(gcoeff(m,2,2),gcoeff(m,1,1)) < 0)
        return 0;
    if (gcmp(gcoeff(m,1,1),gmul2n(gcoeff(m,1,2),1)) < 0)
        return 0;
    if (gsigne(gcoeff(m,1,2)) < 0)
        return 0;
    return 1;
}
int
fail_F2_criterion(GEN L19, GEN tau)
{
   int k;
   for (k = 1; k <= 19; k++)
   {
       GEN T;
       T = test_matrix(k);
       c = gcoeff(tau,2,1);
       d = gcoeff(tau,2,2);
       T = gadd(gmul(c,T),d);
       if (gcmp(gabs(det(T)),gen_1) < 0)
           return k;
   }
   return 0;
}
int
is_F2_reduced(GEN L19, GEN tau)
{
   if (!real_less12(gcoeff(tau,1,1))
     ||!real_less12(gcoeff(tau,1,2))
     ||!real_less12(gcoeff(tau,2,2)))
       return 0;
   if (!is_minkowski_reduced(gimag(tau)))
       return 0;
   return !fail_F2_criterion(L19, tau);
}
/* Minkowski reduction */
static GEN
congr(GEN m, GEN u)
{ return gmul(u, gmul(m, gtrans(u))); }
GEN
minkowski(GEN m)
{
  GEN u;
  pari_sp av = avma;
  int flag = 1;

  u = matid(2);
  do
  {
     GEN t, n;
     if (abs2_less(gcoeff(m,1,2),gcoeff(m,1,1)))
     {
         /* |m[1,1]| <= |m[2,2]| ? */
         if (abs_less(gcoeff(m,1,1),gcoeff(m,2,2)))
         {
             if (gsigne(gcoeff(m,1,2)) <= 0)
             {   
                 static GEN n = mkmat22(gen_1,gen_0,gen_0,gen_m1);
                 m = congr(m,n);
                 u = gmul(n,u);
             }
             flag = 0; //end while loop
         }
         else
         { 
   	      static GEN n = mkmat22(gen_0,gen_1,gen_m1,gem_1);
   	      m = congr(m,n);
   	      u = gmul(n,u);
         }
     }
     t = gneg(ground(dgiv(gcoeff(m,1,2),gcoeff(m,1,1))));
     n = mkmat22(gen_1, gen_0, t, gen_1);
     m = congr(m,n);
     u = gmul(n,u);
  } while (flag);
  return gerepileupto(av, mkvec2(u,m));
}
static GEN
sp4_action(GEN m, GEN tau)
{
    GEN a, b, c, d;
    a = gcoeff(m,1,1);
    b = gcoeff(m,1,2);
    c = gcoeff(m,2,1);
    d = gcoeff(m,2,2);
    a = gadd(gmul(a,tau),b);
    c = gadd(gmul(c,tau),d);
    return gdiv(a,c);
}
/* Compute \gamma and \tau' s.t. tau'\in\F_2 and \tau' = \gamma\tau */
GEN
reduce_to_F2(tau) =
{
   my(gma, U, N, j, a, b, c, d, FORGET);
   gma = matid(2);
   t = 1; \\true
   while (t,
	  \\Minkowski-reduce imaginary part
	  [U, FORGET] = minkowski(imag(tau));
	  N = [U, 0; 0, (U^(-1))~];
	  tau = sp4_action(N, tau);
	  gma = N * gma;
	  \\reduce real part
	  [N, tau] = reduce_real_part(tau);
	  gma = N * gma;
	  \\end loop ?
	  t = 0;
	  if(k = fail_F2_criterion(tau))
      {
        GEN n = test_matrix(k);
		tau = sp4_action(n, tau);
		gma = gmul(n, gma);
      }

	  \\check condition for the 19 test matrices
	  for (j=1, 19,
	       N = F2_test_matrix(j);
	       [c, d] = N[2,];
	       if (abs(matdet(c * tau + d)) < 1,
		   t = 1; \\loop again
		   tau = sp4_action(N, tau);
		   gma = N * gma;
		  );
	      );
	 );
   [gma, tau];
}
#endif
