/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "arb.h"
#include "arb_mat.h"
#include "acb.h"
#include "acb_mat.h"

#include "fmpz_mat_extras.h"
#include "complex_extras.h"
#include "mag_func.h"

#define VBS 0
#define progress(...) if (VBS) flint_printf(__VA_ARGS__);

/******************************************************************************

  data structures

 ******************************************************************************/

typedef struct
{
    /* y^m = prod_{i=1}^d x - roots[i] */
    slong m;             /* degree in y */
    slong n;             /* degree in x */
    acb_ptr roots;       /* branch points */
    slong delta;         /* default = gcd(m, d) */
    slong g;             /* genus, 2g = (m-1)(d-1) - delta + 1 */
}
superelliptic_curve;

typedef superelliptic_curve sec_t;

enum { INT_GC, INT_DE };

typedef struct
{
    arf_t h;       /* step size, exact */
    ulong n;       /* number of points */
    slong prec;    /* precision */
    arb_t factor;  /* lambda*h, lambda=Pi/2 */
    arb_ptr x;     /* tanh(lambda*sinh(k*h)) */
    arb_ptr dx;    /* cosh(k*h)/cosh(lambda*sinh(k*h))^2 */
    arb_ptr ch2m;  /* cosh(lambda*sinh(k*h))^(2/m) */
}
de_integration_struct;
typedef de_integration_struct de_int_t[1];

typedef struct
{
    slong m;
    slong n;
    slong n1;
    acb_ptr u; /* n - 2 + 5 */
    acb_ptr ba2;
    acb_ptr ab;
    acb_ptr c;
    acb_ptr ya;
    acb_ptr yb;
} ydata_struct;
typedef ydata_struct ydata_t[1];

typedef struct
{
    double r;
    slong a;
    slong b;
    /* FIXME: here or below ? */
    ydata_t data;
} edge_t;

typedef struct
{
    slong n;
    slong min;
    double r;
    edge_t * e;
    ydata_struct * data;
} tree_struct;
typedef tree_struct tree_t[1];

/* represents the loop_t(coeff * zeta^shift * tree[index]) */
typedef struct
{
    slong coeff;
    slong shift;
    slong index;
} gamma_t;

typedef struct
{
    slong n; /* length at most 2*g */
    gamma_t * l;
} loop_t;

/* homology basis */
typedef loop_t * homol_t;

/* differential form */

/* i,j s.t. dj >= mi + delta */
inline slong
jmin(slong m, slong n, slong d)
{
    return (1 + (m + d - 1) / n);
}
inline slong
imax(slong j, slong m, slong n, slong d)
{
    slong i = ((j * n - d) / m);
    return (i < n) ? i : n - 1;
}

typedef struct
{
    slong x; /* power of x */
    slong y; /* power of y */
} dform_t;

/* cohomology basis */
typedef dform_t * cohom_t;

typedef struct
{
    /* curve */
    sec_t c;

    /* differentials */
    cohom_t dz;

    /* integration type: Gauss or DE */
    int type;

    /* homology using tree */
    tree_t tree;           /* choice of best loops */

    /* A,B symplectic basis as sums of shifted gamma loops */
    homol_t loop_a;
    homol_t loop_b;

    /* periods  */
    acb_mat_t omega0;
    acb_mat_t omega1;      /* on A,B basis */
    acb_mat_t tau;         /* tau = Omega1^{-1}Omega0 */

    /* for Abel-Jacobi map, choice of a base point */
    slong p0;              /* base point */
    arb_mat_t proj;
}
abel_jacobi_struct;

typedef abel_jacobi_struct abel_jacobi_t[1];

/******************************************************************************

  functions

 ******************************************************************************/
void sec_init(sec_t * c, slong m, acb_srcptr x, slong d);
void sec_clear(sec_t c);

void abel_jacobi_init_roots(abel_jacobi_t aj, slong m, acb_srcptr x, slong d);
void abel_jacobi_init_poly(abel_jacobi_t aj, slong m, acb_poly_t f, slong prec);
void abel_jacobi_compute(abel_jacobi_t aj, slong prec);
void abel_jacobi_clear(abel_jacobi_t aj);

/* parameters for DE integration */
slong de_params_d(double *ph, const cdouble * w, slong len, double r, slong i, slong m, slong prec);
slong de_params(double * h, acb_srcptr u, slong len, double r, slong i, slong m, slong prec);
slong de_params_tree(double * h, const tree_t tree, sec_t c, slong prec);
void de_int_init(de_int_t de, double h, ulong n, ulong m, slong prec);
void de_int_clear(de_int_t de);

/* parameters for GC integration */
slong gc_params_d(const cdouble * w, slong len, double r, slong i, slong prec);
slong gc_params(acb_srcptr u, slong len, double r, slong i, slong prec);
slong gc_params_tree(const tree_t tree, sec_t c, slong prec);

/* compute maximum spanning tree */
void tree_init(tree_t tree, slong d);
void tree_clear(tree_t tree);
void spanning_tree(tree_t tree, acb_srcptr x, slong d, int type);
void shift_info_tree(tree_t tree, cdouble * w, slong d);
void tree_print(const tree_t tree);

/* reduced points of spanning tree */
void ydata_init_edge(ydata_t yab, acb_srcptr x, edge_t e, slong n, slong m, slong prec);
void ydata_clear(ydata_t yab);
void tree_ydata_init(tree_t tree, acb_srcptr x, slong n, slong m, slong prec);
void tree_ydata_clear(tree_t tree);

/* compute local intersections between tree edges */
/* -> (d-1)*(d-1) intersection matrix */
void intersection_tree(fmpz_mat_t c, const tree_t tree, slong n, slong m);

/* find g+g symplectic homology basis from tree */
/* two lists of g loops */
void loop_init(loop_t * l, slong len);
void homol_init(homol_t * cyc, slong len);
void homol_clear(homol_t l, slong len);
void symplectic_reduction(fmpz_mat_t p, fmpz_mat_t m, slong g);
void symplectic_basis(homol_t alpha, homol_t beta, const tree_t tree, sec_t c);

/* find basis of holomorphic differentials */
/* g elementary differentials */
void holomorphic_differentials(cohom_t dz, slong d, slong m);

/* numerically compute d-1 integrals along tree edges */
/* (d-1)*(g-1) matrix, tree edges on lines */
void gc_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, slong g, slong n, slong prec);
void de_integrals_precomp(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c, const cohom_t dz, const de_int_t de, slong prec);
void de_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c, slong prec);
void integrals_edge_factors_gc(acb_ptr res, const acb_t cab, const acb_t ba2, sec_t c, slong prec);
void integrals_edge_factors(acb_ptr res, const acb_t cab, const acb_t ba2, sec_t c, const cohom_t dz, slong prec);
void integrals_tree_de(acb_mat_t integrals, sec_t c, const tree_t tree, const cohom_t dz, slong prec);
void integrals_tree_gc(acb_mat_t integrals, sec_t c, const tree_t tree, slong prec);

/* get all periods on a, b basis */
/* two g*g matrices */
void period_matrix(acb_mat_t omega, const homol_t basis, const cohom_t dz, const acb_mat_t integrals, sec_t c, slong prec);

/* get tau reduced matrix */
void tau_matrix(acb_mat_t tau, const acb_mat_t omega0, const acb_mat_t omega1, slong prec);

/* core functions */
void sqrt_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong prec);
void mth_root_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong m, slong prec);
void mth_root_pol_prod(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong m, slong prec);
void mth_root_pol_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec);

/* vec utilities */
void acb_vec_polynomial_shift(acb_ptr x, slong len, const acb_t c, slong prec);
void acb_vec_mul_geom(acb_ptr x, slong len, acb_t c0, const acb_t c, slong prec);
void acb_vec_add_geom_arb(acb_ptr x, slong len, acb_t c0, const arb_t c, slong prec);
void acb_vec_sub_geom_arb(acb_ptr x, slong len, acb_t c0, const arb_t c, slong prec);

/* tests/bench only */
void acb_vec_set_random(acb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits);
void acb_vec_set_random_u(acb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits, double eps);
void _acb_vec_printd(acb_srcptr u, slong len, slong d, const char * sep);
void _acb_vec_arf_printd(acb_srcptr u, slong len, slong d, const char * sep);
