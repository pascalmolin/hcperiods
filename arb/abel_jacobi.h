/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include <arb.h>
#include <arb_mat.h>
#include <acb.h>
#include <acb_mat.h>
#include <flint/fmpz_poly.h>
#include <arb_fmpz_poly.h>

#include "fmpz_mat_extras.h"
#include "complex_extras.h"
#include "mag_func.h"

#define VBS 0
#define progress(...) if (VBS) flint_printf(__VA_ARGS__);
#define DEBUG 0

/******************************************************************************

  data structures

 ******************************************************************************/

typedef struct
{
    /* y^m = prod_{i=1}^d x - roots[i] */
    slong m;             /* degree in y */
    slong n;             /* degree in x */
    fmpz_poly_t pol;     /* polynomial */
    slong delta;         /* default = gcd(m, d) */
    slong g;             /* genus, 2g = (m-1)(d-1) - delta + 1 */
    slong j1;
    slong nj;
    slong * ni;
}
superelliptic_curve;

typedef superelliptic_curve sec_t;

enum { INT_D2, INT_GC, INT_DE };

typedef struct
{
    arf_t l;       /* lambda, take Pi/2 */
    arf_t h;       /* step size, exact */
    ulong n;       /* number of points */
    slong prec;    /* precision */
    arb_t factor;  /* lambda*h */
    arb_ptr x;     /* tanh(lambda*sinh(k*h)) */
    arb_ptr dx;    /* cosh(k*h)/cosh(lambda*sinh(k*h))^2 */
    arb_ptr ch2m;  /* cosh(lambda*sinh(k*h))^(2/m) */
    mag_t e;       /* error bound */
}
de_integration_struct;
typedef de_integration_struct de_int_t[1];

typedef struct
{
    slong n;
    slong len;
    arb_ptr x;
}
gc_integration_struct;
typedef gc_integration_struct gc_int_t[1];

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

/* differential forms */

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
    /* curve */
    sec_t c;

    /* branch points */
    acb_ptr roots;

    /* integration type: Gauss or DE */
    int type;

    /* homology using tree */
    tree_t tree;           /* choice of best loops */

    /* A,B symplectic basis as sums of shifted gamma loops */
    homol_t loop_a;
    homol_t loop_b;

    /* periods  */
    acb_mat_t integrals;
    acb_mat_t omega0;
    acb_mat_t omega1;      /* on A,B basis */
    acb_mat_t tau;         /* tau = Omega1^{-1}Omega0 */

    /* for Abel-Jacobi map, choice of a base point */
    slong p0;              /* base point */
    arb_mat_t proj;
}
abel_jacobi_struct;

typedef abel_jacobi_struct abel_jacobi_t[1];

/*  flag bits */
enum {
    AJ_USE_DE    = 1 << 0,  /* force DE integration */
    AJ_NO_TAU    = 1 << 1,  /* do not compute tau */
    AJ_NO_AB     = 1 << 2,
    AJ_NO_INT    = 1 << 3,
    AJ_TRIM      = 1 << 4,
    AJ_ROOT_DEF  = 1 << 5,
    AJ_ROOT_TURN = 1 << 6,
    AJ_ROOT_PROD = 1 << 7,
    AJ_DE_SAME   = 1 << 8,
    AJ_VERBOSE   = 5 << 9
};


/******************************************************************************

  functions

 ******************************************************************************/
void sec_init(sec_t * c, slong m, slong n);
void sec_init_poly(sec_t * c, slong m, const fmpz_poly_t pol);
void sec_clear(sec_t c);

void abel_jacobi_init_poly(abel_jacobi_t aj, slong m, const fmpz_poly_t f);
void abel_jacobi_compute(abel_jacobi_t aj, int flag, slong prec);
void abel_jacobi_clear(abel_jacobi_t aj);

/* parameters for DE integration */
void arb_vec_r(arb_t r0, arb_ptr vr, void (*f)(arb_t r, const acb_t u, arb_srcptr l, slong prec), acb_srcptr u, slong len, slong prec);
slong de_params_d(double * h, double * lambda, double * r, const cdouble * w, slong len, slong i, slong m, slong prec);
slong de_params(mag_t e, arf_t h, arf_t l, acb_srcptr u, slong len, double r, slong d, slong j, slong m, slong prec);
slong de_params_tree(mag_t e, arf_t h, arf_t l, const tree_t tree, sec_t c, slong prec);
void de_int_init(de_int_t de, const arf_t h, const arf_t l, ulong n, const mag_t e, ulong m, slong prec);
void de_int_clear(de_int_t de);

/* parameters for GC integration */
slong gc_params(mag_t e, acb_srcptr u, slong len, slong i, slong prec);
slong gc_params_tree(mag_t e, const tree_t tree, sec_t c, slong prec);
void gc_int_init(gc_int_t gc, slong n, slong prec);
void gc_int_clear(gc_int_t gc);

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
slong extraprec_tree(tree_t tree, acb_srcptr x, sec_t c);

/* compute local intersections between tree edges */
/* -> (d-1)*(d-1) intersection matrix */
void intersection_tree(fmpz_mat_t c, const tree_t tree, slong m);

/* find g+g symplectic homology basis from tree */
/* two lists of g loops */
void loop_init(loop_t * l, slong len);
void loop_print(const loop_t l);
void homol_init(homol_t * cyc, slong len);
void homol_clear(homol_t l, slong len);
void symplectic_reduction(fmpz_mat_t p, fmpz_mat_t m, slong g);
void symplectic_basis(homol_t alpha, homol_t beta, const tree_t tree, sec_t c);

/* numerically compute d-1 integrals along tree edges */
/* (d-1)*(g-1) matrix, tree edges on lines */
void gc_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, slong g, slong n, int flag, slong prec);
void de_integrals_precomp(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c, const de_int_t de, int flag, slong prec);
void de_integrals(acb_ptr res, acb_srcptr u, slong d1, slong d, sec_t c, int flag, slong prec);
void integrals_edge_factors_gc(acb_ptr res, const acb_t ba2, const acb_t ab, const acb_t cab, sec_t c, slong prec);
void integrals_edge_factors(acb_ptr res, const acb_t ba2, const acb_t ab, const acb_t cab, sec_t c, slong prec);
void integrals_tree_de(acb_mat_t integrals, sec_t c, const tree_t tree, int flag, slong prec);
void integrals_tree_gc(acb_mat_t integrals, sec_t c, const tree_t tree, int flag, slong prec);
void integral_d2(acb_ptr res, ydata_t ye, sec_t c, slong prec);

/* get all periods on a, b basis */
/* two g*g matrices */
void period_matrix(acb_mat_t omega, const homol_t basis, const acb_mat_t integrals, sec_t c, slong prec);

/* get tau reduced matrix */
void tau_matrix(acb_mat_t tau, const acb_mat_t omega0, const acb_mat_t omega1, slong prec);

/* core functions */
void sqrt_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong prec);
void sqrt_pol_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, slong prec);
void mth_root_pol_def(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec);
void mth_root_pol_prod(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec);
void mth_root_pol_turn(acb_t y, acb_srcptr u, slong d1, slong d, const arb_t x, acb_srcptr z, slong m, slong prec);

/* other */
typedef struct { double r; slong k; } comp_t;
int comp_cmp(const comp_t * x, const comp_t * y);

/* vec utilities */
/*void _acb_vec_scalar_addmul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec);*/
void _acb_vec_add_error_mag(acb_ptr res, slong len, const mag_t e);
void _acb_vec_scalar_addmul_fmpz(acb_ptr res, acb_srcptr vec, slong len, const fmpz_t c, slong prec);
void acb_vec_polynomial_shift(acb_ptr x, slong len, const acb_t c, slong prec);
void acb_vec_mul_geom(acb_ptr x, slong len, acb_t c0, const acb_t c, slong prec);
void acb_vec_add_geom_arb(acb_ptr x, slong len, acb_t c0, const arb_t c, slong prec);
void acb_vec_sub_geom_arb(acb_ptr x, slong len, acb_t c0, const arb_t c, slong prec);
void _acb_vec_sort_lex(acb_ptr vec, slong len);
#define acb_mat_trim(m) _acb_vec_trim(m->entries, m->entries, m->r * m->c)

/* tests/bench only */
void arb_vec_set_random_nonzero(arb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits);
void acb_vec_set_random(acb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits);
void acb_vec_set_random_u(acb_ptr u, slong len, flint_rand_t state, slong prec, slong mag_bits, double eps);
void _acb_vec_printd(acb_srcptr u, slong len, slong d, const char * sep);
void _acb_vec_arf_printd(acb_srcptr u, slong len, slong d, const char * sep);
void acb_mat_print_gp(const acb_mat_t m, slong digits);
void acb_mat_print_error(const acb_mat_t m, slong digits);
