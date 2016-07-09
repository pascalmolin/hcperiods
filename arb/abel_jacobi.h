/******************************************************************************
 
 Copyright (C) 2016 Pascal Molin
 
 ******************************************************************************/

#include "arb.h"
#include "arb_mat.h"
#include "acb.h"
#include "acb_mat.h"

typedef struct
{
    double tau;
    double h;
    double D;
    arb_t factor;
    ulong n;
    arb_ptr * x;
    arb_ptr * dx;
}
de_integration_struct;
typedef struct de_integration_struct de_int_t[1];

typedef struct tree_edge slong[2];
typedef struct spanning_tree tree_edge[];
/* represents the loop coeff * zeta^shift * tree[index] */
typedef struct
{
    slong coeff;
    slong shift;
    slong index;
} gamma_loop;
typedef loop gamma_loop[]; /* length at most 2*g */

/* intersection between two edges of the tree */
typedef struct
{
    slong coeff;
    slong shift;
} int_tree;

/*
typedef intersection slong;
typedef shift slong;
typedef struct edges_sum intersection[];
typedef struct edges_shift intersection[];
*/

typedef struct **slong slong_mat;

typedef struct
{
    /* y^n = prod_{i=1}^d x - roots[i] */
    slong n;             /* degree in y */
    slong d;             /* degree in x, odd */
    slong g;             /* genus, 2 g = (N-1)(d-1) */
    acb_ptr * roots;     /* weierstrass points */

    /* differentials */
    slong * dx;          /* power of x */
    slong * dy;          /* power of x */

    /* homology using tree */
    spanning_tree tree;        /* choice of nice loops */
    int_tree **inter;          /* intersections in tree, size d-1 * d-1 */

    loop ABtoC[];              /* A,B symplectic basis as sums of shifted gamma loops */

    /* periods  */
    acb_mat_t bigperiods;  /* on A,B basis */
    acb_mat_t tau;         /* tau = Omega1^{-1}Omega0 */


    /* for Abel-Jacobi map, choice of a bse point */
    slong p0;            /* base point */
    arb_mat_t proj;

}
super_elliptic_curve_struct;

typedef struct super_elliptic_curve_struct se_curve_t[1];
