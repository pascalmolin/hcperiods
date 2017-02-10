/******************************************************************************

 Copyright (C) Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

#define SKIP 1

typedef slong entry[3];
typedef slong pair[2];

static void
fmpz_mat_set_entries(fmpz_mat_t m, entry * v, slong len)
{
    slong k;
    fmpz_mat_zero(m);
    for (k = 0; k < len; k++)
        fmpz_set_si(fmpz_mat_entry(m,v[k][0],v[k][1]),v[k][2]);
}
static void
acb_vec_set_si(acb_ptr x, slong * v, slong len)
{
    slong k;
    for (k = 0; k < len; k++)
        acb_set_si(x + k,v[k]);
}
static void
acb_vec_set_pair(acb_ptr x, pair * v, slong len)
{
    slong k;
    for (k = 0; k < len; k++)
        acb_set_si_si(x + k, v[k][0], v[k][1]);
}
void
do_example(slong m, slong n, acb_srcptr x, entry * ref, slong size, int flag)
{
    slong len, prec = 64;
    sec_t c;
    tree_t tree;
    fmpz_mat_t ca;
    fmpz_mat_t cm;
    sec_init(&c, m, x, n);
    tree_init(tree,n-1);
    spanning_tree(tree, x, n, INT_DE);
    tree_ydata_init(tree, x, n, m, prec);
    len = (n-1)*(m-1);
    fmpz_mat_init(ca, len, len);
    intersection_tree(ca, tree, n, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_set_entries(cm, ref, size);

    if (!flag && !fmpz_mat_equal(ca,cm))
    {
        flint_printf("invalid intersection matrix\n");
        flint_printf("\nintersection matrix arb:\n");
        fmpz_mat_print_pretty(ca);
        flint_printf("\nintersection matrix magma:\n");
        fmpz_mat_print_pretty(cm);
        abort();
    }

    fmpz_mat_clear(ca);
    fmpz_mat_clear(cm);
    tree_ydata_clear(tree);
    tree_clear(tree);
    sec_clear(c);
}

void
do_example_si(slong m, slong n, slong * coeff, entry * ref, slong size, int flag)
{
    acb_ptr x;
    x = _acb_vec_init(n);
    acb_vec_set_si(x, coeff, n);
    do_example(m, n, x, ref, size, flag);
    _acb_vec_clear(x, n);
}

void
do_example_si_si(slong m, slong n, pair * coeff, entry * ref, slong size, int flag)
{
    acb_ptr x;
    x = _acb_vec_init(n);
    acb_vec_set_pair(x, coeff, n);
    do_example(m, n, x, ref, size, flag);
    _acb_vec_clear(x, n);
}

void
do_example_pol(slong m, slong n, slong * coeff, pair * swap, slong ns, entry * ref, slong size, int flag)
{
    slong k, prec = 64;
    acb_poly_t f;
    acb_ptr x;
    acb_poly_init(f);
    acb_poly_fit_length(f, n + 1);
    _acb_poly_set_length(f, n + 1);
    acb_one(f->coeffs + n);
    acb_vec_set_si(f->coeffs, coeff, n);
    x = _acb_vec_init(n);
    acb_poly_find_roots(x, f, NULL, 0, prec);
    for (k = 0; k < ns; k++)
        acb_swap(x+swap[k][0],x+swap[k][1]);
    do_example(m, n, x, ref, size, flag);
    _acb_vec_clear(x, n);
    acb_poly_clear(f);
}

int main()
{
    slong roots_1[3] = {-3, -1, 2};
    entry mat_1[2] = { {0,1,1}, {1,0,-1} };
    slong roots_2[3] = {-1, 0, 1};
    entry mat_2[16] = { {0,1,1}, {0,5,-1}, {1,0,-1}, {1,2,1}, {1,3,1},
        {2,1,-1}, {2,3,-1}, {2,4,1}, {3,1,-1}, {3,2,1}, {3,4,1}, {4,2,-1},
        {4,3,-1}, {4,5,1}, {5,0,1}, {5,4,-1} };
    pair roots_3[4] = {{-5,-1},{0,0},{0,-1},{0,2}};
    entry mat_3[30] = { {0,1,1}, {0,5,-1}, {0,6,1}, {0,7,-1}, {1,0,-1},
        {1,2,1}, {1,3,1}, {1,7,1}, {1,8,-1}, {2,1,-1}, {2,3,-1}, {2,4,1},
        {2,8,1}, {3,1,-1}, {3,2,1}, {3,4,1}, {4,2,-1}, {4,3,-1}, {4,5,1},
        {5,0,1}, {5,4,-1}, {6,0,-1}, {6,7,1}, {7,0,1}, {7,1,-1}, {7,6,-1},
        {7,8,1}, {8,1,1}, {8,2,-1},{8,7,-1} };
    slong pol_4[4] = {8, 92, -40, 84};
    pair swap_4[2] = { { 1, 2}, {2, 3} };
    entry mat_4[36] = { {0,1,1}, {0,4,-1}, {0,5,1}, {0,8,-1}, {1,0,-1},
        {1,2,1}, {1,5,-1}, {1,6,1}, {2,1,-1}, {2,3,1}, {2,6,-1}, {2,7,1},
        {3,2,-1}, {3,4,1}, {3,7,-1}, {3,8,1}, {4,0,1}, {4,3,-1},
        {4,5,1}, {4,8,-1}, {5,0,-1}, {5,1,1}, {5,4,-1}, {5,6,1}, {6,1,-1},
        {6,2,1}, {6,5,-1}, {6,7,1}, {7,2,-1}, {7,3,1}, {7,6,-1}, {7,8,1},
        {8,0,1}, {8,3,-1}, {8,4,1}, {8,7,-1} };
    slong pol_5[5] = {20, 13, 3, -3, -7};
    pair swap_5[4] = { { 0, 4}, {0, 3}, {0, 2}, {2, 1} };
    entry mat_5[82] = { {0,1,1}, {0,6,-1}, {0,7,1}, {1,0,-1}, {1,2,1},
        {1,7,-1}, {1,8,1}, {2,1,-1}, {2,3,1}, {2,8,-1}, {2,9,1}, {3,2,-1},
        {3,4,1}, {3,9,-1}, {4,3,-1}, {4,5,1}, {5,4,-1}, {5,6,1},
        {5,10,-1}, {6,0,1}, {6,5,-1}, {6,7,1},  {6,10,1}, {6,11,-1}, {7,0,-1},
        {7,1,1}, {7,6,-1}, {7,8,1}, {7,11,1},  {7,12,-1}, {8,1,-1}, {8,2,1},
        {8,7,-1}, {8,9,1}, {8,12,1}, {8,13,1}, {9,2,-1}, {9,3,1}, {9,8,-1},
        {9,13,1}, {9,14,-1}, {10,5,1}, {10,6,-1}, {10,11,1}, {10,16,1},
        {10,17,-1}, {11,6,1}, {11,7,-1}, {11,10,-1}, {11,12,1}, {11,17,1},
        {11,18,-1}, {12,7,1}, {12,8,-1}, {12,11,-1}, {12,13,1}, {12,18,1},
        {12,19,-1}, {13,8,1}, {13,9,-1}, {13,12,-1}, {13,14,1}, {13,19,1},
        {14,9,1}, {14,13,-1}, {14,15,-1}, {15,14,1}, {15,16,1}, {16,10,-1},
        {16,15,-1}, {16,17,1}, {17,10,1}, {17,11,-1}, {17,16,-1}, {17,18,1},
        {18,11,1}, {18,12,-1}, {18,17,-1}, {18,19,1}, {19,12,1}, {19,13,-1},
        {19,18,-1} };
	
    flint_printf("intersection matrix examples..."); fflush(stdout);

    /* Example 1 : y^2 = (x+3)(x+1)(x-2), g = 1 */
    /*do_example(2, 3, {-3, -1, 2}, mat_1, 2);*/
    progress("example 1\n");
    do_example_si(2, 3, roots_1, mat_1, 2, SKIP);

    /* Example 2 : y^4 = (x+1)(x)(x-1), g = 3 */
    /*do_example(4, 3, {-1, 0, 1}, mat_2, 16);*/
    progress("example 2\n");
    do_example_si(4, 3, roots_2, mat_2, 16, SKIP);

    /* Example 3 : y^4 = (x+5+i)(x)(x+i)(x-2i), g = 3 */
    /*do_example(4, 4, {-5-I, 0, -I, 2*I}, mat_3, 30);*/
    progress("example 3\n");
    do_example_si_si(4, 4, roots_3, mat_3, 30, SKIP);

    /* Example 4 : y^4 = x^4 + 84x^3 - 40x^2 + 92x + 8, g = 3 */
    progress("example 4\n");
    do_example_pol(4, 4, pol_4, swap_4, 2, mat_4, 36, SKIP);

    /* Example 5 (s1) : y^6 = x^5 - 7*x^4 - 3*x^3 + 3*x^2 + 13*x + 20, g = 10 */
    progress("example 5\n");
    do_example_pol(6, 5, pol_5, swap_5, 4, mat_5, 82, SKIP);

    flint_cleanup();
    printf("PASS\n");
    return 0;
}
