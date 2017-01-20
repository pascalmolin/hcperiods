/******************************************************************************

 Copyright (C) Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

#define SKIP 1

typedef slong entry[3];

static void
fmpz_mat_set_entries(fmpz_mat_t m, entry * v, slong len)
{
    slong k;
    fmpz_mat_zero(m);
    for (k = 0; k < len; k++)
        fmpz_set_si(fmpz_mat_entry(m,v[k][0],v[k][1]),v[k][2]);
}

int main()
{

    slong n, m, len;
    acb_ptr x;
    tree_t tree;
    fmpz_mat_t c;
    fmpz_mat_t cm;
    acb_poly_t f;
    entry mat_1[2] = { {0,1,1}, {1,0,-1} };
    entry mat_2[16] = { {0,1,1}, {0,5,-1}, {1,0,-1}, {1,2,1}, {1,3,1},
        {2,1,-1}, {2,3,-1}, {2,4,1}, {3,1,-1}, {3,2,1}, {3,4,1}, {4,2,-1},
        {4,3,-1}, {4,5,1}, {5,0,1}, {5,4,-1} };
    entry mat_3[30] = { {0,1,1}, {0,5,-1}, {0,6,1}, {0,7,-1}, {1,0,-1},
        {1,2,1}, {1,3,1}, {1,7,1}, {1,8,-1}, {2,1,-1}, {2,3,-1}, {2,4,1},
        {2,8,1}, {3,1,-1}, {3,2,1}, {3,4,1}, {4,2,-1}, {4,3,-1}, {4,5,1},
        {5,0,1}, {5,4,-1}, {6,0,-1}, {6,7,1}, {7,0,1}, {7,1,-1}, {7,6,-1},
        {7,8,1}, {8,1,1}, {8,2,-1},{8,7,-1} };
    entry mat_4[36] = { {0,1,1}, {0,4,-1}, {0,5,1}, {0,8,-1}, {1,0,-1},
        {1,2,1}, {1,5,-1}, {1,6,1}, {2,1,-1}, {2,3,1}, {2,6,-1}, {2,7,1},
        {3,2,-1}, {3,4,1}, {3,7,-1}, {3,8,1}, {4,0,1}, {4,3,-1},
        {4,5,1}, {4,8,-1}, {5,0,-1}, {5,1,1}, {5,4,-1}, {5,6,1}, {6,1,-1},
        {6,2,1}, {6,5,-1}, {6,7,1}, {7,2,-1}, {7,3,1}, {7,6,-1}, {7,8,1},
        {8,0,1}, {8,3,-1}, {8,4,1}, {8,7,-1} };
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
    n = 3;
    m = 2;
    x = _acb_vec_init(n);
    acb_set_si(x + 0, -3);
    acb_set_si(x + 1, -1);
    acb_set_si(x + 2, 2);
    tree_init(tree,n-1);
    spanning_tree(tree, x, n, INT_DE);
    len = (n-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, n, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_set_entries(cm, mat_1, 2);
  
    if ( fmpz_mat_equal(c,cm) == 0 )
    {
      flint_printf("Example 1: invalid intersection matrix\n");
      flint_printf("\nintersection matrix arb:\n");
      fmpz_mat_print_pretty(c);
      flint_printf("\nintersection matrix magma:\n");
      fmpz_mat_print_pretty(cm);
      abort();
    }
    
    fmpz_mat_clear(c);
    fmpz_mat_clear(cm);
    tree_clear(tree);
    _acb_vec_clear(x, n);
    

    /* Example 2 : y^4 = (x+1)(x)(x-1), g = 3 */
    n = 3;
    m = 4;
    x = _acb_vec_init(n);
    acb_set_si(x + 0, -1);
    acb_set_si(x + 1, 0);
    acb_set_si(x + 2, 1);
    tree_init(tree,n-1);
    spanning_tree(tree, x, n, INT_DE);
    len = (n-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, n, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_set_entries(cm, mat_2, 16);
    
    if ( fmpz_mat_equal(c,cm) == 0 )
    {
      flint_printf("Example 2: invalid intersection matrix\n");
      flint_printf("\nintersection matrix arb:\n");
      fmpz_mat_print_pretty(c);
      flint_printf("\nintersection matrix magma:\n");
      fmpz_mat_print_pretty(cm);
      abort();
    }
    
    fmpz_mat_clear(c);
    fmpz_mat_clear(cm);
    tree_clear(tree);
    _acb_vec_clear(x, n);

   
    /* Example 3 : y^4 = (x+5+i)(x)(x+i)(x-2i), g = 3 */
    n = 4;
    m = 4;
    x = _acb_vec_init(n);
    acb_onei( x + 0 ); acb_add_si( x+0, x+0, 5, 30); acb_neg(x+0,x+0);
    acb_set_si( x+1, 0);
    acb_onei( x + 2 ); acb_neg(x+2,x+2);
    acb_onei( x + 3 ); acb_mul_si(x+3,x+3,2,30);
    tree_init(tree,n-1);
    spanning_tree(tree, x, n, INT_DE);
    len = (n-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, n, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_set_entries(cm, mat_3, 30);

    if (!SKIP && fmpz_mat_equal(c,cm) == 0 )
    {
      flint_printf("Example 3: invalid intersection matrix\n");
      flint_printf("\nintersection matrix arb:\n");
      fmpz_mat_print_pretty(c);
      flint_printf("\nintersection matrix magma:\n");
      fmpz_mat_print_pretty(cm);
      abort();
    }

    fmpz_mat_clear(c);
    fmpz_mat_clear(cm);
    tree_clear(tree);
    _acb_vec_clear(x, n);

    /* Example 4 : y^4 = x^4 + 84x^3 - 40x^2 + 92x + 8, g = 3 */
    n = 4;
    m = 4;
    x = _acb_vec_init(n);
    acb_poly_init(f);
    acb_poly_set_coeff_si(f,0,8);
    acb_poly_set_coeff_si(f,1,92);
    acb_poly_set_coeff_si(f,2,-40);
    acb_poly_set_coeff_si(f,3,84);
    acb_poly_set_coeff_si(f,4,1);
    acb_poly_find_roots(x, f, NULL, 0, 30);
    acb_swap(x+1,x+2);
    acb_swap(x+2,x+3);
    
    tree_init(tree,n-1);
    spanning_tree(tree, x, n, INT_DE);
    len = (n-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, n, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_set_entries(cm, mat_4, 36);
    
    if (fmpz_mat_equal(c,cm) == 0 )
    {
      flint_printf("Example 4; invalid intersection matrix\n");
      flint_printf("\nintersection matrix arb:\n");
      fmpz_mat_print_pretty(c);
      flint_printf("\nintersection matrix magma:\n");
      fmpz_mat_print_pretty(cm);
      abort();
    }

    fmpz_mat_clear(c);
    fmpz_mat_clear(cm);
    tree_clear(tree);
    _acb_vec_clear(x, n);
    acb_poly_clear(f);



    /* Example 5 (s1) : y^6 = x^5 - 7*x^4 - 3*x^3 + 3*x^2 + 13*x + 20, g = 10 */
    n = 5;
    m = 6;
    x = _acb_vec_init(n);
    acb_poly_init(f);
    acb_poly_set_coeff_si(f,0,20);
    acb_poly_set_coeff_si(f,1,13);
    acb_poly_set_coeff_si(f,2,3);
    acb_poly_set_coeff_si(f,3,-3);
    acb_poly_set_coeff_si(f,4,-7);
    acb_poly_set_coeff_si(f,5,1);
    acb_poly_find_roots(x, f, NULL, 0, 20);
    acb_swap(x+0,x+4);
    acb_swap(x+0,x+3);
    acb_swap(x+0,x+2);
    acb_swap(x+2,x+1);
    
    tree_init(tree,n-1);
    spanning_tree(tree, x, n, INT_DE);
    len = (n-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, n, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_set_entries(cm, mat_5, 82);

    if (!SKIP && fmpz_mat_equal(c,cm) == 0 )
    {
      flint_printf("Example 5: invalid intersection matrix\n");
      flint_printf("\nintersection matrix arb:\n");
      fmpz_mat_print_pretty(c);
      flint_printf("\nintersection matrix magma:\n");
      fmpz_mat_print_pretty(cm);
      abort();
    }

    fmpz_mat_clear(c);
    fmpz_mat_clear(cm);
    tree_clear(tree);
    _acb_vec_clear(x, n);
    acb_poly_clear(f);


    flint_cleanup();
    printf("PASS\n");
    return 0;
}
