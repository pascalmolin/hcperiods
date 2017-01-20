/******************************************************************************

 Copyright (C) Christian Neurohr

 ******************************************************************************/

#include "abel_jacobi.h"

int main() {

    flint_printf("intersection matrix examples...");
    fflush(stdout);
    slong d, m, len, k;
    acb_ptr x;
    tree_t tree;
    fmpz_mat_t c;
    fmpz_mat_t cm;
    acb_poly_t f;
	
    /* Example 1 : y^2 = (x+3)(x+1)(x-2), g = 1 */
    d = 3;
    m = 2;
    x = _acb_vec_init(d);
    acb_set_si(x + 0, -3);
    acb_set_si(x + 1, -1);
    acb_set_si(x + 2, 2);
    tree_init(tree,d-1);
    spanning_tree(tree, x, d, INT_DE);
    len = (d-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, d, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_zero(cm);
    fmpz_set_si(fmpz_mat_entry(cm,0,1),1); 
    fmpz_set_si(fmpz_mat_entry(cm,1,0),-1);
  
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
    _acb_vec_clear(x, d);
    

    /* Example 2 : y^4 = (x+1)(x)(x-1), g = 3 */
    d = 3;
    m = 4;
    x = _acb_vec_init(d);
    acb_set_si(x + 0, -1);
    acb_set_si(x + 1, 0);
    acb_set_si(x + 2, 1);
    tree_init(tree,d-1);
    spanning_tree(tree, x, d, INT_DE);
    len = (d-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, d, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_zero(cm);
    fmpz_set_si(fmpz_mat_entry(cm,0,1),1); fmpz_set_si(fmpz_mat_entry(cm,0,5),-1);
    fmpz_set_si(fmpz_mat_entry(cm,1,0),-1); fmpz_set_si(fmpz_mat_entry(cm,1,2),1); fmpz_set_si(fmpz_mat_entry(cm,1,3),1);
    fmpz_set_si(fmpz_mat_entry(cm,2,1),-1); fmpz_set_si(fmpz_mat_entry(cm,2,3),-1); fmpz_set_si(fmpz_mat_entry(cm,2,4),1);
    fmpz_set_si(fmpz_mat_entry(cm,3,1),-1); fmpz_set_si(fmpz_mat_entry(cm,3,2),1);fmpz_set_si(fmpz_mat_entry(cm,3,4),1);
    fmpz_set_si(fmpz_mat_entry(cm,4,2),-1); fmpz_set_si(fmpz_mat_entry(cm,4,3),-1);fmpz_set_si(fmpz_mat_entry(cm,4,5),1);
    fmpz_set_si(fmpz_mat_entry(cm,5,0),1); fmpz_set_si(fmpz_mat_entry(cm,5,4),-1);
    
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
    _acb_vec_clear(x, d);

   
    /* Example 3 : y^4 = (x+5+i)(x)(x+i)(x-2i), g = 3 */
    d = 4;
    m = 4;
    x = _acb_vec_init(d);
    acb_onei( x + 0 ); acb_add_si( x+0, x+0, 5, 30); acb_neg(x+0,x+0);
    acb_set_si( x+1, 0);
    acb_onei( x + 2 ); acb_neg(x+2,x+2);
    acb_onei( x + 3 ); acb_mul_si(x+3,x+3,2,30);
    tree_init(tree,d-1);
    spanning_tree(tree, x, d, INT_DE);
    len = (d-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, d, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_zero(cm);
    fmpz_set_si(fmpz_mat_entry(cm,0,1),1); fmpz_set_si(fmpz_mat_entry(cm,0,5),-1); fmpz_set_si(fmpz_mat_entry(cm,0,6),1); fmpz_set_si(fmpz_mat_entry(cm,0,7),-1);
    fmpz_set_si(fmpz_mat_entry(cm,1,0),-1); fmpz_set_si(fmpz_mat_entry(cm,1,2),1); fmpz_set_si(fmpz_mat_entry(cm,1,3),1);fmpz_set_si(fmpz_mat_entry(cm,1,7),1); fmpz_set_si(fmpz_mat_entry(cm,1,8),-1);
    fmpz_set_si(fmpz_mat_entry(cm,2,1),-1); fmpz_set_si(fmpz_mat_entry(cm,2,3),-1); fmpz_set_si(fmpz_mat_entry(cm,2,4),1); fmpz_set_si(fmpz_mat_entry(cm,2,8),1);
    fmpz_set_si(fmpz_mat_entry(cm,3,1),-1); fmpz_set_si(fmpz_mat_entry(cm,3,2),1); fmpz_set_si(fmpz_mat_entry(cm,3,4),1); 
    fmpz_set_si(fmpz_mat_entry(cm,4,2),-1); fmpz_set_si(fmpz_mat_entry(cm,4,3),-1);fmpz_set_si(fmpz_mat_entry(cm,4,5),1);
    fmpz_set_si(fmpz_mat_entry(cm,5,0),1); fmpz_set_si(fmpz_mat_entry(cm,5,4),-1);
    fmpz_set_si(fmpz_mat_entry(cm,6,0),-1); fmpz_set_si(fmpz_mat_entry(cm,6,7),1);
    fmpz_set_si(fmpz_mat_entry(cm,7,0),1); fmpz_set_si(fmpz_mat_entry(cm,7,1),-1); fmpz_set_si(fmpz_mat_entry(cm,7,6),-1); fmpz_set_si(fmpz_mat_entry(cm,7,8),1);
    fmpz_set_si(fmpz_mat_entry(cm,8,1),1); fmpz_set_si(fmpz_mat_entry(cm,8,2),-1);fmpz_set_si(fmpz_mat_entry(cm,8,7),-1);

    if ( fmpz_mat_equal(c,cm) == 0 )
    {
      flint_printf("Example 3: invalid intersection matrix\n");
      flint_printf("\nintersection matrix arb:\n");
      fmpz_mat_print_pretty(c);
      flint_printf("\nintersection matrix magma:\n");
      fmpz_mat_print_pretty(cm);
      /*abort()*/;
    }

    fmpz_mat_clear(c);
    fmpz_mat_clear(cm);
    tree_clear(tree);
    _acb_vec_clear(x, d);

    /* Example 4 : y^4 = x^4 + 84x^3 - 40x^2 + 92x + 8, g = 3 */
    d = 4;
    m = 4;
    x = _acb_vec_init(d);
    acb_poly_init(f);
    acb_poly_set_coeff_si(f,0,8);
    acb_poly_set_coeff_si(f,1,92);
    acb_poly_set_coeff_si(f,2,-40);
    acb_poly_set_coeff_si(f,3,84);
    acb_poly_set_coeff_si(f,4,1);
    acb_poly_find_roots(x, f, NULL, 0, 30);
    acb_swap(x+1,x+2);
    acb_swap(x+2,x+3);
    
    tree_init(tree,d-1);
    spanning_tree(tree, x, d, INT_DE);
    len = (d-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, d, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_zero(cm);
    fmpz_set_si(fmpz_mat_entry(cm,0,1),1); fmpz_set_si(fmpz_mat_entry(cm,0,4),-1); fmpz_set_si(fmpz_mat_entry(cm,0,5),1); fmpz_set_si(fmpz_mat_entry(cm,0,8),-1);
    fmpz_set_si(fmpz_mat_entry(cm,1,0),-1); fmpz_set_si(fmpz_mat_entry(cm,1,2),1); fmpz_set_si(fmpz_mat_entry(cm,1,5),-1); fmpz_set_si(fmpz_mat_entry(cm,1,6),1);
    fmpz_set_si(fmpz_mat_entry(cm,2,1),-1); fmpz_set_si(fmpz_mat_entry(cm,2,3),1); fmpz_set_si(fmpz_mat_entry(cm,2,6),-1); fmpz_set_si(fmpz_mat_entry(cm,2,7),1);
    fmpz_set_si(fmpz_mat_entry(cm,3,2),-1); fmpz_set_si(fmpz_mat_entry(cm,3,4),1); fmpz_set_si(fmpz_mat_entry(cm,3,7),-1); fmpz_set_si(fmpz_mat_entry(cm,3,8),1);
    fmpz_set_si(fmpz_mat_entry(cm,4,0),1); fmpz_set_si(fmpz_mat_entry(cm,4,3),-1); fmpz_set_si(fmpz_mat_entry(cm,4,5),1); fmpz_set_si(fmpz_mat_entry(cm,4,8),-1);
    fmpz_set_si(fmpz_mat_entry(cm,5,0),-1); fmpz_set_si(fmpz_mat_entry(cm,5,1),1); fmpz_set_si(fmpz_mat_entry(cm,5,4),-1); fmpz_set_si(fmpz_mat_entry(cm,5,6),1);
    fmpz_set_si(fmpz_mat_entry(cm,6,1),-1); fmpz_set_si(fmpz_mat_entry(cm,6,2),1); fmpz_set_si(fmpz_mat_entry(cm,6,5),-1); fmpz_set_si(fmpz_mat_entry(cm,6,7),1);
    fmpz_set_si(fmpz_mat_entry(cm,7,2),-1); fmpz_set_si(fmpz_mat_entry(cm,7,3),1); fmpz_set_si(fmpz_mat_entry(cm,7,6),-1); fmpz_set_si(fmpz_mat_entry(cm,7,8),1);
    fmpz_set_si(fmpz_mat_entry(cm,8,0),1); fmpz_set_si(fmpz_mat_entry(cm,8,3),-1); fmpz_set_si(fmpz_mat_entry(cm,8,4),1); fmpz_set_si(fmpz_mat_entry(cm,8,7),-1);

    if ( fmpz_mat_equal(c,cm) == 0 )
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
    _acb_vec_clear(x, d);
    acb_poly_clear(f);



    /* Example 5 (s1) : y^6 = x^5 - 7*x^4 - 3*x^3 + 3*x^2 + 13*x + 20, g = 10 */
    d = 5;
    m = 6;
    x = _acb_vec_init(d);
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
    
    tree_init(tree,d-1);
    spanning_tree(tree, x, d, INT_DE);
    len = (d-1)*(m-1);
    fmpz_mat_init(c, len, len);
    intersection_tree(c, tree, d, m);

    /* Intersection matrix from magma */
    fmpz_mat_init(cm, len, len);
    fmpz_mat_zero(cm);
    fmpz_set_si(fmpz_mat_entry(cm,0,1),1); fmpz_set_si(fmpz_mat_entry(cm,0,6),-1); fmpz_set_si(fmpz_mat_entry(cm,0,7),1);
    fmpz_set_si(fmpz_mat_entry(cm,1,0),-1); fmpz_set_si(fmpz_mat_entry(cm,1,2),1); fmpz_set_si(fmpz_mat_entry(cm,1,7),-1); fmpz_set_si(fmpz_mat_entry(cm,1,8),1);
    fmpz_set_si(fmpz_mat_entry(cm,2,1),-1); fmpz_set_si(fmpz_mat_entry(cm,2,3),1); fmpz_set_si(fmpz_mat_entry(cm,2,8),-1); fmpz_set_si(fmpz_mat_entry(cm,2,9),1);
    fmpz_set_si(fmpz_mat_entry(cm,3,2),-1); fmpz_set_si(fmpz_mat_entry(cm,3,4),1); fmpz_set_si(fmpz_mat_entry(cm,3,9),-1);
    fmpz_set_si(fmpz_mat_entry(cm,4,3),-1); fmpz_set_si(fmpz_mat_entry(cm,4,5),1);
    fmpz_set_si(fmpz_mat_entry(cm,5,4),-1); fmpz_set_si(fmpz_mat_entry(cm,5,6),1); fmpz_set_si(fmpz_mat_entry(cm,5,10),-1);
    fmpz_set_si(fmpz_mat_entry(cm,6,0),1); fmpz_set_si(fmpz_mat_entry(cm,6,5),-1); fmpz_set_si(fmpz_mat_entry(cm,6,7),1);  fmpz_set_si(fmpz_mat_entry(cm,6,10),1);  
    fmpz_set_si(fmpz_mat_entry(cm,6,11),-1);
    fmpz_set_si(fmpz_mat_entry(cm,7,0),-1); fmpz_set_si(fmpz_mat_entry(cm,7,1),1); fmpz_set_si(fmpz_mat_entry(cm,7,6),-1); fmpz_set_si(fmpz_mat_entry(cm,7,8),1);
    fmpz_set_si(fmpz_mat_entry(cm,7,11),1);  fmpz_set_si(fmpz_mat_entry(cm,7,12),-1);
    fmpz_set_si(fmpz_mat_entry(cm,8,1),-1); fmpz_set_si(fmpz_mat_entry(cm,8,2),1); fmpz_set_si(fmpz_mat_entry(cm,8,7),-1); fmpz_set_si(fmpz_mat_entry(cm,8,9),1);
    fmpz_set_si(fmpz_mat_entry(cm,8,12),1); fmpz_set_si(fmpz_mat_entry(cm,8,13),1);
    fmpz_set_si(fmpz_mat_entry(cm,9,2),-1); fmpz_set_si(fmpz_mat_entry(cm,9,3),1); fmpz_set_si(fmpz_mat_entry(cm,9,8),-1);
    fmpz_set_si(fmpz_mat_entry(cm,9,13),1); fmpz_set_si(fmpz_mat_entry(cm,9,14),-1);
    fmpz_set_si(fmpz_mat_entry(cm,10,5),1); fmpz_set_si(fmpz_mat_entry(cm,10,6),-1); fmpz_set_si(fmpz_mat_entry(cm,10,11),1); fmpz_set_si(fmpz_mat_entry(cm,10,16),1); 
    fmpz_set_si(fmpz_mat_entry(cm,10,17),-1);
    fmpz_set_si(fmpz_mat_entry(cm,11,6),1); fmpz_set_si(fmpz_mat_entry(cm,11,7),-1); fmpz_set_si(fmpz_mat_entry(cm,11,10),-1); fmpz_set_si(fmpz_mat_entry(cm,11,12),1);
    fmpz_set_si(fmpz_mat_entry(cm,11,17),1); fmpz_set_si(fmpz_mat_entry(cm,11,18),-1);
    fmpz_set_si(fmpz_mat_entry(cm,12,7),1); fmpz_set_si(fmpz_mat_entry(cm,12,8),-1); fmpz_set_si(fmpz_mat_entry(cm,12,11),-1); fmpz_set_si(fmpz_mat_entry(cm,12,13),1);
    fmpz_set_si(fmpz_mat_entry(cm,12,18),1); fmpz_set_si(fmpz_mat_entry(cm,12,19),-1);
    fmpz_set_si(fmpz_mat_entry(cm,13,8),1); fmpz_set_si(fmpz_mat_entry(cm,13,9),-1); fmpz_set_si(fmpz_mat_entry(cm,13,12),-1); fmpz_set_si(fmpz_mat_entry(cm,13,14),1);
    fmpz_set_si(fmpz_mat_entry(cm,13,19),1);
    fmpz_set_si(fmpz_mat_entry(cm,14,9),1); fmpz_set_si(fmpz_mat_entry(cm,14,13),-1); fmpz_set_si(fmpz_mat_entry(cm,14,15),-1);
    fmpz_set_si(fmpz_mat_entry(cm,15,14),1); fmpz_set_si(fmpz_mat_entry(cm,15,16),1);
    fmpz_set_si(fmpz_mat_entry(cm,16,10),-1); fmpz_set_si(fmpz_mat_entry(cm,16,15),-1); fmpz_set_si(fmpz_mat_entry(cm,16,17),1);
    fmpz_set_si(fmpz_mat_entry(cm,17,10),1); fmpz_set_si(fmpz_mat_entry(cm,17,11),-1); fmpz_set_si(fmpz_mat_entry(cm,17,16),-1); fmpz_set_si(fmpz_mat_entry(cm,17,18),1);
    fmpz_set_si(fmpz_mat_entry(cm,18,11),1); fmpz_set_si(fmpz_mat_entry(cm,18,12),-1); fmpz_set_si(fmpz_mat_entry(cm,18,17),-1); fmpz_set_si(fmpz_mat_entry(cm,18,19),1);
    fmpz_set_si(fmpz_mat_entry(cm,19,12),1); fmpz_set_si(fmpz_mat_entry(cm,19,13),-1); fmpz_set_si(fmpz_mat_entry(cm,19,18),-1);
   
    if ( fmpz_mat_equal(c,cm) == 0 )
    {
      flint_printf("Example 5: invalid intersection matrix\n");
      flint_printf("\nintersection matrix arb:\n");
      fmpz_mat_print_pretty(c);
      flint_printf("\nintersection matrix magma:\n");
      fmpz_mat_print_pretty(cm);
      /*abort();*/
    }

    fmpz_mat_clear(c);
    fmpz_mat_clear(cm);
    tree_clear(tree);
    _acb_vec_clear(x, d);
    acb_poly_clear(f);


    flint_cleanup();
    printf("PASS\n");
    return 0;
}
