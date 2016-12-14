/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include <fmpz_mat.h>
/* integer matrix */
typedef fmpz_mat_t si_mat_t;

void si_mat_init(si_mat_t m, long nr, long nc);
void si_mat_clear(si_mat_t m, long nr, long nc);

void row_swap(si_mat_t m, long i1, long i2, long len);
void col_swap(si_mat_t m, long j1, long j2, long len);

void col_neg(si_mat_t m, long j, long len);
void row_neg(si_mat_t m, long i, long len);

void row_addmul(si_mat_t m, long i, long k, fmpz_t v, long len);
void col_addmul(si_mat_t m, long j, long k, fmpz_t v, long len);
void row_bezout(si_mat_t m, long i1, long i2, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, long len);
void col_bezout(si_mat_t m, long j1, long j2, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, long len);

int is_symplectic_j(si_mat_t m, long g, long g2);
void si_mat_set_id(si_mat_t p, long len);

void si_mat_print(si_mat_t m, long nr, long nc);
void si_mat_print_gp(si_mat_t m, long nr, long nc);
