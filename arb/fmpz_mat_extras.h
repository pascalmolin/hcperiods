/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include <fmpz.h>
#include <fmpz_mat.h>

void row_swap(fmpz_mat_t m, long i1, long i2);
void col_swap(fmpz_mat_t m, long j1, long j2);

void col_neg(fmpz_mat_t m, long j);
void row_neg(fmpz_mat_t m, long i);

void row_addmul(fmpz_mat_t m, long i, long k, fmpz_t v);
void col_addmul(fmpz_mat_t m, long j, long k, fmpz_t v);
void row_bezout(fmpz_mat_t m, long i1, long i2, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d);
void col_bezout(fmpz_mat_t m, long j1, long j2, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d);

int is_symplectic_j(fmpz_mat_t m, long g, long g2);

/*void fmpz_mat_print_gp(fmpz_mat_t m, long nr, long nc);*/
