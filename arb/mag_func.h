/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include <arb.h>
#define VERBOSE 0

typedef int (* arb_func_t)(arb_t abs, const arb_t, void * params, slong prec);
void arb_bound_func_arf(arb_t m, arb_func_t f, void * params, const arf_t tmin, const arf_t tmax, slong n, slong prec);
void arb_bound_func_arb(arb_t m, arb_func_t f, void * params, const arb_t b, slong n, slong prec);
slong mag_func_arf(mag_t m, arb_func_t f, void * params, const arf_t tmin, const arf_t tmax, slong n, slong prec);
slong mag_func_arb(mag_t m, arb_func_t f, void * params, const arb_t b, slong n, slong prec);
