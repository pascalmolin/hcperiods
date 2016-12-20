/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include <complex.h>
#include <acb.h>

#define PI acos(-1.)
#define PI2 acos(0.)
#define LAMBDA PI2

typedef complex double cdouble;

cdouble acb_get_cdouble(const acb_t z);
