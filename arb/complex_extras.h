/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include <complex.h>
#include <acb.h>

#define PI acos(-1.)
#define PI2 acos(0.)
#define LAMBDA PI2
#define LOG2 log(2)

typedef complex double cdouble;

inline cdouble
acb_get_cdouble(const acb_t z)
{
    return arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR)
    + _Complex_I * arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
}
