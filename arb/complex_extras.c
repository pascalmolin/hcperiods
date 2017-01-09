/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "complex_extras.h"

cdouble
acb_get_cdouble(const acb_t z)
{
    return
        arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR)
    + _Complex_I * arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
}
