/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void tau_matrix(acb_mat_t tau, const acb_mat_t omega0, const acb_mat_t omega1, slong prec)
{
    acb_mat_inv(tau, omega1, prec);
    acb_mat_mul(tau, tau, omega0, prec);
}
