/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
periods_loop(acb_ptr res, const loop_t loop, const cohom_t dz, const acb_mat_t integrals, acb_srcptr z, sec_t c, slong prec)
{
    slong i, k;
    acb_t tmp;

    acb_init(tmp);

    for (k = 0; k < c.g; k++)
        acb_zero(res + k);

    for (i = 0; i < loop.n; i++)
    {
        slong e = loop.l[i].index, l = loop.l[i].shift;
        fmpz * coeff = &loop.l[i].coeff;

        if (l == 0)
        {
            for (k = 0; k < c.g; k++)
                acb_addmul_fmpz(res + k, acb_mat_entry(integrals, e, k), coeff, prec); 
        }
        else
        {
            for (k = 0; k < c.g; k++)
            {
                slong lj;
                lj = (l * dz[k].y) % c.m;
                acb_mul_fmpz(tmp, z + lj, coeff, prec);
                acb_addmul(res + k, acb_mat_entry(integrals, e, k), tmp, prec); 
            }
        }
    }

    acb_clear(tmp);
}

void
period_matrix(acb_mat_t omega, const homol_t basis, const acb_mat_t integrals, sec_t c, slong prec)
{
    return;
}
