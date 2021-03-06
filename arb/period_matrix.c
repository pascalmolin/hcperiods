/******************************************************************************

 Copyright (C) 2016 Pascal Molin

 ******************************************************************************/

#include "abel_jacobi.h"

void
periods_loop(acb_ptr res, const loop_t loop, const acb_mat_t integrals, acb_srcptr z, sec_t c, slong prec)
{
    slong i;
    acb_t tmp;

    acb_init(tmp);
    _acb_vec_zero(res, c.g);

    for (i = 0; i < loop.n; i++)
    {
        slong e = loop.l[i].index, l = loop.l[i].shift;
        fmpz * coeff = &loop.l[i].coeff;
        acb_ptr r = res, ii = integrals->rows[e];

        /* coeff * (int_e^l) */
        if (l == 0)
            _acb_vec_scalar_addmul_fmpz(r, ii, c.g, coeff, prec);
        else
        {
            slong j;
            for (j = 0; j < c.nj; j++)
            {
                slong lj;
                /* z^(-lj) * coeff */
                lj = ((c.m-l) * (c.j1 + j)) % c.m;
                acb_mul_fmpz(tmp, z + lj, coeff, prec);
                /* add on value j */
                _acb_vec_scalar_addmul(r, ii, c.ni[j], tmp, prec);
                r += c.ni[j];
                ii += c.ni[j];
            }
        }
    }

    acb_clear(tmp);
}

void
period_matrix(acb_mat_t omega, const homol_t basis, const acb_mat_t integrals, sec_t c, slong prec)
{
    slong k;
    acb_ptr z;
    z = _acb_vec_init(c.m);
    _acb_vec_unit_roots(z, c.m, c.m, prec);

    for (k = 0; k < c.g; k++)
        periods_loop(omega->rows[k], basis[k], integrals, z, c, prec);
    acb_mat_transpose(omega, omega);

    _acb_vec_clear(z, c.m);
    return;
}
