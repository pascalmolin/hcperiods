#include <acb.h>
#include "mag_func.h"

typedef struct
{
    slong len;
    acb_ptr z;
} params_t;

int
f_1x2(arb_t max, const arb_t t, params_t * p, slong prec)
{
    arb_mul(max, t, t, prec);
    arb_add_si(max, max, 1, prec);
    arb_inv(max, max, prec);
    return 1;
}

int
f_pol(arb_t max, const arb_t t, params_t * p, slong prec)
{
    slong k;
    acb_t z;
    arb_t tmp;

    acb_init(z);
    arb_init(tmp);

    arb_one(max);

    for (k = 0; k < p->len; k++)
    {
        acb_sub_arb(z, p->z + k, t, prec);
        acb_abs(tmp, z, prec);
        if (arb_contains_zero(tmp))
            return 0;
        arb_mul(max, max, tmp, prec);
    }

    arb_inv(max, max, prec);

    arb_clear(tmp);
    acb_clear(z);
    return 1;
}

int
f_thsh(arb_t max, const arb_t t, params_t * p, slong prec)
{
    arb_sinh(max, t, prec);
    arb_tanh(max, max, prec);
    arb_abs(max, max);
    return 1;
}

int
f_thsh_shift(arb_t max, const arb_t t, params_t * p, slong prec)
{
    acb_t z;
    acb_init(z);
    arb_set(acb_realref(z), t);
    arb_set_d(acb_imagref(z), .7);

    acb_sinh(z, z, prec);
    acb_tanh(z, z, prec);
    acb_abs(max, z, prec);

    acb_clear(z);
    return 1;
}

int
f_aj(arb_t m, const arb_t t, params_t * p, slong prec)
{
    slong k;
    acb_t z, zu;
    arb_t abs;

    arb_init(abs);
    acb_init(z);
    acb_init(zu);

    arb_const_pi(abs, prec);
    arb_mul_2exp_si(abs, abs, -2); /* Pi/4 */

    arb_set(acb_realref(z), t);
    arb_set(acb_imagref(z), abs);

    acb_sinh(z, z, prec);
    arb_mul_2exp_si(abs, abs, 1); /* Pi/2 */
    acb_mul_arb(z, z, abs, prec);
    acb_tanh(z, z, prec);

    arb_one(m);
    for (k = 0; k < p->len; k++)
    {
        acb_sub(zu, z, p->z + k, prec);
        if (acb_contains_zero(zu))
        {
            arb_clear(abs);
            acb_clear(zu);
            acb_clear(z);
            return 0;
        }
        acb_abs(abs, zu, prec);
        arb_mul(m, m, abs, prec);
    }
    arb_inv(m, m, prec);

    arb_clear(abs);
    acb_clear(zu);
    acb_clear(z);
    return 1;
}

int main() {

#define ni 5
    slong n, i, f;
    double b[ni][2] = { { 0, 1}, {-1, 1}, {0, 10}, {-20, 3}, {0, 3.14} };
    const slong prec = 40;
#define nf 5
    arb_func_t func[nf] = { (arb_func_t)&f_1x2, (arb_func_t)&f_pol, (arb_func_t)&f_thsh, (arb_func_t)&f_thsh_shift, (arb_func_t)&f_aj };
#if VERBOSE
    char * fn[nf] = { "1/(1+x^2)", "1/p(x)", "th(sh(x))", "th(sh(x+.7I))" , "1/y" };
#endif
    params_t p[nf];
    
    arf_t tmin, tmax;
    mag_t max;

    flint_printf("max_func...");
    fflush(stdout);

    p[1].len = 3;
    p[1].z = _acb_vec_init(3);
    acb_set_d_d(p[1].z + 0, 2, 1);
    acb_set_d_d(p[1].z + 1, 2, .1);
    acb_set_d_d(p[1].z + 2, 1, .1);

    p[4] = p[1];

    arf_init(tmin);
    arf_init(tmax);
    mag_init(max);

    for (i = 0; i < ni; i++)
    {
        arf_set_d(tmin, b[i][0]);
        arf_set_d(tmax, b[i][1]);

        for (f = 0; f < nf; f++)
        { 
            for (n = 5; n < 100; n *= 2)
            {
                slong count;
                count = mag_func_arf(max, func[f], (void *)&p[f], tmin, tmax, n, prec);
#if VERBOSE
                flint_printf("\nmax %s on [%lf, %lf] <= ",fn[f],b[i][0],b[i][1]);
                mag_printd(max,8);
                flint_printf(" [asked %ld, did %ld]", n, count);
#endif
            }
        }

    }

    mag_clear(max);
    arf_clear(tmin);
    arf_clear(tmax);

    printf("PASS\n");
    return 0;
}
