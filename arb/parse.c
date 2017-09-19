#include "parse.h"
#define PARSE_DBG 0

/* parse integer polynomial */

typedef struct {
    const char * str;
} parser_ctx_t;

int
accept_term(fmpz_poly_t poly, parser_ctx_t * ctx)
{
    slong c, s = 1, e = 0;
    char * end;

    /* parse coefficient */
    while (*ctx->str == ' ' || *ctx->str == '+' || *ctx->str == '-')
    {
        if (*ctx->str == '-')
            s = -1;
        ctx->str++;
    }

    c = strtol(ctx->str, &end, 10);
    if (ctx->str == end)
    {
        if (*ctx->str == 'x')
            c = 1;
        else
            return 0;
    }
    ctx->str = end;
    c *= s;

#if PARSE_DBG
    flint_printf(":  c=%i -> '%s'\n", c, ctx->str);
#endif

    while (*ctx->str == ' ' || *ctx->str == '*')
        ctx->str++;

#if PARSE_DBG
    flint_printf(":  [ ] -> '%s'\n", ctx->str);
#endif

    /* variable and power */

    if (*ctx->str == 'x')
    {
        e = 1;
        ctx->str++;

#if PARSE_DBG
        flint_printf(":  x -> '%s'\n", ctx->str);
#endif

        if (*ctx->str == '^') {
            ctx->str++;
#if PARSE_DBG
            flint_printf(":  ^ -> '%s'\n", ctx->str);
#endif

            e = strtol(ctx->str, &end, 10);
            if (ctx->str == end)
                e = 1;
            ctx->str = end;
#if PARSE_DBG
            flint_printf(":  e=%i -> '%s'\n", e, ctx->str);
#endif
        }
    }

    fmpz_poly_set_coeff_si(poly, e, c);
    return 1;
}
    
int
fmpz_poly_parse(fmpz_poly_t poly, const char* str)
{
    parser_ctx_t ctx;
    ctx.str = str;

    while (*ctx.str)
    {
        if (!accept_term(poly, &ctx))
            return 0;
#if PARSE_DBG
        fmpz_poly_print_pretty(poly, "x");
        flint_printf(" -> parse %s\n",ctx.str);
#endif
    }
    return 1;
}


