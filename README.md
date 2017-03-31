# hcperiods

Numerical evaluation of periods of superelliptic curves to
arbitrary precision by numerical integration.

## Structure

- `gp`

  initial projet as gp script, period matrix + Abel-Jacobi map for hyperelliptic curves

- `magma`

  generalization for superelliptic curves with Christian Neurohr,
  complete Abel-Jacobi map.

- `arb`

  Rigorous arb implementation, currently only period matrix

## Example

```
cd arb && make example
```
then
```
build/examples/periods --pol 4 1 0 -2 3 0
```
outputs a period matrix for the curve `y^2 = x^4 -2x^2 + 3x`.

Use options `-m 5` to switch to the curve `y^5 = x^4 -2x^2 + 3x`,
or `--prec 1024` for 1024 bits precision.

Other options:
- `--gp` output fo pari/gp
- `--de` force use of double exponential integration (instead of Gauss, if m = 2).
- `--big` return *big* period matrices instead of reduced matrix tau.
