# hcperiods

Numerical evaluation of periods matrices of hyperelliptic and
superelliptic curves to arbitrary precision by numerical integration.

## Structure

- `gp`

  Initial projet as [gp](https://pari.math.u-bordeaux.fr) script,
  period matrix + Abel-Jacobi map for hyperelliptic curves

- `magma`

  Generalization for superelliptic curves with Christian Neurohr,
  complete Abel-Jacobi map. Written in [magma](http://magma.maths.usyd.edu.au/magma/).

- `arb`

  Rigorous [Arb](http://arblib.org/index.html) implementation,
  currently only period matrices.

- [hcperiods.pdf](hcperiods.pdf) description of the algorithm and proofs

## How to use

See the [magma intructions](magma/README.md) to use the magma package.

The arb version can be used as follows
(needs [arb>=v2.12.0 installed](http://arblib.org/setup.html))
```
cd arb && make example
```
then
```
build/examples/periods -m 2 --pol 4 1 0 -2 3 0
```
outputs a period matrix for the curve `y^2 = x^4 -2x^2 + 3x`.

Use options `-m 5` to switch to the curve `y^5 = x^4 -2x^2 + 3x`,
or `--prec 1024` for 1024 bits precision.

Other options:
- `--gp` output fo pari/gp
- `--de` force use of double exponential integration (instead of Gauss, if m = 2).
- `--big` return *big* period matrices instead of reduced matrix tau.
