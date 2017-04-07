#!/bin/bash
if (( $# < 1 )); then
    echo "$0 <prec> <degree>"
    exit 0
fi
prec=$1;
n=$2;
tmp=$(mktemp run.m.XXXX)
cat <<EOF > $tmp
C<I> := ComplexField($prec:Bits:=true);
K<x> := PolynomialRing(C);
f:=BernoulliPolynomial($n);
fx := Evaluate(f,x);
t := Cputime();
A := AnalyticJacobian(fx);
M := BigPeriodMatrix(A);
t := Cputime(t);
printf "$prec; 2; $n; %o\n", t;
quit;
EOF
magma -b $tmp
rm -f $tmp
