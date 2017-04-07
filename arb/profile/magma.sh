#!/bin/bash
if (( $# < 3 )); then
    echo "$0 <pol> <prec> <degree>"
    exit 0
fi
pol=$1; shift
prec=$1; shift
n=$1; shift
tmp=$(mktemp run.m.XXXX)
cat <<EOF > $tmp
C<I> := ComplexField($prec:Bits:=true);
K<x> := PolynomialRing(C);
EOF
if [ "$pol" == "bern" ]; then
    echo "f:=K!BernoulliPolynomial($n);" >> $tmp
elif [ "$pol" == "bernrev" ]; then 
    echo "f:=Reverse(K!BernoulliPolynomial($n));" >> $tmp
elif [ "$pol" == "exp" ]; then
    echo "f:=K![1 / Factorial(k) : k in [0..$n]];" >> $tmp
elif [ "$pol" == "exprev" ]; then
    echo "f:=K![1 / Factorial($n - k) : k in [0..$n]];" >> $tmp
else
    exit 0
fi
cat <<EOF >> $tmp
t := Cputime();
A := AnalyticJacobian(f);
M := BigPeriodMatrix(A);
t := Cputime(t);
printf "$pol; 2; $n; $prec; %o\n", t;
quit;
EOF
magma -b $tmp
rm -f $tmp
