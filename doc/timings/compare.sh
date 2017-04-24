#!/bin/bash
if (( $# < 4 )); then
    echo "$0 <method> <pol> <m> <n> <prec>"
    exit 0
fi
method=$1; shift
pol=$1; shift
m=$1; shift
n=$1; shift
prec=$1; shift

ARB=$HOME/git/hcperiods/arb

bench_magma() {

    tmp=$(mktemp run.m.XXXX)

    if [ "$method" == "oldmagma" ]; then
        echo "K<x> := PolynomialRing(ComplexField($prec:Bits:=true));" > $tmp
    else
        echo 'AttachSpec("../../magma/spec");' > $tmp
        echo "K<x> := PolynomialRing(Rationals());" >> $tmp
    fi;

    if [ "$pol" == "bern" ]; then
        echo "f:=K!BernoulliPolynomial($n);" >> $tmp
    elif [ "$pol" == "bernrev" ]; then 
        echo "f:=Reverse(K!BernoulliPolynomial($n));" >> $tmp
    elif [ "$pol" == "exp" ]; then
        echo "f:=K![1 / Factorial(k) : k in [0..$n]];" >> $tmp
    elif [ "$pol" == "exprev" ]; then
        echo "f:=K![1 / Factorial($n - k) : k in [0..$n]];" >> $tmp
    else
        echo "$0 <method> <pol> <m> <n> <prec>"
        rm -f $tmp
        exit 0
    fi

    echo "t := Cputime();" >> $tmp

    if [ "$method" == "oldmagma" ]; then
        echo "A := AnalyticJacobian(f);" >> $tmp
        echo "M := BigPeriodMatrix(A);" >> $tmp
    else
        echo "prec := Ceiling($prec/Log(2,10));" >> $tmp
        echo "M := SE_BigPeriodMatrix(f,$m:Prec:=prec);" >> $tmp
    fi
    #printf "$pol; 2; $n; $prec; %o\n", t;
cat <<EOF >> $tmp
t := Cputime(t);
printf "%o\n", t;
quit;
EOF
    magma -b $tmp
    rm -f $tmp
}

if [ "$method" == "oldmagma" ] || [ "$method" == "newmagma" ]; then
        bench_magma
elif [ "$method" == "arb" ]; then
    export LD_LIBRARY_PATH=.:$HOME/install/lib:$ARB;
    dumbbench --raw -- $ARB/build/example/periods --bench 1 --int --$pol $n -m $m --prec $prec
else
    echo "$0 <method> <pol> <m> <n> <prec>"
    exit 0
fi
