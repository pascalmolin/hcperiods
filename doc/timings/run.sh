#!/bin/bash
if (( $# < 2 )); then
    echo "$0 <pol> <m> [-alg arb] [--csv]"
    exit 0
fi
pol="$1"; shift
m="$1"; shift
if [ $m -le 2 ]; then
    METHODARG="oldmagma newmagma arb"
else
    METHODARG="newmagma arb"
fi
while (($# > 0)); do
    if [ "$1" == "-alg" ]; then
        shift; METHODARG="$1"; shift
    elif [ "$1" == "--csv" ]; then
        CSV=1; shift
    fi
done
PREC="128 512 2000 4000 10000"
MARG="2 3 4 9 13"
NARG="8 30 80"
echo "# method; pol; m; n; prec; time"
for n in $NARG
do
    for method in $METHODARG
    do
        if (( $CSV )); then
            echo -n "$pol; $method $m; $n"
        fi
        for prec in $PREC
        do
            if [ "$method" != "oldmagma" ] || [ "$prec" -le "4000" ]; then
                if (( $CSV )); then
                    echo -n "; "
                else
                    echo -n "$pol; $method $m; $n; $prec; "
                fi
                ./compare.sh $method $pol $m $n $prec
            fi
        done
    done
done
