#!/bin/bash
if (( $# < 1 )); then
    echo "$0 [callgrind|bench] [periods arguments]"
    exit 0
else
    run="$1"; shift
fi
export LD_LIBRARY_PATH=.:${HOME}/install/lib;
PREC="128 512 2000 10000"
MARG="2 3 4 9 13"
NARG="3 5 7 13 31"
echo "# extra args $*"
echo "# prec; m; n; time"
for p in $PREC
do
    for m in $MARG
    do
        for n in $NARG
        do
            if [ "$run" == "callgrind" ]; then
                echo "$p $m $n $1"
                valgrind -q --tool=callgrind --callgrind-out-file=profile/callgrind.p$p.m$m.n$n$1.out build/example/periods --bench 1 --int --bern $n -m $m --prec $p $1
            elif [ "$run" == "bench" ]; then
                echo "$p; $m; $n; $*"
                dumbbench --raw -- build/example/periods --bench 1 --bern $n -m $m --prec $p $*
            else
                echo "build/example/periods --bench 1 --bern $n -m $m --prec $p $*"
            fi
        done
    done
done
