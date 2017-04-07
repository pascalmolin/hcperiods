#!/bin/bash
if (( $# < 2 )); then
    echo "$0 [callgrind|bench] pol [other arguments]"
    exit 0
else
    run="$1"; shift
    pol="$1"; shift
fi
export LD_LIBRARY_PATH=.:${HOME}/install/lib;
mkdir -p profile/out
PREC="128 512 2000 4000 10000"
MARG="2 3 4 9 13"
MARG="2 3 4"
NARG="3 5 7 8 13 30"
echo "# extra args $*"
echo "# pol; m; n; prec; time"
for m in $MARG
do
    for n in $NARG
    do
        for p in $PREC
        do
            if [ "$run" == "callgrind" ]; then
                echo "$pol $m $n $p $1"
                valgrind -q --tool=callgrind --callgrind-out-file=profile/out/callgrind.$pol.m$m.n$n.p$p$1.out build/example/periods --bench 1 --int --$pol $n -m $m --prec $p $1
            elif [ "$run" == "bench" ]; then
                echo -n "$pol; $m; $n; $p; $*"
                dumbbench --raw -- build/example/periods --bench 1 --int --$pol $n -m $m --prec $p $*
            else
                echo "build/example/periods --bench 1 --$pol $n -m $m --prec $p $*"
            fi
        done
    done
done
