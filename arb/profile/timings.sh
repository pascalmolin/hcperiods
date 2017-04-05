#!/bin/bash

NARG="3 5 8 13 21 30"
MARG="2 3 4 5 9 13"
PREC="64 256 1024 4096"
echo "# extra arg $1"
echo "# prec; m; n; time"
for p in $PREC
do
    for m in $MARG
    do
        for n in $NARG
        do
            printf "$p; $m; $n; "
            DO="build/example/periods --bench 1 --bern $n -m $m --prec $p"
            if $DO &> /dev/null; then
                dumbbench --raw -- $DO
            else
                echo "FAIL"
            fi
        done
    done
done
