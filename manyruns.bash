#!/bin/sh

if [[ $# -ne 1 ]]; then
    >&2 echo "The command expect a prefix. Don't put a wildcard in there ! (*)"
    exit 1
fi

batchlist=`ls -1 $1*`
echo $batchlist
for i in $batchlist
do  bash $i
    echo $i
    echo "===="
    user=`whoami`
    njob=`squeue -u $user -h | grep -c ""`
    nwait=`squeue -u $user -h -o "%t" | grep -c "PD"`
    echo $njob "initial"
    while  [ $njob -gt 1000 ] && [ $nwait -gt 100 ]
    do 
        njob=`squeue -u $user -h | grep -c ""`
        nwait=`squeue -u $user -h -o "%t" | grep -c "PD"`
        echo $njob $nwait
        sleep 30 
    done
done
echo "End of multiple experiments"
