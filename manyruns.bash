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
    nnew=`grep -c -e "sbatch" $i`
    ntot=$((nnew+njob))
    echo $njob "initial"
    echo $nnew "new"
    while  [ $ntot -ge 1000 ]
    do 
        njob=`squeue -u $user -h | grep -c ""`
        ntot=$((nnew+njob))
        echo $ntot
        sleep 300 
    done
done
echo "End of multiple experiments"
