#!/bin/sh

batchlist=`ls -1 $1*`
echo $batchlist
for i in $batchlist
do  bash $i
    echo $i
    echo "===="
    user=`whoami`
    njob=`bqstat 2>/dev/null | grep -c $user`
    echo $njob "initial"
    while  [ $njob -gt 1000 ]
    do njob=`bqstat 2>/dev/null | grep -c $user`
       echo $njob
       sleep 30 
    done
done
echo "End of multiple experiments"
