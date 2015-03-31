batchlist=`ls -1 $1*`
echo $batchlist
for i in $batchlist
do  bash $i
    echo $i
    echo "===="
    bqmon -u | grep maube > currentjobs.tmp
    read toto que req run toto < currentjobs.tmp
    let njob=que+run
    echo $njob "initial"
    while  [ $njob -gt 750 ]
    do bqmon -u | grep maube > currentjobs.tmp
       read toto que req run toto < currentjobs.tmp
       let njob=que+run
       echo $njob
       sleep 30
    done
done
echo "End of multiple experiments"

