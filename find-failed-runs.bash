#!/bin/bash

dir="."

while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)
    echo "Usage: find-failed-runs.bash [OPTIONS]"
    echo "Used to find runs that haven't completed properly"
    echo "Call this from a directory containing the runs to test"
    echo "OPTIONS"
    echo "    -h,--help  Print this message"
    echo "    -e,--exec  Produce the output needed to reruns the failed runs"
    exit 1
    ;;
    -e|--exec)
    dir=`pwd -P`
    full_out="True"
    ;;
esac
shift
done

find $dir -name wl* | grep -v gridmap | while read dirname
do
    fname=`ls $dirname/*.out 2>/dev/null | grep -v -e "slurm" -e "mie"`
    if [ -z "$fname" ]
    then
	n=0
    else
	n=`grep -c -e "Sky radiance" $fname`
    fi
    if [ $n -eq 0 ]
    then
        if [ -z $full_out ]
        then
            echo $dirname
        else
            echo "cd $dirname"
            echo "sbatch ./execute"
            echo "sleep 0.05"
        fi
    fi
done
