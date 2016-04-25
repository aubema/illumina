#!/bin/bash
#    Copyright (C) 2012  Martin Aube
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: martin.aube@cegepsherbrooke.qc.ca
#
# to be executed when located in a directory containing the gridmap directory
# before executing create a file named LIST with the files created with makeBATCH 
# with a command like "cat $HOME/file1 $HOME/file2 > ./LIST" 
# 
#  usage: 
#  bash extract-output-data.bash experiment_name File_containing_all_batch_run_scripts 
#

while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)
    echo "Usage: extract-output-data.bash experiment_name File_containing_all_batch_run_scripts"
    exit 1
    ;;
    -s|--summary)
    summary="do it"
    echo "Producing simplified output"
    ;;
    *)
    [ -z "$expname" ] && expname=$1 || list=`grep cd $1 | sed 's/cd //g'`
    ;;
esac
shift
done

folder=`pwd`

echo $expname

#
#  copying retres   
#
if [ -z "$summary" ]
then
    cp -f $HOME/hg/illumina/bin/retres .
fi
#
#  creating result directory
#    
mkdir Results
cd Results
if [ -z "$summary" ]
then
    mkdir $expname
    echo "" > $expname"/data.txt"
else
    echo "" > $expname".txt"
fi
cd ..
for i in $list
do 
#
#   archive the files used for the run just in case
#
#    mkdir new_files
#    cp -f *lumlp* new_files
#    cp -f *reflect* new_files
   echo $i | sed 's/\//\n/g' > bidon.tmp
   site=`grep x bidon.tmp | grep y`
   wl=`grep wl bidon.tmp| tail -1`
   az=`grep az bidon.tmp| tail -1`
   el=`grep el bidon.tmp| tail -1`
   ta=`grep ta bidon.tmp| tail -1`
   sd=`grep sd bidon.tmp| tail -1`
   rd=`grep rd bidon.tmp| tail -1` 
   echo $i
   echo "site="$site " wl="$wl " az="$az " el="$el " ta="$ta " sd="$sd " rd="$rd
#
#  bring back here original lumlp files and grid.txt 
#
if [ -z "$summary" ]
then
   cp -f $i/original_files/* .
   cp -f $i/grid.txt .
fi
#
#  copy output files hereDenver-summer.list
#
if [ -z "$summary" ]
then
   cp -f $i/$expname".out" .
   cp -f $i/*recombined.pgm .
   cp -f $i/*pcw.pgm .
   cp -f $i/*pcl.pgm .
#   
#  creating the input file for retres
#
   echo $expname"_lumlp_recombined.pgm" > retres.in
   echo $expname"_pcl.pgm" >> retres.in
   echo $expname"_pcl_new.pgm" >> retres.in
#
# executing retres   
#
   ./retres
#   
#  creating the input file for retres
#   
   echo $expname"_lumlp_recombined.pgm" > retres.in
   echo $expname"_pcw.pgm" >> retres.in
   echo $expname"_pcw_new.pgm" >> retres.in
fi

#
#  executing retres   
#
if [ -z "$summary" ]
then
   ./retres
fi
   o=$site"-"$expname"-"$ta"-"$wl"-"$el"-"$az"-"$rd"-"$sd
   data=`tail -4 $i"/"$expname".out" | head -1 | grep -E "E+|E-"`
if [ -z "$summary" ]
then
   echo $o $data >> "Results/"$expname"/data.txt"
else
   echo $o $data >> "Results/"$expname".txt"
fi

echo $o

if [ -z "$summary" ]
then
   mv ./$expname"_pcl_new.pgm" "./Results/"$expname"/PCL-"$o".pgm"
   mv ./$expname"_pcw_new.pgm" "./Results/"$expname"/PCW-"$o".pgm"
   mv ./$expname".out" "./Results/"$expname"/"$o".out"
fi
done
