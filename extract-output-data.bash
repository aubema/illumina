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
ls
# Experiment name as indicated while running makeBATCH 
expname=$1
list=`grep cd $2 | sed 's/cd //g'`
folder=`pwd`
#
#  copying retres   
#
cp -f $HOME/hg/illumina/bin/retres .
#
#  creating result directory
#    
mkdir Results
cd Results
mkdir $expname
echo "" > "Results/"$expname"/data.txt"
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
   ho=`grep ho bidon.tmp| tail -1`
   ro=`grep ro bidon.tmp| tail -1`
   wl=`grep wl bidon.tmp| tail -1`
   az=`grep az bidon.tmp| tail -1`
   el=`grep el bidon.tmp| tail -1`
   ta=`grep ta bidon.tmp| tail -1`
   sd=`grep sd bidon.tmp| tail -1`
   rd=`grep rd bidon.tmp| tail -1` 
   echo "site="$site " ho="$ho "ro="$ro " wl="$wl " az="$az " el="$el " ta="$ta " sd="$sd " rd="$rd
#
#  bring back here original lumlp files and grid.txt 
#
   cp -f $i/original_files/* .
   cp -f $i/grid.txt .
#
#  copy output files hereDenver-summer.list
#
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

#
#  executing retres   
#
   echo "expname=" $expname   
   ./retres
   o=$site"-"$expname"-"$ta"-"$wl"-"$el"-"$az-"-"$ho"-"$ro"-"$rd"-"$sd
   data=`tail -4 $expname".out" | head -1 | grep -E "E(-|+)"`
   echo $o $data >> "Results/"$expname"/data.txt"
   ls
   pwd
echo $o
   mv ./$expname"_pcl_new.pgm" "./Results/"$expname"/PCL-"$o".pgm"
   mv ./$expname"_pcw_new.pgm" "./Results/"$expname"/PCW-"$o".pgm"
   mv ./$expname".out" "./Results/"$expname"/"$o".out"
done
