#!/bin/bash
#    Copyright (C) 2011  Martin Aube
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
ls
echo "Experiment folder name?"
read folder
echo "Experiment name?"
read expname
pos=`pwd`/$folder
chaine=`echo $pos | sed 's/\// /g'`
list=`grep cd $folder/TortureMammouth | sed 's/cd //g'`

for i in $list

do echo $i | sed 's/\// /g' | sed "s/$chaine/ /g" | sed 's/el//g' | sed 's/az//g' > toto.tmp
   read site ho ro sd rd el az ta wl < toto.tmp
   tail -4 $i/$expname".out" | head -1 | grep "E-" > tata.tmp
   read data < tata.tmp
   if [ $site = "x133y178" ]
   then obs="ORM-"
   fi
   if [ $site = "x268y123" ]
   then obs="OT-"
   fi
   echo $el $az $data >> $obs$folder"-"$ta"-"$wl".txt"
   cp $i/$expname"_pcl.pgm" "PCL-"$obs$folder"-"$ta"-"$wl"-el"$el"-az"$az".pgm"
   cp $i/$expname"_pcw.pgm" "PCW-"$obs$folder"-"$ta"-"$wl"-el"$el"-az"$az".pgm"
done
rm -f toto.tmp
rm -f tata.tmp
