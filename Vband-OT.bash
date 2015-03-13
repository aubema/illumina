#!/bin/bash
# Vmag est un script pour generer les V magnitudes
# a partir des radiances 546 et 569 nm
#
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
list=`ls -1 OT*568*txt`
for i in $list
do j=`echo $i | sed 's/568/546/g'`
   k=`echo $i | sed 's/b_/a_/g' | sed 's/0.025/0.050/g' | sed 's/0.100/0.050/g' | sed 's/0.200/0.050/g'`
   l=`echo $i | sed 's/568/546/g' | sed 's/b_/a_/g' | sed 's/0.025/0.050/g' | sed 's/0.100/0.050/g' | sed 's/0.200/0.050/g'`
   cp -f $i input569.txt
   cp -f $j input546.txt
   cp -f $k after569.txt
   cp -f $l after546.txt
   ./L569-546toVmag-OT
   out=`echo $i | sed 's/wl568/VMag/g'`
   cp -f output.txt $out
   mkAllSkyFitMag.bash $out 1 0
   echo $i $j $k $l
done
