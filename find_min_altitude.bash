#!/bin/bash 
# Usage find_min_altitude.bash mna_pgm_file
# dependancies bc
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
myFile=$1
nligne=`/bin/grep -c "" $myFile`
echo "nligne= " $nligne
ncomment=`/bin/grep -c "#" $myFile`
echo "ncomment= " $ncomment
gain=`/bin/grep "gain" $myFile | sed -e 's/#//g' | sed -e 's/gain//g'`
echo "gain= " $gain
offset=`grep "offset" $myFile | sed 's/#//g' | sed 's/offset//g'`
echo "offset= " $offset
let "nheader=ncomment+2"
echo "nheader= " $nheader
let "ndat=nligne-nheader"
echo "ndat= " $ndat
#  Loop over each observation line
   dmin=100000
   content=`tail -$ndat $myFile`
   for d in $content
   do  if [ $d -lt $dmin ] 
       then echo $dmin $d
            dmin=$d
       fi
   done
echo "dmin= " $dmin
echo "Equation= " $dmin "x" $gain "+" $offset
min=`/bin/echo $dmin "*" $gain "+" $offset | /usr/bin/bc -l`
echo "----------------------------"
echo "Minimum altitude=" $min
echo $min > minaltitude.out
