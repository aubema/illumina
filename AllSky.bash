#!/bin/bash
#    Copyright (C) 2013  Martin Aube
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
#
#  Usage: AllSky.bash file_name_1 calibration_factor wavelength aerosol_optical_depth 
#
here=`pwd`
echo $here

rm -f file1.in
# ouput file name
outname=$1"_wl"$3"_ta"$4
rm -f $outname
if [ ! $1 ]
then echo "Usage: AllSky.bash file_name_1 calibration_factor wavelength aerosol_optical_depth "
fi
# modify input files format
cat $here/$1 | grep "wl"$3 | grep "ta"$4 | sed 's/-/ /g' | sed 's/el/ /g' | sed 's/az/ /g' | sed 's/E /E-/g'> $here/file1.tmp
nligne=`/bin/grep -c " " file1.tmp`
n=0
while [ $n -lt  $nligne ]
do let n=n+1
   head -$n file1.tmp | tail -1 > ligne.tmp
   read bidon bidon bidon bidon el az rad bidon < ligne.tmp
   num=`echo ${rad} | sed 's/E/\*10\^/g' | sed 's/+//g'`
   echo ${rad} | sed 's/E/\*10\^/g' | sed 's/+//g'
   radc=`/bin/echo $num "*" $2 | /usr/bin/bc -l`
   echo $el $az $radc >> $outname
done
len1=`grep -c "" $outname`

#rm -f file1.tmp
#rm -f AllSky.tmp
#rm -f file2.in
