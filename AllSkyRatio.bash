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
#  Usage: AllSkyRatio.bash file_name_1 wavelength_1 aerosol_optical_depth_1 file_name_2 wavelength_2 aerosol_optical_depth_2
# 
#  The program use the aerosol optical depth standard model values
#  the wavelength have to be choosen among the 5 standard wavelength.
#
#  file_name_1 and file_name_2 should correspond to 2 distinct ILLUMINA modeling experiment 
#  e.g. 2 different period (summer/winter, two different years etc) or 2 different relative humidity etc.
#
here=`pwd`
echo $here
wlstd="436 498 546 569 616"
tastd="0.05 0.1 0.2 0.5 1.0"
rm -f file1.in
rm -f file2.in
# ouput file name
outname=$1"_wl"$2"_ta"$3"-over-"$4"_wl"$5"_ta"$6
if [ ! $1 ]
then echo "Usage: AllSkyRatio.bash file_name_1 wavelength_1 aerosol_optical_depth_1 file_name_2 wavelength_2 aerosol_optical_depth_2"
     exit 0
fi

if  echo $wlstd | grep $2 
then echo "wavelength1=" $2
else 
     echo "Bad wavelength: "$2 "("$wlstd")"
     exit 0
fi
if  echo $wlstd | grep $5 
then echo "wavelength2=" $5
else 
     echo "Bad wavelength: "$5 "("$wlstd")"
     exit 0
fi
if  echo $tastd | grep $3 
then echo "aerosol optical depth 1=" $3
else 
     echo "Bad aerosol optical depth 1: "$3 "("$tastd")"
     exit 0
fi
if  echo $tastd | grep $6 
then echo "aerosol optical depth 2=" $6
else 
     echo "Bad aerosol optical depth 2: "$6 "("$tastd")"
     exit 0
fi

# modify input files format
cat $here/$1 | grep "wl"$2 | grep "ta"$3 | sed 's/-/ /g' | sed 's/el/ /g' | sed 's/az/ /g' | sed 's/E /E-/g'> $here/file1.tmp
cat $here/$4 | grep "wl"$5 | grep "ta"$6 | sed 's/-/ /g' | sed 's/el/ /g' | sed 's/az/ /g' | sed 's/E /E-/g' > $here/file2.tmp
pwd
nligne=`/bin/grep -c "" file1.tmp`
n=0
while [ $n -lt  $nligne ]
do let n=n+1
   head -$n file1.tmp | tail -1 > ligne.tmp
   read bidon bidon bidon bidon el az rad < ligne.tmp
   echo $el $az $rad >> file1.in 
done
len1=`grep -c "" file1.in`

nligne=`/bin/grep -c "" file2.tmp`
n=0
while [ $n -lt  $nligne ]
do let n=n+1
   head -$n file2.tmp | tail -1 > ligne.tmp
   read bidon bidon bidon bidon el az rad < ligne.tmp
   echo $el $az $rad >> file2.in 
done
len2=`grep -c "" file2.in`

# prepare input file for AllSkyRatio.f
echo "file1.in" $len1 > AllSkyRatio.tmp
echo "file2.in" $len2 >> AllSkyRatio.tmp
echo $3 >> AllSkyRatio.tmp
echo $outname >> AllSkyRatio.tmp
AllSkyRatio < AllSkyRatio.tmp
#rm -f file1.tmp
#rm -f file2.tmp
#rm -f AllSkyRatio.tmp
#rm -f file1.in
#rm -f file2.in
