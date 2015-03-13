#!/bin/bash 
# Usage find_observer_level.bash altitude_from_sea_level domain_minimum_altitude
# dependancies bc
#
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
if [ ! $1 ] 
then echo "Usage find_observer_level.bash topography_file x y"
     exit 0
fi
if [ ! $2 ] 
then echo "Usage find_observer_level.bash topography_file x y"
     exit 0
fi
if [ ! $3 ] 
then echo "Usage find_observer_level.bash topography_file x y"
     exit 0
fi
n=51
absmin=100000000000
# model levels definition
mod_alt=" 27345.16 22636.31 18738.26 15511.4 12840.16 10628.87 8798.33 7282.98 6028.55 4990.11 4130.47 3418.85 2829.76 2342.1 1938.41 1604.23 1327.59 1098.58 909. 752.06 622.14 514.59 425.56 351.86 290.85 240.35 198.55 163.95 135.31 111.6 91.97 75.72 62.27 51.14 41.93 34.31 28. 22.77 18.44 14.86 11.9 9.45 7.42 5.74 4.35 3.2 2.25 1.46 0.8 0.25 "
# domain min altitude
find_min_altitude.bash $1
read minalt < minaltitude.out
readvalue.bash $1 $2 $3
read obsalt < readvalue.out

# multiplier les valeurs par 100 pour travailler en cm (nombre entier)
/bin/echo "scale=0; "$obsalt"* 100." |/usr/bin/bc -l |sed 's/\./ /g' > obsalt.tmp
read obsalt bidon < obsalt.tmp
/bin/echo "scale=0; "$minalt"* 100." |/usr/bin/bc -l |sed 's/\./ /g' > minalt.tmp
read minalt bidon < minalt.tmp


for h in $mod_alt
do let n=n-1
   /bin/echo "scale=0; "$h"* 100." |/usr/bin/bc -l |sed 's/\./ /g' > h.tmp
   read h bidon < h.tmp
   let min=(obsalt-minalt-h)**2
   if [ $min -lt $absmin ]
   then let absmin=$min
        level=$n
   fi
done
echo "Observer level @ x="$2" y="$3" :" $level
