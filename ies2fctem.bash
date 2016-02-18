#!/bin/bash
# programme permettant de transformer un fichier IES en
# fichier *_fctem_*.dat
# sur 36 angles de 2.5 a 177.5 (intervalles de 5 deg)
# le fichier de sortie est un intrant de ILLUMINA
# la fonction photometrique est moyennee en azimuth
#
#
#    Copyright (C) 2009  Martin Aube
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
x=`basename $0`
if test $# = 0; then
  echo usage: ${x} "IES_file_name tilt_angle(deg) output_file_name"
  exit 1
fi
echo $1
echo $1 > ies2tab.in
echo $1"_tab" >> ies2tab.in
ies2tab.py < ies2tab.in
rm -f ies2tab.in
mv $1"_tab" ./ies2fctem_tab
echo "ies2fctem_tab" > ies2fctem.in
echo $2 >> ies2fctem.in 
echo $3 >> ies2fctem.in
head -1 ies2fctem_tab > tab.tmp
read ntheta nphi <  tab.tmp
ies2fctem.exe

#rm -f ies2fctem.in
#rm -f ies2tab.tmp
