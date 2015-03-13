#!/bin/bash
#  Script de reprise d'un run ILLUMINA interrompu
#  2005
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
folder=`pwd`
echo "Enter the root name of your previous experiment (e.g. if you have a file named test5.out, the root name is test5"
read name
# find the level of completion of the run
nb=`grep CAPTEUR $name.out | grep -c ""`
let "nb=nb-1"

grep -m $nb CAPTEUR $name.out > continue_illumina.1.tmp
# preparer illumina.in
grep "! OBSERVER X POSITION" illumina.in > ligne.tmp
cp illumina.in illumina.in.tmp
read toto1 toto2 toto3 toto4 bidon < ligne.tmp
sed -e "s/${toto1} ${toto2} ${toto3} ${toto4}/${toto1} ${toto2} ${toto3} ${nb}/g" illumina.in.tmp > illumina.in
# lancer illumina
illumina
grep CAPTEUR $name.out >> continue_illumina.2.tmp
cat continue_illumina.1.tmp continue_illumina.2.tmp > 
rm -f *.tmp
# calcul du flux total
echo $name > continue_illumina.in
continue_illumina
