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
#  Usage: rad2lrad.bash file_name
#
here=`pwd`
echo $here
# ouput file name
outname=$1"-log"
if [ ! $1 ]
then echo "Usage: rad2lrad.bash file_name"
fi
cp -f $1 rad2lograd.tmp
grep -c "" $1 > nligne.tmp
echo $outname >>  nligne.tmp
rad2lrad

