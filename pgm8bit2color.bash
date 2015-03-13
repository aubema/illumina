#!/bin/bash
#Usage pgm2color.bash pgm_file_name
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
if [ ! $1 ]
then  echo "Usage: pgm2color.bash pgm_file_name"
    echo "Specify the file name (1st parameter)"
    exit 1
fi
base=`echo $1 | sed 's/\.pgm//g'`
echo $base
ncomment=`/bin/grep -c "#" $1`
echo "ncomment= " $ncomment
let "nheader=ncomment+2"
echo "nheader= " $nheader
vect=(`head -$nheader $1 | tail -1`)
nx=${vect[0]}
ny=${vect[1]}
convert  $1 toto.rgb
grep gain $1 > bidon.tmp 
read toto  toto gain toto < bidon.tmp
grep offset $1 > bidon.tmp 
read toto  toto offset toto < bidon.tmp
ga=`echo $gain |sed 's/ //g'`
of=`echo $offset |sed 's/ //g'`
echo "set term png size 730 480" > toto.gplot
echo "set xrange[0:"$nx"]" >> toto.gplot
echo "set yrange[0:"$ny"]" >> toto.gplot
echo "unset key" >> toto.gplot
echo "set title '" $base "'" >> toto.gplot
echo "plot 'toto.rgb' binary array="$nx"x"$ny" flipy  format='%uchar%uchar%uchar' using (((\$1)*"$ga"+"$of")) with image" >>toto.gplot
gnuplot  < toto.gplot > $base".png"
#display $base".png" &
