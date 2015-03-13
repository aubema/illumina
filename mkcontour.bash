#!/bin/bash
# script pour creer les contour plot
# usage mkcontour filename vmin vmax deltav (el, az, rad) 
# dependancy gri <= 2.12.22 (http://sourceforge.net/projects/gri
#
#
rm -f gri-*.ps
rm -f data.polar
rm -f rad2polar.tmp
cp -f $HOME/svn/illumina/trunk/contourplot.gri ./contourplot.gri.tmp
cp -f $1 rad2polar.tmp

nligne=`grep -c "" $1`
echo $nligne >nligne.tmp
rad2polar
#read a b c bidon < rad2polarscale.tmp
#echo $a $b $c
cat contourplot.gri.tmp | sed -e "s/scale/${2} ${3} ${4}/g" > ./contourplot.gri
gri < ./contourplot.gri
mv gri-00.ps $1.ps
convert -density 288 -alpha remove $1.ps $1-ctr.png
display $1-ctr.png&
rm -f *tmp

