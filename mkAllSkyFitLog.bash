#!/bin/bash
#Usage mkAllSkyFitLog.bash -f txt_file_name -c gain -g Azimuth_mode -s start_val -e end_val
# Azimuth_mode = 1 means geographical azimuth (0=N clockwise)
# Azimuth_mode = 0  means mathematical azimuth (0=E counter clockwise)
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
fflag=  
sflag=
eflag=
gflag=
cflag=
          # fval file name to plot
sval=-13  # start z value to plot
eval=-5   # end z value to plot
gval=1    # set azimuth mode 
cval="1." # calibration factor
while getopts 'f:s:e:g:c:' OPTION
do
  case $OPTION in
    f)  fflag=1
        fval="$OPTARG"  # accepted values at or ss for astronomical twilight of sunset
    ;;
    s)  sflag=1
        sval="$OPTARG"
    ;;
    e)  eflag=1
        eval="$OPTARG"
    ;;
    g)  gflag=1
        gval="$OPTARG"
    ;;
    c)  cflag=1
        cval="$OPTARG"
    ;;
    ?)  printf "Usage: %s: [-f value] [-s value] [-e value] [-g value] [-c value] args\n" $(basename $0) >&2
        exit 2
    ;;
  esac
done
shift $(($OPTIND - 1))
if [ ! $fflag ] 
then echo "Please specify a file name with option: -f file_name"
     exit 1
fi
size=`echo "("$eval"-("$sval"))*10"| bc -l | xargs printf "%1.0f"`
base=`echo $fval | sed 's/\.txt//g'`
echo "set term png size 480 480" > toto.gplot
echo "set title '"$base"'" >>toto.gplot
echo "h(x,y)=a+b*x+c*y+d*x**2+e*y*x+f*y**2+g*x**3+h*x**2*y+i*x*y**2+j*y**3+k*x**4+l*x**2*y**2+m*y**4+n*x**5+o*y**5" >> toto.gplot
calib=$cval
echo "a=0" >> toto.gplot
echo "b=0" >> toto.gplot
echo "c=0" >> toto.gplot
echo "d=0" >> toto.gplot
echo "e=0" >> toto.gplot
echo "f=0" >> toto.gplot
echo "g=0" >> toto.gplot
echo "h=0" >> toto.gplot
echo "i=0" >> toto.gplot
echo "j=0" >> toto.gplot
echo "k=0" >> toto.gplot
echo "l=0" >> toto.gplot
echo "m=0" >> toto.gplot
echo "n=0" >> toto.gplot
echo "o=0" >> toto.gplot
echo "p=0" >> toto.gplot
echo "q=0" >> toto.gplot
echo "r=0" >> toto.gplot
echo "s=0" >> toto.gplot
echo "t=0" >> toto.gplot
echo "u=0" >> toto.gplot
echo "v=0" >> toto.gplot
echo "w=0" >> toto.gplot
echo "FIT_LIMIT = 1e-99" >> toto.gplot
echo "FIT_MAXITER=100" >> toto.gplot
echo "set angles degrees" >> toto.gplot
if [ $gval -eq 0 ]
then echo "fit h(x,y) '"$fval"' using ((90-\$1)*cos(\$2)):((90-\$1)*sin(\$2)):(log10(\$3*"$calib ")):(1) via a,b,c,d,e,f,g,h,i,j,k,l,m,n,o" >> toto.gplot
else echo "fit h(x,y) '"$fval"' using ((90-\$1)*cos(90-\$2)):((90-\$1)*sin(90-\$2)):(log10(\$3*"$calib ")):(1) via a,b,c,d,e,f,g,h,i,j,k,l,m,n,o" >> toto.gplot
fi
#echo "set palette defined (  0 '"#000000"' ,  1 '"#0000FF"', 2 '"#0CFFFF"', 3 '"#00FF00"', 4 '"#FFFF00"' , 5 '"#FF0000"' ,6  '"#FF02EA7"',7 '"#FFFFFF"' )"  >> toto.gplot
echo "set palette defined (  0 '"#000000"' ,  1 '"#0000FF"', 2 '"#0CFFFF"', 3 '"#00FF00"', 4 '"#FFFF00"' , 5 '"#FF8000"' , 6 '"#FF0000"' ,7  '"#FF02EA7"',8 '"#FFFFFF"' )"  >> toto.gplot

echo "set pal maxcolors "$size >> toto.gplot
echo "set xrange [-90:90]" >> toto.gplot
echo "set yrange [-90:90]" >> toto.gplot
echo "set cbrange ["$sval":"$eval"]" >> toto.gplot
echo "set size square 1.25,1.25" >> toto.gplot
echo "set origin -0.125,-0.1" >> toto.gplot
echo "set sample 308" >> toto.gplot
echo "set isosample 308" >> toto.gplot
echo "set pm3d map" >> toto.gplot
echo "unset border" >> toto.gplot
echo "unset key" >> toto.gplot
echo "set zeroaxis lt 2" >> toto.gplot
echo "set xtics axis textcolor lt 2" >> toto.gplot
echo "set ytics axis textcolor lt 2" >> toto.gplot
echo "set border linetype 2" >> toto.gplot
echo "set style line 2604 linetype -1" >> toto.gplot
echo "set colorbox border 2604" >> toto.gplot
echo "set colorbox horizontal" >> toto.gplot
echo "set colorbox user origin 0.1,0.04 size 0.8,0.02" >> toto.gplot
echo "set cbtics textcolor lt -1" >> toto.gplot
echo "set cbtics scale 1 1" >> toto.gplot
echo "set arrow from 0,0 to 64,64 nohead lt 2 front" >> toto.gplot
echo "set arrow from 0,0 to 64,-64 nohead lt 2 front" >> toto.gplot
echo "set arrow from 0,0 to -64,-64 nohead lt 2 front" >> toto.gplot
echo "set arrow from 0,0 to -64,64 nohead lt 2 front" >> toto.gplot
echo "set xtics nomirror" >> toto.gplot
echo "set ytics nomirror" >> toto.gplot
echo "ciel0(x,y) = (h(x,y)>0) ? 0 : h(x,y)" >> toto.gplot
echo "ciel(x,y) =(x**2+y**2>=8100) ? 0 : ciel0(x,y)" >> toto.gplot
# cercle horizon
#1er cadran
echo "set arrow from 90,0 to 88,16 nohead lt -1 front" >> toto.gplot
echo "set arrow from 88,16 to 85,31 nohead lt -1 front" >> toto.gplot
echo "set arrow from 85,31 to 78,45 nohead lt -1 front" >> toto.gplot
echo "set arrow from 78,45 to 69,58 nohead lt -1 front" >> toto.gplot
echo "set arrow from 69,58 to 58,69 nohead lt -1 front" >> toto.gplot
echo "set arrow from 58,69 to 45,78 nohead lt -1 front" >> toto.gplot
echo "set arrow from 45,78 to 31,85 nohead lt -1 front" >> toto.gplot
echo "set arrow from 31,85 to 16,88 nohead lt -1 front" >> toto.gplot
echo "set arrow from 16,88 to 0,90 nohead lt -1 front" >> toto.gplot
#2e cadran
echo "set arrow from -90,0 to -88,16 nohead lt -1 front" >> toto.gplot
echo "set arrow from -88,16 to -85,31 nohead lt -1 front" >> toto.gplot
echo "set arrow from -85,31 to -78,45 nohead lt -1 front" >> toto.gplot
echo "set arrow from -78,45 to -69,58 nohead lt -1 front" >> toto.gplot
echo "set arrow from -69,58 to -58,69 nohead lt -1 front" >> toto.gplot
echo "set arrow from -58,69 to -45,78 nohead lt -1 front" >> toto.gplot
echo "set arrow from -45,78 to -31,85 nohead lt -1 front" >> toto.gplot
echo "set arrow from -31,85 to -16,88 nohead lt -1 front" >> toto.gplot
echo "set arrow from -16,88 to 0,90 nohead lt -1 front" >> toto.gplot
#3e cadran
echo "set arrow from -90,0 to -88,-16 nohead lt -1 front" >> toto.gplot
echo "set arrow from -88,-16 to -85,-31 nohead lt -1 front" >> toto.gplot
echo "set arrow from -85,-31 to -78,-45 nohead lt -1 front" >> toto.gplot
echo "set arrow from -78,-45 to -69,-58 nohead lt -1 front" >> toto.gplot
echo "set arrow from -69,-58 to -58,-69 nohead lt -1 front" >> toto.gplot
echo "set arrow from -58,-69 to -45,-78 nohead lt -1 front" >> toto.gplot
echo "set arrow from -45,-78 to -31,-85 nohead lt -1 front" >> toto.gplot
echo "set arrow from -31,-85 to -16,-88 nohead lt -1 front" >> toto.gplot
echo "set arrow from -16,-88 to 0,-90 nohead lt -1 front" >> toto.gplot
#4e cadran
echo "set arrow from 90,0 to 88,-16 nohead lt -1 front" >> toto.gplot
echo "set arrow from 88,-16 to 85,-31 nohead lt -1 front" >> toto.gplot
echo "set arrow from 85,-31 to 78,-45 nohead lt -1 front" >> toto.gplot
echo "set arrow from 78,-45 to 69,-58 nohead lt -1 front" >> toto.gplot
echo "set arrow from 69,-58 to 58,-69 nohead lt -1 front" >> toto.gplot
echo "set arrow from 58,-69 to 45,-78 nohead lt -1 front" >> toto.gplot
echo "set arrow from 45,-78 to 31,-85 nohead lt -1 front" >> toto.gplot
echo "set arrow from 31,-85 to 16,-88 nohead lt -1 front" >> toto.gplot
echo "set arrow from 16,-88 to 0,-90 nohead lt -1 front" >> toto.gplot
# cercle 60
#1er cadran
echo "set arrow from 60,0 to 59,11 nohead lt -1 front" >> toto.gplot
echo "set arrow from 59,11 to 57,21 nohead lt -1 front" >> toto.gplot
echo "set arrow from 57,21 to 52,30 nohead lt -1 front" >> toto.gplot
echo "set arrow from 52,30 to 46,39 nohead lt -1 front" >> toto.gplot
echo "set arrow from 46,39 to 39,46 nohead lt -1 front" >> toto.gplot
echo "set arrow from 39,46 to 30,52 nohead lt -1 front" >> toto.gplot
echo "set arrow from 30,52 to 21,57 nohead lt -1 front" >> toto.gplot
echo "set arrow from 21,57 to 11,59 nohead lt -1 front" >> toto.gplot
echo "set arrow from 11,59 to 0,60 nohead lt -1 front" >> toto.gplot
#2e cadran
echo "set arrow from -60,0 to -59,11 nohead lt -1 front" >> toto.gplot
echo "set arrow from -59,11 to -57,21 nohead lt -1 front" >> toto.gplot
echo "set arrow from -57,21 to -52,30 nohead lt -1 front" >> toto.gplot
echo "set arrow from -52,30 to -46,39 nohead lt -1 front" >> toto.gplot
echo "set arrow from -46,39 to -39,46 nohead lt -1 front" >> toto.gplot
echo "set arrow from -39,46 to -30,52 nohead lt -1 front" >> toto.gplot
echo "set arrow from -30,52 to -21,57 nohead lt -1 front" >> toto.gplot
echo "set arrow from -21,57 to -11,59 nohead lt -1 front" >> toto.gplot
echo "set arrow from -11,59 to 0,60 nohead lt -1 front" >> toto.gplot
#3e cadran
echo "set arrow from -60,0 to -59,-11 nohead lt -1 front" >> toto.gplot
echo "set arrow from -59,-11 to -57,-21 nohead lt -1 front" >> toto.gplot
echo "set arrow from -57,-21 to -52,-30 nohead lt -1 front" >> toto.gplot
echo "set arrow from -52,-30 to -46,-39 nohead lt -1 front" >> toto.gplot
echo "set arrow from -46,-39 to -39,-46 nohead lt -1 front" >> toto.gplot
echo "set arrow from -39,-46 to -30,-52 nohead lt -1 front" >> toto.gplot
echo "set arrow from -30,-52 to -21,-57 nohead lt -1 front" >> toto.gplot
echo "set arrow from -21,-57 to -11,-59 nohead lt -1 front" >> toto.gplot
echo "set arrow from -11,-59 to 0,-60 nohead lt -1 front" >> toto.gplot
#4e cadran
echo "set arrow from 60,0 to 59,-11 nohead lt -1 front" >> toto.gplot
echo "set arrow from 59,-11 to 57,-21 nohead lt -1 front" >> toto.gplot
echo "set arrow from 57,-21 to 52,-30 nohead lt -1 front" >> toto.gplot
echo "set arrow from 52,-30 to 46,-39 nohead lt -1 front" >> toto.gplot
echo "set arrow from 46,-39 to 39,-46 nohead lt -1 front" >> toto.gplot
echo "set arrow from 39,-46 to 30,-52 nohead lt -1 front" >> toto.gplot
echo "set arrow from 30,-52 to 21,-57 nohead lt -1 front" >> toto.gplot
echo "set arrow from 21,-57 to 11,-59 nohead lt -1 front" >> toto.gplot
echo "set arrow from 11,-59 to 0,-60 nohead lt -1 front" >> toto.gplot

# cercle 30
#1er cadran
echo "set arrow from 30,0 to 30,5 nohead lt -1 front" >> toto.gplot
echo "set arrow from 30,5 to 28,10 nohead lt -1 front" >> toto.gplot
echo "set arrow from 28,10 to 26,15 nohead lt -1 front" >> toto.gplot
echo "set arrow from 26,15 to 23,19 nohead lt -1 front" >> toto.gplot
echo "set arrow from 23,19 to 19,23 nohead lt -1 front" >> toto.gplot
echo "set arrow from 19,23 to 15,26 nohead lt -1 front" >> toto.gplot
echo "set arrow from 15,26 to 10,28 nohead lt -1 front" >> toto.gplot
echo "set arrow from 10,28 to 5,29 nohead lt -1 front" >> toto.gplot
echo "set arrow from 5,29 to 0,30 nohead lt -1 front" >> toto.gplot
#2e cadran
echo "set arrow from -30,0 to -29,5 nohead lt -1 front" >> toto.gplot
echo "set arrow from -29,5 to -28,10 nohead lt -1 front" >> toto.gplot
echo "set arrow from -28,10 to -26,15 nohead lt -1 front" >> toto.gplot
echo "set arrow from -26,15 to -23,19 nohead lt -1 front" >> toto.gplot
echo "set arrow from -23,19 to -19,23 nohead lt -1 front" >> toto.gplot
echo "set arrow from -19,23 to -15,26 nohead lt -1 front" >> toto.gplot
echo "set arrow from -15,26 to -10,28 nohead lt -1 front" >> toto.gplot
echo "set arrow from -10,28 to -5,29 nohead lt -1 front" >> toto.gplot
echo "set arrow from -5,29 to 0,30 nohead lt -1 front" >> toto.gplot
#3e cadran
echo "set arrow from -30,0 to -29,-5 nohead lt -1 front" >> toto.gplot
echo "set arrow from -29,-5 to -28,-10 nohead lt -1 front" >> toto.gplot
echo "set arrow from -28,-10 to -26,-15 nohead lt -1 front" >> toto.gplot
echo "set arrow from -26,-15 to -23,-19 nohead lt -1 front" >> toto.gplot
echo "set arrow from -23,-19 to -19,-23 nohead lt -1 front" >> toto.gplot
echo "set arrow from -19,-23 to -15,-26 nohead lt -1 front" >> toto.gplot
echo "set arrow from -15,-26 to -10,-28 nohead lt -1 front" >> toto.gplot
echo "set arrow from -10,-28 to -5,-29 nohead lt -1 front" >> toto.gplot
echo "set arrow from -5,-29 to 0,-30 nohead lt -1 front" >> toto.gplot
#4e cadran
echo "set arrow from 30,0 to 29,-5 nohead lt -1 front" >> toto.gplot
echo "set arrow from 29,-5 to 28,-10 nohead lt -1 front" >> toto.gplot
echo "set arrow from 28,-10 to 26,-15 nohead lt -1 front" >> toto.gplot
echo "set arrow from 26,-15 to 23,-19 nohead lt -1 front" >> toto.gplot
echo "set arrow from 23,-19 to 19,-23 nohead lt -1 front" >> toto.gplot
echo "set arrow from 19,-23 to 15,-26 nohead lt -1 front" >> toto.gplot
echo "set arrow from 15,-26 to 10,-28 nohead lt -1 front" >> toto.gplot
echo "set arrow from 10,-28 to 5,-29 nohead lt -1 front" >> toto.gplot
echo "set arrow from 5,-29 to 0,-30 nohead lt -1 front" >> toto.gplot
echo "set label 1 at -1, 94" >> toto.gplot
echo "set label 1 'N'" >> toto.gplot
echo "set label 2 at 93, 0" >> toto.gplot
echo "set label 2 'E'" >> toto.gplot
echo "set label 3 at -95, 0" >> toto.gplot
echo "set label 3 'W'" >> toto.gplot
echo "set label 4 at -1, -95" >> toto.gplot
echo "set label 4 'S'" >> toto.gplot
echo "set label 5 at 64, 68" >> toto.gplot
echo "set label 5 'NE'"  >> toto.gplot
echo "set label 6 at 64, -68" >> toto.gplot
echo "set label 6 'SE'"  >> toto.gplot
echo "set label 7 at -70, -68" >> toto.gplot
echo "set label 7 'SW'"  >> toto.gplot
echo "set label 8 at -70, 68" >> toto.gplot
echo "set label 8 'NW'" >> toto.gplot
echo "set label 9 at 10, -95" >> toto.gplot
echo "set label 9 'Log10(radiance [W/m^2/sr])'" >> toto.gplot
echo "splot ciel(x,y)" >> toto.gplot
gnuplot  < toto.gplot > $base".png"
# display $base".png" &
