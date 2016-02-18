#!/bin/sh
gnuplot << EOF

set t png size 1280,960 font arial 24

nb = $1

set output "Lamps/zone1_lamp.png"
plot "Lamps/zone1_lamp.bin" binary matrix with image
set xr [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
set yr [GPVAL_DATA_Y_MAX:GPVAL_DATA_Y_MIN]

set xlabel "Wavelength (nm)"
set ylabel "Zenith angle (deg)"

do for [i=1:nb] {
	outfile = sprintf("Lamps/zone%d_lamp.png",i)
	set output outfile
	set title "Zone ".i
	plot "Lamps/zone".i."_lamp.bin" binary matrix with image
}

EOF
