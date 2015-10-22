set t png

do for [i=1:20] {
	outfile = sprintf("Lamps/zone%d_lamp.png",i)
	set output outfile
	plot "Lamps/zone".i."_lamp.dat" matrix with image title "Zone ".i
}

