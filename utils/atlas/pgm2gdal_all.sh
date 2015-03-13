#!/usr/bin/env bash

function lpa_make_geotiff ()
{
	while read fname; do
		if [ ${fname##*.} = "pgm" ]; then
			lpa_pgm2tif.py $fname
		fi
	done
}


ls -1 . | lpa_make_geotiff
