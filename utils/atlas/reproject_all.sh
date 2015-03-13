#!/usr/bin/env bash

if  [ -z $1 ]
then
    echo "Proj4 string must be passed as first argument"
    exit 3
fi


function lpa_reproject ()
{
        while read fname; do
                if [ ${fname##*.} = "tif" ]; then
  			gdalwarp -s_srs $1 -t_srs "+init=epsg:900913"   $fname  ${fname/.tif/_epsg900913.tif}
                fi
        done
}


ls -1 . | lpa_reproject
