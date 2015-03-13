#! /bin/bash

# Set South West corner for OMM domain

for PGM in `ls -1 *.pgm`; do
   mv $PGM _$PGM
   sed -e "2c # lat0 45" -e "3c # lon0 -75" _$PGM > $PGM
   rm _$PGM
done
