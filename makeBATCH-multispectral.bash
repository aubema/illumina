#!/bin/bash
#  Script to prepare a serie of modelling experiment with illumina on MS-II supercomputer
#
#
#    Copyright (C) 2010  Martin Aube
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# usage: makeBATCH-multispectral.bash parameter_file_path_and_name
#
#
# ===================================
#  reading variables for all experiments from makeBATCH.in file
if [ ! $1 ]
then echo "ERROR -- Please provide the path_to_makeBATCH.in"
     echo "usage: makeBATCH-multispectral.bash parameter_file_path_and_name"
     exit 0
fi
if [ ! -f $1 ]; then
   echo "ERROR -- No parameter file found at " $1
   exit 2
fi
read a outscpt <<< $(grep ^batch_file_name $1 | cut -f1 -d"#")
if [ $2 ]; then
	outscpt=$2
fi
read a pixsize <<< $(grep ^pixel_size $1 | cut -f1 -d"#")
read a latitu <<< $(grep ^latitudeN $1 | cut -f1 -d"#")
#read a nbx <<< $(grep ^dimension_WE $1 | cut -f1 -d"#")
#read a nby <<< $(grep ^dimension_SN $1 | cut -f1 -d"#")
read a exp_name <<< $(grep ^experiment_name $1 | cut -f1 -d"#")
read a pressure <<< $(grep ^pressure $1 | cut -f1 -d"#")
read a est_time <<< $(grep ^estimated_computing_time $1 | cut -f1 -d"#")
read a mna_file <<< $(grep ^terrain_elevation_file $1 | cut -f1 -d"#")
read a rh <<< $(grep ^relative_humidity $1 | cut -f1 -d"#")
read a aero_mod <<< $(grep ^aerosol_model $1 | cut -f1 -d"#")
read a cloud <<< $(grep ^cloud_model $1 | cut -f1 -d"#")
read a dmin <<< $(grep ^nearest_source_distance $1 | cut -f1 -d"#")
read a l1 <<< $(grep ^1_radius $1 | cut -f1 -d"#")
read a l3 <<< $(grep ^3_radius $1 | cut -f1 -d"#")
read a stoplim <<< $(grep ^stop_limit $1 | cut -f1 -d"#")
read a sx_sites <<< $(grep ^x_positions $1 | cut -f1 -d"#")
read a sy_sites <<< $(grep ^y_positions $1 | cut -f1 -d"#")
read a sz_sites <<< $(grep ^z_positions $1 | cut -f1 -d"#")
read a ssautdif <<< $(grep ^scattering_skip $1 | cut -f1 -d"#")
read a sr_dif <<< $(grep ^scattering_radius $1 | cut -f1 -d"#")
read a selevation  <<< $(grep ^elevation_angles $1 | cut -f1 -d"#")
read a sazimut <<< $(grep ^azimuth_angles $1 | cut -f1 -d"#")
read a stau <<< $(grep ^aerosol_optical_depth $1 | cut -f1 -d"#")
read a salpha <<< $(grep ^angstrom_coefficients $1 | cut -f1 -d"#")
x_sites=( $sx_sites )
y_sites=( $sy_sites )
z_sites=( $sz_sites )
saut_dif=( $ssautdif )
r_dif=( $sr_dif )
elevation=( $selevation )
azimut=( $sazimut )
tau=( $stau )
alpha=( $salpha )
#
# =========================
# Vectors determined from local files with .lst extensions
#
i=0; for line in $(<wav.lst); do wav[i]="$line";let i=i+1;done
i=0; for line in $(<zon.lst); do lamp_l[i]="$line";let i=i+1;done
n=0; while [ $n -lt ${#wav[*]} ] ; do wav[$n]="$(echo -e "${wav[$n]}" | tr -d ' ')" ; let n=n+1; done
n=0; while [ $n -lt ${#lamp_l[*]} ] ; do lamp_l[$n]="$(echo -e "${lamp_l[$n]}" | tr -d ' ')" ; let n=n+1; done
#
# ===================================
# begin of script
#
# find domain size and resolution
#grep -v "#" $mna_file | head -2 |tail -1 > size.tmp
#read nbx nby bidon < size.tmp
#grep "pixsiz" $mna_file > size.tmp
#read bidon bidon pixsiz bidon < size.tmp
od -j4 -N8 -tdI -An $mna_file > size.tmp
read nbx nby < size.tmp
#
# find max height
find_max_height
read maxh bidon < maxheight.tmp

nmam=0
folder=`pwd`
griddir=$folder"/gridmap"
rm -f $HOME/$outscpt*
rm -fr $griddir
mkdir $griddir
echo $griddir
#
# accelerate computation from a subset of initial lumlp file data
ns=0
while [ $ns -lt ${#x_sites[*]} ]
do xsit=${x_sites[$ns]}
   let xdist=nbx-xsit
   if [ ${x_sites[$ns]} -lt $xdist ]
   then let xdist=xsit
   fi
   ysit=${y_sites[$ns]}
   let ydist=nby-ysit
   if [ $ysit -lt $ydist ]
   then let ydist=ysit
   fi
   let dist=xdist
   if [ $ydist -lt $xdist ]
   then let dist=ydist
   fi
   min_elev=`/bin/echo "180.0/3.14159*a("$maxh"/("$pixsize"*" $dist"))" | /usr/bin/bc -l`
   min_elev=`/bin/echo "("$min_elev"+0.5)/1" | /usr/bin/bc `
   echo "Calculations made at elevation angles lower than " $min_elev " deg may be less accurate in some cases."
#   let min_elev=min_elev/2
   let min_elev=1
   mkdir $griddir/"x"${x_sites[$ns]}"y"${y_sites[$ns]}
   echo $griddir
   echo "Site: x"${x_sites[$ns]}"y"${y_sites[$ns]}
   wl=0
   while [ $wl -lt ${#wav[*]} ]
   do echo "Wavelength: " ${wav[$wl]}
             mkdir $griddir/"x"${x_sites[$ns]}"y"${y_sites[$ns]}/"wl"${wav[$wl]}


             cd $griddir/"x"${x_sites[$ns]}"y"${y_sites[$ns]}/"wl"${wav[$wl]}/
             mkdir original_files
             nl=0
             while [ $nl -lt ${#lamp_l[*]} ]
             do let nolp=nl+1
                cp $folder"/"$exp_name"_"${wav[$wl]}"_lumlp_"${lamp_l[$nl]}".bin" "./"$exp_name"_lumlp_"${lamp_l[$nl]}".bin"
                let nl=nl+1
             done
             cp -f *lumlp* original_files
             cp -f $folder"/modis_"${wav[$wl]}".bin" "./"$exp_name"_reflect.bin"
             cp -f *_reflect* original_files



             # combine lumlp files
echo "Combining lumlp files..."
             listlp=`ls -1 *lumlp*.bin`
             ls *lumlp*.bin | grep -c "" > comblp.in
             echo $exp_name"_lumlp_recombined.bin" >> comblp.in
             for nlp in $listlp
             do echo $nlp >> comblp.in
             done
             cp -f $HOME/hg/illumina/bin/combine .
             ./combine < comblp.in
echo "Making variable resolution grid lumlp files..."
             nl=0
             while [ $nl -lt ${#lamp_l[*]} ]
             do echo $exp_name"_lumlp_"${lamp_l[$nl]}".bin"  > varres.in
                echo $exp_name"_reflect.bin" >> varres.in
                echo $exp_name"_lumlp_"${lamp_l[$nl]}"_new.bin"  >> varres.in
                echo $exp_name"_reflect_new.bin"  >> varres.in
                echo ${x_sites[$ns]} >> varres.in
                echo ${y_sites[$ns]} >> varres.in
                echo $l1 $l3 >> varres.in
                cp -f $HOME/hg/illumina/bin/varres .
                ./varres
                let nl=nl+1
             done

# creating the variable resolution reflectance
# this reflectance is a weighted average with the lumlp files on the variable res grid
                echo $exp_name"_lumlp_recombined.bin"  > varres.in
                echo $exp_name"_reflect.bin" >> varres.in
                echo $exp_name"_lumlp_"${lamp_l[$nl]}"_new.bin"  >> varres.in
                echo $exp_name"_reflect_new.bin"  >> varres.in
                echo ${x_sites[$ns]} >> varres.in
                echo ${y_sites[$ns]} >> varres.in
                echo $l1 $l3 >> varres.in
                cp -f $HOME/hg/illumina/bin/varres .
                ./varres

             listbin=`ls -1 *_new.bin`
             for ibin in $listbin
             do obin=`echo $ibin | sed 's/_new//g'`
                mv -f $ibin $obin
             done
       let wl=wl+1
   done
   let ns=ns+1
done
#
#
#
cd $folder
echo "Starting from "$folder
      ncas=0
      ns=0
      while [ $ns -lt ${#x_sites[*]} ]
      do mkdir "x"${x_sites[$ns]}"y"${y_sites[$ns]}
# echo "Observer position " $ns "("${x_sites[$ns]}"," ${y_sites[$ns]}"," ${z_sites[$ns]}")"

         cd "x"${x_sites[$ns]}"y"${y_sites[$ns]}
         here=`pwd`
         echo "Entering "$here
               nsd=0
               while [ $nsd -lt ${#saut_dif[*]} ]
               do mkdir "sd"${saut_dif[$nsd]}
# 2nd scat. acceleration factor
                  cd "sd"${saut_dif[$nsd]}
                  here=`pwd`
                  echo "Entering "$here
                  nrd=0
                  while [ $nrd -lt ${#r_dif[*]} ]
                  do mkdir "rd"${r_dif[$nrd]}
# 2nd scat. radius
                     cd "rd"${r_dif[$nrd]}
                     here=`pwd`
                     echo "Entering "$here
                     ne=0
                     while [ $ne -lt ${#elevation[*]} ]
                     do if [ ${elevation[$ne]} -ge $min_elev ]
                        then mkdir "el"${elevation[$ne]}
# Elevation viewing angle
                             cd "el"${elevation[$ne]}
                             here=`pwd`
                             echo "Entering "$here
                             if [ ${elevation[$ne]} -eq 90 ]
                             then natot=1
# Pointing toward zenith
                             else natot=${#azimut[*]}
                             fi
                             na=0
                             while [ $na -lt $natot ]
                             do mkdir "az"${azimut[$na]}
# Azimuth viewing angle
                                cd "az"${azimut[$na]}
                                here=`pwd`
                                echo "Entering "$here
                                nt=0
                                while [ $nt -lt ${#tau[*]} ]
                                do mkdir "ta"${tau[$nt]}
# Aerosol Optical Depth
                                   cd "ta"${tau[$nt]}
                                   here=`pwd`
                                   echo "Entering "$here
                                   wl=0
                                   while [ $wl -lt ${#wav[*]} ]
                                   do mkdir "wl"${wav[$wl]}
# wavelength
                                      cd "wl"${wav[$wl]}
                                      here=`pwd`
                                      echo "Entering "$here
                                      ln -s $folder"/"*"_RH"*${wav[$wl]}*".mie.out" ./
                                      nlamp=${#lamp_l[*]}
#
# creating illumina.in
                                      echo "                ! Input file for ILLUMINA" > illumina.in
                                      echo  $exp_name "              ! ROOT FILE NAME (every usefull files have to begin with this)" >> illumina.in
#                                      echo $nbx $nby  "         ! Modeling domain size W-E and S-N in pixel " >> illumina.in
                                      echo $pixsize $pixsize "                ! CELL SIZE [m]" >> illumina.in
#                                      echo $latitu "   ! central latitude north of the domain" >> illumina.in




#################
                                      echo $aero_mod"_RH"$rh"_0."${wav[$wl]}"0um.mie.out               ! AEROSOL OPTICAL CROSS SECTIONS FILE [-]" >> illumina.in
                                      echo "                !" >> illumina.in
                                      echo ${r_dif[$nrd]} ${saut_dif[$nsd]} "   ! DOUBLE SCATTERING RADIUS [m] ; SCATTERING STEP [-]" >> illumina.in
                                      echo "                 ! (1=complete, 2= 2 times faster, ...)" >> illumina.in
                                      echo ${wav[$wl]} "  ! WAVELENGTH [nm]" >> illumina.in
                                      echo $pressure "   ! GROUND LEVEL PRESSURE [kPa]"  >> illumina.in
                                      echo ${tau[$nt]} ${alpha[$nt]} "   ! 500nm AEROSOL OPTICAL DEPTH ; angstrom exponent [-]"  >> illumina.in
                                      echo ${#lamp_l[*]} "    ! NUMBER OF SOURCE TYPES [-]"  >> illumina.in
                                      echo $stoplim "                ! Contribution threshold : Stop computation when the new voxel contribution is less than 1/threshold of the cumulated flux (suggested value = 5000.)" >> illumina.in
                                      echo "                !" >> illumina.in
                                      echo ${x_sites[$ns]} ${y_sites[$ns]} ${z_sites[$ns]}  " 1   ! OBSERVER X POSITION [cell unit] ; OBS Y POS [-] ; OBS elevation above the local ground level [m] ; BEGINNING CELL ALONG THE LINE OF SIGHT (1=complete)" >> illumina.in
                                      echo "    ! (usefull in case of computer crash) (x=1 et y=1 is the south-west cell)" >> illumina.in
                                      echo ${elevation[$ne]} ${azimut[$na]} "   ! ELEVATION VIEWING ANGLE [deg] ; AZIMUTAL VIEWING ANGLE [deg] (0=east, 90=north, 180=west, 270=south)" >> illumina.in
                                      echo "   !   " >> illumina.in
                                      echo ".0001 .00001 .35 .07    ! SLIT WIDTH [m] ; PIXEL SIZE [m] ; FOCAL LENGTH [m] ; APPERTURE DIAMETER [m]" >> illumina.in
                                      echo "                        ! to get the radiance in W/str/m**2, SW*PS*pi*AD**2/4/FL**2 should equal 1" >> illumina.in
                                      echo "                        ! 1. 1. 1. 1.1283791671 is doing exactly that">> illumina.in
                                      echo $cloud "                 ! Cloud model selection 0=clear, 1=Thin Cirrus/Cirrostratus, 2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus, 4=Cumulus/Cumulonimbus, 5=Stratocumulus">> illumina.in
                                      echo $dmin "          ! Minimal distance to the nearest light source (m)">> illumina.in
# copie de illumina dans ./
                                      ln -s $HOME/hg/illumina/bin/illumina .
#
# creation du script pour qsub dans $folder
                                      ici=`pwd`
# create batch run files for 300 jobs blocks
                                      let nmam=nmam+1
                                      let jobno=nmam/300+1
                                      outscript=$outscpt"_"$jobno
                                      echo "cd " $ici >> $HOME/$outscript
                                      #echo "qsub -W umask=0011 -q qwork@ms  ./execute" >> $HOME/$outscript
				      echo "sbatch ./execute" >> $HOME/$outscript
                                      echo "sleep 0.05"  >> $HOME/$outscript
                                      echo "#!/bin/sh" > $ici/execute
                                      #echo "#PBS -l cput="$est_time":00:00" >> $ici/execute
                                      #echo "#PBS -l walltime="$est_time":00:00" >> $ici/execute
                                      #echo "#PBS -m bea" >> $ici/execute
                                      #echo "#PBS -M aubema@gmail.com" >> $ici/execute
                                      echo "#SBATCH --job-name=Illumina" >> $ici/execute
                                      echo "#SBATCH --time="$est_time":00:00" >> $ici/execute
                                      echo "#SBATCH --mem-per-cpu=1920" >> $ici/execute
                                      echo "cd " $ici >> $ici/execute
				      echo "umask 0011" >> $ici/execute
                                      echo "./illumina" >> $ici/execute
                                      chmod u+x $ici/execute
                                      #echo "cd " $ici >> $HOME/$outscpt
                                      #echo "bqsub -P \"command=./illumina >illumina.out\" -q qwork@ms -l walltime="$est_time":00:00" >> $HOME/$outscpt
                                      nl=0
                                      while [ $nl -lt ${#lamp_l[*]} ]
                                      do let nolp=nl+1
                                         ln -s $folder"/fctem_wl_"${wav[$wl]}"_zon_"${lamp_l[$nl]}".dat" "./"$exp_name"_fctem_"${lamp_l[$nl]}".dat"
                                         let nl=nl+1
                                      done
                                      ln -s $folder"/"$exp_name"_altlp.bin" "./"$exp_name"_altlp.bin"
                                      ln -s $folder"/"$exp_name"_obsth.bin" "./"$exp_name"_obsth.bin"
                                      ln -s $folder"/"$exp_name"_obstd.bin" "./"$exp_name"_obstd.bin"
                                      ln -s $folder"/"$exp_name"_obstf.bin" "./"$exp_name"_obstf.bin"
                                      ln -s $griddir/"x"${x_sites[$ns]}"y"${y_sites[$ns]}/"wl"${wav[$wl]}/* .
# copy other bin files
                                      ln -s $folder/$mna_file "./"$exp_name"_topogra.bin"
                                      cd ..
                                      let ncas=ncas+1
                                      let wl=wl+1
                                   done
                                   cd ..
                                   let nt=nt+1
                                done
                                cd ..
                                let na=na+1
                             done
                             cd ..
                        else     echo ${elevation[$ne]} "is lower than the minimal elevation angle allowed for your modelling geometry (min="$min_elev")"

                        fi     #end of the lowest allowed angle verification
                             let ne=ne+1

                     done
                     cd ..
                     let nrd=nrd+1
                  done
                  cd ..
                  let nsd=nsd+1
               done
               cd ..
               let ns=ns+1
      done

echo "Total number of runs:" $ncas
let timetot=est_time*ncas
echo "Total estimated CPU time required:" $timetot "Hour"
cd $folder
