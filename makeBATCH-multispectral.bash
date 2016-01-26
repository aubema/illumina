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
# usage makeBATCH output_script
#
#
# ===================================
#  Fixed variables for all experiments
#

if [ ! $1 ]
then outscpt="TortureMammouth"
else outscpt=$1
fi
pixsize=1000                                                   # Size of the pixel
exp_name="Hawaii"                                              # Base name of the experiment
pressure="101.3"                                               # Atm pressure of the lowest point
est_time="110"                                                 # Estimated number of hour per individual run
mna_file="srtm.pgm"                                            # Numrerical terrain model file name
rh="70"                                                        # Relative humidity used for mie.out files
aero_mod="maritime"                                            # aerosol model for the .mie.out files (urban, rural or maritime)
cloud=0                                                        # Cloud model selection 0=clear, 1=Thin Cirrus/Cirrostratus, 
                                                               # 2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus, 
                                                               # 4=Cumulus/Cumulonimbus, 5=Stratocumulus
dmin="150"                                                     # Minimal distance to the nearest light source (m)
#
# ===================================
#  Vectors
#
# sites in order:  ORM OT
x_sites=( 269 )                                                # list of x observer positions 
y_sites=( 245 )                                                # list of y observer positions (with respect to x_sites)
z_sites=(  38  )                                                # list of z observer positions (with respect to x_sites)
                                                               # 0   0   0   0 0.04 .68 5.9  37 164 1141 2000   94  0   equivalent elevations (m)
d_reflect=( 25 )                                               # list of Mean light free path toward ground 
h_obstacle=( 7 )                                               # list of subgrid obstacle height
saut_dif=( 71 )                                                # list of 2nd scattering computing acceleration factor (ideally a prime number
r_dif=( 4000 )                                                 # list of 2nd scattering radius
elevation=( 90 )                              # list of elevation viewing angles  
azimut=( 0  )                         # list of azimut viewing angles
tau=( 0.11 )                                                # list of AOD values at 500 nm
alpha=( 0.7 )                                       # list of angstrom exponents values
#
# =========================
# Vectors determined from local files with .lst extensions
#
i=0; for line in $(<wav.lst); do wav[i]="$line";let i=i+1;done
i=0; for line in $(<zon.lst); do lamp_l[i]="$line";let i=i+1;done
# i=0; for line in $(<altlp.lst); do lamp_h[i]="$line";let i=i+1;done
# removing unwanted white spaces in arrays
n=0; while [ $n -lt ${#wav[*]} ] ; do wav[$n]="$(echo -e "${wav[$n]}" | tr -d ' ')" ; let n=n+1; done
n=0; while [ $n -lt ${#lamp_l[*]} ] ; do lamp_l[$n]="$(echo -e "${lamp_l[$n]}" | tr -d ' ')" ; let n=n+1; done
# n=0; while [ $n -lt ${#lamp_h[*]} ] ; do lamp_h[$n]="$(echo -e "${lamp_h[$n]}" | tr -d ' ')" ; let n=n+1; done
#
# ===================================
# begin of script
#
# find domain size and resolution
grep -v "#" $mna_file | head -2 |tail -1 > size.tmp
read nbx nby bidon < size.tmp
grep "pixsiz" $mna_file > size.tmp
read bidon bidon pixsiz bidon < size.tmp
#
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
   min_elev=`/bin/echo "180.0/3.14159*a(27345.16/("$pixsiz"*" $dist"))" | /usr/bin/bc -l`
   min_elev=`/bin/echo "("$min_elev"+0.5)/1" | /usr/bin/bc `
   echo "min_elev=" $min_elev
   let min_elev=min_elev/2
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
                cp $folder"/"$exp_name"_"${wav[$wl]}"_lumlp_"${lamp_l[$nl]}".pgm" "./"$exp_name"_lumlp_"${lamp_l[$nl]}".pgm"
                let nl=nl+1
             done
             cp -f *lumlp* original_files
             cp -f $folder"/modis_"${wav[$wl]}".pgm" "./"$exp_name"_reflect.pgm"
             cp -f *_reflect* original_files



             # combine lumlp files
echo "Combining lumlp files..."
             listlp=`ls -1 *lumlp*.pgm`
             ls *lumlp*.pgm | grep -c "" > comblp.in
             echo $exp_name"_lumlp_recombined.pgm" >> comblp.in
             for nlp in $listlp
             do echo $nlp >> comblp.in
             done             
             cp -f $HOME/hg/illumina/bin/pgmcombine16bit .
             ./pgmcombine16bit < comblp.in
echo "Making variable resolution grid lumlp files..."
             nl=0
             while [ $nl -lt ${#lamp_l[*]} ]
             do echo $exp_name"_lumlp_"${lamp_l[$nl]}".pgm"  > varres.in
                echo $exp_name"_reflect.pgm" >> varres.in
                echo $exp_name"_lumlp_"${lamp_l[$nl]}"_new.pgm"  >> varres.in
                echo $exp_name"_reflect_new.pgm"  >> varres.in
                echo ${x_sites[$ns]} >> varres.in
                echo ${y_sites[$ns]} >> varres.in 
                cp -f $HOME/hg/illumina/bin/varres .
                ./varres 
                let nl=nl+1
             done

# creating the variable resolution reflectance
# this reflectance is a weighted average with the lumlp files on the variable res grid
                echo $exp_name"_lumlp_recombined.pgm"  > varres.in
                echo $exp_name"_reflect.pgm" >> varres.in
                echo $exp_name"_lumlp_"${lamp_l[$nl]}"_new.pgm"  >> varres.in
                echo $exp_name"_reflect_new.pgm"  >> varres.in
                echo ${x_sites[$ns]} >> varres.in
                echo ${y_sites[$ns]} >> varres.in 
                cp -f $HOME/hg/illumina/bin/varres .
                ./varres 
             
             listpgm=`ls -1 *_new.pgm`
             for ipgm in $listpgm
             do opgm=`echo $ipgm | sed 's/_new//g'`
                mv -f $ipgm $opgm
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
         nh=0
         while [ $nh -lt ${#h_obstacle[*]} ]
         do mkdir "ho"${h_obstacle[$nh]}  
# Obstacle height
            cd  "ho"${h_obstacle[$nh]} 
            here=`pwd`
            echo "Entering "$here
            ndr=0
            while [ $ndr -lt ${#d_reflect[*]} ]
            do mkdir "ro"${d_reflect[$ndr]}
# Mean free path
               cd "ro"${d_reflect[$ndr]}
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
                                      echo  $exp_name "     ! ROOT FILE NAME (every usefull files have to begin with this)" >> illumina.in
                                      echo "                !" >> illumina.in
                                      echo $pixsize $pixsize "                ! CELL SIZE [m]" >> illumina.in
                                      echo "                !" >> illumina.in




#################
                                      echo $aero_mod"_RH"$rh"_0."${wav[$wl]}"0um.mie.out               ! AEROSOL OPTICAL CROSS SECTIONS FILE [-]" >> illumina.in
                                      echo "                !" >> illumina.in
                                      echo ${r_dif[$nrd]} ${saut_dif[$nsd]} "   ! DOUBLE SCATTERING RADIUS [m] ; SCATTERING STEP [-]" >> illumina.in
                                      echo "                 ! (1=complete, 2= 2 times faster, ...)" >> illumina.in 
                                      echo ${wav[$wl]} "  ! WAVELENGTH [nm]" >> illumina.in
                                      echo $pressure "   ! GROUND LEVEL PRESSURE [kPa]"  >> illumina.in
                                      echo ${tau[$nt]} ${alpha[$nt]} "   ! 500nm AEROSOL OPTICAL DEPTH ; angstrom exponent [-]"  >> illumina.in
                                      echo ${#lamp_l[*]} "    ! NUMBER OF SOURCE TYPES [-]"  >> illumina.in
                                      echo "                !" >> illumina.in
                                      echo "                !" >> illumina.in
                                      echo ${x_sites[$ns]} ${y_sites[$ns]} ${z_sites[$ns]}  " 1   ! OBSERVER X POSITION [cell unit] ; OBS Y POS [-] ; OBS Z POS [-] ; BEGINNING CELL ALONG THE LINE OF SIGHT (1=complete)" >> illumina.in
                                      echo "    ! (usefull in case of computer crash) (x=1 et y=1 is the south-west cell)" >> illumina.in
                                      echo ${elevation[$ne]} ${azimut[$na]} "   ! ELEVATION VIEWING ANGLE [deg] ; AZIMUTAL VIEWING ANGLE [deg]" >> illumina.in
                                      echo "   ! (0=east, 90=north, 180=west, 270=south)" >> illumina.in
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
# create batch run files for 500 jobs blocks
                                      let nmam=nmam+1
                                      let jobno=nmam/500+1
                                      outscript=$outscpt"_"$jobno
                                      echo "cd " $ici >> $HOME/$outscript
                                      echo "qsub  -q qwork@ms  ./execute" >> $HOME/$outscript
                                      echo "sleep 0.05"  >> $HOME/$outscript
                                      echo "#!/bin/csh" > $ici/execute
                                      echo "#PBS -l cput="$est_time":00:00" >> $ici/execute
                                      echo "#PBS -l walltime="$est_time":00:00" >> $ici/execute
                                      echo "#PBS -m bea" >> $ici/execute
#                                 echo "#PBS -M aubema@gmail.com" >> $ici/execute
                                      echo "cd " $ici >> $ici/execute
                                      echo "./illumina" >> $ici/execute
                                      chmod u+x $ici/execute
#                                 echo "cd " $ici >> $HOME/$outscpt
#                                 echo "bqsub -P \"command=./illumina >illumina.out\" -q qwork@ms -l walltime="$est_time":00:00" >> $HOME/$outscpt
                                      nl=0        
                                      while [ $nl -lt ${#lamp_l[*]} ]
                                      do let nolp=nl+1
#                                         if [ $nolp -lt 1000 ]
#                                         then numlp=$nolp
#                                         fi
#                                         if [ $nolp -lt 100 ]
#                                         then numlp="0"$nolp
#                                         fi
#                                         if [ $nolp -lt 10 ]
#                                         then numlp="00"$nolp
#                                         fi                    
#                                         ln -s $folder"/"${lamp_h[$nl]} "./"$exp_name"_altlp_"${lamp_l[$nl]}".pgm"
                                         ln -s $folder"/"$exp_name"_altlp_"${lamp_l[$nl]}".pgm" "./"$exp_name"_altlp_"${lamp_l[$nl]}".pgm"

                                         ln -s $folder"/fctem_wl_"${wav[$wl]}"_zon_"${lamp_l[$nl]}".dat" "./"$exp_name"_fctem_"${lamp_l[$nl]}".dat"

                                         ln -s $folder"/"$exp_name"_obsth_"${lamp_l[$nl]}".pgm" "./"$exp_name"_obsth_"${lamp_l[$nl]}".pgm"
                                         ln -s $folder"/"$exp_name"_obstd_"${lamp_l[$nl]}".pgm" "./"$exp_name"_obstd_"${lamp_l[$nl]}".pgm"
                                         let nl=nl+1
                                      done
                                      ln -s $griddir/"x"${x_sites[$ns]}"y"${y_sites[$ns]}/"wl"${wav[$wl]}/* .
# copie des autres fichiers d intrants
                                      ln -s $folder/$mna_file "./"$exp_name"_topogra.pgm"
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
               let ndr=ndr+1
            done
            cd ..
            let nh=nh+1
         done
         cd ..
         let ns=ns+1
      done 

echo "Total number of runs:" $ncas
let timetot=est_time*ncas
echo "Total estimated CPU time required:" $timetot "Hour"
cd $folder
