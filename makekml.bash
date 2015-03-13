#!/bin/bash
# program to run in the directory containing PCL or PCW png files
# this program will organize the files for the web site and create
# the kml files
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
rm -f toto
pp="../"
list=`du -a | grep png`
for i in $list
do if [ -f $i ] 
   then f='IAC-2010/'
        folder='aubema\/recherches\/data\/IAC-2010\/'
        if [ ! -d $pp$f ]
        then mkdir $pp$f
        fi
        if [ `echo $i | grep ORM ` ]
        then folder=$folder"ORM\/"
             f=$f"ORM/"
        elif [ `echo $i | grep OT ` ]
        then folder=$folder"OT\/"
        f=$f"OT/"
        fi
        if [ ! -d $pp$f ]
        then        mkdir $pp$f
        fi
        
        if [ `echo $i | grep b_midnight ` ]
        then folder=$folder"Before_midnight\/"
        f=$f"Before_midnight/"
        elif [ `echo $i | grep a_midnight ` ]
        then folder=$folder"After_midnight\/"
        f=$f"After_midnight/"
        fi
        if [ ! -d $pp$f ]
        then        mkdir $pp$f
        fi
        
        if [ `echo $i | grep PCL ` ]
        then folder=$folder"RCM\/"
        f=$f"RCM/"
        elif [ `echo $i | grep PCW ` ]
        then folder=$folder"RSM\/"
        f=$f"RSM/"
        fi
        if [ ! -d $pp$f ]
        then        mkdir $pp$f
        fi

        if [ `echo $i | grep ta0.025 ` ]
        then folder=$folder"AOD0.025\/"
        f=$f"AOD0.025/"
        elif [ `echo $i | grep ta0.050 ` ] 
        then folder=$folder"AOD0.050\/"
        f=$f"AOD0.050/"
        elif [ `echo $i | grep ta0.100 ` ] 
        then folder=$folder"AOD0.100\/"
        f=$f"AOD0.100/"
        elif [ `echo $i | grep ta0.200 ` ]
        then folder=$folder"AOD0.200\/"
          f=$f"AOD0.200/"
        fi
        if [ ! -d $pp$f ]
        then        mkdir $pp$f
        fi
        if [ `echo $i | grep wl435 ` ]
        then folder=$folder"435nm\/"
        f=$f"435nm/"
        elif [ `echo $i | grep wl498 ` ] 
        then folder=$folder"498nm\/"
        f=$f"498nm/"
        elif [ `echo $i | grep wl546 ` ] 
        then folder=$folder"546nm\/"
        f=$f"546nm/"
        elif [ `echo $i | grep wl568 ` ]
        then folder=$folder"569nm\/"
        f=$f"569nm/"
        elif [ `echo $i | grep wl615 ` ]
        then folder=$folder"615nm\/"
         f=$f"615nm/"
        fi   
        if [ ! -d $pp$f ]
        then        mkdir $pp$f
        fi
        
        if [ `echo $i | grep PCL-ORM ` ] 
        then echo $i | sed -e 's/PCL-ORM/ ORM-RCM/g' | sed -e 's/_noneg//g' | sed 's/midnight/00/g' > toto 
        read a b < toto
        cp $i $pp$f$b
        c=`echo $b | sed 's/png/kml/g'`
        d=$folder$b

        
        cat generic.kml | sed -e "s/fichier/${d}/g" > $pp$f$c    
        fi
        if [ `echo $i | grep PCW-ORM ` ] 
        then echo $i | sed -e 's/PCW-ORM/ ORM-RSM/g' | sed -e 's/_noneg//g' | sed 's/midnight/00/g' > toto 
        read a b < toto
        cp $i $pp$f$b
        c=`echo $b | sed 's/png/kml/g'` 
        d=$folder$b
        
        cat generic.kml | sed -e "s/fichier/${d}/g" > $pp$f/$c
        fi
        if [ `echo $i | grep PCL-OT ` ] 
        then echo $i | sed -e 's/PCL-OT/ OT-RCM/g' | sed -e 's/_noneg//g' | sed 's/midnight/00/g' > toto 
        read a b < toto
        cp $i $pp$f$b
        c=`echo $b | sed 's/png/kml/g'`  
                d=$folder$b
        cat generic.kml | sed -e "s/fichier/${d}/g" > $pp$f$c
        fi
        if [ `echo $i | grep PCW-OT ` ] 
        then echo $i | sed -e 's/PCW-OT/ OT-RSM/g' | sed -e 's/_noneg//g' | sed 's/midnight/00/g' > toto 
        read a b < toto
        cp $i $pp$f$b
        c=`echo $b | sed 's/png/kml/g'`
                d=$folder$b
        cat generic.kml | sed -e "s/fichier/${d}/g" > $pp$f$c
        fi
   fi 
done
