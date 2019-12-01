c    programme pour ecrire des valeurs de radiance du ciel a en unite de ratio 
c    sur la radiance naturelle a partir d'un fichier texte yxz
c    
c    les axes sont definis comme suit
c    x = direction nord = degres
c    y = direction est = degres
c
c
c    Copyright (C) 2014  Martin Aube
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c    Contact: martin.aube@cegepsherbrooke.qc.ca
c -----------------
c   identification des variables 
c
c
c --------------------
c
c   programme principal
c
      program writeaod
      integer width                                                       ! Matrix dimension in Length/width
      parameter (width=256)
c
c ----------
c
c   declaration des variables
c
      real pixsiz,lat0,lon0,gain,offset,dat
      real datxy(width,width),dmax,s0,magmax,magmin
      real En
      character*72 bidon,datfile,nomxyz
      character*73 nombin
      character*12 nom
      integer i,j,lennom,nx,ny,xmin,ymin,n
      integer maxi,hcnt,valmax
c
c initialisation variables
c
      do i=1,1000
        do j=1,1000
           datxy(i,j)=0.
        enddo
      enddo

c
c -----------
c
c   choix du nom de la racine de fichiers
c

      open(unit=11,file='writevalue.in',status='old')
        read(11,*) nomxyz
        read(11,*) En                                             ! natural radiance in illumina units
        read(11,*) S0                                             ! S0 of the relation Ss=S0 - 2.5 log10(Em+En) where Em is the sky radiance as given by illumina,Ss is the measured surf. brightness in mag/sq arcsec and En is the natural sky radiance. You need a relative large number of data point to find S0 and En 
        read(11,*) pixsiz
      close(unit=11) 
c -----------------
c
c   Lecture des donnees
c
         cmax=0.
         xmin=width
         xmax=0
         ymin=width
         ymax=0

         open(unit=1,file=nomxyz,status='old')
         read(1,*) nsite
         do n=1,nsite
            read(1,*) i,j,dat
            if (dat.gt.dmax) dmax=dat
            if (i.gt.xmax) xmax=i
            if (j.gt.ymax) ymax=j
            if (i.lt.xmin) xmin=i
            if (j.lt.ymin) ymin=j    
         enddo
         rewind 1
         nx=xmax-xmin+1
         ny=ymax-ymin+1
         magmax=0.
         magmin=50.
         read(1,*) nsite
         do n=1,nsite
            read(1,*) i,j,dat

            datxy(i-xmin+1,j-ymin+1)=S0-2.5*log10(dat+En)
            if (datxy(i-xmin+1,j-ymin+1).gt.magmax) magmax=
     +      datxy(i-xmin+1,j-ymin+1)
            if (datxy(i-xmin+1,j-ymin+1).lt.magmin) magmin=
     +      datxy(i-xmin+1,j-ymin+1)

         enddo
        close(unit=1)
        print*,magmax,magmin
         
c
c ----------------
c
c   Ecriture de l image de sortie
c
            nombin='data.bin'
            nom='data   '
            valmax=65535
            offset=magmax
            gain=-1.*(magmax-magmin)/real(valmax)
            print*,magmin,magmax,gain,offset
            call extrants2d (nombin,datxy,nom,lat0,lon0,pixsiz,
     +      gain,offset,nx,ny,valmax)
       end
