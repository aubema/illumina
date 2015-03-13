c programme pour convertir un fichier pgm provenant directement 
c de viirs en valeur de 
c luminosite installee
c
       program viirs2lum
c
c L_viirs=F *(cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))
c L_viirs is the luminance detected by VIIRS
c F is the total lamp flux
c refl is the ground reflectance
c lop is the unnormalized light output pattern
c d_omeg is the incremental solid angle =  2.*pi*sin(z) dz
c 
c Equation above can be inverted to derive F_lamp
c F = L_viirs / (cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))
c
c We subdivise F_lamp into sodium and mercury lamps using the percentage from .zon file
c F_specific_lamp=percent/100.* F
c    Le % a pour but de
c    permettre la prise en compte d'une diminution de la luminance
c    lors de couvre feu.
c    Par ex si la luminosite est reduite de 50% apres
c    minuit on mettra une valeur de 50 car la luminance ols est detectee entre
c    19h30 et 21h30. On peut aussi utiliser ce % pour moduler le pourcentage de
c    la luminosite observee provenant d'un type de lampe responsable d'une raie 
c    en particulier.
c    La premiere ligne contient le nombre 
c    de zones, la valeur par defaut de down pour la ville, le % de la luminance
c    en reference a la luminance ols pour la ville et le fichier lop liee a 
c    la valeur par defaut pour la ville, le seuil ols pour la ville, la valeur 
c    par defaut de down pour la campagne, le % de la luminance
c    en reference a la luminance ols pour la campagne et le fichier lop liee a 
c    la valeur par defaut pour la campagne
c    La position est donnee en pixel. Le pixel en bas a gauche
c    a la coordonnee 1,1 et les valeur x augmentent vers la droite 
c    alors que les y augmentent vers le haut. le rayon est aussi donne 
c    pixel. Les zones ne sont pas additives une
c    nouvelle zone vient ecraser une precedente s il y a superposition
c    partielle ou totale des aires
c    l extension du fichier zone doit etre .zon et a pour nom le nom de 
c    base. Chaque zone resulte en un nouveau type de lampe comme 
c    ILLUMINA est limite a 99 type de lampes, cela implique qu'on ne 
c    peut definir plus de 97 zone car les zones 1 et 2 par defaut fait partie du
c    nombre. Par contre l'implementation du 2e digit pour le no de 
c    lampe ouvre la porte a lever cette limitation le moment venu
c   
c    Copyright (C) 2010  Martin Aube
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
c
c
       integer width
       parameter (width=1024)
       real refl(width,width),viirs(width,width),lum(width,width,99),down
       real pixsiz,gain,offset,r,x,y,dat(width,width)
       real percent,lop,pi,tot,f,theta,seuil,xcell0,ycell0,cityt
       real citysat,citytot,citydown,citycone,counsat,countot,coundown
       real councone,cone,sat
       integer nbx,nby,valmax,i,j,nzon,n
       character*72 reflex_file,ols_file,basename,zonfile
       character*72 outfile,city_lop,coun_lop
       character*72 lop_file
       character*12 nom
       character*3 lampno
       pi=3.14159
       valmax=65535
       gain=0.7
       offset=0.
c le seuil poru le masque water (0=land,1=water). En dessous de seuil 
c on pose la luminosite a 0
       seuil=0.5
       do i=1,width
       do j=1,width
       do n=1,9
         lum(i,j,n)=0.
       enddo
       enddo
       enddo
       cityt=0.
       print*,'Be shure to have an viirs illumina pgm file, a zonal '
       print*,'file, and a illumina pgm reflectance file.'
       print*,'Output root name of the experiment?'
       read*,basename
       print*,'viirs-dnb file name (format pgm illumina)?'
       read*,ols_file
       print*,'modis reflectance file name (format pgm illumina)?'
       read*,reflex_file
       print*,'zonal file name?'
       read*,zonfile
       lenbase=index(basename,' ')-1  
       nom='reflexion '
       call intrants2d(reflex_file,refl,nom,xcell0,ycell0,pixsiz,
     + nbx,nby)
       nom='VIIRS-DNB '
       call intrants2d(ols_file,viirs,nom,xcell0,ycell0,pixsiz,nbx,nby)
       open(unit=1,file=zonfile,status='old')
         read(1,*) nzon,citypercent,city_lop
         read(1,*) cityt,counpercent,coun_lop              
c cityt est le seuil du fichier viirs pour choisir entre ville et campagne
         open(unit=2,file=city_lop,status='old')
c reading city light output pattern - default value cities
           citytot=0.
           citydown=0.
           citycone=0.
           citysat=0.
           do i=1,181
              read(2,*) f,theta
              citytot=citytot+f*2.*pi*sin(pi*theta/180.)*pi/180.
              if (i.gt.90) then                       
                citydown=citydown+f*2.*pi*sin(pi*theta/180.)*pi/180.
              endif
              if (i.le.56) then
                citysat=citysat+f*2.*pi*sin(pi*theta/180.)*pi/180.          ! acceptance angle of VIIRS
                citycone=citycone+2.*pi*sin(pi*theta/180.)*pi/180.
              endif
           enddo
        print*,'tot=',citytot,'down=',citydown,'cone=',citycone,'sat='
     +  ,citysat
         close(unit=2)
            do i=1,nbx
             do j=1,nby
             if (viirs(i,j).ge.cityt) then
c F = L_viirs / (cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))
                lum(i,j,1)=citypercent/100.*viirs(i,j)/(refl(i,j)/pi*
     +          citydown/citytot+citysat/(citycone*citytot))
             endif
             enddo
            enddo
         open(unit=2,file=coun_lop,status='old')
c reading light output pattern - default value countryside
           countot=0.
           coundown=0.
           councone=0.
           counsat=0.
           do i=1,181
              read(2,*) f,theta
              countot=countot+f*2.*pi*sin(pi*theta/180.)*pi/180.
              if (i.gt.90) then                       
                coundown=coundown+f*2.*pi*sin(pi*theta/180.)*pi/180.
              endif
              if (i.le.56) then
                counsat=counsat+f*2.*pi*sin(pi*theta/180.)*pi/180.          ! acceptance angle of VIIRS
                councone=councone+2.*pi*sin(pi*theta/180.)*pi/180.
              endif
           enddo
         close(unit=2)
            do i=1,nbx
             do j=1,nby
             if (viirs(i,j).lt.cityt) then
c F = L_viirs / (cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))
                lum(i,j,2)=counpercent/100.*viirs(i,j)/(refl(i,j)/pi*
     +          coundown/countot+counsat/(councone*countot))
             endif
             enddo
            enddo
c other zones (circular zones)
         do n=1,nzon
            print*,'Zone #',n,'/',nzon
            read(1,*) x,y,r,percent,lop_file
         open(unit=2,file=lop_file,status='old')
c reading light output pattern - other zones
           tot=0.
           down=0.
           cone=0.
           sat=0.
           do i=1,181
              read(2,*) f,theta
              tot=tot+f*2.*pi*sin(pi*theta/180.)*pi/180.
              if (i.gt.90) then                       
                down=down+f*2.*pi*sin(pi*theta/180.)*pi/180.
              endif
              if (i.le.56) then
                sat=sat+f*2.*pi*sin(pi*theta/180.)*pi/180.                ! acceptance angle of VIIRS
                cone=cone+2.*pi*sin(pi*theta/180.)*pi/180.
              endif
           enddo
         close(unit=2)
      print*,'Down zone #',n,': ',down
         close(unit=2)
            do i=1,nbx
             do j=1,nby
              dist=sqrt((real(i)-x)**2.+(real(j)-y)**2.)
              if (dist.le.r) then
c F = L_viirs / (cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))
                lum(i,j,n+2)=percent/100.*viirs(i,j)/(refl(i,j)/pi*
     +          down/tot+sat/(cone*tot))
               do k=1,n+1
                 lum(i,j,k)=0.                                            ! Reset previous values inside the new zone
               enddo                                                      ! This allow the new zone to get full priority
              endif                                                       ! and overwrite previous zones 
             enddo
            enddo
         enddo
c write output  
       do n=1,nzon+2
        write(lampno, '(I3.3)' ) n
        outfile=basename(1:lenbase)//'_lumlp_'//lampno//'.pgm'
        nom='Luminosity'
        do i=1,nbx
         do j=1,nby
           dat(i,j)=lum(i,j,n)
         enddo
        enddo
        call extrants2d (outfile,dat,nom,xcell0,ycell0,pixsiz,
     +  gain,offset,nbx,nby,valmax)
       enddo
       stop
       end
