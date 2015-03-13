c programme pour convertir un fichier pgm provenant directement 
c de DMSP-OLS en luminance non calibrees de 0-63 en valeur de 
c luminosite installee
c
       program ols2lum
c le principe est de tenir compte de la reflectance (refl), de la fonction
c d'emission de la lampe pour z=0 (lop(0)) et du  
c pourcentage de uplight (pup) pour passer d une intensite upward (lup)
c a une luminosite de la lampe (lum). lup=lum percent/100. (1./pi*(1.-pup)*cos(0)*refl+lop(0))
c Une fois inversee on obtiens 
c  lum=percent/100.*lup/(1./pi*(1.-pup)*refl+lop(0))
c    le pup et lop sont definis par zones. les zones sont lues dans un fichier
c    contenant a chaque ligne la position centrale de la zone (circulaire)
c    son rayon, la valeur de pup, le % de la luminance en reference  
c    a ols et le fichier lop. Le % a pour but de
c    permettre la prise en compte d'une diminution de la luminance
c    lors de couvre feu.
c    Par ex si la luminosite est reduite de 50% apres
c    minuit on mettra une valeur de 50 car la luminance ols est detectee entre
c    19h30 et 21h30. On peut aussi utiliser ce % pour moduler le pourcentage de
c    la luminosite observee provenant d'un type de lampe responsable d'une raie 
c    en particulier.
c    La premiere ligne contient le nombre 
c    de zones, la valeur par defaut de pup pour la ville, le % de la luminance
c    en reference a la luminance ols pour la ville et le fichier lop liee a 
c    la valeur par defaut pour la ville, le seuil ols pour la ville, la valeur 
c    par defaut de pup pour la campagne, le % de la luminance
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
       real refl(width,width),lup(width,width),lum(width,width,99),pup
       real pixsiz,gain,offset,r,x,y,dat(width,width)
       real water(width,width)
       real percent,lop,pi,tlop,f,theta,seuil,xcell0,ycell0,cityt
       real citylop,citytlop,counlop,countlop,citypup,counpup
       integer nbx,nby,valmax,i,j,nzon,n
       character*72 reflex_file,ols_file,basename,zonfile
       character*72 outfile,city_lop,coun_lop,eau_file
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
       print*,'Be shure to have an ols illumina pgm file, a zonal %'
       print*,'uplight, and a illumina pgm reflectance file.'
       print*,'Output root name of the experiment?'
       read*,basename
       print*,'Fichier ols (format pgm illumina)?'
       read*,ols_file
       print*,'Fichier reflectance (format pgm illumina)?'
       read*,reflex_file
       print*,'Fichier masque eau (format pgm illumina)?'
       read*,eau_file
       print*,'Fichier zone?'
       read*,zonfile
       lenbase=index(basename,' ')-1  
c       reflex_file=basename(1:lenbase)//'_reflect.pgm'

c       zonfile=basename(1:lenbase)//'.zon'
       nom='reflexion '
       call intrants2d(reflex_file,refl,nom,xcell0,ycell0,pixsiz,
     + nbx,nby)
       nom='DMSP-OLS  '
       call intrants2d(ols_file,lup,nom,xcell0,ycell0,pixsiz,nbx,nby)
       nom='WaterMask  '
       call intrants2d(eau_file,water,nom,xcell0,ycell0,pixsiz,nbx,nby)

       open(unit=1,file=zonfile,status='old')
         read(1,*) nzon,citypercent,city_lop
         read(1,*) cityt,counpercent,coun_lop              
c cityt est le seuil du fichier ols pour choisir entre ville et campagne
         open(unit=2,file=city_lop,status='old')
c reading and normalizing city light output pattern - default value cities
           read(2,*) citylop
           rewind 2
           citytlop=0.
           citypup=0.
           do i=1,181
              read(2,*) f,theta
              citytlop=citytlop+f*2.*pi*sin(pi*theta/180.)*pi/180.
              if (i.le.90) then
                citypup=citypup+f*2.*pi*sin(pi*theta/180.)*pi/180.
              endif
           enddo
           citylop=citylop/citytlop
           citypup=citypup/citytlop
      print*,'city up:',citypup
         close(unit=2)
            do i=1,nbx
             do j=1,nby
             lum(i,j,1)=0.
c filter out water surface where there is not suppose to be any lights
c this feature have been added to correct the overspill problem of ols data
         if (lup(i,j).ge.cityt) then
              if (water(i,j).gt.seuil) then
                lum(i,j,1)=0.
              else
                lum(i,j,1)=citypercent/100.*lup(i,j)/(1./pi*(1.-
     +          citypup)*refl(i,j)+citylop)
              endif
         endif
             enddo
            enddo
         open(unit=2,file=coun_lop,status='old')
c reading and normalizing light output pattern - default value countryside
           read(2,*) counlop
           rewind 2
           countlop=0.
           counpup=0.
           do i=1,181
              read(2,*) f,theta
              countlop=countlop+f*2.*pi*sin(pi*theta/180.)*pi/180.
              if (i.le.90) then
                counpup=counpup+f*2.*pi*sin(pi*theta/180.)*pi/180.
              endif
           enddo
           counlop=counlop/countlop
           counpup=counpup/countlop
      print*,'country up:',counpup
         close(unit=2)
            do i=1,nbx
             do j=1,nby
              lum(i,j,2)=0.
c filter out water surface where there is not suppose to be any lights
c this feature have been added to correct the overspill problem of ols data
         if (lup(i,j).lt.cityt) then                                    ! countryside case
              if (water(i,j).gt.seuil) then
                lum(i,j,2)=0.
              else
                lum(i,j,2)=counpercent/100.*lup(i,j)/(1./pi*(1.-
     +          counpup)*refl(i,j)+counlop)
              endif
         endif
             enddo
            enddo


         do n=1,nzon
            print*,'Zone #',n,'/',nzon

            read(1,*) x,y,r,percent,lop_file
         open(unit=2,file=lop_file,status='old')
c reading and normalizing light output pattern - other zones
           read(2,*) lop
           rewind 2
           tlop=0.
           pup=0.
           do i=1,181
              read(2,*) f,theta
              tlop=tlop+f*2.*pi*sin(pi*theta/180.)*pi/180.
              if (i.le.90) then
                pup=pup+f*2.*pi*sin(pi*theta/180.)*pi/180.
              endif
           enddo
           lop=lop/tlop
           pup=pup/tlop
      print*,'Up zone #',n,': ',pup
         close(unit=2)
            do i=1,nbx
             do j=1,nby
              dist=sqrt((real(i)-x)**2.+(real(j)-y)**2.)
              if (dist.le.r) then
c filter out ocean surface where there is not suppose to be any lights
c this feature have been added to correct the overspill problem of ols data
              if (water(i,j).gt.seuil) then
                lum(i,j,n+2)=0.
              else
               lum(i,j,n+2)=percent/100.*lup(i,j)/(1./pi*(1.-pup)*
     +         refl(i,j)+lop)
              endif
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
