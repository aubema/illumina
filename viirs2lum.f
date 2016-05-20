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
      real viirs(width,width),lum(width,width,99)
      real pixsiz,gain,offset,dat(width,width)
      real lop,pi,f,theta,seuil,xcell0,ycell0,cityt
      real citysat,citytot,citydown,citycone,counsat
      real councone,cone,sat,tot,down,countot,coundown
      real x(99),y(99),r(99)
      integer nbx,nby,valmax,i,j,nzon,n,ii,jj
      character*72 reflex_file,viirs_file,basename,zonfile
      character*72 outfile,city_lop,coun_lop
      character*72 lop_file,viirs_resp,lambdas_file
      character*12 nom
      character*3 lambno,lambda
c       
      CHARACTER(LEN=1) :: junk
      INTEGER :: n_lambda,ios,n_int,n_modis,wl,int_ind,refl_ind
      REAL :: intp_factor,reflex,maxim,viirs_sum,mid_wl
      REAL :: modis_out(width,width)
      REAL, DIMENSION(:), ALLOCATABLE :: viirs_data,lambdas,
     +modis_lambdas,thetas
      REAL, DIMENSION(:,:), ALLOCATABLE :: int_limits
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: refl,lamp_spct
c
      pi=3.14159
      valmax=65535
      gain=0.7
      offset=0.
      n_angles=181
c le seuil pour le masque water (0=land,1=water). En dessous de seuil 
c on pose la luminosite a 0
      seuil=0.5
      do i=1,width
        do j=1,width
          do n=1,99
            lum(i,j,n)=0.
          enddo
        enddo
      enddo
      cityt=0.
      print*,'Be sure to have an viirs illumina pgm file, a zonal '
      print*,'file, and a illumina pgm reflectance file.'
      print*,'Output root name of the experiment [this name will be use
     +d for all the subsequent files]?'
      read "(a)",basename
c      basename='Hawaii'
      print*,'viirs-dnb file name? [e.g. stable-light.pgm]'
      read "(a)",viirs_file
c      viirs_file='stable_lights_lowcut.pgm'
      print*,'modis reflectance file list file name? [e.g. modis.dat]'
      read "(a)",reflex_file
c      reflex_file='modis.dat'
      print*,'viirs response file name? [e.g. viirs.dat]'
      read "(a)",viirs_resp
c      viirs_resp='viirs.dat'
      print*,'zonal file name? [e.g. Hawaii.zon]'
      read "(a)",zonfile
c      zonfile='Hawaii.zon'
      print*,'spectral integration limits file name ? [e.g. integration
     +_limits.dat]'
      read "(a)",lambdas_file
c      lambdas_file='integration_limits.dat'
c
c Ajouts Alexandre Simoneau
c
c Lecture de la sensibilité du satellite VIIRS
c
c     compter le nombre de ligne du fichier de sensibilite viirs
      n_lambda = 0
      OPEN(UNIT=1,FILE=viirs_resp,STATUS='OLD')
        DO
          READ(1,*,IOSTAT=ios) junk
          IF (ios /= 0) EXIT
          n_lambda = n_lambda + 1
        ENDDO
        REWIND(1)
c       lecture des donnees de sensibilite du viirs
        n_lambda = n_lambda - 1
        ALLOCATE(viirs_data(n_lambda))
        ALLOCATE(lambdas(n_lambda))
        READ(1,*) junk
        viirs_sum = 0.
        DO i=1,n_lambda
          READ(1,*) lambdas(i),viirs_data(i)
          viirs_sum = viirs_sum + viirs_data(i)
        ENDDO
        viirs_sum = viirs_sum * (lambdas(2)-lambdas(1))
        viirs_data = viirs_data / viirs_sum
      CLOSE(1)
c
c Lecture des bornes d'integrations
c        
      OPEN(UNIT=1,FILE=lambdas_file,STATUS='OLD')
        READ(1,*) n_int
        ALLOCATE(int_limits(2,n_int))
        DO i=1,n_int
          READ(1,*) int_limits(:,i)
        ENDDO
      CLOSE(1)
c
c Lecture des fichiers de reflectance modis
c
      OPEN(UNIT=42,FILE=reflex_file,STATUS='OLD')
        READ(42,*) n_modis
        ALLOCATE(modis_lambdas(n_modis))
        ALLOCATE(refl(width,width,n_modis))
        DO i=1,n_modis
          READ(42,*) modis_lambdas(i),lop_file
          CALL intrants2d(lop_file,refl(:,:,i),xcell0,ycell0,
     +    pixsiz,nbx,nby)
        ENDDO
      CLOSE(42)
c be sure to replace negative reflectances by 0
      do i=1,n_modis
        do ii=1,nbx
          do jj=1,nby
            if (refl(ii,jj,i).lt.0.) then
              print*,'WARNIG - negative reflectance replaced by 0'
              refl(ii,jj,i)=0.
            endif
          enddo
        enddo
      enddo

c EST-CE BIEN N_ANGLES QUI VA LA? EST-CE QUE CETTE MATRICE NE DEVRAIT PAS ETRE EN 2D? SOIT WL ET N_ANGLES?
c Seul le spectre pour l'angle actuel est conservé en mémoire, donc 1D suffit.
      ALLOCATE(lamp_spct(n_lambda,n_angles,99))
      ALLOCATE(thetas(n_angles))

c Retour code Martin
      lum = 0.
      lenbase=index(basename,' ')-1  
      call intrants2d(viirs_file,viirs,xcell0,ycell0,pixsiz,nbx,nby)
        open(unit=1,file=zonfile,status='old')
          read(1,*) nzon,city_lop
          read(1,*) cityt,coun_lop              
c cityt est le seuil du fichier viirs pour choisir entre ville et campagne
          open(unit=2,file=city_lop,status='old')
            do i=1,n_angles
              read(2,*) thetas(i), lamp_spct(:,i,1)
            enddo
          close(unit=2)
          open(unit=2,file=coun_lop,status='old')
            do i=1,n_angles
              read(2,*) thetas(i), lamp_spct(:,i,2)
            enddo
          close(unit=2)
          do n=1,nzon
            read(1,*) x(n),y(n),r(n),lop_file
            open(unit=2,file=lop_file,status='old')
              do i=1,n_angles
                read(2,*) thetas(i), lamp_spct(:,i,n+2)
              enddo
            close(unit=2)
          enddo
        close(1)
c CITY
          do int_ind=1,n_int
            do wl=1,n_lambda
              if ((lambdas(wl).lt.int_limits(1,int_ind)) .or. 
     +        (lambdas(wl).gt.int_limits(2,int_ind))) then
                cycle
              endif
              if (lambdas(wl).lt.modis_lambdas(1)) then
                intp_factor = 1
                refl_ind = 0
              else if (lambdas(wl).gt.modis_lambdas(n_modis)) then
                intp_factor = 0
                refl_ind = 0
              else
                do refl_ind=1,n_modis
                  if (lambdas(wl).ge.modis_lambdas(refl_ind)) exit
                enddo
                intp_factor = (lambdas(wl)-modis_lambdas(refl_ind))/
     +          (modis_lambdas(refl_ind+1)-modis_lambdas(refl_ind))
              endif
c reading city light output pattern - default value cities
              citytot=0.
              citydown=0.
              citycone=0.
              citysat=0.
              do i=1,n_angles
                theta=thetas(i)
C REVOIR LA LECTURE DE LAMP_SPCT (IDEM PARTOUT AILLEURS DANS LE CODE)
C IL ME SEMBLE AUSSI QUE CES VARIABLES INTEGREES DEVRAIT INCLURE LA REPONSE DU VIIRS MAIS PEUT-ETRE QUE CET ELEMENT ES FAIT AILLEURS?
                citytot=citytot+lamp_spct(wl,i,1)*2.*pi*
     +                  sin(pi*theta/180.)*pi/180.
                if (i.gt.90) then                       
                  citydown=citydown+lamp_spct(wl,i,1)*2.*pi*sin(pi*
     +            theta/180.)*pi/180.
                endif
                if (i.le.56) then
                  citysat=citysat+lamp_spct(wl,i,1)*2.*pi*sin(pi*theta
     +            /180.)*pi/180.          ! acceptance angle of VIIRS
                  citycone=citycone+2.*pi*sin(pi*theta/180.)*pi/180.
                endif
              enddo
              print*,'tot=',citytot,'down=',citydown,'cone=',citycone,
     +        'sat=',citysat,lamp_spct(wl,i-1,1)
              if (citytot.eq.0.) citytot=1.
              do i=1,nbx
                do j=1,nby
                  if (viirs(i,j).ge.cityt) then
c F = L_viirs / (cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))
                    if (refl_ind.eq.0) then
                      if (intp_factor.eq.0) then
                        reflex = refl(i,j,n_modis)
                      else
                        reflex = refl(i,j,1)
                      endif
                    else
                      reflex = refl(i,j,refl_ind) + intp_factor*
     +                (refl(i,j,refl_ind+1)-refl(i,j,refl_ind))
                    endif
                      lum(i,j,1)=lum(i,j,1)+(viirs_data(wl)*(reflex/
     +                pi*citydown/citytot+citysat/(citycone*citytot)))
                  endif
                enddo
              enddo
            enddo
            do i=1,nbx
              do j=1,nby
                if (lum(i,j,1).NE.0.) then
                  lum(i,j,1) = viirs(i,j) / ( lum(i,j,1)
     +                                      * (lambdas(2)-lambdas(1)) )
                endif
              enddo
            enddo
c COUNTRY
            do wl=1,n_lambda
              if ((lambdas(wl).lt.int_limits(1,int_ind)) .or. 
     +        (lambdas(wl).gt.int_limits(2,int_ind))) then
                cycle
              endif
              if (lambdas(wl).lt.modis_lambdas(1)) then
                intp_factor = 1
                refl_ind = 0
              else if (lambdas(wl).gt.modis_lambdas(n_modis)) then
                intp_factor = 0
                refl_ind = 0
              else
                do refl_ind=1,n_modis
                  if (lambdas(wl).ge.modis_lambdas(refl_ind)) exit
                enddo
                intp_factor = (lambdas(wl)-modis_lambdas(refl_ind))/
     +          (modis_lambdas(refl_ind+1)-modis_lambdas(refl_ind))
              endif
c reading light output pattern - default value countryside
              countot=0.
              coundown=0.
              councone=0.
              counsat=0.
              do i=1,n_angles
                theta=thetas(i)
                countot=countot+lamp_spct(wl,i,2)*2.*pi*
     +                  sin(pi*theta/180.)*pi/180.
                if (i.gt.90) then                       
                  coundown=coundown+lamp_spct(wl,i,2)*2.*pi*sin(pi*
     +            theta/180.)*pi/180.
                endif
                if (i.le.56) then
                  counsat=counsat+lamp_spct(wl,i,2)*2.*pi*sin(pi*theta
     +            /180.)*pi/180.          ! acceptance angle of VIIRS
                  councone=councone+2.*pi*sin(pi*theta/180.)*pi/180.
                endif
              enddo
              if (countot.eq.0.) countot=1.
              do i=1,nbx
                do j=1,nby
                  if (viirs(i,j).lt.cityt) then
c F = L_viirs / (cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))
                    if (refl_ind.eq.0) then
                      if (intp_factor.eq.0) then
                        reflex = refl(i,j,n_modis)
                      else
                        reflex = refl(i,j,1)
                      endif
                    else
                      reflex = refl(i,j,refl_ind) + intp_factor*
     +                (refl(i,j,refl_ind+1)-refl(i,j,refl_ind))
                    endif
                    lum(i,j,2)=lum(i,j,2)+(viirs_data(wl)*(reflex/pi*
     +              coundown/countot+counsat/(councone*countot)))
                  endif
                enddo
              enddo
            enddo
            do i=1,nbx
              do j=1,nby
                if (lum(i,j,2).NE.0.) then
                  lum(i,j,2) = viirs(i,j) / ( lum(i,j,2)
     +                                      * (lambdas(2)-lambdas(1)) )
                endif
              enddo
            enddo
c other zones (circular zones)
            do n=1,nzon
              print*,'Zone #',n,'/',nzon
              do wl=1,n_lambda
                if ((lambdas(wl).lt.int_limits(1,int_ind)) .or.
     +            (lambdas(wl).gt.int_limits(2,int_ind))) then
                  cycle
                endif
                if (lambdas(wl).lt.modis_lambdas(1)) then
                  intp_factor = 1
                  refl_ind = 0
                else if (lambdas(wl).gt.modis_lambdas(n_modis)) then
                  intp_factor = 0
                  refl_ind = 0
                else
                  do refl_ind=1,n_modis
                    if (lambdas(wl).ge.modis_lambdas(refl_ind)) exit
                  enddo
                  intp_factor = (lambdas(wl)-modis_lambdas(refl_ind))/
     +            (modis_lambdas(refl_ind+1)-modis_lambdas(refl_ind))
                endif
c reading light output pattern - other zones
                tot=0.
                down=0.
                cone=0.
                sat=0.
                do i=1,n_angles
                  theta=thetas(i)
                  tot=tot+lamp_spct(wl,i,n+2)*2.*pi*sin(pi*theta/180.)
     +                *pi/180.
                  if (i.gt.90) then                       
                    down=down+lamp_spct(wl,i,n+2)*2.*pi*
     +                   sin(pi*theta/180.)*pi/180.
                  endif
                  if (i.le.56) then
                    sat=sat+lamp_spct(wl,i,n+2)*2.*pi*sin(pi*theta                ! acceptance angle of VIIRS
     +              /180.)*pi/180.
                    cone=cone+2.*pi*sin(pi*theta/180.)*pi/180.
                  endif
                enddo
                if (tot.eq.0.) tot=1.
c               print*,'Down zone #',n,': ',down
                do i=1,nbx
                  do j=1,nby
                    dist=sqrt((real(i)-x(n))**2.+(real(j)-y(n))**2.)
                    if (dist.le.r(n)) then
c F = L_viirs / (cos(0)*refl/pi*int_90-180[lop*d_omeg]/int_0-180[lop*d_omeg]+ int_0-56[lop*d_omeg]/(int_0-56[d_omeg]*int_0-180[lop*d_omeg]))

c Martin
c F = L_viirs / int_273-899.5[(RVIIRS_lamb * SP_lamb * (refl/pi*int_90-180[lop*d_omeg]+ int_0-56[sin(z)*lop*d_omeg/(cos(56)-1)])))]

                      if (refl_ind.eq.0) then
                        if (intp_factor.eq.0) then
                          reflex = refl(i,j,n_modis)
                        else
                          reflex = refl(i,j,1)
                        endif
                      else
                        reflex = refl(i,j,refl_ind) + intp_factor*
     +                  (refl(i,j,refl_ind+1)-refl(i,j,refl_ind))
                      endif
                      lum(i,j,n+2)=lum(i,j,n+2)+(viirs_data(wl)*
     +                (reflex/pi*down/tot+sat/(cone*tot)))
                      do k=1,n+1
                        lum(i,j,k)=0.                                           ! Reset previous values inside the new zone
                      enddo                                                      ! This allow the new zone to get full priority
                    endif                                                       ! and overwrite previous zones 
                  enddo
                enddo
              enddo
              do i=1,nbx
                do j=1,nby
                  if (lum(i,j,n+2).NE.0.) then
                    lum(i,j,n+2) = viirs(i,j)/( lum(i,j,n+2)
     +                                        *(lambdas(2)-lambdas(1)) )
                  endif
                enddo
              enddo
            enddo
           
c write output  
            mid_wl = ( int_limits(1,int_ind)+int_limits(2,int_ind) ) /2.
            do n=1,nzon+2
              write(lambno, '(I3.3)' ) n
              write(lambda, '(I3.3)' ) int(mid_wl)
              outfile=basename(1:lenbase)//'_'//lambda//'_lumlp_'//
     +                lambno//'.pgm'
              nom='Luminosity'
              maxim=0.
              do i=1,nbx
                do j=1,nby
                  dat(i,j)=lum(i,j,n)
                  if (maxim.lt.lum(i,j,n)) then
                    maxim = lum(i,j,n)
                  endif
c                 maxim = MAX(lum(i,j,n),maxim)
                enddo
              enddo
              gain = maxim/real(valmax)
              call extrants2d (outfile,dat,nom,xcell0,ycell0,pixsiz,
     +        gain,offset,nbx,nby,valmax)
              print*,n,maxim
            enddo
            do wl=1,n_lambda
              if (lambdas(wl).ge.mid_wl) exit
            enddo
            if ((lambdas(wl)-mid_wl).gt.(mid_wl-lambdas(wl-1))) wl=wl-1
            if (lambdas(wl).lt.modis_lambdas(1)) then
              intp_factor = 1
              refl_ind = 0
            else if (lambdas(wl).gt.modis_lambdas(n_modis)) then
              intp_factor = 0
              refl_ind = 0
            else
              do refl_ind=1,n_modis
                if (lambdas(wl).ge.modis_lambdas(refl_ind)) exit
              enddo
              intp_factor = (lambdas(wl)-modis_lambdas(refl_ind))/
     +        (modis_lambdas(refl_ind+1)-modis_lambdas(refl_ind))
            endif
            do i=1,nbx
              do j=1,nby
                if (refl_ind.eq.0) then
                  if (intp_factor.eq.0) then
                    modis_out(i,j) = refl(i,j,n_modis)
                  else
                    modis_out(i,j) = refl(i,j,1)
                  endif
                else
                  modis_out(i,j) = refl(i,j,refl_ind) + intp_factor*
     +                    (refl(i,j,refl_ind+1)-refl(i,j,refl_ind))
                endif
              enddo
            enddo
            nom="Reflectance"
            maxim=0.
            do i=1,nbx
              do j=1,nby
                dat(i,j)=modis_out(i,j)
                if (maxim.lt.modis_out(i,j)) then
                  maxim = modis_out(i,j)
                endif
c               maxim = MAX(modis_out(i,j),maxim)
              enddo
            enddo
            write(lambda, '(I3.3)' ) int(mid_wl)
            outfile='modis_'//lambda//'.pgm'
            gain = maxim/real(valmax)
            call extrants2d (outfile,dat,nom,xcell0,ycell0,pixsiz,
     +      gain,offset,nbx,nby,valmax)
          enddo

      DEALLOCATE(viirs_data)
      DEALLOCATE(lambdas)
      DEALLOCATE(int_limits)
      DEALLOCATE(modis_lambdas)
      DEALLOCATE(refl)
      DEALLOCATE(lamp_spct)
      DEALLOCATE(thetas)
      stop
      end
