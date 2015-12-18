c programme pour convertir un fichier pgm provenant directement 
c de viirs en valeur de 
c
c To compile:
c gfortran viirs2lum-multispectral.f extrants2d.f intrants2d.f interp_modis.f -o viirs2lum
c
c possible option if the array sizes are too big: -mcmodel=large 
c   
c    Copyright (C) 2015  Martin Aube
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
c    along with this program. If not, see <http://www.gnu.org/licensces/> 
c    Contact: martin.aube@cegepsherbrooke.qc.ca
c
      program viirs2lum
      integer wid,nwa,nag,nzo,nbd
      parameter (wid=1024,nwa=1254,nag=181,nzo=99,nbd=150)
c wid is the max size of the domain
c nwa is the number of wavelength in lamp spectrum
c nang is the number of zenith angles in the lop
c nzo is the number of zones
c ndb is the number of spectral bands
      integer na,i,j,n
      integer zone(wid,wid)
      integer valmax,nbx,nby,nw,nb,n_bands,lenbase
      real dist,spct(nwa,nag)
      real dnb(wid,wid),sumwav,Phi_c(wid,wid),wavel(nwa),viirs_sens(nwa)
      real thetas(nag),obsth(wid,wid),obstd(wid,wid),lamph(wid,wid)
c     obsth=obstacle height, obstd=mean free path to the ground, 
c     lamph=lamp height above the ground
      character*72 Gn,zonfile,viirs_resp,bands_file,fctfile
      character*72 lumfile,outfile,viirs_file,rfile
c     Gn=barG file name
      character*12 nom,basename,junk
      real x,y,r,hobst,dobst,hlamp,dat(wid,wid),datf(wid,wid)
      real pi,maxim,gain,offset,xcell0,ycell0,pixsiz,rho(wid,wid),dwav
      real sp(nbd,nag),lum(wid,wid,nbd),rad,bands(nbd,2)
      real G_moy(nwa,nzo),Fdown(nwa,nzo),fctem(nzo,nbd),avgwav(nbd)
      real tlamb,nwi,dbnunits
c nzo nombre max de bandes
      character*3 zonenu,waven,lambda
c default values
c converting to nanoW/cm2/sr en W/m2/sr DNB=DNB*dnbunits
      dnbunits=1.-5
      pi=3.14159
      rad=pi/180.
      valmax=65535
      offset=0.
      n_angles=nag
      pixsiz=1000.
      ycell0=0.
      xcell0=0.
      print*,'Be sure to have in the execution folder:'
      print*,'1- viirs-dnb pgm file'
      print*,'2- zonal file [e.g. Hawaii.zon]'
      print*,'3- modis pgm files for 3 modis bands'
      print*,'4- modis wavelength and name file [modis.dat]'
      print*,'5- spectral bands file [integration_limits.dat]'
      print*,'6- viirs-dnb spectral response [viirs.dat]'
      print*,'7- spectral/angular bar_G functions for each zone'
      print*,'   files 7- are created with make_inputs.py'
      print*,' '
      print*,'================'
      print*,'Output root name of the experiment ?'
c      read*,basename
 
      basename='Hawaii'
      lenbase=index(basename,' ')-1 
      print*,'viirs-dnb file name? [e.g. stable-light.pgm]'
c      read*,viirs_file
      viirs_file='stable_lights_lowcut.pgm'
      print*,'zonal file name? [e.g. Hawaii.zon]'
c      read*,zonfile
      zonfile='Hawaii.zon'
      bands_file='integration_limits.dat'
      viirs_resp='viirs.dat'
c Reading the spectral bands file
      print*,'Reading spectral bands...'
      OPEN(UNIT=42,FILE=bands_file,STATUS='OLD')
        READ(42,*) n_bands
        DO nb=1,n_bands
          READ(42,*) bands(nb,1),bands(nb,2)
          avgwav(nb)=(bands(nb,2)+bands(nb,1))/2.
        ENDDO
      close(unit=42)
c Reading VIIRS-dnb spectral response
      print*,'Reading VIIRS-dnb spectral response...'
      OPEN(UNIT=1,FILE=viirs_resp,STATUS='OLD')
        READ(1,*) junk
        viirs_sum = 0.
        DO nw=1,nwa
          READ(1,*) wavel(nw),viirs_sens(nw)
          viirs_sum = viirs_sum + viirs_sens(nw)
        ENDDO
        dwav=(wavel(2)-wavel(1))
        viirs_sum=viirs_sum*dwav
c Normalizing VIIRS-dnb spectral response
        print*,'Normalizing VIIRS-dnb spectral response...'
        do nw=1,nwa
          viirs_sens(nw)=viirs_sens(nw)/viirs_sum
        enddo
      CLOSE(unit=1)
c lecture de l'image viirs dnb
        print*,'Reading VIIRS-dnb image...'
        nom='VIIRS-DNB '
        call intrants2d(viirs_file,dnb,nom,xcell0,ycell0,pixsiz,nbx,nby)
c interpolate modis images to modelling bands wavelength
c defined in file integration_limits.dat
        print*,'Interpolation MODIS images...'
        call interp_modis()
c reading modis reflectance at ~700nm to filter out the water pixels 
c rho<0.05 . Water show a very low reflectance in that part of the 
c spectrum compared to soil and vegetation
c the problem with low reflectance of water is that in many cases we 
c enconter large cities next to water and then the scattered light 
c over the water surface behave like an important direct source over a
c very dark surface. To overcome cleanly this problem we should make a
c correction for the estimated scattered light. It is not done yet.
        tlamb=700.
        do nb=1,n_bands
          if (avgwav(nb).ge.tlamb) then
c reading the right modis file
             write(lambda, '(I3.3)' ) int(avgwav(nb))
             rfile='modis_'//lambda//'.pgm'
             nom='MODIS ref.'
             call intrants2d(rfile,rho,nom,xcell0,ycell0,pixsiz,nbx,nby)
             cycle
          endif
        enddo
        do i=1,nbx
          do j=1,nbx
              dnb(i,j)*dnb(i,j)*dnbunits
              if (rho(i,j).lt.0.05) then
                 dnb(i,j)=0.
              endif
          enddo
        enddo



c initialisation
        do i=1,nbx
          do j=1,nby
            datf(i,j)=0.
            obsth(i,j)=0.
            obstd(i,j)=0.
            lamph(i,j)=0.
          enddo
        enddo
c reading zones properties
      print*,'Reading zones properties...'
      open(unit=11,file=zonfile,status='old')
        read(11,*) nzon
c
c ==============
c debut boucle sur les zones
c
        do n=1,nzon
c
c
c initialisation
          print*,'Initializing arrays...'
          do i=1,nbx
            do j=1,nby
              Phi_c(i,j)=0.
              obsth(i,j)=0.
              obstd(i,j)=0.
              lamph(i,j)=0.
              do nb=1,n_bands
                lum(i,j,nb)=0.
              enddo
            enddo
          enddo
          write(zonenu, '(I3.3)' ) n
          read(11,*) x,y,r,Gn,hobst,dobst,hlamp
c   reading the barG
          print*,'Reading the barG...',Gn
          open(unit=2,file=Gn,status='old')
            do na=1,nag
              read(2,*) thetas(na), spct(:,na)
            enddo
          close(unit=2)
c     loop over the modelling domain
            do i=1,nbx
              do j=1,nby
                dist=sqrt((real(i)-x)**2.+(real(j)-y)**2.)
                if (dist.le.r) then
c       determining the zone for each pixel
                   zone(i,j)=n
c       defining the mean obstacle height for each pixel
                   obsth(i,j)=hobst
c       defining the mean free path to the ground for each pixel
                   obstd(i,j)=dobst
c       defining the lamp height for each pixel
                   lamph(i,j)=hlamp
                endif
              enddo
            enddo
c     writing obtacle height pgm files
            print*,'Writing obtacle height pgm files...'
            outfile=basename(1:lenbase)//'_obsth_'//zonenu//'.pgm'
            nom='obstacle h'
            maxim=0.
            do i=1,nbx
              do j=1,nby
                if (maxim.lt.obsth(i,j)) then
                  maxim = obsth(i,j)
                endif
              enddo
            enddo
            nom='ObstacleH'
            gain = maxim/real(valmax)
            call extrants2d (outfile,obsth,nom,xcell0,ycell0,pixsiz,
     +      gain,offset,nbx,nby,valmax)
c     writing obtacle mean free path pgm files
            print*,'Writing obtacle mean free path pgm files...'
            outfile=basename(1:lenbase)//'_obstd_'//
     +                zonenu//'.pgm'
            nom='obstacle d'
            maxim=0.
            do i=1,nbx
              do j=1,nby
                if (maxim.lt.obstd(i,j)) then
                  maxim = obstd(i,j)
                endif
              enddo
            enddo
            nom='MeanFreePath'
            gain = maxim/real(valmax)
            call extrants2d (outfile,obstd,nom,xcell0,ycell0,pixsiz,
     +      gain,offset,nbx,nby,valmax)
c     writing light fixture height relative to the ground pgm files
            print*,'Writing light fixture height file...'
            outfile=basename(1:lenbase)//'_altlp_'//
     +                zonenu//'.pgm'
            nom='LampHeight'
            maxim=0.
            do i=1,nbx
              do j=1,nby
                if (maxim.lt.lamph(i,j)) then
                  maxim = lamph(i,j)
                endif
              enddo
            enddo
            gain = maxim/real(valmax)
            call extrants2d (outfile,lamph,nom,xcell0,ycell0,pixsiz,
     +      gain,offset,nbx,nby,valmax)
c calcul du G_moy pour chaque zone et chaque lambda
            print*,'Computing G_moy...'
            do nw=1,nwa 
              G_moy(nw,n)=0.  
              do na=1,57
                G_moy(nw,n)=G_moy(nw,n)+spct(nw,na)*
     +          sin(thetas(na)*rad)*1.*rad/
     +          (1.-cos(56.*rad))
              enddo
            enddo
c calcul du Fdown pour chaque zone et chaque lambda
            print*,'Computing Fdown...'
            do nw=1,nwa 
              Fdown(nw,n)=0.  
              do na=92,181
                Fdown(nw,n)=Fdown(nw,n)+2.*pi*
     +          spct(nw,na)*sin(thetas(na)*rad)*1.*rad
              enddo
            enddo
c determination de la reflectance la plus proche en lambda
        do nb=1,n_bands
          if (avgwav(nb).ge.wavel(nw)) then
c reading the right modis file
             write(lambda, '(I3.3)' ) int(avgwav(nb))
             rfile='modis_'//lambda//'.pgm'
             nom='MODIS ref.'
             call intrants2d(rfile,rho,nom,xcell0,ycell0,pixsiz,nbx,nby)
            cycle
          endif
        enddo
c calcul Phi_c pour chaque pixel
        do i=1,nbx
          do j=1,nby
            sumwav=0.
            do nw=1,nwa
c integrale au denominateur
c sumwav=R(lambda)*(1/pi*rho(lambda)*Fdown(lambda)+G_moy(lambda))
c *dlambda
                sumwav=sumwav+viirs_sens(nw)*(1./pi*rho(i,j)*Fdown(nw,
     +          n)+G_moy(nw,n))*dwav
            enddo
c Phi_c=dnb/sum_lambda
            if (zone(i,j).eq.n) then
              if (sumwav.ne.0.) then
                Phi_c(i,j)=dnb(i,j)/sumwav
              else
                Phi_c(i,j)=0.
              endif
            else
              Phi_c(i,j)=0.
            endif
          enddo
        enddo
c calcul du lumlp pour chaque zone
c integrer le spectre sur chaque bande pour chaque pixel
c et faire le ratio de chaque bande sur le total
          do nb=1,n_bands
            do na=1,n_angles
              sp(nb,na)=0.
              nwi=0.
              do nw=1,nwa
                if ((wavel(nw).ge.bands(nb,1)).and.(wavel(nw).lt.
     +          bands(nb,2))) then
                  sp(nb,na)=sp(nb,na)+spct(nw,na)
                  nwi=nwi+1.
                endif
              enddo
              if (nwi.ne.0.) then
                sp(nb,na)=sp(nb,na)/nwi
              else
                sp(nb,na)=0.
              endif
            enddo
          enddo

c debug pour valider le spectre
c          if (n.eq.11) then
c             na=181
c             do i=1,nwa
c               print*,wavel(i),spct(i,na)
c             enddo
c             print*,'==========='
c             do i=1,n_bands
c               print*,avgwav(i),sp(i,na)
c             enddo
c           stop
c           endif




c extraire le fctem pour chaque zone et chaque bande
        print*,'Extracting LOP...'
        do nb=1,n_bands
          write(waven, '(I3.3)' )  int(avgwav(nb))
            write(zonenu, '(I3.3)' ) n
            fctfile='fctem_wl_'//waven//'_zon_'//zonenu//'.dat'
            open (unit=25,file=fctfile,status='unknown')
              do na=1,n_angles
                write(25,*) sp(nb,na),thetas(na)
              enddo
            close(unit=25)
        enddo
c et multiplier par Phi_c avec la fraction de bande
c pour obtenir le lumlp de chaque zone et chaque bande
        do nb=1,n_bands 
          do i=1,nbx
            do j=1,nby 
              do na=1,nag  
                lum(i,j,nb)=lum(i,j,nb)+Phi_c(i,j)
     +          *sp(nb,na)*2.*pi*sin(thetas(na)*rad)*1.*rad
              enddo
            enddo
          enddo
        enddo
        do nb=1,n_bands
          write(waven, '(I3.3)' ) int(avgwav(nb))
            write(zonenu, '(I3.3)' ) n
            lumfile=basename(1:lenbase)//'_'//waven//'_lumlp_'//
     +                zonenu//'.pgm'
            nom='Luminosity'
            maxim=0.
            do i=1,nbx
              do j=1,nby
                dat(i,j)=lum(i,j,nb)
                  if (maxim.lt.dat(i,j)) then
                    maxim=dat(i,j)
                  endif
                enddo
              enddo
              gain=maxim/real(valmax)
              if (gain.eq.0.) gain=1.
              call extrants2d(lumfile,dat,nom,xcell0,ycell0,pixsiz,
     +        gain,offset,nbx,nby,valmax)
        enddo
c
c
c fin bouche zones
          enddo
c
c ===================
c 
c creating combined lumlp for each band

        do nb=1,n_bands
          do n=1,nzon
          write(waven, '(I3.3)' ) int(avgwav(nb))
            write(zonenu, '(I3.3)' ) n
            lumfile=basename(1:lenbase)//'_'//waven//'_lumlp_'//
     +                zonenu//'.pgm'
            nom='Luminosity'
       call intrants2d(lumfile,dat,nom,xcell0,ycell0,pixsiz,nbx,nby)
            
            do i=1,nbx
              do j=1,nby
                datf(i,j)=datf(i,j)+dat(i,j)
                  if (maxim.lt.datf(i,j)) then
                    maxim=datf(i,j)
                  endif
                enddo
              enddo


              enddo
            lumfile=basename(1:lenbase)//'_'//waven//'_lumlp.pgm'
            nom='Luminosity'
              gain=maxim/real(valmax)
              if (gain.eq.0.) gain=1.
              call extrants2d(lumfile,datf,nom,xcell0,ycell0,pixsiz,
     +        gain,offset,nbx,nby,valmax)
        enddo
        close(unit=11)
      stop
      end
