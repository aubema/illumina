c programme pour convertir un fichier bin provenant directement 
c de viirs en valeur de 
c
c To compile:
c gfortran viirs2lum-multispectral.f 2dout.f 2din.f interp_modis.f -o viirs2lum
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
      integer nbx,nby,nw,nb,n_bands,lenbase
      real dist,spct(nwa,nag),zone(wid,wid)
      real dnb(wid,wid),sumwav,Phi_c(wid,wid),wavel(nwa),viirs_sens(nwa)
      real thetas(nag),obsth(wid,wid),obstd(wid,wid),lamph(wid,wid)
      real obstf(wid,wid),fobst(nzo)
c     obsth=obstacle height, obstd=mean free path to the ground, 
c     lamph=lamp height above the ground
      character*72 Gn(nzo),zonfile,viirs_resp,bands_file,fctfile
      character*72 lumfile,outfile,viirs_file,rfile
c     Gn=barG file name
      character*12 basename,junk
      real x(nzo),y(nzo),r(nzo),hobst(nzo),dobst(nzo),hlamp(nzo)
      real dat(wid,wid),datf(wid,wid)
      real pi,pixsiz,rho(wid,wid),dwav
      real sp(nbd,nag),lum(wid,wid,nbd),rad,bands(nbd,2)
      real G_moy(nwa,nzo),Fdown(nwa,nzo),fctem(nzo,nbd),avgwav(nbd)
      real tlamb,nwi,dbnunits
c nzo nombre max de bandes
      character*3 zonenu,waven,lambda
c default values
c converting to nanoW/cm2/sr en W/m2/sr DNB=DNB*dnbunits
      dnbunits=1.E-5
      pi=3.14159
      rad=pi/180.
      valmax=65535
      offset=0.
      n_angles=181
      pixsiz=1000.
      ycell0=0.
      xcell0=0.
      print*,'Be sure to have in the execution folder:'
      print*,'1- viirs-dnb bin file'
      print*,'2- zonal file [e.g. Hawaii.zon]'
      print*,'3- modis bin files for 4 modis bands'
      print*,'4- modis wavelength and name file [modis.dat]'
      print*,'5- spectral bands file [integration_limits.dat]'
      print*,'6- viirs-dnb spectral response [viirs.dat]'
      print*,'7- spectral/angular bar_G functions for each zone'
      print*,'   files 7- are created with make_inputs.py'
      print*,' '
      print*,'================'
      print*,'Output root name of the experiment ?'
      read*,basename
      lenbase=index(basename,' ')-1 
      print*,'viirs-dnb file name? [e.g. stable-light.bin]'
      read*,viirs_file
      print*,'zonal file name? [e.g. Hawaii.zon]'
      read*,zonfile
      bands_file='integration_limits.dat'
      viirs_resp='viirs.dat'
c Reading the spectral bands file
      print*,'Reading spectral bands...'
      OPEN(UNIT=42,FILE=bands_file,STATUS='OLD')
        READ(42,*) n_bands
        DO nb=1,n_bands+1
          READ(42,*) bands(nb,1)
        ENDDO
        do nb=1,n_bands
          bands(nb,2)=bands(nb+1,1)
          avgwav(nb)=(bands(nb,2)+bands(nb,1))/2.
        enddo
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
        call 2din(nbx,nby,viirs_file,dnb)
c interpolate modis images to modelling bands wavelength
c defined in file integration_limits.dat
        print*,'Interpolation MODIS images...'
        call interp_modis()
c reading modis reflectance at 700nm to filter out the water pixels 
c rho<0.01 . Water show a very low reflectance in that part of the 
c spectrum compared to soil and vegetation
c the problem with low reflectance of water is that in many cases we 
c enconter large cities next to water and then the scattered light 
c over the water surface behave like an important direct source over a
c very dark surface. To overcome cleanly this problem we should make a
c correction for the estimated scattered light. It is not done yet.
c the 0.01 threshold is not really validated
        rfile='modis_700.bin'
        call 2din(nbx,nby,rfile,rho)
        do i=1,nbx
          do j=1,nbx
              dnb(i,j)=dnb(i,j)*dnbunits*pixsiz*pixsiz
              if (rho(i,j).lt.0.01) then
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
          read(11,*) x(n),y(n),r(n),Gn(n),hobst(n),dobst(n),fobst(n),
     +    hlamp(n)
c       loop over the modelling domain
            do i=1,nbx
              do j=1,nby
                dist=sqrt((real(i)-x(n))**2.+(real(j)-y(n))**2.)
                if (dist.le.r(n)) then
c       determining the zone for each pixel
                   zone(i,j)=real(n)
c       defining the mean obstacle height for each pixel
                   obsth(i,j)=hobst(n)
c       defining the mean free path to the ground for each pixel
                   obstd(i,j)=dobst(n)
c       defining obstacle filling factor for each pixel
                   obstf(i,j)=fobst(n)
c       defining the lamp height for each pixel
                   lamph(i,j)=hlamp(n)
                endif
              enddo
            enddo          
        enddo
        do n=1,nzon
c
c initialisation
          print*,'Initializing arrays...'
          do i=1,nbx
            do j=1,nby
              Phi_c(i,j)=0.
              do nb=1,n_bands
                lum(i,j,nb)=0.
              enddo
            enddo
          enddo
          write(zonenu, '(I3.3)' ) n

c   reading the barG
          print*,'Reading the barG...',Gn(n)
          open(unit=2,file=Gn(n),status='old')
            do na=1,nag
              read(2,*) thetas(na), spct(:,na)
            enddo
          close(unit=2)


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
             rfile='modis_'//lambda//'.bin'
             call 2din(nbx,nby,rfile,rho)
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
            if (zone(i,j).eq.real(n)) then
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
     +    zonenu//'.bin'
          call 2dout(nbx,nby,lumfile,dat)
        enddo
c
c
c fin bouche zones
          enddo
c     writing light fixture height relative to the ground bin files
            print*,'Writing light fixture height file...'
            outfile=basename(1:lenbase)//'_altlp.bin'
            call 2dout(nbx,nby,outfile,lamph)
c
c
c     writing obtacle height bin file
            print*,'Writing obtacle height bin file...'
            outfile=basename(1:lenbase)//'_obsth.bin'
            call 2dout(nbx,nby,outfile,obsth)
c     writing obtacle mean free path bin file
            print*,'Writing obtacle mean free path bin file...'
            outfile=basename(1:lenbase)//'_obstd.bin'
            call 2dout(nbx,nby,outfile,obstd)
c     writing obtacle filling factor bin file
            print*,'Writing obtacle filling factor bin file...'
            outfile=basename(1:lenbase)//'_obstf.bin'
            call 2dout(nbx,nby,outfile,obstf)
c     writing zone definition bin file
            print*,'Writing zone definition bin file...'
            outfile=basename(1:lenbase)//'_zone.bin'
            call 2dout(nbx,nby,outfile,zone)
c
c ===================
c 
c creating combined lumlp for each band

        do nb=1,n_bands
          do n=1,nzon
            write(waven, '(I3.3)' ) int(avgwav(nb))
            write(zonenu, '(I3.3)' ) n
            lumfile=basename(1:lenbase)//'_'//waven//'_lumlp_'//
     +                zonenu//'.bin'
            call 2din(nbx,nby,lumfile,dat)
            do i=1,nbx
              do j=1,nby
                datf(i,j)=datf(i,j)+dat(i,j)
              enddo
            enddo
          enddo
          lumfile=basename(1:lenbase)//'_'//waven//'_lumlp.bin'
          call 2dout(nbx,nby,lumfile,datf)
        enddo
        close(unit=11)
      stop
      end
