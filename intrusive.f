c intrusive light calculation program
c
c gfortran intrusive.f intrants2D.f extrants2d.f -o intrusive
c 
      program intrusive
      integer width                                                       ! Matrix dimension in Length/width and height
      parameter (width=1024)

      real Irad(width,width)                                              ! intrusive radiance
      real h_w                                                            ! height of the center of a window relative to the ground
      real pvalno(181),pval(181)                                          ! Light output pattern
      real srei(width,width)                                              ! Modis reflectance
      real d_o(width,width)                                               ! avg distance between obstacles
      real h_o(width,width)                                               ! avg height of obstacles relative to the ground
      real h_l(width,width)                                               ! avg lamp height relative to the ground
      real intlu(width,width)                                             ! lamp flux
      real filfac                                                         ! filling factor of the facades in front of the window e.g. if the distance between the house is twice the width of the house the filling factor will be 0.33
      real val2d(width,width)
      real z,pi,dz
      real inteo,integ,pvalto
      real Iradmax,gain,offset,xcell0,ycell0,pixsiz
      real ratmoy,rat1
      integer stype,i,j,k,iw,nwav,nw,nzon,nz,lenbase,valmax
      integer nbx,nby
      integer nmoy,n1
      character*3 zon(120),wav(400)
      character*72 basenm,ohfile,odfile,lhfile,rfile,lfile,lopf,intrufi
      character*72 pafile
      character*12 nom
      character*4 suffix
      pi=3.14159
      dz=pi/180.
      ratmoy=0.
      rat1=0.
      nmoy=0
      n1=0
      print*,'Name of the experiment?'
      read*,basenm
      print*,'Height to the center of a window from the ground (m)?'
      read*,h_w
      print*,'Surface filling factor of the facades (0-1)?'
      read*,filfac
      if (h_w.ge.100.) then
        write(suffix, '(F4.0)' ) h_w
      elseif (h_w.ge.10.) then
        write(suffix, '(F4.1)' ) h_w
      else
        write(suffix, '(F4.2)' ) h_w
      endif
c
c  determine the Length of basenm
c 
      lenbase=index(basenm,' ')-1  
c   reading wav.lst
c
      nwav=0
      open(unit=1,file='wav.lst',status='old')
          do k=1,400
             read(1,*,end=10) wav(k)
             nwav=nwav+1
          enddo
 10   close(unit=1)
c   reading zon.lst
c
      nzon=0
      open(unit=1,file='zon.lst',status='old')
          do k=1,120
             read(1,*,end=20) zon(k)
             nzon=nzon+1
          enddo
 20   close(unit=1)
c reading obstacle heights
        ohfile=basenm(1:lenbase)//'_obsth.pgm'
        call intrants2d(ohfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
         do i=1,nbx                                                       ! beginning of the loop over all cells along x.
           do j=1,nby                                                     ! beginning of the loop over all cells along y.
             h_o(i,j)=val2d(i,j)                                          ! Filling of the array
           enddo                                                          ! end of the loop over all cells along y.
         enddo                                                            ! end of the loop over all cells along x.
c reading obstacle distances
        odfile=basenm(1:lenbase)//'_obstd.pgm'
        call intrants2d(odfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
         do i=1,nbx                                                       ! beginning of the loop over all cells along x.
           do j=1,nby                                                     ! beginning of the loop over all cells along y.
             d_o(i,j)=val2d(i,j)                                          ! Filling of the array
           enddo                                                          ! end of the loop over all cells along y.
         enddo                                                            ! end of the loop over all cells along x.
c readind lamp heights
        lhfile=basenm(1:lenbase)//'_altlp.pgm'
        call intrants2d(lhfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
         do i=1,nbx                                                       ! beginning of the loop over all cells along x.
           do j=1,nby                                                     ! beginning of the loop over all cells along y.
             h_l(i,j)=val2d(i,j)                                          ! Filling of the array 
           enddo                                                          ! end of the loop over all cells along y.
         enddo                                                            ! end of the loop over all cells along x.
c
c 
      do nw=1,nwav
        do i=1,nbx
          do j=1,nby
            Irad(i,j)=0.
          enddo
        enddo
        Iradmax=0.
c reading modis reflectances
        rfile='modis_'//wav(nw)//'.pgm'
        call intrants2d(rfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
          do j=1,nby                                                      ! beginning of the loop over all cells along y.
            srei(i,j)=val2d(i,j)                                          ! Filling of the array 
          enddo                                                           ! end of the loop over all cells along y.
        enddo      

        do nz=1,nzon
c reading LOPs
          lopf='fctem_wl_'//wav(nw)//'_zon_'//zon(nz)//'.dat'
          pvalto=0.
          open(UNIT=7, FILE=lopf,status='OLD')                            ! opening file pa#.dat, angular photometry.
            do k=1,181                                                    ! beginning of the loop for the 181 data points
              read(7,*) pval(k)                                           ! reading of the data in the array pval.
              pvalto=pvalto+pval(k)*2.*pi*                                ! Sum of the values of the  photometric function 
     a        sin(real(k-1)*dz)*dz                                        ! (pvaleur x 2pi x sin z x z) (ou z egale 
c                                                                         ! (i-1) x 1 degrees).
            enddo                                                         ! end of the loop over the 181 donnees du fichier pa#.dat.
          close(7)                                                        ! closing file pa#.dat, angular photometry.
          do k=1,181
            if (pvalto.ne.0.) pvalno(k)=pval(k)/pvalto                    ! Normalisation of the photometric function.
          enddo   

c reading du fluxes
          lfile=basenm(1:lenbase)//'_'//wav(nw)//'_lumlp_'//zon(nz)//
     +'.pgm'
          call intrants2d(lfile,val2d,xcell0,ycell0,pixsiz,nbx,nby)
          do i=1,nbx                                                      ! beginning of the loop over all cells along x.
            do j=1,nby                                                    ! beginning of the loop over all cells along y.
              intlu(i,j)=val2d(i,j)                                       ! Filling of the array 
            enddo                                                         ! end of the loop over all cells along y.
          enddo      

          do i=1,nbx
            do j=1,nby
              if ((d_o(i,j).gt.0.).and.(intlu(i,j).gt.0.)) then
c calculate the basic angles of the geometry
                z_o=pi/2.-atan((h_o(i,j)-h_l(i,j))/(d_o(i,j)/2.))
                z_g=pi-atan((d_o(i,j)/2.)/h_l(i,j))
                z_w=pi/2.-atan((h_w-h_l(i,j))/(d_o(i,j)/2.))
c           print*,z_g,z_o,z_w

c integrate the LOP from obstacle base to obstacle top (inteo)
c and LOP from nadir to obstacle base (integ)
                inteo=0.
                integ=0.
                do k=1,181
                  z=real(k-1)*pi/180.
                  if ((z.ge.z_o).and.(z.lt.z_g)) then
                   inteo=inteo+pvalno(k)*sin(z)*dz
                  endif
                  if ((z.ge.z_g).and.(z.lt.pi)) then
                   integ=integ+pvalno(k)*sin(z)*dz
                  endif
                  if (abs(z-z_w).lt.dz/2.) then 
                     iw=k
                  endif
 
                enddo

                Irad(i,j)=Irad(i,j)+intlu(i,j)*(pvalno(iw)+2.*srei(i,j)*
     +          integ+0.25*srei(i,j)*inteo*filfac)
       if (((0.25*srei(i,j)*inteo*filfac).ne.0.).and.
     + ((2.*srei(i,j)*integ).ne.0.)) then
c       print*,pvalno(iw),2.*srei(i,j)*integ,0.25*srei(i,j)*inteo*filfac,
c     +srei(i,j),Irad(i,j)
        rat1=rat1+pvalno(iw)/(2.*srei(i,j)*integ)
        ratmoy=ratmoy+(2.*srei(i,j)*integ)/(0.25*srei(i,j)*inteo*filfac)
        nmoy=nmoy+1
        n1=n1+1
       endif
                if (Irad(i,j).gt.Iradmax) then
                 Iradmax=Irad(i,j)
                endif
              endif
            enddo
          enddo
        enddo                                                             ! end of loop over zones
c writing the instrusive light map for each wavelength
c
        intrufi=basenm(1:lenbase)//'_'//suffix//'_'//wav(nw)//
     +  '_intrus.pgm'
        nom='Intrusive '
        valmax=65535       
        gain=Iradmax/real(valmax)
        offset=0.
        call extrants2d (intrufi,Irad,nom,xcell0,ycell0,pixsiz,
     +  gain,offset,nbx,nby,valmax)
        print*,rat1/real(n1),ratmoy/real(nmoy)
      enddo                                                               ! end of loop over wavelength
        print*,'final',rat1/real(n1),ratmoy/real(nmoy)
      stop
      end

