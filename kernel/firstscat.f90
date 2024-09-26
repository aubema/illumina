!-----------------------------------------------------------------------
!
!=======================================================================
! Routine firstscat
! Calculation of the first scattering 
! part of Illumina v3
!-----------------------------------------------------------------------
!   
!    Copyright (C) 2024  Martin Aube
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    Contact: aubema@gmail.com
!
!
      subroutine firstscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,iz,lamplu,ofill,srefl,drefle, &
      reflsiz,obsH,altsol,inclix,incliy,pvalno,stype,tranam,tabs,tranaa,tranal,secdif,secdil, &
      fdifan,fdifl,haer,hlay,dx,dy,cloudt,cloudbase,omefov,scal,portio,idif,icloud)
      integer width,nzon
      real*8 pi                                                 
      parameter (width=512,nzon=256)
      parameter (pi=3.141592654)
      integer rho,x_s,y_s,x_c,y_c,x_sr,y_sr,stype,is3a,jda,ida,anglez,na,naz,cloudt
      integer x_obs,y_obs
      real*8 z_s,z_c,z_sr,lamplu(width,width,nzon)
      real*8 haer,hlay,dx,dy,cloudbase,ofill(width,width),drefle(width,width),obsH(width,width)
      real*8 altsol(width,width),fdifan(181),fdifl(181),pvalno(181,nzon)
      real*8 un,zenith,rx_difa,ry_difa,z_difa,dxp,dyp
      real*8 rx_sr,ry_sr,rx_s,ry_s,rx_c,ry_c
      real*8 ds1,ds3
      real*8 angdif,inclix(width,width),incliy(width,width),epsilx,epsily,reflsiz
      real*8 ouvang,iz,rcloud,doc2,dsc2,azcl1,azcl2
      real*8 rx_obs,ry_obs,z_obs,nbang
      real*8 tranam,tranaa,tranal,secdif,secdil
      real*8 flux,icloud,idif,flrefl,irefl,P_indir,omefov,transl,P_dif,srefl,scal,portio
      real*8 omega,omeg2,omemax,hh,ff1,ff2,distd,transm,tabs,transa
      real*8 pdif
      real*8 xc,yc,zc,xn,yn,zn                                            ! Position (meter) of the elements (starting point, final point) for the calculation of the solid angle.
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! Components of the vectors used in the solid angle calculation routine.    
      un=1.
      ! omemax: exclude calculations too close (<10m) this is a sustended angle of 1 deg.
      ! the calculated flux is highly sensitive to that number for a very high
      ! pixel resolution (a few 10th of meters). We assume anyway that somebody
      ! observing the sky will never lies closer than that distance to a
      ! light fixture. This number is however somehow subjective and that means
      ! that the value of sky brightness near sources will be affected by this
      ! choice
      omemax=1./((10.D0)**2.)
      rx_sr=real(x_sr)*dx
      ry_sr=real(y_sr)*dy
      rx_s=real(x_s)*dx
      ry_s=real(y_s)*dy
      rx_c=real(x_c)*dx
      ry_c=real(y_c)*dy
      rx_obs=real(x_obs)*dx
      ry_obs=real(y_obs)*dy 
      if (rho.eq.0) then ! from source
        call anglezenithal(rx_s,ry_s,z_s,rx_c,ry_c,z_c,zenith)
        call blocking(x_s,y_s,z_s,x_c,y_c,z_c,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)
        ! computation of the transmittance between the reflection surface and the 3rd scattering voxel
        distd=sqrt((rx_c-rx_s)**2.+(ry_c-ry_s)**2.+(z_c-z_s)**2.)
        call transmitm(zenith,z_s,z_c,distd,transm,tranam,tabs)
        call transmita(zenith,z_s,z_c,distd,haer,transa,tranaa)
        call transmita(zenith,z_s,z_c,distd,hlay,transl,tranal)
        ! computation of the solid angle of the 3rd scattering voxel seen from the reflecting surface
        omega=1./distd**2.
        if (omega.gt.omemax) omega=0.D0
        anglez=idnint(180.D0*zenith/pi)+1
        P_dif=pvalno(anglez,stype)
        flux=lamplu(x_s,y_s,stype)*P_dif*omega*transm*transa*transl*(1.-ff1)*(1.-ff2)*hh ! flux crossing line of sight voxel
        !flux=lamplu(x_s,y_s,stype)*P_dif*omega*transm*transa*transl*(1.-ff1)*hh
        ! computing the scattering probability toward 2nd scat voxel
        if (omega.ne.0.D0) then
          call angle3points(rx_s,ry_s,z_s,rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,angdif) ! scattering angle.
          call diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil,fdifan,fdifl,haer,hlay,pdif,z_c) ! scattering probability pdif of the direct light.
        else
          pdif=0.D0
        endif
        idif=flux*pdif*scal*portio        
        if (cloudt.ne.0) then ! line of sight voxel = cloud
          if (cloudbase-z_c.le.iz*scal) then ! this is the cloud base interface
            call anglezenithal(rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,azcl1) ! zenith angle from cloud to observer
            call anglezenithal(rx_c,ry_c,z_c,rx_s,ry_s,z_s,azcl2) ! zenith angle from source to cloud
            doc2=(rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.+(z_c-z_obs)**2.
            dsc2=(rx_s-rx_c)**2.+(ry_s-ry_c)**2.+(z_s-z_c)**2.
            call cloudreflectance(zenith,cloudt,rcloud) ! cloud intensity from direct illum
            ! computing the cloud intensity toward the observer
            icloud=flux/omega*rcloud*doc2*omefov*abs(cos(azcl2)/cos(azcl1))/dsc2/pi ! return to main
          endif
        endif
        
        
      else ! from ground
        call anglezenithal(rx_s,ry_s,z_s,rx_sr,ry_sr,z_sr,zenith)
        ! computation of the transmittance between the source and the surface
        distd=sqrt((rx_s-rx_sr)**2.+(ry_s-ry_sr)**2.+(z_s-z_sr)**2.)
        call transmitm(zenith,z_s,z_sr,distd,transm,tranam,tabs)
        call transmita(zenith,z_s,z_sr,distd,haer,transa,tranaa)
        call transmita(zenith,z_s,z_sr,distd,hlay,transl,tranal)
        ! computation of the solid angle of the reflecting cell seen from the source
        xc=dble(x_sr)*dble(dx) ! Position in meters of the observer voxel (longitude).
        yc=dble(y_sr)*dble(dy) ! Position in meters of the observer voxel (latitu).
        zc=dble(z_sr) ! Position in meters of the observer voxel (altitude).
        xn=dble(x_s)*dble(dx) ! Position in meters of the source (longitude).
        yn=dble(y_s)*dble(dy) ! Position in meters of the source (latitu).
        zn=dble(z_s) ! Position in meters of the source (altitude).
        epsilx=inclix(x_sr,y_sr) ! tilt along x of the ground reflectance
        epsily=incliy(x_sr,y_sr) ! tilt along x of the ground reflectance
        if (dx.gt.reflsiz) then ! use a sub-grid surface when the reflectance radius is smaller than the cell size
          if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then
            dxp=reflsiz
          else
            dxp=dx
          endif
        else
          dxp=dx
        endif
        if (dy.gt.reflsiz) then
          if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then
            dyp=reflsiz
          else
            dyp=dy
          endif
        else
          dyp=dy
        endif
        r1x=xc-dble(dxp)/2.-xn ! computation of the composante along x of the first vector.
        r1y=yc+dble(dyp)/2.-yn ! computation of the composante along y of the first vector.
        r1z=zc-tan(dble(epsilx))*dble(dxp)/2.+tan(dble(epsily))*dble(dyp)/2.-zn ! computation of the composante en z of the first vector.
        r2x=xc+dble(dxp)/2.-xn ! computation of the composante along x of the second vector.
        r2y=yc+dble(dyp)/2.-yn ! computation of the composante along y of the second vector.
        r2z=zc+tan(dble(epsilx))*dble(dxp)/2.+tan(dble(epsily))*dble(dyp)/2.-zn ! computation of the composante en z of the second vector.
        r3x=xc-dble(dxp)/2.-xn ! computation of the composante along x of the third vector.
        r3y=yc-dble(dyp)/2.-yn ! computation of the composante along y of the third vector.
        r3z=zc-tan(dble(epsilx))*dble(dxp)/2.-tan(dble(epsily))*dble(dyp)/2.-zn ! computation of the composante en z of the third vector.
        r4x=xc+dble(dxp)/2.-xn ! computation of the composante along x of the fourth vector.
        r4y=yc-dble(dyp)/2.-yn ! computation of the composante along y of the fourth vector.
        r4z=zc+tan(dble(epsilx))*dble(dxp)/2.-tan(dble(epsily))*dble(dyp)/2.-zn ! computation of the composante en z of the fourth vector.
        call anglesolide(omeg2,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z) ! Call of the routine anglesolide to compute the angle solide.
        if (omeg2.lt.0.D0) then
          print*,'ERROR: Solid angle of the reflecting surface < 0.'
          stop
        endif
        ! estimation of the half of the underlying angle of the solid angle ! this angle servira a obtenir un meilleur isime (moyenne) of
        ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
        ouvang=sqrt(omeg2/pi) ! Angle in radian.
        ouvang=ouvang*180.D0/pi ! Angle in degrees.
        ! computation of the photometric function of the light fixture toward the reflection surface
        anglez=idnint(180.D0*zenith/pi)
        if (anglez.lt.0) anglez=-anglez
        if (anglez.gt.180) anglez=360-anglez
        anglez=anglez+1 ! Transform the angle in integer degree into the position in the array.
        ! average +- ouvang
        naz=0
        nbang=0.D0
        P_indir=0.D0
        do na=-idnint(ouvang),idnint(ouvang)
          naz=anglez+na
          if (naz.lt.0) naz=-naz
          if (naz.gt.181) naz=362-naz ! symetric function
          if (naz.eq.0) naz=1
          P_indir=P_indir+pvalno(naz,stype)*abs(sin(pi*real(naz)/180.D0))/2.
          nbang=nbang+1.*abs(sin(pi*real(naz)/180.D0))/2.
        enddo
        P_indir=P_indir/nbang 
        flrefl=lamplu(x_s,y_s,stype)*P_indir*omeg2*transm*transa*transl
        ! computation of the reflected intensity leaving the ground surface
        irefl=flrefl*srefl/pi ! The factor 1/pi comes from the normalisation of the fonction 
        ! computing the scattering probability toward 2nd scat voxel
        call anglezenithal(rx_sr,ry_sr,z_sr,rx_c,ry_c,z_c,zenith)
        call blocking(x_sr,y_sr,z_sr,x_c,y_c,z_c,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)
        ! computation of the transmittance between the reflection surface and the 3rd scattering voxel
        distd=sqrt((rx_c-rx_sr)**2.+(ry_c-ry_sr)**2.+(z_c-z_sr)**2.)
        call transmitm(zenith,z_sr,z_c,distd,transm,tranam,tabs)
        call transmita(zenith,z_sr,z_c,distd,haer,transa,tranaa)
        call transmita(zenith,z_sr,z_c,distd,hlay,transl,tranal)
        ! computation of the solid angle of the 2nd scattering voxel seen from the reflecting surface
        omega=1./distd**2.
        if (omega.gt.omemax) omega=0.D0   
        flux=irefl*omega*transm*transa*transl*(1.-ff1)*(1.-ff2)*hh
        !flux=irefl*omega*transm*transa*transl*(1.-ff1)*hh              
        if (omega.ne.0.D0) then
          call angle3points(rx_sr,ry_sr,z_sr,rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,angdif) ! scattering angle.
          call diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil,fdifan,fdifl,haer,hlay,pdif,z_c) ! scattering probability of the direct light.
        else
          pdif=0.D0
        endif
        idif=flux*pdif*scal*portio        
 
        
        
        
        if (cloudt.ne.0) then ! line of sight voxel = cloud
          if (cloudbase-z_c.le.iz*scal) then ! this is the cloud base interface
            call anglezenithal(rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,azcl1) ! zenith angle from cloud to observer
            call anglezenithal(rx_c,ry_c,z_c,rx_sr,ry_sr,z_sr,azcl2) ! zenith angle from source to cloud
            doc2=(rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.+(z_c-z_obs)**2.
            dsc2=(rx_sr-rx_c)**2.+(ry_sr-ry_c)**2.+(z_sr-z_c)**2.
            call cloudreflectance(zenith,cloudt,rcloud) ! cloud intensity from direct illum
            ! computing the cloud intensity toward the observer
            icloud=flux/omega*rcloud*doc2*omefov*abs(cos(azcl2)/cos(azcl1))/dsc2/pi ! return to main
        endif
        endif                 
      endif ! end of source vs ground cell cases
      
      

      
      
      return
      end
