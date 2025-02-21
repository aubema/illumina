!-----------------------------------------------------------------------
!
!=======================================================================
! Routine thirdscat
! Calculation of the third scattering 
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
      subroutine thirdscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,iz,lamplu,ofill,srefl,drefle, &
      reflsiz,obsH,altsol,inclix,incliy,pvalno,stype,zondif,siz,volu,ndif,tranam,tabs,tranaa,tranal,secdif,secdil, &
      fdifan,fdifl,haer,hlay,dx,dy,nbx,nby,cloudt,cloudbase,omefov,scal,portio,idift,icloud)
      integer width,nzon
      real*8 pi                                                 
      parameter (width=512,nzon=256)
      parameter (pi=3.141592654D0)
      integer rho,x_s,y_s,x_c,y_c,x_sr,y_sr,ndif,is3a,is3b,jda,jdb,ida,idb,stype,anglez,na,naz,cloudt
      integer x_obs,y_obs,nbx,nby
      real*8 z_s,z_c,z_sr,lamplu(width,width,nzon),zondif(3000000,3),siz,tranam,tranaa,tranal,secdif,secdil
      real*8 haer,hlay,dx,dy,cloudbase,volu,icloud,ofill(width,width),drefle(width,width),obsH(width,width)
      real*8 altsol(width,width),fdifan(181),fdifl(181),pvalno(181,nzon)
      real*8 un,pdif,omega,omeg2,omemax,hh,ff1,ff2,zenith,flux,rx_difa,ry_difa,rx_difb,ry_difb,z_difa,z_difb,dxp,dyp
      real*8 idif,idift,rx_sr,ry_sr,rx_s,ry_s,rx_c,ry_c,dss,ds1,ds2,ds3,ds4,distd,transm,tabs,transa
      real*8 transl,P_dif,angdif,inclix(width,width),incliy(width,width),epsilx,epsily,reflsiz
      real*8 ouvang,P_indir,flrefl,irefl,iz,rcloud,doc2,dsc2,omefov,azcl1,azcl2
      real*8 rx_obs,ry_obs,z_obs,srefl,nbang,scal,portio
      real*8 xc,yc,zc,xn,yn,zn                                            ! Position (meter) of the elements (starting point, final point) for the calculation of the solid angle.
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! Components of the vectors used in the solid angle calculation routine.    
      un=1.D0
      ! omemax: exclude calculations too close (<10m) this is a sustended angle of 1 deg.
      ! the calculated flux is highly sensitive to that number for a very high
      ! pixel resolution (a few 10th of meters). We assume anyway that somebody
      ! observing the sky will never lies closer than that distance to a
      ! light fixture. This number is however somehow subjective and that means
      ! that the value of sky brightness near sources will be affected by this
      ! choice
      omemax=1.D0/((10.D0)**2.D0)
      rx_sr=dble(x_sr)*dx
      ry_sr=dble(y_sr)*dy
      rx_s=dble(x_s)*dx
      ry_s=dble(y_s)*dy
      rx_c=dble(x_c)*dx
      ry_c=dble(y_c)*dy
      rx_obs=dble(x_obs)*dx
      ry_obs=dble(y_obs)*dy         
      dss=1.D0*siz/2.
      !dss=1000.
      icloud=0.D0
      idif=0.D0
      idift=0.D0
      do is3a=1,ndif ! beginning of the loop over the scattering voxels.
        rx_difa=zondif(is3a,1)
        ry_difa=zondif(is3a,2)
        ida=idnint(rx_difa/dx)
        jda=idnint(ry_difa/dy)
        if (ida.gt.width) ida=width
        if (ida.lt.1) ida=1
        if (jda.gt.width) jda=width
        if (jda.lt.1) jda=1
        z_difa=zondif(is3a,3)                                     
        do is3b=1,ndif
          rx_difb=zondif(is3b,1)
          ry_difb=zondif(is3b,2)
          idb=idnint(rx_difb/dx)          
          jdb=idnint(ry_difb/dy)
          if (idb.gt.width) idb=width
          if (idb.lt.1) idb=1
          if (jdb.gt.width) jdb=width
          if (jdb.lt.1) jdb=1
          z_difb=zondif(is3b,3)
          if ((ida.ne.idb).and.(jda.ne.jdb).and.(z_difa.ne.z_difb)) then
            if (((z_difa-siz/2..le.altsol(ida,jda)).or.(z_difa.gt.35000.).or.(z_difa.gt.cloudbase)).or.  &
            ((z_difb-siz/2..le.altsol(idb,jdb)).or.(z_difb.gt.35000.).or.(z_difb.gt.cloudbase)))  then
            else        
              ds1=dsqrt((rx_sr-rx_difa)**2.+(ry_sr-ry_difa)**2.+(z_sr-z_difa)**2.)
              ds2=dsqrt((rx_difb-rx_difa)**2.+(ry_difb-ry_difa)**2.+(z_difb-z_difa)**2.)
              ds3=dsqrt((rx_s-rx_difa)**2.+(ry_s-ry_difa)**2.+(z_s-z_difa)**2.)
              ds4=dsqrt((rx_c-rx_difb)**2.+(ry_c-ry_difb)**2.+(z_c-z_difb)**2.)
              if (rho.eq.0) then ! from source
                if ((ds2.lt.dss).or.(ds3.lt.dss).or.(ds2.lt.dss)) then
                  idif=0.D0
                else
                  call anglezenithal(rx_s,ry_s,z_s,rx_difa,ry_difa,z_difa,zenith)
                  call blocking(x_s,y_s,z_s,ida,jda,z_difa,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)
                  ! computation of the transmittance between the reflection surface and the 3rd scattering voxel
                  distd=dsqrt((rx_difa-rx_s)**2.+(ry_difa-ry_s)**2.+(z_difa-z_s)**2.)
                  call transmitm(zenith,z_s,z_difa,distd,transm,tranam,tabs)
                  call transmita(zenith,z_s,z_difa,distd,haer,transa,tranaa)
                  call transmita(zenith,z_s,z_difa,distd,hlay,transl,tranal)
                  ! computation of the solid angle of the 3rd scattering voxel seen from the reflecting surface
                  omega=1./distd**2.
                  if (omega.gt.omemax) omega=0.
                  anglez=idnint(180.*zenith/pi)+1
                  P_dif=pvalno(anglez,stype)
                  flux=lamplu(x_s,y_s,stype)*P_dif*omega*transm*transa*transl*(1.-ff1)*(1.-ff2)*hh
                  !flux=lamplu(x_s,y_s,stype)*P_dif*omega*transm*transa*transl*(1.-ff1)*hh
                  ! computing the scattering probability toward 2nd scat voxel
                  if (omega.ne.0.) then
                    call angle3points(rx_s,ry_s,z_s,rx_difa,ry_difa,z_difa,rx_difb,ry_difb,z_difb,angdif) ! scattering angle.
                    call diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil,fdifan,fdifl,haer,hlay,pdif,z_difa) ! scattering probability pdif of the direct light.
                  else
                    pdif=0.
                  endif
                  idif=flux*pdif*volu
                  ! AJOUTER LES NUAGES ICI
                  
                  
                  
                  
                endif
              else ! from ground
                if ((ds1.lt.dss).or.(ds2.lt.dss).or.(ds3.lt.dss)) then
                  idif=0.
                else
                  call anglezenithal(rx_sr,ry_sr,z_sr,rx_difa,ry_difa,z_difa,zenith)
                  !call blocking(x_sr,y_sr,z_sr,ida,jda,z_difa,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)
                  ! computation of the transmittance between the reflection surface and the 3rd scattering voxel
                  distd=dsqrt((rx_difa-rx_sr)**2.+(ry_difa-ry_sr)**2.+(z_difa-z_sr)**2.)
                  call transmitm(zenith,z_sr,z_difa,distd,transm,tranam,tabs)
                  call transmita(zenith,z_sr,z_difa,distd,haer,transa,tranaa)
                  call transmita(zenith,z_sr,z_difa,distd,hlay,transl,tranal)
                  ! computation of the solid angle of the 3rd scattering voxel seen from the reflecting surface
                  omega=1./distd**2.
                  if (omega.gt.omemax) omega=0.
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
                  r1x=xc-dxp/2.-xn ! computation of the composante along x of the first vector.
                  r1y=yc+dyp/2.-yn ! computation of the composante along y of the first vector.
                  r1z=zc-dtan(epsilx)*dxp/2.+dtan(epsily)*dyp/2.-zn ! computation of the composante en z of the first vector.
                  r2x=xc+dxp/2.-xn ! computation of the composante along x of the second vector.
                  r2y=yc+dyp/2.-yn ! computation of the composante along y of the second vector.
                  r2z=zc+dtan(epsilx)*dxp/2.+dtan(epsily)*dyp/2.-zn ! computation of the composante en z of the second vector.
                  r3x=xc-dxp/2.-xn ! computation of the composante along x of the third vector.
                  r3y=yc-dyp/2.-yn ! computation of the composante along y of the third vector.
                  r3z=zc-dtan(epsilx)*dxp/2.-dtan(epsily)*dyp/2.-zn ! computation of the composante en z of the third vector.
                  r4x=xc+dxp/2.-xn ! computation of the composante along x of the fourth vector.
                  r4y=yc-dyp/2.-yn ! computation of the composante along y of the fourth vector.
                  r4z=zc+dtan(epsilx)*dxp/2.-dtan(epsily)*dyp/2.-zn ! computation of the composante en z of the fourth vector.
                  call anglesolide(omeg2,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z) ! Call of the routine anglesolide to compute the angle solide.
                  if (omeg2.lt.0.) then
                    print*,'ERROR: Solid angle of the reflecting surface < 0.'
                    stop
                  endif
                  ! estimation of the half of the underlying angle of the solid angle ! this angle servira a obtenir un meilleur isime (moyenne) of
                  ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                  ouvang=dsqrt(omeg2/pi) ! Angle in radian.
                  ouvang=ouvang*180./pi ! Angle in degrees.
                  ! computation of the photometric function of the light fixture toward the reflection surface
                  anglez=idnint(180.*zenith/pi)+1
                  if (anglez.lt.0) anglez=-anglez
                  if (anglez.gt.180) anglez=360-anglez
                  ! average +- ouvang
                  naz=0
                  nbang=0.
                  P_indir=0.
                  do na=-idnint(ouvang),idnint(ouvang)
                    naz=anglez+na
                    if (naz.lt.0) naz=-naz
                    if (naz.gt.181) naz=362-naz ! symetric function
                    if (naz.eq.0) naz=1
                    P_indir=P_indir+pvalno(naz,stype)*dabs(dsin(pi*dble(naz)/180.))/2.
                    nbang=nbang+1.*dabs(dsin(pi*dble(naz)/180.))/2.
                  enddo
                  P_indir=P_indir/nbang 
                  flrefl=lamplu(x_s,y_s,stype)*P_indir*omeg2*transm*transa*transl
                  ! computation of the reflected intensity leaving the ground surface
                  irefl=flrefl*srefl/pi ! The factor 1/pi comes from the normalisation of the fonction 
                  ! computing the scattering probability toward 2nd scat voxel
                  call anglezenithal(rx_sr,ry_sr,z_sr,rx_difa,ry_difa,z_difa,zenith)
                  call blocking(x_sr,y_sr,z_sr,ida,jda,z_difa,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)
                  ! computation of the transmittance between the reflection surface and the 3rd scattering voxel
                  distd=dsqrt((rx_difa-rx_sr)**2.+(ry_difa-ry_sr)**2.+(z_difa-z_sr)**2.)
                  call transmitm(zenith,z_sr,z_difa,distd,transm,tranam,tabs)
                  call transmita(zenith,z_sr,z_difa,distd,haer,transa,tranaa)
                  call transmita(zenith,z_sr,z_difa,distd,hlay,transl,tranal)
                  ! computation of the solid angle of the 3rd scattering voxel seen from the reflecting surface
                  omega=1./distd**2.
                  if (omega.gt.omemax) omega=0.   
                  flux=irefl*omega*transm*transa*transl*(1.-ff1)*(1.-ff2)*hh
                  !flux=irefl*omega*transm*transa*transl*(1.-ff1)*hh
                  if (omega.ne.0.) then
                    call angle3points(rx_sr,ry_sr,z_sr,rx_difa,ry_difa,z_difa,rx_difb,ry_difb,z_difb,angdif) ! scattering angle.
                    call diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil,fdifan,fdifl,haer,hlay,pdif,z_difa) ! scattering probability of the direct light.
                  else
                    pdif=0.
                  endif
                  idif=flux*pdif*volu
                endif
              endif ! end of source vs ground cell cases
              ! computing scattered intensity at the 2nd scat voxel toward teh 3rd scat voxel
              call anglezenithal(rx_difa,ry_difa,z_difa,rx_difb,ry_difb,z_difb,zenith)
              call blocking(ida,jda,z_difa,idb,jdb,z_difb,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)            
              ! computing transmittance between the 3rd scattering voxel 2nd
              distd=dsqrt((rx_difa-rx_difb)**2.+(ry_difa-ry_difb)**2.+(z_difa-z_difb)**2.)
              call transmitm(zenith,z_difa,z_difb,distd,transm,tranam,tabs)
              call transmita(zenith,z_difa,z_difb,distd,haer,transa,tranaa)
              call transmita(zenith,z_difa,z_difb,distd,hlay,transl,tranal)
              ! computing the solid angle of the 2nd scat voxel as seen from the 3rdscattering voxel
              omega=1./distd**2.
              if (omega.gt.omemax) omega=0.
              ! computation of the scattered flux reaching 2nd scattering voxel
              flux=idif*omega*transm*transa*transl*(1.-ff1)*(1.-ff2)*hh
              !flux=idif*omega*transm*transa*transl*(1.-ff1)*hh
              ! computing the scattering probability toward line of sight voxel
              if (omega.ne.0.) then
                call angle3points (rx_difa,ry_difa,z_difa,rx_difb,ry_difb,z_difb,rx_c,ry_c,z_c,angdif) ! scattering angle
                call diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil,fdifan,fdifl,haer,hlay,pdif,z_difb) ! scattering probability 
              else
                pdif=0.
              endif
              ! computing scattered intensity at the 3rd scat voxel toward the line of sight voxel
              idif=flux*pdif*volu
              call anglezenithal(rx_difb,ry_difb,z_difb,rx_c,ry_c,z_c,zenith)
              call blocking(idb,jdb,z_difb,x_c,y_c,z_c,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)            
              ! computing transmittance between the 2nd scattering voxel line of sight
              distd=dsqrt((rx_difb-rx_c)**2.+(ry_difb-ry_c)**2.+(z_difb-z_c)**2.)
              call transmitm(zenith,z_difb,z_c,distd,transm,tranam,tabs)
              call transmita(zenith,z_difb,z_c,distd,haer,transa,tranaa)
              call transmita(zenith,z_difb,z_c,distd,hlay,transl,tranal)
              ! computing the solid angle of the 2nd scat voxel as seen from the 3rdscattering voxel
              omega=1./distd**2.
              if (omega.gt.omemax) omega=0.
              ! computation of the scattered flux reaching 2nd scattering voxel
              flux=idif*omega*transm*transa*transl*(1.-ff1)*(1.-ff2)*hh    
              !flux=idif*omega*transm*transa*transl*(1.-ff1)*hh        
              ! computation of the scattering probability of the scattered light toward line of sight voxel
              if (omega.ne.0.) then
                call angle3points(rx_difb,ry_difb,z_difb,rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,angdif) ! scattering angle.
                call diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil,fdifan,fdifl,haer,hlay,pdif,z_c) ! scattering probability of the direct light.
              else
                pdif=0.
              endif
              ! computing scattered intensity at the line of sight voxel toward the observer 
              idift=idift+flux*pdif*scal*portio


              if (cloudt.ne.0) then ! line of sight voxel = cloud
                if (cloudbase-z_c.le.iz*scal) then ! this is the cloud base interface
                  call anglezenithal(rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,azcl1) ! zenith angle from cloud to observer
                  call anglezenithal(rx_c,ry_c,z_c,rx_difb,ry_difb,z_difb,azcl2) ! zenith angle from source to cloud
                  doc2=(rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.+(z_c-z_obs)**2.
                  dsc2=(rx_difb-rx_c)**2.+(ry_difb-ry_c)**2.+(z_difb-z_c)**2.
                  call cloudreflectance(zenith,cloudt,rcloud) ! cloud intensity from direct illum
                  ! computing the cloud intensity toward the observer
                  icloud=icloud+flux/omega*rcloud*doc2*omefov*dabs(dcos(azcl2)/dcos(azcl1))/dsc2/pi ! return to main
                endif
              endif
              
              
            endif ! above ground below clouds and below 35km
            

            
            
          endif ! 2nd and 3rd scat cells not equal
        enddo ! end scattering volume b
      enddo ! end scattering volume a
      return
      end
