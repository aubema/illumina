!     **********************************************************************************************************************
!     ** Illumina - health in Fortran 77                                                                                  **
!     ** Programmers in decreasing order of contribution  :                                                               **
!     **                            Martin Aube                                                                           **
!     **                                                                                                                  **
!     ** Illumina can be downloaded via:                                                                                  **
!     **    git clone https://github.com/aubema/illumina.git                                                              **
!     **    git checkout illum-health                                                                                     **
!     ** To compile:                                                                                                      **
!     **    cd git/illumina                                                                                               **
!     **    bash bin/makeILLUMINA                                                                                         **
!     **                                                                                                                  **
!     **  Current version features/limitations :                                                                          **
!     **                                                                                                                  **
!     **    - Lambertian reflexion on the ground                                                                          **
!     **    - Terrain slope considered (apparent surface and shadows)                                                     **
!     **    - Angular photometry of a lamp is considered assuming ies file orientation aligned with nearest street        **
!     **    - Sub-grid obstacles considered (with the mean free path of light toward ground, mean obstacle height, and    **
!     **      obstacles transparency (filling factor) The typical mean free path is the distance between two streets      **
!     **    - Accounting for heterogeneity luminaires number, luminaires heights, luminaire spectrum,                     **
!     **      angular photometry, obstacle properties                                                                     **
!     **    - Wavelength dependant                                                                                        **
!     **    - Consider reflectance of streets and building facades and other.                                             **
!     **    - Support direct observation of a source                                                                      **
!     **    - Direct observation of the ground and facades are implemented                                                **
!     **                                                                                                                  **
!     **********************************************************************************************************************
!
!     Copyright (C) 2022 Martin Aube PhD
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Publi! License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Publi! License for more details.
!
!     You should have received a copy of the GNU General Publi! License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!     Contact: martin.aube@cegepsherbrooke.qc.ca
!
!
!
program illumhealth                           ! Beginning
  implicit none
  !
  !     Variables declaration
  !
  INTEGER width,ntype                               ! ntype Nb of light source types (different combination of spectrum and uplight)
  !     Matrix dimensions
  ! parameter (width=1024,ntype=36)
  real pi
  INTEGER verbose                                   ! verbose = 1 to have more print out, 0 for silent
  parameter (pi=3.141592654)
  character*72 mnaf                                 ! Terrain elevation file
  character*72 outfile                              ! Results file
  character*72 basenm                               ! Base name of files
  INTEGER lenbase                                   ! Length of the Base name of the experiment
  real lambda,pressi                                ! Wavelength (nanometer), atmospheri! pressure (kPa)
  real pvalto                                       !
  
  real, dimension(:,:,:), allocatable :: pval       ! pval(181,181,ntype)
  real, dimension(:,:,:), allocatable :: pvalno     ! pvalno(181,181,ntype) values of the angular photometry functions 1sr_axis=azim, 2nd=zenith 
  real, dimension(:,:,:), allocatable :: lamplu     ! lamplu(width,width,ntype)  Source fluxes
  real, dimension(:,:), allocatable :: obsD         ! obsD(width,width)  mean free path to the ground (meter)
  real, dimension(:,:), allocatable ::  val2d       ! val2d(width,width)  Temporary input array 2d
  real, dimension(:,:), allocatable ::  altsol      ! altsol(width,width)  Ground elevation (meter)
  real, dimension(:,:), allocatable ::  altsob      ! altsob(width,width)  Elevation including buildings
  real, dimension(:,:), allocatable ::  reflec      ! reflec(width,width)  reflectance  
  real, dimension(:,:), allocatable ::  lampal      ! lampal(width,width) Height of the light sources relative to the ground (meter)
  real, dimension(:,:), allocatable ::  inclix      ! inclix(width,width)  tilt of the ground pixel along x (radian)
  real, dimension(:,:), allocatable ::  incliy      ! incliy(width,width)  tilt of the ground pixel along y (radian)  
  real, dimension(:,:), allocatable ::  obsH        ! obsH(width,width)  averaged height of the sub-grid obstacles
  real, dimension(:,:), allocatable ::  ofill       ! ofill(width,width)  fill factor giving the probability to hit an obstacle when pointing in its direction real 0-1
  real, dimension(:,:), allocatable ::  irradi      ! irradi(width,width)  direct irradiance on a surface normal to the line of sight (no scattering)
  INTEGER, dimension(:,:), allocatable ::  gndty    ! gndty(width,width)  ground type flag  0=observer, 1=building facade, 2=street, 3=other
  real, dimension(:,:), allocatable ::  azims       ! azims(width,width)  azimuth toward nearest point on a street
  INTEGER, dimension(:), allocatable ::  imin       ! imin(ntype)
  INTEGER, dimension(:), allocatable :: imax        ! imax(ntype)
  INTEGER, dimension(:), allocatable :: jmin        ! jmin(ntype)
  INTEGER, dimension(:), allocatable :: jmax        ! jmax(ntype)     x and y limits containing a type of lamp
  real, dimension(:), allocatable :: totlu          ! totlu(ntype)   total flux of a source type
  
  
  real angmin                                 ! minimum angle
  real reflsiz                                ! Size of the reflecting surface
  real largx                                  ! Width (x axis) of the modeling domain (meter)
  real largy                                  ! Length (y axis) of the modeling domain (meter)
  INTEGER nbx,nby                             ! Number of pixels in the modeling domain

  real srefl,brefl,orefl                      ! Street reflectance , Building facades reflectance , Other reflectance
  INTEGER stype                               ! Source type or zone index
  character*72 pafile,lufile,alfile,ohfile,odfile,offile,azfile 
  ! Files related to light sources and obstacles (photometri
  ! function of the sources (sr-1), flux (W), height (m), obstacles           
  ! height (m), obstacle distance (m), obstacle filling factor (0-1).

  real dtheta               ! Angle increment of the photometri! function of the sources
  real dx,dy,dxp,dyp        ! Width of the voxel (meter)
  INTEGER boxx,boxy         ! reflection window size (pixels)
  INTEGER x_obs,y_obs       ! Position of the observer (INTEGER)
  real rx_obs,ry_obs
  real z_o                  ! observer height relative to the ground (meter)
  real z_obs                ! Height of the observer (meter) to the vertical grid scale
  INTEGER x_s,y_s,x_sr,y_sr ! Positions of the source, the reflecting surface
  real z_s,z_sr,z_dif       ! Heights of the source, the reflecting surface, and the scattering voxel (metre).
  real rx_s,ry_s,rx_sr,ry_sr
  real angzen,ouvang        ! Zenithal angle between two voxels (radians) and opening angle of the solid angle in degrees.
  INTEGER anglez            ! Emitting zenithal angle from the luminaire.
  INTEGER anglep            ! Emitting azimuth angle from the luminaire
  real P_dir,P_indir        ! photometri! function of the light sources (direct,reflected)
  real*8 xc,yc,zc,xn,yn,zn  ! Position (meter) of the elements (starting point, final point) for the calculation of the solid angle.
  real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z ! Components of the vectors used in the solid angle calculation routine.
  real omega                ! Solid angles
  real haut                 ! Haut (negative indicate that the surface is lighted from inside the ground. I.e. not considered in the calculation
  real epsilx,epsily        ! tilt of the ground pixel
  real flrefl               ! flux reaching a reflecting surface (watts).
  real irefl,irefl1         ! intensity leaving a reflecting surface toward the line of sight voxel.
  real angvis,azim          ! viewing angles of the sensor.
  real nbang                ! for the averaging of the photometri! function


  INTEGER naz,na,nap,np
  character*3 lampno        ! lamp number string

  real angazi               ! azimuth angle between two points in rad, max dist for the horizon determination
  real latitu               ! approximate latitude of the domain center
  INTEGER xsrmi,xsrma,ysrmi,ysrma ! limits of the loop valeur for the reflecting surfaces

  real ff,ff2,hh            ! temporary obstacle filling factor and horizon blocking factor
  real distd                ! distance to compute the scattering probability
  real angvi1,angaz1,angze1 ! viewing angles in radian
  real ix,iy,iz             ! base vector of the viewing (length=1)
  real azcl1,azcl2          ! zenith angle from the (source, refl surface, or scattering voxel) to line of path and observer to line p.
  real dh,dho               ! distance of the horizon limit
  INTEGER n2nd              ! desired number of voxel in the calculation of the 2nd scattering
  INTEGER step              ! skiping 2nd scat on 1 dim
  INTEGER i,j,k,id,jd
  real tranam,tranaa,tranal ! TOA atmospheri! transmittancess of a path (molecular, aerosol, layer)
  real transa,transl,transm
  real zhoriz               ! zenith angle of the horizon

  real dang                 ! Angle between the line of sight and the direction of a source
  real dzen                 ! zenith angle of the source-observer line
  real ddir_obs             ! distance between the source and the observer
  real rx,ry,rz             ! driving vector for the calculation of the projection angle for direct radiance. It is 20km long

  character*72 gifile       ! name of the ground type flag file
  real dh0,dhmax            ! horizontal distance along the line of sight and maximum distance before beeing blocked by topography
  real layaod               ! 500 nm aod of the particle layer
  real hlay                 ! exponential vertical scale height of the particle layer
  real haer                 ! exponential vertical scale height of the background aerosol layer
  real bandw                ! bandwidth of the spectral bin
  real tabs                 ! TOA transmittance related to molecule absorption
  INTEGER oi,oj             ! scanning position of observers
  real taua,alpha           ! aerosol layer optical properties
  verbose=1                 ! Very little printout=0, Many printout = 1, even more=2
  ff=0.
  ff2=0.
  step=1
  if (verbose.ge.1) then
     print*,'Starting ILLUMINA computations...'
  endif
  !     reading of the input file (illum_health.in)
  print*,'Reading illum-health.in input file'
  open(unit=1,file='illum-health.in',status='old')
  read(1,*)
  read(1,*) basenm
  read(1,*) dx,dy
  read(1,*) ntype
  read(1,*) width
  read(1,*) taua,alpha,haer
  read(1,*) lambda,bandw
  read(1,*) srefl,brefl,orefl
  read(1,*) pressi
  read(1,*) z_o
  close(1)
  ! allocate arrays 
  allocate ( pval(181,181,ntype) )
  allocate ( pvalno(181,181,ntype) )
  allocate ( lamplu(width,width,ntype) )
  allocate ( obsD(width,width) )
  allocate ( val2d(width,width) )
  allocate ( altsol(width,width) )
  allocate ( altsob(width,width) )
  allocate ( reflec(width,width) )
  allocate ( lampal(width,width) )
  allocate ( inclix(width,width) )
  allocate ( incliy(width,width) )
  allocate ( obsH(width,width) )
  allocate ( ofill(width,width) )
  allocate ( irradi(width,width) )
  allocate ( gndty(width,width) )
  allocate ( azims(width,width) )
  allocate ( imin(ntype) )
  allocate ( imax(ntype) )
  allocate ( jmin(ntype) )
  allocate ( jmax(ntype) )
  allocate ( totlu(ntype) )
  
  
  
  
  !     windows pointing toward horizon
  angvis=0.
  angvi1 = (pi*angvis)/180.
  angze1 = pi/2.-angvi1
  !     max distance to consider the reflexion in meter
  reflsiz=30.
  boxx=nint(reflsiz/dx)     ! Number of column to consider left/right of the source for the reflection.
  boxy=nint(reflsiz/dy)     ! Number of column to consider up/down of the source for the reflection.
  if (verbose.gt.0) then
     print*,'Pixel size = ',dx,' x ',dy
  endif
  !     computing the actual AOD at the wavelength lambda
  if (verbose.ge.1) print*,'500nm AOD=',taua,'500nm angstrom coeff.=',alpha
  taua=taua*(lambda/500.)**(-1.*alpha)
  print*,'Wavelength (nm):',lambda,' Aerosol optical depth:',taua
  print*,'Elevation angle:',angvis,' azim angle (counterclockwise from east)',azim
  !     Initialisation of arrays and variables
  if (verbose.ge.1) print*,'Initializing variables...'
  do i=1,width
     do j=1,width
        val2d(i,j)=0.
        altsol(i,j)=0.
        altsob(i,j)=0.
        obsH(i,j)=0.
        gndty(i,j)=0
        ofill(i,j)=0.
        inclix(i,j)=0.
        incliy(i,j)=0.
        azims(i,j)=0.
        lampal(i,j)=0.
        do k=1,ntype
           lamplu(i,j,k)=0.
        enddo
     enddo
  enddo
  do i=1,181
     do j=1,181
        do k=1,ntype
           pval(i,j,k)=0.
           pvalno(i,j,k)=0.
        enddo
     enddo
  enddo
  angmin=0.
  flrefl=0.
  irefl=0.
  irefl1=0.
  hh=1.
  !     reading of the environment variables
  !     determine the Length of basenm
  lenbase=index(basenm,' ')-1
  outfile=basenm(1:lenbase)//'.bin'
  mnaf=basenm(1:lenbase)//'_topogra.bin' ! determine the names of input and output files
  ohfile=basenm(1:lenbase)//'_obsth.bin'
  odfile=basenm(1:lenbase)//'_obstd.bin'
  alfile=basenm(1:lenbase)//'_altlp.bin' ! setting the file name of height of the sources lumineuse.
  offile=basenm(1:lenbase)//'_obstf.bin'
  azfile=basenm(1:lenbase)//'_azimu.bin'
  gifile='groundtype.bin'
  dtheta=.017453293
  !     reading of the elevation file
  call twodin(nbx,nby,mnaf,altsol,width)
  !     reading lamp heights
  call twodin(nbx,nby,alfile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        lampal(i,j)=val2d(i,j) ! filling of the array for the lamp stype
     enddo                  ! end of the loop over all cells along y.
  enddo                     ! end of the loop over all cells along x.
  !     reading subgrid obstacles average height
  call twodin(nbx,nby,ohfile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        obsH(i,j)=val2d(i,j) ! filling of the array
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     reading subgrid obstacles average distance
  call twodin(nbx,nby,odfile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        obsD(i,j)=val2d(i,j)/2.
        if (obsD(i,j).eq.0.) obsD(i,j)=dx ! when outside a zone, block to the size of the cell (typically 1km)
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     reading subgrid obstacles filling factor
  call twodin(nbx,nby,offile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        ofill(i,j)=val2d(i,j) ! Filling of the array 0-1
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     reading ground type flag
  call twodin(nbx,nby,gifile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        gndty(i,j)=nint(val2d(i,j)) ! ground type flag flag array 0 or 1
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     reading azimuth
  call twodin(nbx,nby,azfile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        azims(i,j)=val2d(i,j) ! Filling of the array 0-1
        !     conversion of the geographical viewing angles toward the cartesian
        !     angle we assume that the angle in the file illumina.in
        !     is consistent with the geographical definition
        !     geographical, azim=0 toward north, 90 toward east, 180 toward south
        !     etc
        !     cartesian, azim=0 toward east, 90 toward north, 180 toward west etc
        azims(i,j)=90.-azims(i,j)
        if (azims(i,j).lt.0.) azims(i,j)=azims(i,j)+360.
        if (azims(i,j).ge.360.) azims(i,j)=azims(i,j)-360.
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     Some preliminary tasks
  do stype=1,ntype          ! beginning of the loop over ntype types of sources.
     imin(stype)=nbx
     jmin(stype)=nby
     imax(stype)=1
     jmax(stype)=1
     pvalto=0.
     write(lampno, '(I3.3)' ) stype ! support of ntype different sources (3 digits)
     pafile=basenm(1:lenbase)//'_fctem_'//lampno//'.dat' ! setting the file name of angular photometry.
     lufile=basenm(1:lenbase)//'_lumlp_'//lampno//'.bin' ! setting the file name of the luminosite of the cases.
     !     reading photometry files
     open(UNIT=1, FILE=pafile,status='OLD') ! opening file pa#.dat, angular photometry.
     do i=1,181
        do j=1,181          ! beginning of the loop for the 181 data points
           read(1,*) pval(i,j,stype) ! reading of the data in the array pval.
           pvalto=pvalto+pval(i,j,stype)*2.*sin(real(j-1)*dtheta)*dtheta**2.  ! units of 1/sr.
        enddo
     enddo                  ! end of the loop over the 181 donnees of the fichier pa#.dat.
     close(1)               ! closing file pa#.dat, angular photometry.
     do i=1,181
        do j=1,181
           if (pvalto.ne.0.) pvalno(i,j,stype)=pval(i,j,stype)/pvalto ! Normalisation of the photometri! function.
        enddo
     enddo
     !     reading luminosity files
     call twodin(nbx,nby,lufile,val2d,width)
     do i=1,nbx             ! beginning of the loop over all cells along x.
        do j=1,nby          ! beginning of the loop over all cells along y.
           if (val2d(i,j).lt.0.) then ! searching of negative fluxes
              print*,'***Negative lamp flux!, stopping execution'
              stop
           endif
        enddo               ! end of the loop over all cells along y.
     enddo
     do i=1,nbx             ! searching of the smallest rectangle containing the zone
        do j=1,nby          ! of non-null luminosity to speedup the calculation
           if (val2d(i,j).ne.0.) then
              if (i-1.lt.imin(stype)) imin(stype)=i-2
              if (imin(stype).lt.1) imin(stype)=1
              goto 333
           endif
        enddo
     enddo
     imin(stype)=1
333  do i=nbx,1,-1
        do j=1,nby
           if (val2d(i,j).ne.0.) then
              if (i+1.gt.imax(stype)) imax(stype)=i+2
              if (imax(stype).gt.nbx) imax(stype)=nbx
              goto 334
           endif
        enddo
     enddo
     imax(stype)=1
334  do j=1,nby
        do i=1,nbx
           if (val2d(i,j).ne.0.) then
              if (j-1.lt.jmin(stype)) jmin(stype)=j-2
              if (jmin(stype).lt.1) jmin(stype)=1
              goto 335
           endif
        enddo
     enddo
     jmin(stype)=1
335  do j=nby,1,-1
        do i=1,nbx
           if (val2d(i,j).ne.0.) then
              if (j+1.gt.jmax(stype)) jmax(stype)=j+2
              if (jmax(stype).gt.nby) jmax(stype)=nby
              goto 336
           endif
        enddo
     enddo
     jmax(stype)=1
336  do i=1,nbx             ! beginning of the loop over all cells along x.
        do j=1,nby          ! beginning of the loop over all cells along y.
           lamplu(i,j,stype)=val2d(i,j) ! remplir the array of the lamp type: stype
           totlu(stype)=totlu(stype)+lamplu(i,j,stype) ! the total lamp flux should be non-null to proceed to the calculations
        enddo               ! end of the loop over all cells along y.
     enddo                  ! end of the loop over all cells along x.
  enddo  ! enddo stype
  !     distribute reflectance values and adding fake buildings
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        if (gndty(i,j).eq.2) then
           reflec(i,j)=srefl
        elseif (gndty(i,j).eq.1) then
           reflec(i,j)=srefl
           altsob(i,j)=altsol(i,j)+obsH(i,j)
        elseif ((gndty(i,j).eq.0).or.(gndty(i,j).eq.3)) then
           reflec(i,j)=brefl
           altsob(i,j)=altsol(i,j)+obsH(i,j)/2.
        else
           reflec(i,j)=orefl
           altsob(i,j)=altsol(i,j)
        endif
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     computation of the basi! tilt of the pixels along x and along y
  !
  !     0=building front
  !     1=building top               111
  !     2=street                     111
  !     3=building rear             01113
  !     4=other                     01113
  !     fake building profile   222201113444444444
  do i=1,nbx                ! beginning of the loop over the column (longitude) of the domain.
     do j=1,nby             ! beginning of the loop over the rows (latitu) of the domain.
        if (gndty(i,j).eq.2) then
           if (i.eq.1) then ! specifi! case close to the border of the domain (vertical side left).
              inclix(i,j)=atan((altsol(i+1,j)-altsol(i,j))/real(dx)) ! computation of the tilt along x of the surface.
           elseif (i.eq.nbx) then ! specifi! case close to the border of the domain (vertical side right).
              inclix(i,j)=atan((altsol(i-1,j)-altsol(i,j))/(real(dx))) ! computation of the tilt along x of the surface.
           else
              inclix(i,j)=atan((altsol(i+1,j)-altsol(i-1,j))/(2.*real(dx)))
           endif
           if (j.eq.1) then ! specifi! case close to the border of the domain (horizontal side down).
              incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j))/(real(dy))) ! computation of the tilt along y of the surface.
           elseif (j.eq.nby) then ! specifi! case close to the border of the domain (horizontal side up).
              incliy(i,j)=atan((altsol(i,j-1)-altsol(i,j))/(real(dy))) ! computation of the tilt along y of the surface.
           else
              incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j-1))/(2.*real(dy)))
           endif
        else
           if (i.eq.1) then
              inclix(i,j)=atan((altsob(i+1,j)-altsob(i,j))/real(dx))
           elseif (i.eq.nbx) then
              inclix(i,j)=atan((altsob(i-1,j)-altsob(i,j))/(real(dx)))
           else
              inclix(i,j)=atan((altsob(i+1,j)-altsob(i-1,j))/(2.*real(dx)))
           endif
           if (j.eq.1) then
              incliy(i,j)=atan((altsob(i,j+1)-altsob(i,j))/(real(dy)))
           elseif (j.eq.nby) then
              incliy(i,j)=atan((altsob(i,j-1)-altsob(i,j))/(real(dy)))
           else
              incliy(i,j)=atan((altsob(i,j+1)-altsob(i,j-1))/(2.*real(dy)))
           endif
        endif
     enddo                  ! end of the loop over the rows (latitu) of the domain
  enddo                     ! end of the loop over the column (longitude) of the domain
  largx=dx*real(nbx)        ! computation of the Width along x of the case.
  largy=dy*real(nby)        ! computation of the Width along y of the case.
  !     remove 2nd aerosol layer
  layaod=0.
  hlay=2000.
  !     loop over potential observers
  do oi=1,nbx
     do oj=1,nby
        !     only calculate if it is an observer position
        if (gndty(oi,oj).eq.0) then
           !     convert to radians
           angaz1 = (pi*azims(oi,oj))/180.
           ix = ( sin((pi/2.)-angvi1) ) * (cos(angaz1)) ! viewing vector components
           iy = ( sin((pi/2.)-angvi1) ) * (sin(angaz1))
           iz = (sin(angvi1))
           x_obs=oi
           y_obs=oj
           z_obs=z_o+altsol(x_obs,y_obs)
           rx_obs=real(x_obs)*dx
           ry_obs=real(y_obs)*dy
           if (z_obs.eq.0.) z_obs=0.001
           !     determination of the vertical atmospheri! transmittance
           call transtoa(lambda,bandw,taua,layaod,pressi,tranam,tranaa,tranal,tabs)
           irradi(oi,oj)=0. ! initialize the total direct irradiance from sources to observer
           !     ==============================
           !     Calculation of the direct radiances
           !
           if (verbose.ge.1) print*,' Calculating obtrusive light...'
           !     loop over source types
           do stype=1,ntype
              !     Is any flux in that lamp type?
              if (totlu(stype).ne.0.) then
                 if (verbose.ge.1) print*,' Turning on lamps',stype
                 do x_s=imin(stype),imax(stype) ! beginning of the loop over the column (longitude the) of the domain.
                    do y_s=jmin(stype),jmax(stype) ! beginning of the loop over the rows (latitud) of the domain.
                       rx_s=real(x_s)*dx
                       ry_s=real(y_s)*dy
                       if (lamplu(x_s,y_s,stype).ne.0.) then ! if the luminosite of the case is null, the program ignore this case.
                          z_s=(altsol(x_s,y_s)+lampal(x_s,y_s)) ! Definition of the position (metre) vertical of the source.
                          !
                          !     *********************************************************************************************************
                          !     calculation of the direct radiance of sources falling on a surface perpendicular
                          !     to the viewing angle Units of W/nm/m2/sr
                          !     *********************************************************************************************************
                          !     defining a line of sight wit a point located far away
                          rx=rx_obs+20000.*ix
                          ry=ry_obs+20000.*iy
                          rz=z_obs+20000.*iz
                          dho=sqrt((rx_obs-rx_s)**2.+(ry_obs-ry_s)**2.)
                          if ((dho.gt.0.).and.(z_s.ne.z_obs)) then  ! if same horizontal position we cannot have same z elevation
                             call anglezenithal(rx_obs,ry_obs,z_obs,rx_s,ry_s,z_s,dzen)
                             call angleazimutal(rx_obs,ry_obs,rx_s,ry_s,angazi)
                             if (dzen.gt.pi/4.) then ! 45deg. it is unlikely to have a 1km high mountain less than 1
                                call horizon(width,x_obs,y_obs,z_obs,dx,dy,altsol,angazi,zhoriz,dh)
                                if (dh.le.dho) then
                                   if (dzen-zhoriz.lt.0.00001) then ! shadow the path line of sight-source is not below the horizon => we compute
                                      hh=1.
                                   else
                                      hh=0.
                                   endif
                                else
                                   hh=1.
                                endif
                             else
                                hh=1.
                             endif
                             ff=0.
                             !     sub-grid obstacles
                             if (dho.gt.obsD(x_obs,y_obs)+obsD(x_s,y_s)) then ! light path to observer larger than the mean free path -> subgrid obstacles
                                angmin=pi/2.-atan2((altsol(x_obs,y_obs)+obsH(x_obs,y_obs)-z_obs),obsD(x_obs,y_obs))
                                if (dzen.lt.angmin) then ! condition sub-grid obstacles direct.
                                   ff=0.
                                else
                                   ff=ofill(x_obs,y_obs)
                                endif
                             endif ! end light path to the observer larger than mean free path
                             call anglezenithal(rx_s,ry_s,z_s,rx_obs,ry_obs,z_obs,dzen)
                             ff2=0.
                             if (dho.gt.obsD(x_s,y_s)) then ! light path from source larger than the mean free path -> subgrid obstacles
                                angmin=pi/2.-atan2((altsol(x_s,y_s)+obsH(x_s,y_s)-z_s),obsD(x_s,y_s))
                                if (dzen.lt.angmin) then ! condition sub-grid obstacles direct.
                                   ff2=0.
                                else
                                   ff2=ofill(x_s,y_s)
                                endif
                             endif ! end light path to the observer larger than mean free path
                             call anglezenithal(rx_obs,ry_obs,z_obs,rx_s,ry_s,z_s,dzen)
                             !     projection angle of line to the lamp and the viewing angle
                             call angle3points (rx_s,ry_s,z_s,rx_obs,ry_obs,z_obs,rx,ry,rz,dang)
                             dang=pi-dang
                             !     computation of the solid angle of the line of sight voxel seen from the source
                             anglez=nint(180.*(pi-dzen)/pi)+1
                             anglep=nint(180.*(pi-abs(angazi-angaz1))/pi)+1
                             P_dir=pvalno(anglep,anglez,stype)
                             !     computation of the flux direct reaching the line of sight voxel
                             if ((cos(dang).gt.0.).and.(dang.lt.pi/2.)) then
                                ddir_obs=sqrt((rx_obs-rx_s)**2.+(ry_obs-ry_s)**2.+(z_obs-z_s)**2.)
                                !     computation of the solid angle 1m^2 at the observer as seen from the source
                                omega=1.*abs(cos(dang))/ddir_obs**2.
                                call transmitm(dzen,z_obs,z_s,ddir_obs,transm,tranam,tabs)
                                call transmita(dzen,z_obs,z_s,ddir_obs,haer,transa,tranaa)
                                call transmitl(dzen,z_obs,z_s,ddir_obs,hlay,transl,tranal)
                                irradi(oi,oj)=irradi(oi,oj)+lamplu(x_s,y_s,stype)*transa*transm*transl*P_dir*omega &
                                *(1.-ff)*(1.-ff2)*hh
                             endif
                             !
                             !     **********************************************************************************
                             !     * computation of the direct light toward the observer by the ground reflection   *
                             !     **********************************************************************************
                             !
                             xsrmi=x_s-boxx
                             if (xsrmi.lt.1) xsrmi=1
                             xsrma=x_s+boxx
                              if (xsrma.gt.nbx) xsrma=nbx
                             ysrmi=y_s-boxy
                             if (ysrmi.lt.1) ysrmi=1
                             ysrma=y_s+boxy
                             if (ysrma.gt.nby) ysrma=nby
                             do x_sr=xsrmi,xsrma ! beginning of the loop over the column (longitude) reflecting.
                                rx_sr=real(x_sr)*dx
                                do y_sr=ysrmi,ysrma ! beginning of the loop over the rows (latitu) reflecting.
                                   ry_sr=real(y_sr)*dy
                                   irefl=0.
                                   z_sr=altsol(x_sr,y_sr)
                                   if((x_sr.gt.nbx).or.(x_sr.lt.1).or.(y_sr.gt.nby).or.(y_sr.lt.1)) then
                                      if (verbose.eq.2) then
                                         print*,'Ground cell out of borders'
                                      endif
                                   else
                                      if((x_s.eq.x_sr).and.(y_s.eq.y_sr).and.(z_s.eq.z_sr)) then
                                         if (verbose.eq.2) then
                                            print*,'Source pos = Ground cell'
                                         endif
                                      else
                                         haut=-(rx_s-rx_sr)*tan(inclix(x_sr,y_sr))-(ry_s-ry_sr)*tan(incliy(x_sr,y_sr))+z_s-z_sr
                                         if (haut .gt. 0.) then ! Condition: the ground cell is lighted from above
                                            !     computation of the zenithal angle between the source and the surface reflectance
                                            call anglezenithal(rx_s,ry_s,z_s,rx_sr,ry_sr,z_sr,angzen)
                                            call angleazimutal(rx_s,ry_s,rx_sr,ry_sr,angazi)
                                            !     computation of the transmittance between the source and the ground surface
                                            distd=sqrt((rx_s-rx_sr)**2.+(ry_s-ry_sr)**2.+(z_s-z_sr)**2.)
                                            call transmitm(angzen,z_s,z_sr,distd,transm,tranam,tabs)
                                            call transmita(angzen,z_s,z_sr,distd,haer,transa,tranaa)
                                            call transmitl(angzen,z_s,z_sr,distd,hlay,transl,tranal)
                                            !     computation of the solid angle of the reflecting cell seen from the source
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
                                            r1z=zc-tan(dble(epsilx))*dble(dxp)/2.+tan(dble(epsily))*dble(dyp)/2.-zn
                                            r2x=xc+dble(dxp)/2.-xn ! computation of the composante along x of the second vector.
                                            r2y=yc+dble(dyp)/2.-yn ! computation of the composante along y of the second vector.
                                            r2z=zc+tan(dble(epsilx))*dble(dxp)/2.+tan(dble(epsily))*dble(dyp)/2.-zn
                                            r3x=xc-dble(dxp)/2.-xn ! computation of the composante along x of the third vector.
                                            r3y=yc-dble(dyp)/2.-yn ! computation of the composante along y of the third vector.
                                            r3z=zc-tan(dble(epsilx))*dble(dxp)/2.-tan(dble(epsily))*dble(dyp)/2.-zn
                                            r4x=xc+dble(dxp)/2.-xn ! computation of the composante along x of the fourth vector.
                                            r4y=yc-dble(dyp)/2.-yn ! computation of the composante along y of the fourth vector.
                                            r4z=zc+tan(dble(epsilx))*dble(dxp)/2.-tan(dble(epsily))*dble(dyp)/2.-zn
                                            call anglesolide(omega,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                                            if (omega.lt.0.) then
                                               print*,'ERROR: Solid angle of the reflecting surface < 0.'
                                               stop
                                            endif
                                            !  estimation of the half of the underlying angle of the solid angle
                                            ! this angle servira a obtenir un meilleur isime (moyenne) of
                                            ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                                            ouvang=sqrt(omega/pi) ! Angle in radian.
                                            ouvang=ouvang*180./pi ! Angle in degrees.
                                            !     computation of the photometri! function of the light fixture toward the reflection surface
                                            anglez=nint(180.*angzen/pi)
                                            anglep=nint(180.*(pi-abs(angazi-angaz1))/pi)+1
                                            if (anglez.lt.0) anglez=-anglez
                                            if (anglez.gt.180) anglez=360-anglez
                                            anglez=anglez+1 ! Transform the angle in integer degree into the position in the array.
                                            if (anglep.lt.0) anglep=-anglep
                                            if (anglep.gt.180) anglep=360-anglep
                                            anglep=anglep+1 ! Transform the angle in integer degree into the position in the array.
                                            !     average +- ouvang
                                            naz=0
                                            nap=0
                                            nbang=0.
                                            P_indir=0.
                                            do na=-nint(ouvang),nint(ouvang)
                                               naz=anglez+na
                                               do np=-nint(ouvang),nint(ouvang)
                                                  nap=anglep+np
                                                  if (naz.lt.0) naz=-naz
                                                  if (naz.gt.181) naz=362-naz ! symetri! function
                                                  if (naz.eq.0) naz=1
                                                  if (nap.lt.0) nap=-nap
                                                  if (nap.gt.181) nap=362-nap ! symetri! function
                                                  if (nap.eq.0) nap=1
                                                  P_indir=P_indir+pvalno(nap,naz,stype)*abs(sin(pi*real(naz)/180.))/2.
                                                  nbang=nbang+1.*abs(sin(pi*real(naz)/180.))/2.
                                               enddo
                                            enddo
                                            P_indir=P_indir/nbang
                                            !     computation of the flux reaching the reflecting surface
                                            flrefl=lamplu(x_s,y_s,stype)*P_indir*omega*transm*transa*transl
                                            !     computation of the reflected intensity leaving the ground surface
                                            irefl1=flrefl*reflec(i,j)/pi ! The factor 1/pi comes from the normalisation of the fonction
                                            dho=sqrt((rx_obs-rx_sr)**2.+(ry_obs-ry_sr)**2.)
                                            if ((dho.gt.0.).and.(z_s.ne.z_obs)) then
                                               call anglezenithal(rx_obs,ry_obs,z_obs,rx_sr,ry_sr,z_sr,dzen)
                                               call angleazimutal(rx_obs,ry_obs,rx_sr,ry_sr,angazi)
                                               if (dzen.gt.pi/4.) then ! 45deg. it is unlikely to have a 1km high mountain less than 1
                                                  call horizon(width,x_obs,y_obs,z_obs,dx,dy,altsol,angazi,zhoriz,dh)
                                                  if (dh.le.dho) then
                                                     if (dzen-zhoriz.lt.0.00001) then ! shadow the path line of sight-source is not below the horizon => we compute
                                                        hh=1.
                                                     else
                                                        hh=0.
                                                     endif
                                                  else
                                                     hh=1.
                                                  endif
                                               else
                                                  hh=1.
                                               endif
                                               !     sub-grid obstacles
                                               ff=0.
                                               if (dho.gt.obsD(x_obs,y_obs)+obsD(x_sr,y_sr)) then ! light path to observer larger than the mean free path -> subgrid obstacles
                                                  angmin=pi/2.-atan2((altsol(x_obs,y_obs)+obsH(x_obs,y_obs)-z_obs), &
                                                  obsD(x_obs,y_obs))
                                                  if (dzen.lt.angmin) then ! condition sub-grid obstacles direct.
                                                     ff=0.
                                                  else
                                                     ff=ofill(x_obs,y_obs)
                                                  endif
                                               endif ! end light path to the observer larger than mean free path
                                               call anglezenithal(rx_sr,ry_sr,z_sr,rx_obs,ry_obs,z_obs,dzen)
                                               ff2=0.
                                               if (dho.gt.obsD(x_sr,y_sr)) then ! light path from reflecting surface larger than the mean free path -> subgrid obstacles
                                                  angmin=pi/2.-atan2((altsol(x_sr,y_sr)+obsH(x_sr,y_sr)-z_sr),obsD(x_sr,y_sr))
                                                  if (dzen.lt.angmin) then ! condition sub-grid obstacles direct.
                                                     ff2=0.
                                                  else
                                                     ff2=ofill(x_sr,y_sr)
                                                  endif
                                               endif ! end light path to the observer larger than mean free path
                                               call anglezenithal(rx_obs,ry_obs,z_obs,rx_sr,ry_sr,z_sr,dzen)
                                               !     projection angle of line to the lamp and the viewing angle
                                               call angle3points (rx_sr,ry_sr,z_sr,rx_obs,ry_obs,z_obs,rx,ry,rz,dang)
                                               dang=pi-dang
                                               !     computation of the flux direct reaching the line of sight voxel
                                               if ((cos(dang).gt.0.).and.(dang.lt.pi/2.)) then
                                                  ddir_obs=sqrt((rx_obs-rx_sr)**2.+(ry_obs-ry_sr)**2.+(z_obs-z_sr)**2.)
                                                  !     computation of the solid angle of the line of sight voxel seen from the source
                                                  omega=1.*abs(cos(dang))/ddir_obs**2.
                                                  call transmitm(dzen,z_obs,z_sr,ddir_obs,transm,tranam,tabs)
                                                  call transmita(dzen,z_obs,z_sr,ddir_obs,haer,transa,tranaa)
                                                  call transmitl(dzen,z_obs,z_sr,ddir_obs,hlay,transl,tranal)
                                                  irradi(oi,oj)=irradi(oi,oj)+irefl1*omega*transa*transm*transl*hh*(1.-ff)*(1.-ff2)
                                               endif
                                            endif
                                         endif
                                      endif
                                   endif !     end if inside borders
                                enddo  !     loops x_sr y_sr
                             enddo
                          endif ! source and observer not at the same exact position
                       endif   !     end if lamplu ne 0
                    enddo ! loop over source of stype
                 enddo  ! loop over source of stype
              endif !     end if totlu ne 0
           enddo !     end loop styp
        endif !     end if for gntty=0
           if (verbose.ge.1) print*,'Direct irradiance from sources (W/m**2/nm) at (',oi,',',oj,'):',irradi(oi,oj)
     enddo !     end of loop over every observer
  enddo
  !     save the irradiance file
  open(unit=2,form='unformatted',file=outfile,action='write')
     write(2) nbx,nby
     do j=nby,1,-1
        do i=1,nbx
           write(2) irradi(i,j)
        enddo
     enddo
  close(unit=1)
  deallocate ( pval )
  deallocate ( pvalno )
  deallocate ( lamplu )
  deallocate ( obsD )
  deallocate ( val2d )
  deallocate ( altsol )
  deallocate ( altsob )
  deallocate ( reflec )
  deallocate ( lampal )
  deallocate ( inclix )
  deallocate ( incliy )
  deallocate ( obsH )
  deallocate ( ofill )
  deallocate ( irradi )
  deallocate ( gndty )
  deallocate ( azims )
  deallocate ( imin )
  deallocate ( imax )
  deallocate ( jmin )
  deallocate ( jmax )
  deallocate ( totlu )
  stop
end program illumhealth
