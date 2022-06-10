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
!     **    - Buildings considered (based on distances to the nearest street and averaged building heights per zones)     **
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
  real, dimension(:,:,:), allocatable :: pvalno     ! pvalno(181,181,ntype) values of the angular photometry functions 1sr(x)_axis=zenith, 2nd(y)_axis=azimuth
  real, dimension(:,:), allocatable :: lamplu       ! lamplu(width,width)  Source fluxes
  real, dimension(:,:), allocatable ::  val2d       ! val2d(width,width)  Temporary input array 2d
  real, dimension(:,:), allocatable ::  altsol      ! altsol(width,width)  Ground elevation (meter)
  real, dimension(:,:), allocatable ::  altsob      ! altsob(width,width)  Elevation including buildings
  real, dimension(:,:), allocatable ::  reflec      ! reflec(width,width)  reflectance  
  real, dimension(:,:), allocatable ::  lampal      ! lampal(width,width) Height of the light sources relative to the ground (meter)
  real, dimension(:,:), allocatable ::  inclix      ! inclix(width,width)  tilt of the ground pixel along x (radian)
  real, dimension(:,:), allocatable ::  incliy      ! incliy(width,width)  tilt of the ground pixel along y (radian)  
  real, dimension(:,:), allocatable ::  obsH        ! obsH(width,width)  averaged height of the buildings
  real, dimension(:,:), allocatable ::  ofill       ! ofill(width,width)  fill factor giving the probability to hit an obstacle when pointing in its direction real 0-1
  real, dimension(:,:), allocatable ::  irradi      ! irradi(width,width)  direct irradiance on a surface normal to the line of sight (no scattering)
  INTEGER, dimension(:,:), allocatable ::  gndty    ! gndty(width,width)  ground type flag  0=street 1=observer/building facade, 2=building top , 3=rear facace, 4=other
  INTEGER, dimension(:,:), allocatable ::  lmpty    ! lmpty(width,width)  5 lamp types according to the 5 spectral classes from Alejandro
  real, dimension(:,:), allocatable ::  azims       ! azims(width,width)  azimuth toward nearest point on a street
  integer imin,imax,jmin,jmax                       ! x and y limits containing lamps
  real totlu                                        ! total lamp power over the domain
  real angmin                                 ! minimum angle
  real reflsiz                                ! Size of the reflecting surface
  real largx                                  ! Width (x axis) of the modeling domain (meter)
  real largy                                  ! Length (y axis) of the modeling domain (meter)
  INTEGER nbx,nby                             ! Number of pixels in the modeling domain

  real srefl,brefl,orefl                      ! Street reflectance , Building facades reflectance , Other reflectance
  INTEGER stype                               ! Source type or zone index
  character*72 pafile,lufile,alfile,ohfile,odfile,offile,azfile,lifile,gifile 
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
  real angzen,ouvang        ! Zenithal angle between points (radians) and opening angle of the solid angle in degrees.
  INTEGER anglez            ! Emitting zenithal angle from the luminaire.
  INTEGER anglea            ! Emitting azimuth angle from the luminaire
  real P_dir,P_indir        ! photometri! function of the light sources (direct,reflected)
  real*8 xc,yc,zc,xn,yn,zn  ! Position (meter) of the elements (starting point, final point) for the calculation of the solid angle.
  real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z ! Components of the vectors used in the solid angle calculation routine.
  real omega                ! Solid angles
  real haut                 ! Haut (negative indicate that the surface is lighted from inside the ground. I.e. not considered in the calculation
  real epsilx,epsily        ! tilt of the ground pixel
  real flrefl               ! flux reaching a reflecting surface (watts).
  real irefl,irefl1         ! intensity leaving a reflecting surface toward the line of sight voxel.
  real nbang                ! for the averaging of the photometri! function

  INTEGER nzp,nz,nap,na
  character*3 lampno        ! lamp number string

  real angazi               ! azimuth angle between two points in rad, max dist for the horizon determination
  real latitu               ! approximate latitude of the domain center
  INTEGER xsrmi,xsrma,ysrmi,ysrma ! limits of the loop valeur for the reflecting surfaces

  real ff,ff2,hh            ! temporary obstacle filling factor and horizon blocking factor
  real distd                ! distance to compute the scattering probability
  real angazor              ! viewing azimuth angles in radian
  real angzeor              ! horizontal zenith angle (normal to a window)
  real angazsr              ! azimuth angle between a source and the nearest road
  real angzeos              ! zenith angle between an observer and a source
  real angazos              ! azimith angle between an observer and a source
  real angazssr             ! azimuth angle between source and surface
  real angazosr             ! azimuth angle between observer and surface
  real angzessr             ! zenith angle between source and surface
  real angzeosr             ! zenith angle between observer and surface
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
  real ddir_obs             ! distance between the source and the observer
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
  ntype=5
  angazor=90.                ! the window is pointing horizontally
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
  allocate ( lamplu(width,width) )  
  allocate ( val2d(width,width) )
  allocate ( altsol(width,width) )
  allocate ( altsob(width,width) )
  allocate ( reflec(width,width) )
  allocate ( lampal(width,width) )
  allocate ( inclix(width,width) )
  allocate ( incliy(width,width) )
  allocate ( obsH(width,width) )
  allocate ( irradi(width,width) )
  allocate ( gndty(width,width) )
  allocate ( azims(width,width) )
  allocate ( lmpty(width,width) )

  !     windows pointing toward horizon
  angzeor = pi/2.
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
  !     Initialisation of arrays and variables
  if (verbose.ge.1) print*,'Initializing variables...'
  do i=1,width
     do j=1,width
        val2d(i,j)=0.
        altsol(i,j)=0.
        altsob(i,j)=0.
        obsH(i,j)=0.
        gndty(i,j)=0
        lmpty(i,j)=0
        inclix(i,j)=0.
        incliy(i,j)=0.
        azims(i,j)=0.
        lampal(i,j)=0.
        lamplu(i,j)=0.
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
  alfile=basenm(1:lenbase)//'_altlp.bin' ! setting the file name of height of the sources lumineuse.
  azfile=basenm(1:lenbase)//'_azimu.bin'
  lifile=basenm(1:lenbase)//'lmpty.bin' 
  gifile=basenm(1:lenbase)//'gndty.bin'
  lufile=basenm(1:lenbase)//'_lumlp_.bin' ! setting the file name of the luminosite of the cases.
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
  !     reading buildings average height
  call twodin(nbx,nby,ohfile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        obsH(i,j)=val2d(i,j) ! filling of the array
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     reading ground type flag
  call twodin(nbx,nby,gifile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        gndty(i,j)=nint(val2d(i,j)) ! ground type flag flag array 0 or 1
     enddo                  ! end of the loop over all cells along y.
  enddo
  !     reading lamp type flag
  call twodin(nbx,nby,lifile,val2d,width)
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        lmpty(i,j)=nint(val2d(i,j)) ! ground type flag flag array 0 or 1
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
  !     reading luminosity files
  imin=nbx
  jmin=nby
  imax=1
  jmax=1
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
           if (i.lt.imin) imin=i
           if (i.gt.imax) imax=i  
           if (j.lt.jmin) jmin=j
           if (j.gt.jmax) jmax=i           
        endif
     enddo
  enddo
  if (jmin.lt.1) jmin=1
  if (imin.lt.1) imin=1
  if (jmax.gt.nbx) jmax=nby
  if (imax.gt.nbx) imax=nbx
  do i=1,nbx             ! beginning of the loop over all cells along x.
     do j=1,nby          ! beginning of the loop over all cells along y.
        lamplu(i,j)=val2d(i,j) ! filling lamp power array
        totlu=totlu+lamplu(i,j) ! the total lamp flux should be non-null to proceed to the calculations
     enddo               ! end of the loop over all cells along y.
  enddo                  ! end of the loop over all cells along x.  
  
  
  
  
  ! loading photometry files
  do stype=1,ntype          ! beginning of the loop over ntype types of sources.
     pvalto=0.
     write(lampno, '(I1.1)' ) stype ! support of ntype different sources (1 digit)
     pafile=basenm(1:lenbase)//'_fctem_'//lampno//'.dat' ! setting the file name of angular photometry.

     !     reading photometry files
     open(UNIT=1, FILE=pafile,status='OLD') ! opening file pa#.dat, angular photometry.
     ! each line (j) is a new azimuth starting at 0 deg (perpendicular to the street) and ending at 180 deg (behind)
     ! on each line we have the zenith angle (i) beginning at 0 (zenith) end ending at 180 (nadir)
     do i=1,181
        do j=1,181          ! beginning of the loop for the 181 data points
           read(1,*) pval(i,j,stype) ! reading of the data in the array pval.
           pvalto=pvalto+pval(i,j,stype)*2.*sin(real(i-1)*dtheta)*dtheta**2.  ! units of 1/sr.
        enddo
     enddo                  ! end of the loop over the 181 donnees of the fichier pa#.dat.
     close(1)               ! closing file pa#.dat, angular photometry.
     do i=1,181
        do j=1,181
           if (pvalto.ne.0.) pvalno(i,j,stype)=pval(i,j,stype)/pvalto ! Normalisation of the photometri! function.
        enddo
     enddo
  enddo  ! enddo stype
  deallocate ( pval )
  deallocate ( val2d )
  !     distribute reflectance values and adding fake buildings
  do i=1,nbx                ! beginning of the loop over all cells along x.
     do j=1,nby             ! beginning of the loop over all cells along y.
        if (gndty(i,j).eq.0) then
           reflec(i,j)=srefl
        elseif (gndty(i,j).eq.2) then
           reflec(i,j)=srefl
           altsob(i,j)=altsol(i,j)+obsH(i,j)
        elseif ((gndty(i,j).eq.1).or.(gndty(i,j).eq.3)) then
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
  !     0=street 
  !     1=building front and observer               222
  !     2=building top                              222
  !     3=building rear                            12223
  !     4=other                                    12223
  !     fake building profile                  000012223444444444
  do i=1,nbx                ! beginning of the loop over the column (longitude) of the domain.
     do j=1,nby             ! beginning of the loop over the rows (latitu) of the domain.
        if (gndty(i,j).eq.0) then  ! on street keep the inclinaison without the buildings    
           if (i.eq.1) then ! specific case close to the border of the domain (vertical side left).
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
        else  ! all other case we include the buildings to calculate the inclinaison in particular to reproduce the facades surfaces for the reflexion
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
  deallocate ( altsob )
  !     determination of the vertical atmospheric transmittance
  call transtoa(lambda,bandw,taua,layaod,pressi,tranam,tranaa,tranal,tabs)
  largx=dx*real(nbx)        ! computation of the Width along x of the case.
  largy=dy*real(nby)        ! computation of the Width along y of the case.
  !     remove 2nd aerosol layer
  layaod=0.
  hlay=2000.
  !     loop over potential observers
  do oi=1,nbx
     do oj=1,nby
        !     only calculate if it is an observer position
        if (gndty(oi,oj).eq.1) then
           angazor = (pi*azims(oi,oj))/180.         ! angle to the road  from an observer (gndty=1)
           ix = sin(angzeor)*cos(angazor)           ! viewing vector components
           iy = sin(angzeor)*sin(angazor)           ! looking horizontally toward street from the observer position
           iz = cos(angzeor)
           x_obs=oi
           y_obs=oj
           z_obs=z_o+altsol(x_obs,y_obs)
           rx_obs=real(x_obs)*dx
           ry_obs=real(y_obs)*dy
           if (z_obs.eq.0.) z_obs=0.001

           irradi(oi,oj)=0. ! initialize the total direct irradiance from sources to observer
           !     ==============================
           !     Calculation of the direct radiances
           !
           if (verbose.ge.1) print*,' Calculating obtrusive light...'
           !     loop over source types
              !     Is any flux in that lamp type?
              if (totlu.ne.0.) then
                 do x_s=imin,imax ! beginning of the loop over the column (longitude the) of the domain.
                    do y_s=jmin,jmax ! beginning of the loop over the rows (latitud) of the domain.
                       rx_s=real(x_s)*dx
                       ry_s=real(y_s)*dy
                       if (lamplu(x_s,y_s).ne.0.) then ! if the luminosite of the case is null, the program ignore this case.
                          z_s=(altsol(x_s,y_s)+lampal(x_s,y_s)) ! Definition of the position (metre) vertical of the source.
                          !
                          !     *********************************************************************************************************
                          !     calculation of the direct radiance of sources falling on a surface perpendicular
                          !     to the viewing angle Units of W/nm/m2/sr
                          !     *********************************************************************************************************
                          dho=sqrt((rx_obs-rx_s)**2.+(ry_obs-ry_s)**2.)
                          if ((dho.gt.0.).and.(z_s.ne.z_obs)) then  ! if same horizontal position we cannot have same z elevation
                             ! azimuth source to road
                             angazsr=azims(x_s,y_s)
                             call angleazimutal(rx_s,ry_obs,rx_s,ry_s,angazos)                             
                             ! azimuth and zenith from obsever to source
                             call anglezenithal(rx_obs,ry_obs,z_obs,rx_s,ry_s,z_s,angzeos)
                             call angleazimutal(rx_obs,ry_obs,rx_s,ry_s,angazos)
                             if (angzeos.gt.pi/4.) then ! 45deg. it is unlikely to have a 1km high mountain less than 1
                                call horizon(width,x_obs,y_obs,z_obs,dx,dy,altsol,angazos,zhoriz,dh)
                                if (dh.le.dho) then
                                   if (angzeos-zhoriz.lt.0.00001) then ! shadow the path line of sight-source is not below the horizon => we compute
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
                             !     projection angle of line to the lamp and the viewing angle
                             call angle3points (rx_s,ry_s,z_s,rx_obs,ry_obs,z_obs,ix,iy,iz,dang)
                             dang=pi-dang
                             anglez=nint(180.*(pi-angzeos)/pi)+1           ! index of zenith angle from source to observer
                             anglea=nint(180.*(pi-angazos+angazsr)/pi)+1   ! relative azimuth between the obs-source and the direction of the nearest street from the source
                             if (anglez.lt.0) anglez=-anglez
                             if (anglez.gt.181) anglez=362-anglez ! symetric function
                             if (anglez.eq.0) anglez=1
                             if (anglea.lt.0) anglea=-anglea
                             if (anglea.gt.181) anglea=362-anglea ! symetric function
                             if (anglea.eq.0) anglea=1
                             P_dir=pvalno(anglez,anglea,lmpty(x_s,y_s))
                             if ((cos(dang).gt.0.).and.(dang.lt.pi/2.)) then
                                ddir_obs=sqrt((rx_obs-rx_s)**2.+(ry_obs-ry_s)**2.+(z_obs-z_s)**2.)
                                !     computation of the solid angle 1m^2 at the observer as seen from the source
                                omega=1.*abs(cos(dang))/ddir_obs**2.
                                call transmitm(angzeos,z_obs,z_s,ddir_obs,transm,tranam,tabs)
                                call transmita(angzeos,z_obs,z_s,ddir_obs,haer,transa,tranaa)
                                call transmitl(angzeos,z_obs,z_s,ddir_obs,hlay,transl,tranal)
                                irradi(oi,oj)=irradi(oi,oj)+lamplu(x_s,y_s)*transa*transm*transl*P_dir*omega &
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
                                            call anglezenithal(rx_s,ry_s,z_s,rx_sr,ry_sr,z_sr,angzessr)
                                            call angleazimutal(rx_s,ry_s,rx_sr,ry_sr,angazssr)
                                            !     computation of the transmittance between the source and the ground surface
                                            distd=sqrt((rx_s-rx_sr)**2.+(ry_s-ry_sr)**2.+(z_s-z_sr)**2.)
                                            call transmitm(angzessr,z_s,z_sr,distd,transm,tranam,tabs)
                                            call transmita(angzessr,z_s,z_sr,distd,haer,transa,tranaa)
                                            call transmitl(angzessr,z_s,z_sr,distd,hlay,transl,tranal)
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
                                            ! estimation of the half of the underlying angle of the solid angle
                                            ! this angle servira a obtenir un meilleur isime (moyenne) of
                                            ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                                            ouvang=sqrt(omega/pi) ! Angle in radian.
                                            ouvang=ouvang*180./pi ! Angle in degrees.
                                            !     computation of the photometric function of the light fixture toward the reflection surface
                                            anglez=nint(180.*angzessr/pi)
                                            anglea=nint(180.*(pi-(2.*pi-angazssr)+angazsr)/pi)+1
                                            if (anglez.lt.0) anglez=-anglez
                                            if (anglez.gt.180) anglez=360-anglez
                                            anglez=anglez+1 ! Transform the angle in integer degree into the position in the array.
                                            if (anglea.lt.0) anglea=-anglea
                                            if (anglea.gt.180) anglea=360-anglea
                                            anglea=anglea+1 ! Transform the angle in integer degree into the position in the array.
                                            !     average +- ouvang
                                            nzp=0
                                            nap=0
                                            nbang=0.
                                            P_indir=0.
                                            do nz=-nint(ouvang),nint(ouvang)
                                               nzp=anglez+nz
                                               do na=-nint(ouvang),nint(ouvang)
                                                  nap=anglea+na
                                                  if (nzp.lt.0) nzp=-nzp
                                                  if (nzp.gt.181) nzp=362-nzp ! symetri! function
                                                  if (nzp.eq.0) nzp=1
                                                  if (nap.lt.0) nap=-nap
                                                  if (nap.gt.181) nap=362-nap ! symetri! function
                                                  if (nap.eq.0) nap=1
                                                  P_indir=P_indir+pvalno(nzp,nap,lmpty(x_s,y_s))*abs(sin(pi*real(nzp)/180.))/2.
                                                  nbang=nbang+1.*abs(sin(pi*real(nzp)/180.))/2.
                                               enddo
                                            enddo
                                            P_indir=P_indir/nbang
                                            !     computation of the flux reaching the reflecting surface
                                            flrefl=lamplu(x_s,y_s)*P_indir*omega*transm*transa*transl
                                            !     computation of the reflected intensity leaving the ground surface
                                            irefl1=flrefl*reflec(i,j)/pi ! The factor 1/pi comes from the normalisation of the fonction
                                            dho=sqrt((rx_obs-rx_sr)**2.+(ry_obs-ry_sr)**2.)
                                            if ((dho.gt.0.).and.(z_s.ne.z_obs)) then
                                               call anglezenithal(rx_obs,ry_obs,z_obs,rx_sr,ry_sr,z_sr,angzeosr)
                                               call angleazimutal(rx_obs,ry_obs,rx_sr,ry_sr,angazosr)
                                               if (angzen.gt.pi/4.) then ! 45deg. it is unlikely to have a 1km high mountain less than 1
                                                  call horizon(width,x_obs,y_obs,z_obs,dx,dy,altsol,angazosr,zhoriz,dh)
                                                  if (dh.le.dho) then
                                                     if (angzeosr-zhoriz.lt.0.00001) then ! shadow the path line of sight-source is not below the horizon => we compute
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
                                               !     projection angle of line to the lamp and the viewing angle
                                               call angle3points (rx_sr,ry_sr,z_sr,rx_obs,ry_obs,z_obs,ix,iy,iz,dang)
                                               dang=pi-dang
                                               !     computation of the flux direct reaching the line of sight voxel
                                               if ((cos(dang).gt.0.).and.(dang.lt.pi/2.)) then
                                                  ddir_obs=sqrt((rx_obs-rx_sr)**2.+(ry_obs-ry_sr)**2.+(z_obs-z_sr)**2.)
                                                  !     computation of the solid angle of the line of sight voxel seen from the source
                                                  omega=1.*abs(cos(dang))/ddir_obs**2.
                                                  call transmitm(angzeosr,z_obs,z_sr,ddir_obs,transm,tranam,tabs)
                                                  call transmita(angzeosr,z_obs,z_sr,ddir_obs,haer,transa,tranaa)
                                                  call transmitl(angzeosr,z_obs,z_sr,ddir_obs,hlay,transl,tranal)
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
                    enddo ! loop over sources
                 enddo  ! loop over sources
              endif !     end if totlu ne 0
        endif !     end if for gntty=0
           if (verbose.ge.1) print*,'Direct irradiance from sources (W/m**2/nm) at (',oi,',',oj,'):',irradi(oi,oj)
     enddo !     end of loop over every observer
  enddo
  deallocate ( pvalno )  
  deallocate ( lamplu )
  deallocate ( altsol )
  deallocate ( reflec )
  deallocate ( lampal )
  deallocate ( inclix )
  deallocate ( incliy )
  deallocate ( obsH )
  deallocate ( gndty )
  deallocate ( lmpty )
  deallocate ( azims )

  !     save the irradiance file
  open(unit=2,form='unformatted',file=outfile,action='write')
     write(2) nbx,nby
     do j=nby,1,-1
        do i=1,nbx
           write(2) irradi(i,j)
        enddo
     enddo
  close(unit=1)
  deallocate ( irradi )
  stop
end program illumhealth
