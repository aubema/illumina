!                         *                          *                       iii                   *     *
!                                                                           iiiii
!  IIIIII    lLLLL    *    lLLLL         UUU    UUU      MMMMM      MMMMM    iii        NNNN     NN          AAAA
!   IIII     LLLL          LLLL   *     UUU      UUU     MMMMMMM  MMMMMMM          *    NNNNN    NN        AAAaaAAA
!   IIII     LLLL          LLLL        UUU *      UUU    MMM MMMMMMMM MMM    iii        NNNNNN   NN       AAA    AAA
!   IIII     LLLL   *      LLLL        UUU        UUU    MMM *        MMM  iii          NNN  NNN NN     AAAAAAAAAAAAAA
!   IIII     LLLl          LLLl        UUUu      uUUU    MMM          MMM  iiii    ii   NNN   NNNNN    AAAa        aAAA
!   IIII    LLLLLLLLLL    LLLLLLLLLL    UUUUUuuUUUUU     MMM          MMM   iiiiiiiii   NNN    NNNN   aAAA    *     AAAa
!  IIIIII   LLLLLLLLLLL   LLLLLLLLLLL     UUUUUUUU      mMMMm        mMMMm   iiiiiii   nNNNn    NNNn  aAAA          AAAa
!
! **********************************************************************************************************************
! ** Illumina VERSION 2 - in Fortran 77                                                                               **
! ** Programmers in decreasing order of contribution  :                                                               **
! **                            Martin Aube                                                                           **
! **              Still having very few traces of their contributions :                                               **
! **                            Loic Franchomme-Fosse,  Mathieu Provencher, Andre Morin                               **
! **                            Alex Neron, Etienne Rousseau                                                          **
! **                            William Desroches, Maxime Girardin, Tom Neron                                         **
! **                                                                                                                  **
! ** Illumina can be downloaded via:   git clone https://github.com/aubema/illumina.git                               **
! ** To compile:                                                                                                      **
! **    cd hg/illumina/bin                                                                                            **
! **    bash makeILLUMINA                                                                                             **
! **                                                                                                                  **
! **  Current version features/limitations :                                                                          **
! **                                                                                                                  **
! **    - Calculation of artificial sky radiance in a given line of sight                                             **
! **    - Calculation of the atmospheric transmittance and 1st and 2nd 3rd orders of scattering                       **
! **    - Lambertian reflexion on the ground                                                                          **
! **    - Terrain slope considered (apparent surface and shadows)                                                     **
! **    - Angular photometry of a lamp is considered uniform along the azimuth                                        **
! **    - Sub-grid obstacles considered (with the mean free path of light toward ground, mean obstacle height, and    **
! **      obstacles transparency (filling factor)                                                                     **
! **    - Molecules and aerosol optics (phase function, scattering probability, aerosol absorption)                   **
! **    - Exponential concentrations vertical profile                                                                 **
! **    - Accounting for heterogeneity luminaires number, luminaires heights, luminaire spectrum,                     **
! **      angular photometry, obstacle properties                                                                     **
! **    - Wavelength dependant                                                                                        **
! **    - Cloud models (type and cloud base height) only the overhead clouds are considered with cloud fraction       **
! **    - Support direct observation of a source                                                                      **
! **    - Direct observation of the ground is implemented                                                             **
! **                                                                                                                  **
! **********************************************************************************************************************
!
!  Copyright (C) 2024 Martin Aube PhD
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Contact: martin.aube@cegepsherbrooke.qc.ca
!
!
!
      program illumina ! Beginning
      implicit none
      ! Variables declaration
      integer width,nzon                                           ! Arrays dimension in Length/width and zones
      parameter (width=512,nzon=256) 
      real*8 pi,pix4
      real*8 zero,un                                                 ! value of 0. and 1.
      integer verbose                                              ! verbose = 1 to have more print out, 0 for silent
      parameter (pi=3.141592654)
      parameter (pix4=4.*pi)
      character*72 mnaf                                            ! Terrain elevation file
      character*72 diffil                                          ! Aerosol file
      character*72 outfile                                         ! Results file
      character*72 pclimg
      character*72 basenm                                          ! Base name of files
      integer lenbase                                              ! Length of the Base name of the experiment
      real*8 lambda,pressi,drefle(width,width)                       ! Wavelength (nanometer), atmospheric pressure (kPa), mean free path to the ground (meter).
      real*8 reflsiz                                                 ! Size of the reflecting surface
      integer ntype                                                ! Number of light source types or zones considered
      real*8 largx                                                   ! Width (x axis) of the modeling domain (meter)
      real*8 largy                                                   ! Length (y axis) of the modeling domain (meter)
      integer nbx,nby                                              ! Number of pixels in the modeling domain
      real*8 val2d(width,width)                                      ! Temporary input array 2d
      real*8 altsol(width,width)                                     ! Ground elevation (meter)
      real*8 srefl                                                   ! Ground reflectance
      integer stype                                                ! Source type or zone index
      character*72 pafile,lufile,alfile,ohfile,odfile,offile       ! Files related to light sources and obstacles (photometric function of the sources (sr-1), flux (W), height (m), obstacles c                                                                  ! height (m), obstacle distance (m), obstacle filling factor (0-1).
      real*8 lamplu(width,width,nzon)                                ! Source fluxes
      real*8 lampal(width,width)                                     ! Height of the light sources relative to the ground (meter)
      real*8 pval(181,nzon),pvalto,pvalno(181,nzon)                  ! Values of the angular photometry functions (unnormalized, integral, normalized)
      real*8 dtheta                                                  ! Angle increment of the photometric function of the sources
      real*8 dx,dy,dxp,dyp                                           ! Width of the voxel (meter)
      integer boxx,boxy                                            ! reflection window size (pixels)
      real*8 fdifa(181),fdifan(181)                                  ! Aerosol scattering functions (unnormalized and normalized)
      real*8 extinc,scatte,anglea(181)                               ! Aerosol cross sections (extinction and scattering), scattering angle (degree)
      real*8 secdif                                                  ! Contribution of the scattering to the extinction
      real*8 inclix(width,width)                                     ! tilt of the ground pixel along x (radian)
      real*8 incliy(width,width)                                     ! tilt of the ground pixel along y (radian)
      integer x_obs,y_obs                                          ! Position of the observer (INTEGER)
      real*8 rx_obs,ry_obs 
      real*8 z_o                                                     ! observer height relative to the ground (meter)
      real*8 z_obs                                                   ! Height of the observer (meter) to the vertical grid scale
      integer ncible,icible                                        ! Number of line of sight voxels, number loops over the voxels
      integer x_c,y_c                                              ! Position of the line of sight voxel (INTEGER)
      real*8 rx_c,ry_c 
      real*8 z_c                                                     ! Height of the line of sight voxel (meter)
      integer x_s,y_s,x_sr,y_sr                                    ! Positions of the source, the reflecting surface, and the scattering voxels
      real*8 z_s,z_sr,z_dif2,z_dif3                                  ! Heights of the source, the reflecting surface, and the scattering voxel (metre).
      real*8 rx_s,ry_s,rx_sr,ry_sr,rx_dif2,ry_dif2,rx_dif3,ry_dif3
      real*8 angzen,ouvang                                           ! Zenithal angle between two voxels (radians) and opening angle of the solid angle in degrees.
      integer anglez                                               ! Emitting zenithal angle from the luminaire.
      real*8 zenith                                                   ! zenith angle used for blocking
      real*8 P_indir,P_dif1                                    ! photometric function of the light sources (direct,indirect,scattered)
      real*8 P_dir
      real*8 transa,transm,transl                                    ! Transmittance between two voxels (aerosols,molecules,particle layer).
      real*8 taua                                                    ! Aerosol optical depth @ 500nm.
      real*8 alpha                                                   ! Angstrom coefficient of aerosol AOD
      real*8 xc,yc,zc,xn,yn,zn                                     ! Position (meter) of the elements (starting point, final point) for the calculation of the solid angle.
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z       ! Components of the vectors used in the solid angle calculation routine.
      real*8 omega                                                   ! Solid angles
      real*8 idif3                                                   ! 3rd scat  intensity
      real*8 itodif                                                  ! Total contribution of the source to the scattered intensity toward the sensor.
      real*8 fcapt                                                   ! Flux reaching the observer voxel from all FOV voxels in a given model level
      real*8 ftocap                                                  ! Total flux reaching the observer voxel
      real*8 haut                                                    ! Haut (negative indicate that the surface is lighted from inside the ground. I.e. not considered in the calculation
      real*8 epsilx,epsily                                           ! tilt of the ground pixel
      real*8 flrefl                                                  ! flux reaching a reflecting surface (watts).
      real*8 irefl,irefl1                                            ! intensity leaving a reflecting surface toward the line of sight voxel.
      real*8 effdif2,effdif3                                         ! Distance around the source voxel and line of sight voxel considered to compute the 2nd order of scattering.
      real*8 zondi2(3000000,3),zondi3(3000000,3)                   ! Array for the scattering voxels positions
      integer ndiff2                                               ! Number of 2nd scattering voxels, counter of the loop over the 2nd scattering voxels
      integer ndiff3                                               ! Number of 3rd scattering voxels, counter of the loop over the 3rd scattering voxels
      integer ssswit                                               ! activate single double and triple scattering (0= direct only, 1= 1st scat max, 2= 2nd scat max, 3= 3rd scat max)
      real*8 idif1                                                   ! intensity toward a line of sight voxel from a scattering voxel (without and with reflexion).
      real*8 idif2
      real*8 portio                                                  ! ratio of voxel surface to the solid angle of the sensor field of view.
      real*8 dis_obs                                                 ! Distance between the line of sight and the observer.
      real*8 ometif                                                  ! Solid angle of the telescope objective as seen from the line of sight voxel
      real*8 omefov                                                  ! Solid angle of the spectrometer slit.
      real*8 angvis,azim                                             ! viewing angles of the sensor.
                                                                   ! Useful for the calculation of the lambertian reflectance.
      real*8 nbang                                                   ! for the averaging of the photometric function
      real*8 obsH(width,width)                                       ! averaged height of the sub-grid obstacles, minimum angle under wich
                                                                   ! a light ray cannot propagate because it is blocked by a sub-grid obstable
      real*8 ofill(width,width)                                      ! fill factor giving the probability to hit an obstacle when pointing in its direction real 0-1
      integer naz,na
      real*8 ITT(width,width,nzon)                                   ! total intensity per type of lamp
      real*8 ITC(width,width)                                        ! total intensity per line of sight voxel
      real*8 FTC(width,width)                                        ! fraction of the total flux at the sensor level
      real*8 FCA(width,width)                                        ! sensor flux array
      real*8 lpluto(width,width)                                     ! total luminosity of the ground cell for all lamps
      character*3 lampno                                           ! lamp number string
      integer imin(nzon),imax(nzon),jmin(nzon),jmax(nzon)          ! x and y limits containing a type of lamp
      real*8 latitu                                                  ! approximate latitude of the domain center
      integer prmaps                                               ! flag to enable the tracking of contribution and sensitivity maps
      integer cloudt                                               ! cloud type 0=clear, 1=Thin Cirrus/Cirrostratus, 2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus,
                                                                   ! 4=Stratocumulus/stratus, 5=Cumulus/Cumulonimbus
      real*8 cloudslope                                              ! slope of the radiance dependency on the cloud fraction (in percentage) According to
                                                                   ! Sciezoras 2020 the slope vary depending on the level of LP and how it is distributed.
                                                                   ! We decided instead to simplify this by using an average slope of -0.013.
                                                                   ! Rad=Rad_100 * 10**(0.4*(100-cloudfrac)*cloudslope) this equation is derived from
                                                                   ! Tomasz Sciezor , The impact of clouds on the brightness of the night sky, Journal of
                                                                   ! Quantitative Spectroscopy & Radiative Transfer (2020),
                                                                   ! doi: https://doi.org/10.1016/j.jqsrt.2020.106962
      real*8 cloudfrac                                               ! cloud fraction in percentage
      integer xsrmi,xsrma,ysrmi,ysrma                              ! limits of the loop valeur for the reflecting surfaces
      real*8 icloud                                                  ! cloud reflected intensity
      real*8 fccld                                                   ! correction for the FOV to the flux reaching the intrument from the cloud voxel
      real*8 fctcld                                                  ! total flux from cloud at the sensor level
      real*8 totlu(nzon)                                             ! total flux of a source type
      real*8 stoplim                                                 ! Stop computation when the new voxel contribution is less than 1/stoplim of the cumulated flux
      real*8 ff1,ff2,hh                                              ! temporary obstacle filling factor and horizon blocking factor
      real*8 cloudbase,cloudtop                                      ! cloud base and top altitude (m)
      real*8 distd                                                   ! distance to compute the scattering probability
      real*8 volu2,volu3                                             ! volume of a voxel
      real*8 scal                                                    ! stepping along the line of sight
      real*8 scalo                                                   ! previous value of scal
      real*8 siz2,siz3                                               ! resolution of the 2nd scat grid in meter
      real*8 angvi1,angaz1,angze1                                    ! viewing angles in radian
      real*8 ix,iy,iz                                                ! base vector of the viewing (length=1)
      real*8 omemax                                                  ! max solid angle allowed
      real*8 flcld(width,width)                                      ! flux crossing a low cloud
      real*8 ds1,ds2,ds3,ds4,dss2,dss3                               ! double scattering distances
      integer ndi2,ndi3                                            ! number of cell under ground
      real*8 diamobj                                                 ! instrument objective diameter
      integer i,j,k
      integer dirck                                                       ! Test for the position of the source (case source=line of sight voxel)
      real*8 tranam,tranaa                                           ! atmospheric transmittancess of a path (molecular, aerosol)
      real*8 zhoriz                                                  ! zenith angle of the horizon
      real*8 direct                                                  ! direct radiance from sources on a surface normal to the line of sight (no scattering)
      real*8 rdirect                                                 ! direct radiance from a reflecting surface on a surface normal to the line of sight (no scattering)
      real*8 irdirect                                                ! direct irradiance from sources on a surface normal to the line of sight (no scattering)
      real*8 irrdirect                                               ! direct irradiance from a reflecting surface on a surface normal to the line of sight (no scattering)
      real*8 dang                                                    ! Angle between the line of sight and the direction of a source
      real*8 ddir_obs                                                ! distance between the source and the observer
      real*8 rx,ry,rz                                                ! driving vector for the calculation of the projection angle for direct radiance. It is 20km long
      real*8 dfov                                                    ! field of view in degrees for the calculation of the direct radiance this number will be a kind of smoothing
                                                                   ! effect. The angular grid resolution to create a direct radiance panorama should be finer than that number
      real*8 Fo                                                      ! flux correction factor for obstacles
      real*8 thetali                                                 ! limit angle for the obstacles blocking of viirs
      integer viirs(width,width)                                   ! viirs flag 1=yes 0=no
      character*72 vifile                                          ! name of the viirs flag file
      real*8 dh0,dhmax                                               ! horizontal distance along the line of sight and maximum distance before beeing blocked by topography
      character*72 layfile                                         ! filename of the optical properties of the particle layer
      real*8 layaod                                                  ! 500 nm aod of the particle layer
      real*8 layalp                                                  ! spectral exponent of the aod for the particle layer
      real*8 hlay                                                    ! exponential vertical scale height of the particle layer
      real*8 secdil                                                  ! scattering/extinction ratio for the particle layer
      real*8 fdifl(181)                                              ! scattering phase function of the particle layer
      real*8 tranal                                                  ! top of atmos transmission of the particle layer
      real*8 haer                                                    ! exponential vertical scale height of the background aerosol layer
      real*8 bandw                                                   ! bandwidth of the spectral bin
      real*8 tabs                                                    ! TOA transmittance related to molecule absorption

      real*8 itoclou                                                 ! cloud intensity
      real*8 itodif1                                                 ! First scattering intensity
      real*8 itodif2                                                 ! Second scattering intensity
      real*8 itodif3                                                 ! total 3rd order intentity
      real*8 itot
      real*8 tot1,tot2,tot3
      integer rho                                                  ! switch between source (rho=0) or ground pixel (rho=1)
      verbose=1                                                    ! Very little printout=0, Many printout = 1, even more=2
      diamobj=1.D0                                                   ! A dummy value for the diameter of the objective of the instrument used by the observer.
      volu2=0.D0
      volu3=0.D0
      zero=0.D0
      un=1.D0
      ff1=0.D0
      ff2=0.D0
      ncible=1024
      cloudslope=-0.013D0
      cloudfrac=100.D0
      effdif2=5000.D0
      effdif3=10000.D0
      siz2=1000.D0
      !siz3=20000.D0
      siz3=(4.*pi*effdif3**3./3/sqrt(1000000.))**(1./3.)
      volu2=siz2**3.D0
      volu3=siz3**3.D0
      if (verbose.ge.1) then
        print*,'Starting ILLUMINA computations...'
      endif
      ! reading of input file (illumina.in)
      print*,'Reading illumina.in input file'
      open(unit=1,file='illumina_new.in',status='old')
        read(1,*)
        read(1,*) basenm
        read(1,*) dx,dy
        read(1,*) diffil
        read(1,*) layfile, layaod, layalp, hlay
        read(1,*) ssswit
        read(1,*) 
        read(1,*) lambda,bandw
        read(1,*) srefl
        read(1,*) pressi
        read(1,*) taua,alpha,haer
        read(1,*) ntype
        read(1,*) stoplim
        read(1,*)
        read(1,*) x_obs,y_obs,z_o
        read(1,*) 
        read(1,*) angvis,azim
        read(1,*) dfov
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) reflsiz
        read(1,*) cloudt, cloudbase, cloudfrac
        read(1,*)
      close(1)
      if (angvis.gt.90.) then
        print*,'Error: elevation angle larger than 90 deg'
        stop
      endif
      if (angvis.lt.-90.) then
        print*,'Error: elevation angle smaller than -90 deg'
        stop
      endif
      ! conversion of the geographical viewing angles toward the cartesian
      ! angle we assume that the angle in the file illumina.in
      ! is consistent with the geographical definition
      ! geographical, azim=0 toward north, 90 toward east, 180 toward south
      ! etc
      ! cartesian, azim=0 toward east, 90 toward north, 180 toward west etc
      azim=90.-azim
      if (azim.lt.0.) azim=azim+360.
      if (azim.ge.360.) azim=azim-360.
      angvi1=(pi*angvis)/180.
      angze1=pi/2.-angvi1
      angaz1=(pi*azim)/180.
      ix=(sin((pi/2.)-angvi1))*(cos(angaz1))                        ! viewing vector components
      iy=(sin((pi/2.)-angvi1))*(sin(angaz1))
      iz=(sin(angvi1))
      dfov=(dfov*pi/180.)/2.

      scal=19.D0
      scalo=scal
      boxx=idnint(reflsiz/dx) ! Number of column to consider left/right of the source for the reflection.
      boxy=idnint(reflsiz/dy) ! Number of column to consider up/down of the source for the reflection.
      ! omemax: exclude calculations too close (<10m) this is a sustended angle of 1 deg.
      ! the calculated flux is highly sensitive to that number for a very high
      ! pixel resolution (a few 10th of meters). We assume anyway that somebody
      ! observing the sky will never lies closer than that distance to a
      ! light fixture. This number is however somehow subjective and that means
      ! that the value of sky brightness near sources will be affected by this
      ! choice
      omemax=1./((10.D0)**2.D0)
      if (verbose.gt.0) then
        print*,'2nd order scattering grid = ',siz2,'m'
        print*,'2nd order scattering radius=',effdif2,'m'
        print*,'3rd order scattering grid = ',siz3,'m'
        print*,'3rd order scattering radius=',effdif3,'m'        
        print*,'Pixel size = ',dx,' x ',dy
        print*,'Maximum radius for reflection = ',reflsiz
      endif
      ! computing the actual AOD at the wavelength lambda
      if (verbose.ge.1) then
        print*,'500nm AOD=',taua,'500nm angstrom coeff.=',alpha
      endif
      taua=taua*(lambda/500.D0)**(-1.*alpha)
      layaod=layaod*(lambda/500.D0)**(-1.*layalp)
      !  determine the Length of basenm
      lenbase=index(basenm,' ')-1
      mnaf=basenm(1:lenbase)//'_topogra.bin' ! determine the names of input and output files
      outfile=basenm(1:lenbase)//'.out'
      pclimg=basenm(1:lenbase)//'_pcl.bin'
      ! opening output file
      open(unit=2,file=outfile,status='unknown')
        write(2,*) "ILLUMINA version __version__"
        write(2,*) 'FILE USED:'
        write(2,*) mnaf,diffil
        print*,'Wavelength (nm):',lambda,' Aerosol optical depth:',taua
        write(2,*) 'Wavelength (nm):',lambda,' Aerosol optical depth:',taua
        write(2,*) '2nd order scattering radius:',effdif2,' m'
        print*,'2nd order scattering radius:',effdif2,' m'
         write(2,*) '3rd order scattering radius:',effdif3,' m'
        print*,'3rd order scattering radius:',effdif3,' m'       
        write(2,*) 'Observer position (x,y,z)',x_obs,y_obs,z_o
        print*,'Observer position (x,y,z)',x_obs,y_obs,z_o
        write(2,*) 'Elevation angle:',angvis,' azim angle (counterclockwise from east)',azim
        print*,'Elevation angle:',angvis,' azim angle (counterclockwise from east)',azim
        ! Initialisation of the arrays and variables
        if (verbose.ge.1) print*,'Initializing variables...'
        if (cloudt.eq.0) then
          cloudbase=1000000000.
        endif
        prmaps=1
        icloud=0.
        do i=1,width
          do j=1,width
            val2d(i,j)=0.
            altsol(i,j)=0.
            obsH(i,j)=0.
            drefle(i,j)=0.
            viirs(i,j)=0
            ofill(i,j)=0.
            inclix(i,j)=0.
            incliy(i,j)=0.
            lpluto(i,j)=0.
            ITC(i,j)=0.
            FTC(i,j)=0.
            FCA(i,j)=0.
            flcld(i,j)=0.
            lampal(i,j)=0.
            do k=1,nzon
              lamplu(i,j,k)=0.
              ITT(i,j,k)=0.
            enddo
          enddo
        enddo
        do i=1,181
          fdifa(i)=0.
          fdifan(i)=0.
          fdifl(i)=0.
          anglea(i)=0.
          do j=1,nzon
            pval(i,j)=0.
            pvalno(i,j)=0.
          enddo
        enddo
        do i=1,3000000
          do j=1,3
            zondi2(i,j)=1.
            zondi3(i,j)=1.
          enddo
        enddo
        do i=1,nzon
          totlu(i)=0.
        enddo
        idif1=0.
        idif2=0.
        idif3=0.
        flrefl=0.
        irefl=0.
        irefl1=0.
        portio=0.
        fccld=0.
        fctcld=0.
        ometif=0.
        omefov=0.
        hh=1.
        itodif1=0.
        itodif2=0.
        itodif3=0.
        ! determine the 2nd scattering zone
        if (ssswit.gt.1) then
          call zone_diffusion(effdif2,zondi2,ndiff2,siz2)
          dss2=1.D0*siz2
          if (verbose.gt.0) then
            print*,'2nd order scattering grid points =',ndiff2
            print*,'2nd order scattering smoothing radius =',dss2,'m'
          endif
        endif
 
        ! determine the 3rd scattering zones
        if (ssswit.gt.2) then
          call zone_diffusion(effdif3,zondi3,ndiff3,siz3)
          dss3=1.D0*siz3
          if (verbose.gt.0) then
            print*,'3rd order scattering grid points =',ndiff3*ndiff3
            print*,'3rd order scattering smoothing radius =',dss3,'m'
          endif
        endif   
        ! determination of the vertical atmospheric transmittance
        call transtoa(lambda,bandw,taua,layaod,pressi,tranam,tranaa,tranal,tabs) ! tranam and tranaa are the top of atmosphere transmittance (molecules and aerosols)
        ! reading of the environment variables
        ! reading of the elevation file
        call twodin(nbx,nby,mnaf,altsol)
        ! computation of the tilt of the pixels along x and along y
        do i=1,nbx ! beginning of the loop over the column (longitude) of the domain.
          do j=1,nby ! beginning of the loop over the rows (latitu) of the domain.
            if (i.eq.1) then ! specific case close to the border of the domain (vertical side left).
              inclix(i,j)=atan((altsol(i+1,j)-altsol(i,j))/real(dx)) ! computation of the tilt along x of the surface.
            elseif (i.eq.nbx) then ! specific case close to the border of the domain (vertical side right).
              inclix(i,j)=atan((altsol(i-1,j)-altsol(i,j))/(real(dx))) ! computation of the tilt along x of the surface.
            else
              inclix(i,j)=atan((altsol(i+1,j)-altsol(i-1,j))/(2.*real(dx))) ! computation of the tilt along x of the surface.
            endif
            if (j.eq.1) then ! specific case close to the border of the domain (horizontal side down).
              incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j))/(real(dy))) ! computation of the tilt along y of the surface.
            elseif (j.eq.nby) then ! specific case close to the border of the domain (horizontal side up).
              incliy(i,j)=atan((altsol(i,j-1)-altsol(i,j))/(real(dy))) ! computation of the tilt along y of the surface.
            else
              incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j-1))/(2.*real(dy))) ! computation of the tilt along y of the surface
            endif
          enddo ! end of the loop over the rows (latitu) of the domain
        enddo ! end of the loop over the column (longitude) of the domain
        ! reading of the values of P(theta), height, luminosities and positions
        ! of the sources, obstacle height and distance
        ohfile=basenm(1:lenbase)//'_obsth.bin'
        odfile=basenm(1:lenbase)//'_obstd.bin'
        alfile=basenm(1:lenbase)//'_altlp.bin' ! setting the file name of height of the sources lumineuse.
        offile=basenm(1:lenbase)//'_obstf.bin'
        vifile='origin.bin'
        dtheta=.017453293 ! one degree
        ! reading lamp heights
        call twodin(nbx,nby,alfile,val2d)
        do i=1,nbx ! beginning of the loop over all cells along x.
          do j=1,nby ! beginning of the loop over all cells along y.
            lampal(i,j)=val2d(i,j) ! filling of the array for the lamp stype
          enddo ! end of the loop over all cells along y.
        enddo ! end of the loop over all cells along x.
        ! reading subgrid obstacles average height
        call twodin(nbx,nby,ohfile,val2d)
        do i=1,nbx ! beginning of the loop over all cells along x.
          do j=1,nby ! beginning of the loop over all cells along y.
            obsH(i,j)=val2d(i,j) ! filling of the array
          enddo ! end of the loop over all cells along y.
        enddo
        ! reading subgrid obstacles average distance
        call twodin(nbx,nby,odfile,val2d)
        do i=1,nbx ! beginning of the loop over all cells along x.
          do j=1,nby ! beginning of the loop over all cells along y.
            drefle(i,j)=val2d(i,j)/2.
            if (drefle(i,j).eq.0.) drefle(i,j)=dx ! when outside a zone, block to the size of the cell (typically 1km)
          enddo ! end of the loop over all cells along y.
        enddo
        ! reading subgrid obstacles filling factor
        call twodin(nbx,nby,offile,val2d)
        do i=1,nbx ! beginning of the loop over all cells along x.
          do j=1,nby ! beginning of the loop over all cells along y.
            ofill(i,j)=val2d(i,j) ! Filling of the array 0-1
          enddo ! end of the loop over all cells along y.
        enddo
        ! reading viirs flag
        call twodin(nbx,nby,vifile,val2d)
        do i=1,nbx ! beginning of the loop over all cells along x.
          do j=1,nby ! beginning of the loop over all cells along y.
            viirs(i,j)=idnint(val2d(i,j)) ! viirs flag array 0 or 1
          enddo ! end of the loop over all cells along y.
        enddo
        ! reading of the scattering parameters for background aerosols
        open(unit = 1, file = diffil,status= 'old') ! opening file containing the scattering parameters
          read(1,*)  secdif ! the scattering / extinction ratio
          read(1,*)
          do i=1,181
            read(1,*) anglea(i), fdifa(i) ! reading of the scattering functions
            fdifan(i)=fdifa(i)/pix4 ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle
          enddo
        close(1)
        ! reading scattering parameters of particle layer
        open(unit = 1, file = layfile,status= 'old') ! opening file containing the scattering parameters
          read(1,*)  secdil ! the scattering / extinction ratio of particle layer
          read(1,*)
          do i=1,181
            read(1,*) anglea(i), fdifl(i) ! reading of the scattering functions of the particle layer
            fdifl(i)=fdifl(i)/pix4 ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle
          enddo
        close(1)
        ! Some preliminary tasks
        do stype=1,ntype ! beginning of the loop 1 for the nzon types of sources.
          pvalto=0.
          write(lampno, '(I3.3)' ) stype ! support of nzon different sources (3 digits)
          pafile=basenm(1:lenbase)//'_fctem_'//lampno//'.dat' ! setting the file name of angular photometry.
          lufile=basenm(1:lenbase)//'_lumlp_'//lampno//'.bin' ! setting the file name of the luminosite of the cases.
          ! reading photometry files
          open(UNIT=1, FILE=pafile,status='OLD') ! opening file pa#.dat, angular photometry.
            do i=1,181 ! beginning of the loop for the 181 data points
              read(1,*) pval(i,stype) ! reading of the data in the array pval.
              pvalto=pvalto+pval(i,stype)*2.*pi*sin(real(i-1)*dtheta)*dtheta ! Sum of the values of the  photometric function
              ! (pvaleur x 2pi x sin theta x dtheta) (ou theta egale (i-1) x 1 degrees).
            enddo ! end of the loop over the 181 donnees of the fichier pa#.dat.
          close(1) ! closing file pa#.dat, angular photometry.
          do i=1,181
            if (pvalto.ne.0.) pvalno(i,stype)=pval(i,stype)/pvalto ! Normalisation of the photometric function.
          enddo
          
          
          
          
          ! reading luminosity files
          call twodin(nbx,nby,lufile,val2d)
          imin(stype)=nbx
          imax(stype)=1
          jmin(stype)=nby
          jmax(stype)=1
 336      do i=1,nbx ! beginning of the loop over all cells along x.
            do j=1,nby ! beginning of the loop over all cells along y.
              lamplu(i,j,stype)=val2d(i,j) ! remplir the array of the lamp type: stype
              if (lamplu(i,j,stype).lt.0.) then ! searching of negative fluxes
                print*,'***Negative lamp flux!, stopping execution'
                stop
              endif              
              if (lamplu(i,j,stype).ne.0.) then ! searching of the smallest rectangle containing the zone of non-null luminosity to speedup the calculation
                if (i.lt.imin(stype)) imin(stype)=i
                if (j.lt.jmin(stype)) jmin(stype)=j 
                if (i.gt.imax(stype)) imax(stype)=i
                if (j.gt.jmax(stype)) jmax(stype)=j                 
              endif
              
           
              
              
              
              
              
              ! Atmospheric correction and obstacles masking corrections to the lamp
              ! flux arrays (lumlp)
              if (viirs(i,j).eq.1) then
                lamplu(i,j,stype)=lamplu(i,j,stype)/(tranam*tranaa*tranal)
                thetali=atan2(drefle(i,j),obsH(i,j))
                if (thetali .lt. 70.*pi/180.) then
                  Fo=(1.-cos(70.*pi/180.))/(1.-ofill(i,j)*cos(thetali)+(ofill(i,j)-1.)*cos(70.*pi/180.))
                  lamplu(i,j,stype)=lamplu(i,j,stype)*Fo
                else
                  Fo=1.
                endif
              endif
              totlu(stype)=totlu(stype)+lamplu(i,j,stype) ! the total lamp flux should be non-null to proceed to the calculations

            enddo ! end of the loop over all cells along y.
          enddo ! end of the loop over all cells along x.
        enddo ! end of the loop 1 over the nzon types of sources.
        dy=dx
        omefov=0.00000001 ! solid angle of the spectrometer slit on the sky. Here we only need a small value
        z_obs=z_o+altsol(x_obs,y_obs) ! z_obs = the local observer elevation plus the height of observation above ground (z_o)
        rx_obs=real(x_obs)*dx
        ry_obs=real(y_obs)*dy
        if (z_obs.eq.0.) z_obs=0.001
        largx=dx*real(nbx) ! computation of the Width along x of the case.
        largy=dy*real(nby) ! computation of the Width along y of the case.
        write(2,*) 'Width of the domain [NS](m):',largx,'#cases:',nbx
        write(2,*) 'Width of the domain [EO](m):',largy,'#cases:',nby
        write(2,*) 'Size of a cell (m):',dx,' X ',dy
        write(2,*) 'latitu center:',latitu
        direct=0. ! initialize the total direct radiance from sources to observer
        rdirect=0. ! initialize the total reflected radiance from surface to observer
        irdirect=0. ! initialize the total direct irradiance from sources to observer
        irrdirect=0. ! initialize the total reflected irradiance from surface to observer
! Calculation of the direct radiances and irradiances
        if (verbose.ge.1) print*,' Calculating obtrusive light...'
        do stype=1,ntype ! beginning of the loop over the source types
          if (totlu(stype).ne.0.) then ! check if there are any flux in that source type otherwise skip this lamp
            if (verbose.ge.1) print*,' Turning on lamps zone',stype
            if (verbose.ge.1) write(2,*) ' Turning on lamps zone',stype
            do x_s=imin(stype),imax(stype) ! beginning of the loop over the source in x
              do y_s=jmin(stype),jmax(stype) ! beginning of the loop over source in y
                rx_s=real(x_s)*dx
                ry_s=real(y_s)*dy
                if (lamplu(x_s,y_s,stype).ne.0.) then ! ignore null luminosity
                  z_s=(altsol(x_s,y_s)+lampal(x_s,y_s)) ! Definition of the position (metre) vertical of the source.
                  rx=rx_obs+20000.*ix
                  ry=ry_obs+20000.*iy
                  rz=z_obs+20000.*iz
                  call anglezenithal(rx_s,ry_s,z_s,rx_obs,ry_obs,z_obs,zenith)
                  call blocking(x_s,y_s,z_s,x_obs,y_obs,z_obs,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2)
                  ! projection angle of line to the lamp and the viewing angle
                  call angle3points(rx_s,ry_s,z_s,rx_obs,ry_obs,z_obs,rx,ry,rz,dang) ! scattering angle.
                  dang=pi-dang
                  ! computation of the solid angle of the line of sight voxel seen from the source
                  anglez=idnint(180.*zenith/pi)+1
                  P_dir=pvalno(anglez,stype)
                  ! computation of the flux direct reaching the line of sight voxel
                  if ((cos(dang).gt.0.).and.(dang.lt.pi/2.)) then
                    ddir_obs=sqrt((rx_obs-rx_s)**2.+(ry_obs-ry_s)**2.+(z_obs-z_s)**2.) ! distance direct sight between source and observer
                    ! computation of the solid angle 1m^2 at the observer as seen from the source
                    omega=1.*abs(cos(dang))/ddir_obs**2.
                    call transmitm(zenith,z_obs,z_s,ddir_obs,transm,tranam,tabs)
                    call transmita(zenith,z_obs,z_s,ddir_obs,haer,transa,tranaa)
                    call transmita(zenith,z_obs,z_s,ddir_obs,hlay,transl,tranal)
                    if (dang.lt.dfov) then ! check if the reflecting surface enter the field of view of the observer
                      direct=direct+lamplu(x_s,y_s,stype)*transa*transm*transl*P_dir*omega*(1.-ff1)*(1.-ff2)*hh/(pi*dfov**2.) ! correction for obstacle filling factor
                    endif
                    irdirect=irdirect+lamplu(x_s,y_s,stype)*transa*transm*transl*P_dir*omega*(1.-ff1)*(1.-ff2)*hh ! correction for obstacle filling factor
                  endif
                  xsrmi=x_s-boxx
                  if (xsrmi.lt.1) xsrmi=1
                  xsrma=x_s+boxx
                  if (xsrma.gt.nbx) xsrma=nbx
                  ysrmi=y_s-boxy
                  if (ysrmi.lt.1) ysrmi=1
                  ysrma=y_s+boxy
                  if (ysrma.gt.nby) ysrma=nby
                  do x_sr=xsrmi,xsrma ! beginning of the loop reflection x
                    rx_sr=real(x_sr)*dx
                    do y_sr=ysrmi,ysrma ! beginning of the loop reflection y
                      ry_sr=real(y_sr)*dy
                      irefl=0.
                      z_sr=altsol(x_sr,y_sr)
                      if((x_sr.gt.nbx).or.(x_sr.lt.1).or.(y_sr.gt.nby).or.(y_sr.lt.1)) then ! condition out of borders
                        if (verbose.eq.2) then
                          print*,'Ground cell out of borders'
                        endif
                      else
                        if((x_s.eq.x_sr).and.(y_s.eq.y_sr).and.(z_s.eq.z_sr)) then
                          if (verbose.eq.2) then
                            print*,'Source pos = Ground cell'
                          endif
                        else
                          ! if haut is negative, the ground cell is lighted from below
                          haut=-(rx_s-rx_sr)*tan(inclix(x_sr,y_sr))-(ry_s-ry_sr)*tan(incliy(x_sr,y_sr))+z_s-z_sr
                          if (haut.gt.0.) then ! Condition: the ground cell is illuminated from above
                            call anglezenithal(rx_s,ry_s,z_s,rx_sr,ry_sr,z_sr,angzen) ! computation of the zenithal angle between the source ground surface
                            distd=sqrt((rx_s-rx_sr)**2.+(ry_s-ry_sr)**2.+(z_s-z_sr)**2.)
                            call transmitm(angzen,z_s,z_sr,distd,transm,tranam,tabs) ! computation of the transmittance between the source and the ground surface
                            call transmita(angzen,z_s,z_sr,distd,haer,transa,tranaa)
                            call transmita(angzen,z_s,z_sr,distd,hlay,transl,tranal)
                            ! computation of the solid angle of the groud surface seen from the source
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
                            call anglesolide(omega,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z) ! Call of the routine anglesolide to compute the angle solide.
                            if (omega.lt.0.) then
                              print*,'ERROR: Solid angle of the reflecting surface < 0.'
                              stop
                            endif
                            ! estimation of the half of the underlying angle of the solid angle       
                            ! this angle allow an estimate of the average emission function over the ground surface                                                                      
                            ouvang=sqrt(omega/pi) ! Angle in radian.
                            ouvang=ouvang*180./pi ! Angle in degrees.
                            ! computation of the photometric function of the light fixture toward the ground surface
                            anglez=idnint(180.*angzen/pi)
                            if (anglez.lt.0)anglez=-anglez
                            if (anglez.gt.180) anglez=360-anglez
                            anglez=anglez+1 ! Transform the angle in integer degree into the position in the array.
                            ! average +- ouvang
                            naz=0
                            nbang=0.
                            P_indir=0.
                            do na=-idnint(ouvang),idnint(ouvang)
                              naz=anglez+na
                              if (naz.lt.0) naz=-naz
                              if (naz.gt.181) naz=362-naz ! symetric function
                              if (naz.eq.0) naz=1
                              P_indir=P_indir+pvalno(naz,stype)*abs(sin(pi*real(naz)/180.))/2.
                              nbang=nbang+1.*abs(sin(pi*real(naz)/180.))/2.
                            enddo
                            P_indir=P_indir/nbang
                            ! computation of the flux reaching the reflecting ground surface
                            flrefl=lamplu(x_s,y_s,stype)*P_indir*omega*transm*transa*transl
                            ! computation of the reflected intensity leaving the ground surface
                            irefl1=flrefl*srefl/pi ! The factor 1/pi comes from the normalisation of the fonction
                            call anglezenithal(rx_sr,ry_sr,z_sr,rx_obs,ry_obs,z_obs,zenith)
                            call blocking(x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,dx,dy,nbx,nby,altsol,drefle,ofill,obsH,hh,ff1,ff2) 
                            ! projection angle of line to the lamp and the viewing angle
                            call angle3points(rx_sr,ry_sr,z_sr,rx_obs,ry_obs,z_obs,rx,ry,rz,dang)
                            dang=pi-dang
                            if ((cos(dang).gt.0.).and.(dang.lt.pi/2.)) then
                              ddir_obs=sqrt((rx_obs-rx_sr)**2.+(ry_obs-ry_sr)**2.+(z_obs-z_sr)**2.) ! distance direct sight between source and observer
                              omega=1.*abs(cos(dang))/ddir_obs**2.
                              call transmitm(zenith,z_obs,z_sr,ddir_obs,transm,tranam,tabs)
                              call transmita(zenith,z_obs,z_sr,ddir_obs,haer,transa,tranaa)
                              call transmita(zenith,z_obs,z_sr,ddir_obs,hlay,transl,tranal)
                              if (dang.lt.dfov) then ! check if the reflecting surface enter the field of view of the observer
                                rdirect=rdirect+irefl1*omega*transa*transm*transl*hh*(1.-ff1)*(1.-ff2)/(pi*dfov**2.)
                              endif
                              irrdirect=irrdirect+irefl1*omega*transa*transm*transl*hh*(1.-ff1)*(1.-ff2)
                            endif ! Condition: the ground cell is illuminated from above
                          endif ! end haut
                        endif
                      endif ! end of out of borders
                    enddo ! end of the loop reflection y
                  enddo ! end of the loop reflection x
                endif ! ignore null luminosity
              enddo ! end of the loop over source in y
            enddo ! end of the loop over the source in x
          endif ! check if there are any flux in that source type otherwise skip this lamp
        enddo ! end of the loop over the source types.
        ! end of direct calculations
          
          

        itot=0. 
        itoclou=0.
        fctcld=0.
        fcapt=0.
        ftocap=0.   
        tot1=0.
        tot2=0.
        tot3=0.     
        if (ssswit.gt.0) then ! if any need for scattering calculations
! Calculation of the scattered radiances - 1st, 2nd, and 3rd scattering
          cloudtop=100000.
          if ((z_obs.ge.cloudbase).and.(z_obs.le.cloudtop)) then
            print*,'The observer is inside the cloud! Abort computing.',z_obs,cloudbase
            stop
          endif
          call horizon(x_obs,y_obs,z_obs,dx,dy,altsol,angaz1,zhoriz,dhmax) ! calculating the distance before the line of sight beeing blocked by topography
          rx_c=real(x_obs)*dx-ix*scal/2.
          ry_c=real(y_obs)*dx-iy*scal/2.
          z_c=z_obs-iz*scal/2.
          do icible=1,ncible ! beginning of the loop over the line of sight voxels
            itodif1=0.
            itodif2=0.
            itodif3=0.
            rx_c=rx_c+ix*(scalo/2.+scal/2.)
            ry_c=ry_c+iy*(scalo/2.+scal/2.)
            dh0=sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2)
            if ((dh0.le.dhmax).or.((dh0.gt.dhmax).and.(angze1-zhoriz.lt.0.00001))) then ! the line of sight is not yet blocked by the topography
              x_c=idnint(rx_c/dx)
              if (x_c.lt.1) x_c=1
              if (x_c.gt.width) x_c=width
              y_c=idnint(ry_c/dy)
              if (y_c.lt.1) y_c=1
              if (y_c.gt.width) y_c=width
              z_c=z_c+iz*(scalo/2.+scal/2.)
              if (z_c.gt.altsol(x_c,y_c)) then
                if ((fcapt.ge.ftocap/stoplim).and.(z_c.lt.cloudbase).and.(z_c.lt.35000.)) then 
                  ! stop the calculation of the viewing line when the increment is lower than 1/stoplim
                  ! or when hitting a cloud or when z>35km (scattering probability =0 (given precision)
                  fcapt=0.
                  do i=1,nbx
                    do j=1,nby
                      FCA(i,j)=0.
                    enddo
                  enddo
                  ! Calculate the solid angle of the line of sight voxel unit voxel
                  ! (1 m^3) given the fixed FOV of the observer.
                  ! For line of sight voxel near the observer
                  ! we need to calculate the scattering on a part of the voxel. For far
                  ! voxels we may be needed to increase the solid angle since the FOV can
                  ! encompass more than the voxel size. This correction is done with the
                  ! portio parameter calculated as the ratio of the solid angle of the
                  ! observer FOV over the line of sight voxel solid angle as seen from the
                  ! observer.
                  distd=sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.+(z_c-z_obs)**2.)
                  ! computation of the Solid angle of the line of sight voxel seen from the observer
                  omega=1./distd**2.
                  if (omega.gt.omemax) then
                    omega=0.
                    portio=0.
                  else
                    portio=(omefov/omega)
                  endif

                  if( (rx_c.gt.real(nbx*dx)).or.(rx_c.lt.dx).or.(ry_c.gt.(nby*dy)).or.(ry_c.lt.dy)) then ! Condition line of sight inside the modelling domain
                  else
                    if (verbose.ge.1) then 
                      print*,'================================================'
                      print*,' Progression along the line of sight :',icible
                      print*,' Horizontal dist. line of sight =',idnint(sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.)),' m'
                      print*,' Vertical dist. line of sight =',idnint(abs(z_c-z_obs)),' m'
                      write(2,*) '============================================='
                      write(2,*) ' Progression along the line of sight :',icible
                      write(2,*) ' Horizontal dist. line of sight =',idnint(sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.)),' m'
                      write(2,*) ' Vertical dist. line of sight =',idnint(abs(z_c-z_obs)),' m'
                    endif
                    dis_obs=sqrt((z_c-z_obs)**2.+(ry_c-ry_obs)**2.+(rx_c-rx_obs)**2.)
                    if (dis_obs.eq.0.) then
                      print*,'ERROR problem with dis_obs',dis_obs
                      stop
                    endif
                    ometif=pi*(diamobj/2.)**2./dis_obs**2.
  
                    do stype=1,ntype ! beginning of the loop over the source types.
                      if (totlu(stype).ne.0.) then ! check if there are any flux in that source type otherwise skip this lamp
                        if (verbose.ge.1) print*,' Turning on lamps zone',stype
                        if (verbose.ge.1) write(2,*) ' Turning on lamps zone',stype

                        do x_s=1,nbx ! the intensity toward the sensor by a line of sight voxel.
                          do y_s=1,nby
                            ITT(x_s,y_s,stype)=0.
                          enddo
                        enddo
                        do x_s=imin(stype),imax(stype) ! beginning of the loop over the column (longitude the) of the domain.
                          do y_s=jmin(stype),jmax(stype)       
! 1st scattering from source and ground                              
                            rho=0 ! from source                   
                            icloud=0.
                            if (lamplu(x_s,y_s,stype).ne.0.) then ! if the luminosite of the case is null, the program ignore this case.
                              z_s=(altsol(x_s,y_s)+lampal(x_s,y_s)) ! Definition of the position (metre) vertical of the source
                              dirck=0
                              if ((x_s.eq.x_c).and.(y_s.eq.y_c).and.(z_s.eq.z_c)) then ! if the source is not at the line of sight voxel
                                dirck=1
                                if (verbose.ge.1) then
                                  print*,'Source = line of sight'
                                endif
                              endif   
                              if (dirck.ne.1) then ! the source is not at the line of sight voxel position
                                idif1=0.                             
                                call firstscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,iz,lamplu,ofill, &
                                srefl,drefle,reflsiz,obsH,altsol,inclix,incliy,pvalno,stype,tranam,tabs,tranaa,tranal,secdif, &
                                secdil,fdifan,fdifl,haer,hlay,dx,dy,cloudt,cloudbase,omefov,scal,portio,idif1,icloud)
                                itoclou=itoclou+icloud
                                itodif1=itodif1+idif1 ! Correct the result for the skipping of scattering voxels to accelerate the calculation 
                                ! outputs are idif1 and icloud
                                rho=1 ! from ground
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
                                    icloud=0.
                                    z_sr=altsol(x_sr,y_sr)
                                    if ((x_sr.eq.x_s).and.(y_sr.eq.y_s).and.(z_sr.eq.z_s)) then 
                                    else
                                      haut=-(rx_s-rx_sr)*tan(inclix(x_sr,y_sr))-(ry_s-ry_sr)*tan(incliy(x_sr,y_sr))+z_s-z_sr ! if haut < 0, ground cell is lighted from below
                                      if (haut.gt.0.) then ! Condition: the ground cell is lighted from above  
                                        if((x_sr.le.nbx).and.(x_sr.ge.1).and.(y_sr.le.nby).and.(y_sr.ge.1)) then
                                          if((x_c.ne.x_sr).and.(y_c.ne.y_sr).and.(z_c.ne.z_sr)) then       ! if the ground is not at the line of sight voxel  
                                            idif1=0.                         
                                            call firstscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,iz, &
                                            lamplu,ofill,srefl,drefle,reflsiz,obsH,altsol,inclix,incliy,pvalno,stype,tranam,tabs, &
                                            tranaa,tranal,secdif,secdil,fdifan,fdifl,haer,hlay,dx,dy,cloudt,cloudbase,omefov, &
                                            scal,portio,idif1,icloud)
                                            itoclou=itoclou+icloud
                                            itodif1=itodif1+idif1
                                          endif
                                        endif
                                      endif
                                    endif
                                  enddo
                                enddo  
                           
                                           
! 2nd scattering from source and ground
                                if (ssswit.gt.1) then
                                  rho=0 ! from source
                                  icloud=0.
                                  idif2=0.
                                  call secondscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,iz,lamplu, &
                                  ofill,srefl,drefle,reflsiz,obsH,altsol,inclix,incliy,pvalno,stype,zondi2,siz2,volu2,ndiff2, &
                                  tranam,tabs,tranaa,tranal,secdif,secdil,fdifan,fdifl,haer,hlay,dx,dy,cloudt,cloudbase,omefov, &
                                  scal,portio,idif2,icloud)
                                  itoclou=itoclou+icloud
                                  itodif2=itodif2+idif2 ! Correct the result for the skipping of scattering voxels to accelerate the calculation
                                  ! outputs are idif2 and icloud
                                  rho=1 ! from ground
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
                                      icloud=0.
                                      z_sr=altsol(x_sr,y_sr)
                                      if ((x_sr.eq.x_s).and.(y_sr.eq.y_s).and.(z_sr.eq.z_s)) then ! if the ground is not at the line of sight voxel
                                      else
                                        haut=-(rx_s-rx_sr)*tan(inclix(x_sr,y_sr))-(ry_s-ry_sr)*tan(incliy(x_sr,y_sr))+z_s-z_sr ! if haut < 0, ground cell is lighted from below
                                        if (haut.gt.0.) then ! Condition: the ground cell is lighted from above                                      
                                          if((x_sr.le.nbx).and.(x_sr.ge.1).and.(y_sr.le.nby).and.(y_sr.ge.1)) then
                                            if((x_c.ne.x_sr).and.(y_c.ne.y_sr).and.(z_c.ne.z_sr)) then 
                                              idif2=0.
                                              call secondscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,iz, &
                                              lamplu,ofill,srefl,drefle,reflsiz,obsH,altsol,inclix,incliy,pvalno,stype,zondi2, &
                                              siz2,volu2,ndiff2,tranam,tabs,tranaa,tranal,secdif,secdil,fdifan,fdifl,haer,hlay, &
                                              dx,dy,cloudt,cloudbase,omefov,scal,portio,idif2,icloud)
                                              itoclou=itoclou+icloud
                                              itodif2=itodif2+idif2 ! Correct for the skipping of scat voxels to accelerate the calculation                                
                                            endif
                                          endif
                                        endif
                                      endif                                  
                                    enddo
                                  enddo        
                                endif ! end 2nd scat
! 3rd scattering from source and ground                             
                                if (ssswit.gt.2) then
                                  rho=0 ! from source
                                  icloud=0.                                 
                                  call thirdscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs,iz,lamplu, &
                                  ofill,srefl,drefle,reflsiz,obsH,altsol,inclix,incliy,pvalno,stype,zondi3,siz3,volu3, &
                                  ndiff3,tranam,tabs,tranaa,tranal,secdif,secdil,fdifan,fdifl,haer,hlay,dx,dy,cloudt,cloudbase, &
                                  omefov,scal,portio,idif3,icloud)
                                  itodif3=itodif3+idif3
                                  itoclou=itoclou+icloud
                                  ! outputs are idif3 and icloud
                                  rho=1 ! from ground
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
                                      icloud=0.
                                      z_sr=altsol(x_sr,y_sr)
                                      if ((x_sr.eq.x_s).and.(y_sr.eq.y_s).and.(z_sr.eq.z_s)) then 
                                      else
                                        haut=-(rx_s-rx_sr)*tan(inclix(x_sr,y_sr))-(ry_s-ry_sr)*tan(incliy(x_sr,y_sr))+z_s-z_sr ! if haut < 0, ground cell is lighted from below
                                        if (haut.gt.0.) then ! Condition: the ground cell is lighted from above
                                          if((x_sr.le.nbx).and.(x_sr.ge.1).and.(y_sr.le.nby).and.(y_sr.ge.1)) then
                                            if((x_c.ne.x_sr).and.(y_c.ne.y_sr).and.(z_c.ne.z_sr)) then    ! if the ground is not at the line of sight voxel                              
                                              call thirdscat(rho,x_s,y_s,z_s,x_c,y_c,z_c,x_sr,y_sr,z_sr,x_obs,y_obs,z_obs, &
                                              iz,lamplu,ofill,srefl,drefle,reflsiz,obsH,altsol,inclix,incliy,pvalno,stype, &
                                              zondi3,siz3,volu3,ndiff3,tranam,tabs,tranaa,tranal,secdif,secdil,fdifan,fdifl,haer, &
                                              hlay,dx,dy,cloudt,cloudbase,omefov,scal,portio,idif3,icloud)
                                              itodif3=itodif3+idif3
                                              itoclou=itoclou+icloud    
                                            endif
                                          endif
                                        endif
                                      endif
                                    enddo
                                  enddo
                                endif ! end 3rd scat
                                ! computation of the total light coming from a source
                                ! In the order 1st scat; 2nd scat; 3rd scat
                              else
                                itodif1=0.
                                itodif2=0.
                                itodif3=0.
                              endif ! end the source is not at the line of sight voxel position
                              if (verbose.eq.2) then
                                print*,' Total intensity per component for type ',ntype,':'
                                print*,' First scattering=',itodif1
                                print*,' Second scattering=',itodif2
                                print*,' Third scattering=',itodif3
                                if ((itodif1.lt.0.).or.(itodif2.lt.0.).and.(itodif3.lt.0.)) then
                                  print*,'PROBLEM! Negative intensity.'
                                  stop
                                endif
                              endif

                            
                                                        


                            

                            !ITT(x_s,y_s,stype)=ITT(x_s,y_s,stype)+isourc ! ITT stores itotty in a matrix
                            endif
                          enddo ! end the loop over the lines (latitude) of the domain (y_s).
                        enddo ! end the loop over the column (longitude) of the domain (x_s).

                        
                           !if (icible.eq.9)   stop
                        
                        

                        do x_s=imin(stype),imax(stype)
                          do y_s=jmin(stype),jmax(stype)
                            !ITC(x_s,y_s)=ITC(x_s,y_s)+ITT(x_s,y_s,stype)
                          enddo
                        enddo
                        ! calculate total lamp flux matrix for all lamp types
                        do x_s=1,nbx
                          do y_s=1,nby
                            !lpluto(x_s,y_s)=lpluto(x_s,y_s)+lamplu(x_s,y_s,stype)
                          enddo
                        enddo
                      endif ! end of condition if there are any flux in that source type
                      
                    enddo ! end of the loop over the types of sources (stype).
                    itot=itodif1+itodif2+itodif3
                    tot1=tot1+itodif1
                    tot2=tot2+itodif2
                    tot3=tot3+itodif3

                    
                    ! computation of the luminous flux reaching the observer
                    ! computation of the zenithal angle between the observer and the line of sight voxel
                    call anglezenithal(rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,angzen) ! computation of the zenithal angle between the line of sight voxel and the observer.
                    ! computation of the transmittance between the line of sight voxel and the observer
                    distd=sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.+(z_c-z_obs)**2.)
                    call transmitm(angzen,z_c,z_obs,distd,transm,tranam,tabs)
                    call transmita(angzen,z_c,z_obs,distd,haer,transa,tranaa)
                    call transmita(angzen,z_c,z_obs,distd,hlay,transl,tranal)
                    ! computation of the flux reaching the objective of the telescope from the line of sight voxel
                    fcapt=itot*ometif*transa*transm*transl ! computation of the flux reaching the intrument from the line of sight voxel
                    
                    

                    
                    
                    do x_s=1,nbx
                      do y_s=1,nby
                        !FCA(x_s,y_s)=ITC(x_s,y_s)*ometif*transa*transm*transl
                      enddo
                    enddo
                    if (cos(pi-angzen).eq.0.) then
                      print*,'ERROR perfectly horizontal sight is forbidden'
                      stop
                    endif
                    ! end of the computation of the flux reaching the observer voxel from the line of sight voxel
                    ftocap=ftocap+fcapt ! flux for all source all type all line of sight element
                    
                    
                    
                    !print*,'fcatp',icible,fcapt,ftocap,itotci,ometif,transa,transm,transl
                    !if (icible.eq.10) stop
                    
                    
                    do x_s=1,nbx
                      do y_s=1,nby
                        !FTC(x_s,y_s)=FTC(x_s,y_s)+FCA(x_s,y_s) 
                        ! FTC is the array of the flux total at the sensor to identify
                        ! the contribution of each ground pixel to the total flux at the observer level
                        ! The % is simply given by the ratio FTC/ftocap
                      enddo
                    enddo
                    ! correction for the FOV to the flux reaching the intrument from the cloud voxel
                    if (cloudt.ne.0) then
                      ! computation of the flux reaching the intrument from the cloud voxel
                      fccld=itoclou*ometif*transa*transm*transl
                      fctcld=fctcld+fccld ! cloud flux for all source all type all line of sight element
                    endif
                    
                    
                                             
                    
                    
                    
                    
                    if (verbose.ge.1) print*,'Added radiance =',fcapt/omefov/(pi*(diamobj/2.)**2.)
                    if (verbose.ge.1) print*,'Radiance accumulated =',ftocap/omefov/(pi*(diamobj/2.)**2.)
                    if (verbose.ge.1) write(2,*) 'Added radiance =',fcapt/omefov/(pi*(diamobj/2.)**2.)
                    if (verbose.ge.1) write(2,*) 'Radiance accumulated =',ftocap/omefov/(pi*(diamobj/2.)**2.)
                  endif ! end of the condition line of sight voxel inside the modelling domain
                endif ! end condition line of sight voxel 1/stoplim
                ! accelerate the computation as we get away from the sources
                scalo=scal
                if (scal.le.3000.)  scal=scal*1.12
              endif
            endif ! end the line of sight is not yet blocked by the topography
          enddo ! end of the loop over the line of sight voxels.
      
                          print*,tot1,tot2,tot3
          print*,'Ratio 2nd/1st scat=',tot2/tot1
          print*,'Ratio 3rd/1st scat=',tot3/tot1
        
      
          fctcld=fctcld*10**(0.4*(100.-cloudfrac)*cloudslope) ! correction for the cloud fraction (defined from 0 to 100)
          if (prmaps.eq.1) then
            do x_s=1,nbx
              do y_s=1,nby
                FTC(x_s,y_s)=FTC(x_s,y_s)/ftocap ! Here FTC becomes the flux fraction of each pixel. The sum of FTC values over all pixels give the total flux
              enddo
            enddo
            if (verbose.eq.2) then
              print*,'Writing normalized contribution array'
              print*,'Warning Cloud contrib. excluded from that array.'
            endif
            call twodout(nbx,nby,pclimg,FTC)
          endif ! end of condition for creating contrib and sensit maps
        endif ! end of scattered light
        ! End of calculation of the scattered radiances
        if (verbose.ge.1) print*,'====================================================='
        print*,'         Direct irradiance from sources (W/m**2/nm)'
        write(*,2001)  irdirect
        print*,'       Direct irradiance from reflexion (W/m**2/nm)'
        write(*,2001)  irrdirect
        print*,'         Direct radiance from sources (W/str/m**2/nm)'
        write(*,2001)  direct
        print*,'         Direct radiance from reflexion (W/str/m**2/nm)'
        write(*,2001)  rdirect
        print*,'             Cloud radiance (W/str/m**2/nm)'
        write(*,2001) fctcld/omefov/(pi*(diamobj/2.)**2.)
        print*,'            Diffuse radiance (W/str/m**2/nm) including clouds'
        write(*,2001) (ftocap+fctcld)/omefov/(pi*(diamobj/2.)**2.)
        if (verbose.ge.1) write(2,*) '==================================================='
        write(2,*) '     Direct irradiance from sources (W/m**2/nm)'
        write(2,2001)  irdirect
        write(2,*) '     Direct irradiance from reflexion (W/m**2/nm)'
        write(2,2001)  irrdirect
        write(2,*) '     Direct radiance from sources (W/str/m**2/nm)'
        write(2,2001)  direct
        write(2,*) '     Direct radiance from reflexion (W/str/m**2/nm)'
        write(2,2001)  rdirect
        write(2,*) '           Cloud radiance (W/str/m**2/nm)         '
        write(2,2001) fctcld/omefov/(pi*(diamobj/2.)**2.)
        write(2,*) '         Diffuse radiance (W/str/m**2/nm) including clouds       '
        write(2,2001) (ftocap+fctcld)/omefov/(pi*(diamobj/2.)**2.)
      close(2)
2001  format('                   ',E10.3E2)
      stop
      end
! ***********************************************************************
! *                                                                      
! *                                         end                          
! *                                                                      
! ***********************************************************************
