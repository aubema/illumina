c                         *                          *                       iii                   *     *
c                                                                           iiiii
c  IIIIII    lLLLL    *    lLLLL         UUU    UUU      MMMMM      MMMMM    iii        NNNN     NN          AAAA
c   IIII     LLLL          LLLL   *     UUU      UUU     MMMMMMM  MMMMMMM          *    NNNNN    NN        AAAaaAAA
c   IIII     LLLL          LLLL        UUU *      UUU    MMM MMMMMMMM MMM    iii        NNNNNN   NN       AAA    AAA
c   IIII     LLLL   *      LLLL        UUU        UUU    MMM *        MMM  iii          NNN  NNN NN     AAAAAAAAAAAAAA
c   IIII     LLLl          LLLl        UUUu      uUUU    MMM          MMM  iiii    ii   NNN   NNNNN    AAAa        aAAA
c   IIII    LLLLLLLLLL    LLLLLLLLLL    UUUUUuuUUUUU     MMM          MMM   iiiiiiiii   NNN    NNNN   aAAA    *     AAAa
c  IIIIII   LLLLLLLLLLL   LLLLLLLLLLL     UUUUUUUU      mMMMm        mMMMm   iiiiiii   nNNNn    NNNn  aAAA          AAAa
c
c **********************************************************************************************************************
c ** Illumina VERSION 2 - in Fortran 77                                                                               **
c ** Programmers in decreasing order of contribution  :                                                               **
c **                            Martin Aube                                                                           **
c **              Still having very few traces of their contributions :                                               **
c **                            Loic Franchomme-Fosse,  Mathieu Provencher, Andre Morin                               **
c **                            Alex Neron, Etienne Rousseau                                                          **
c **                            William Desroches, Maxime Girardin, Tom Neron                                         **
c **                                                                                                                  **
c ** Illumina can be downloaded via:   git clone https://github.com/aubema/illumina.git                               **
c ** To compile:                                                                                                      **
c **    cd hg/illumina                                                                                                **
c **    mkdir bin                                                                                                     **
c **    bash makeILLUMINA                                                                                             **
c **                                                                                                                  **
c **  Current version features/limitations :                                                                          **
c **                                                                                                                  **
c **    - Calculation of artificial sky radiance in a given line of sight                                             **
c **    - Calculation of the atmospheric transmittance and 1st and 2nd order of scattering                            **
c **    - Lambertian reflexion on the ground                                                                          **
c **    - Terrain slope considered (apparent surface and shadows)                                                     **
c **    - Angular photometry of a lamp is considered uniform along the azimuth                                        **
c **    - Sub-grid obstacles considered (with the mean free path of light toward ground, mean obstacle height, and    **
c **      obstacles transparency (filling factor)                                                                     **
c **    - Molecules and aerosol optics (phase function, scattering probability, aerosol absorption)                   **
c **    - Exponential concentrations vertical profile                                                                 **
c **    - Accounting for heterogeneity luminaires number, luminaires heights, luminaire spectrum,                     **
c **      angular photometry, obstacle properties                                                                     **
c **    - Wavelength dependant                                                                                        **
c **    - Cloud models (type and cloud base height) only the overhead clouds are considered with cloud fraction       **
c **    - Support direct observation of a source                                                                      **
c **    - Direct observation of the ground is implemented                                                             **
c **                                                                                                                  **
c **********************************************************************************************************************
c
c  Copyright (C) 2021 Martin Aube PhD
c
c  This program is free software: you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation, either version 3 of the License, or
c  (at your option) any later version.
c
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c
c  You should have received a copy of the GNU General Public License
c  along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c  Contact: martin.aube@cegepsherbrooke.qc.ca
c
c
c
      program illumina                                                    ! Beginning
      implicit none
c
c=======================================================================
c     Variables declaration
c=======================================================================
c
      integer width,nzon                                                  ! Matrix dimension in Length/width and height
      parameter (width=512,nzon=256)
      integer iun,ideux
      real pi,pix4
      real zero,un                                                        ! value of 0. and 1.
      integer verbose                                                     ! verbose = 1 to have more print out, 0 for silent
      parameter (pi=3.141592654)
      parameter (pix4=4.*pi)
      character*72 mnaf                                                   ! Terrain elevation file
      character*72 diffil                                                 ! Aerosol file
      character*72 outfile                                                ! Results file
      character*72 pclf,pclgp                                             ! Files containing contribution and sensitivity maps
      character*72 pclimg,pcwimg
      character*72 basenm                                                 ! Base name of files
      integer lenbase                                                     ! Length of the Base name of the experiment
      real lambda,pressi,drefle(width,width)                              ! Wavelength (nanometer), atmospheric pressure (kPa), mean free path to the ground (meter).
      real reflsiz                                                        ! Size of the reflecting surface
      integer ntype                                                       ! Number of light source types or zones considered
      real largx                                                          ! Width (x axis) of the modeling domain (meter)
      real largy                                                          ! Length (y axis) of the modeling domain (meter)
      integer nbx,nby                                                     ! Number of pixels in the modeling domain
      real val2d(width,width)                                             ! Temporary input array 2d
      real altsol(width,width)                                            ! Ground elevation (meter)
      real srefl                                                          ! Ground reflectance
      integer stype                                                       ! Source type or zone index
      character*72 pafile,lufile,alfile,ohfile,odfile,offile              ! Files related to light sources and obstacles (photometric function of the sources (sr-1), flux (W), height (m), obstacles c                                                               ! height (m), obstacle distance (m), obstacle filling factor (0-1).
      real lamplu(width,width,nzon)                                       ! Source fluxes
      real lampal(width,width)                                            ! Height of the light sources relative to the ground (meter)
      real pval(181,nzon),pvalto,pvalno(181,nzon)                         ! Values of the angular photometry functions (unnormalized, integral, normalized)
      real dtheta                                                         ! Angle increment of the photometric function of the sources
      real dx,dy,dxp,dyp                                                  ! Width of the voxel (meter)
      integer boxx,boxy                                                   ! reflection window size (pixels)
      real fdifa(181),fdifan(181)                                         ! Aerosol scattering functions (unnormalized and normalized)
      real extinc,scatte,anglea(181)                                      ! Aerosol cross sections (extinction and scattering), scattering angle (degree)
      real secdif                                                         ! Contribution of the scattering to the extinction
      real inclix(width,width)                                            ! tilt of the ground pixel along x (radian)
      real incliy(width,width)                                            ! tilt of the ground pixel along y (radian)
      integer x_obs,y_obs                                                 ! Position of the observer (INTEGER)
      real rx_obs,ry_obs
      real z_o                                                            ! observer height relative to the ground (meter)
      real z_obs                                                          ! Height of the observer (meter) to the vertical grid scale
      integer ncible,icible                                               ! Number of line of sight voxels, number loops over the voxels
      integer x_c,y_c                                                     ! Position of the line of sight voxel (INTEGER)
      real rx_c,ry_c
      real z_c                                                            ! Height of the line of sight voxel (meter)
      integer dirck                                                       ! Test for the position of the source (case source=line of sight voxel)
      integer x_s,y_s,x_sr,y_sr,x_dif,y_dif,zceldi                        ! Positions of the source, the reflecting surface, and the scattering voxels
      real z_s,z_sr,z_dif                                                 ! Heights of the source, the reflecting surface, and the scattering voxel (metre).
      real rx_s,ry_s,rx_sr,ry_sr,rx_dif,ry_dif
      real angzen,ouvang                                                  ! Zenithal angle between two voxels (radians) and opening angle of the solid angle in degrees.
      integer anglez                                                      ! Emitting zenithal angle from the luminaire.
      real P_dir,P_indir,P_dif1                                           ! photometric function of the light sources (direct,indirect,scattered)
      real transa,transm,transl                                           ! Transmittance between two voxels (aerosols,molecules,particle layer).
      real tran1a,tran1m                                                  ! Transmittance of the voxel (aerosols,molecules).
      real taua                                                           ! Aerosol optical depth @ 500nm.
      real alpha                                                          ! Angstrom coefficient of aerosol AOD
      real*8 xc,yc,zc,xn,yn,zn                                            ! Position (meter) of the elements (starting point, final point) for the calculation of the solid angle.
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! Components of the vectors used in the solid angle calculation routine.
      real omega,omega1                                                   ! Solid angles
      real fldir                                                          ! Flux coming from a source (watt).
      real flindi                                                         ! Flux coming from a reflecting ground element (watt).
      real fldiff                                                         ! Flux coming from a scattering voxel (watt).
      real zidif,zfdif                                                    ! initial and final limits of a scattering path.
      real angdif                                                         ! scattering angle.
      real pdifdi,pdifin,pdifd1,pdifd2                                    ! scattering probability (direct,indirect,1st and 2nd order of scattering
      real intdir                                                         ! Direct intensity toward the sensor from a scattering voxel.
      real intind                                                         ! Contribution of the reflecting cell to the reflected intensity toward the sensor.
      real itotind                                                        ! Total contribution of the source to the reflected intensity toward the sensor.
      real idiff2                                                         ! Contribution of the scattering voxel to the scattered intensity toward the sensor.
      real itodif                                                         ! Total contribution of the source to the scattered intensity toward the sensor.
      real isourc                                                         ! Total contribution of the source to the intensity from a line of sight voxel toward the sensor.
      real itotty                                                         ! Total contribution of a source type to the intensity coming from a line of sight voxel toward the sensor.
      real itotci                                                         ! total intensity from a line of sight voxel toward the sensor.
      real itotrd                                                         ! total intensity a voxel toward the sensor after reflexion and double scattering.
      real flcib                                                          ! Flux reaching the observer voxel from a line of sight voxel.
      real fcapt                                                          ! Flux reaching the observer voxel from all FOV voxels in a given model level
      real ftocap                                                         ! Total flux reaching the observer voxel
      real haut                                                           ! Haut (negative indicate that the surface is lighted from inside the ground. I.e. not considered in the calculation
      real epsilx,epsily                                                  ! tilt of the ground pixel
      real flrefl                                                         ! flux reaching a reflecting surface (watts).
      real irefl,irefl1                                                   ! intensity leaving a reflecting surface toward the line of sight voxel.
      real effdif                                                         ! Distance around the source voxel and line of sight voxel considered to compute the 2nd order of scattering.
      real zondif(3000000,3)                                              ! Array for the scattering voxels positions
      integer ndiff,idi                                                   ! Number of scattering voxels, counter of the loop over the scattering voxels
      integer stepdi                                                      ! scattering step to speedup the calculation e.g. if =2 one computation over two will be done
      integer ssswit                                                      ! activate double scattering (1=yes, 0=no)
      integer fsswit                                                      ! activate first scattering (1=yes, 0=no)
      integer nvis0                                                       ! starting value for the calculation along of the viewing line.
c                                                                         ! by default the value is 1 but it can be larger
c                                                                         ! when we resume a previous interrupted calculation.
      real fldif1,fldif2                                                  ! flux reaching a scattering voxel.
      real fdif2                                                          ! flux reaching the line of sight voxel after reflexion > scattering
      real idif1,idif2,idif2p                                             ! intensity toward a line of sight voxel from a scattering voxel (without and with reflexion).
      real portio                                                         ! ratio of voxel surface to the solid angle of the sensor field of view.
      real dis_obs                                                        ! Distance between the line of sight and the observer.
      real ometif                                                         ! Solid angle of the telescope objective as seen from the line of sight voxel
      real omefov                                                         ! Solid angle of the spectrometer slit.
      real angvis,azim                                                    ! viewing angles of the sensor.
c                                                                         ! Useful for the calculation of the lambertian reflectance.
      real nbang                                                          ! for the averaging of the photometric function
      real obsH(width,width),angmin                                       ! averaged height of the sub-grid obstacles, minimum angle under wich
c                                                                         ! a light ray cannot propagate because it is blocked by a sub-grid obstable
      real ofill(width,width)                                             ! fill factor giving the probability to hit an obstacle when pointing in its direction real 0-1
      integer naz,na
      real ITT(width,width,nzon)                                          ! total intensity per type of lamp
      real ITC(width,width)                                               ! total intensity per line of sight voxel
      real FTC(width,width)                                               ! fraction of the total flux at the sensor level
      real FCA(width,width)                                               ! sensor flux array
      real lpluto(width,width)                                            ! total luminosity of the ground cell for all lamps
      character*3 lampno                                                  ! lamp number string
      integer imin(nzon),imax(nzon),jmin(nzon),jmax(nzon)                 ! x and y limits containing a type of lamp
      real angazi                                                         ! azimuth angle between two points in rad, max dist for the horizon determination
      real latitu                                                         ! approximate latitude of the domain center
      integer prmaps                                                      ! flag to enable the tracking of contribution and sensitivity maps
      integer cloudt                                                      ! cloud type 0=clear, 1=Thin Cirrus/Cirrostratus, 2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus,
                                                                          ! 4=Stratocumulus/stratus, 5=Cumulus/Cumulonimbus
      real cloudslope                                                     ! slope of the radiance dependency on the cloud fraction (in percentage) According to
                                                                          ! Sciezoras 2020 the slope vary depending on the level of LP and how it is distributed.
                                                                          ! We decided instead to simplify this by using an average slope of -0.013.
                                                                          ! Rad=Rad_100 * 10**(0.4*(100-cloudfrac)*cloudslope) this equation is derived from
                                                                          ! Tomasz Sciezor , The impact of clouds on the brightness of the night sky, Journal of
                                                                          ! Quantitative Spectroscopy & Radiative Transfer (2020),
                                                                          ! doi: https://doi.org/10.1016/j.jqsrt.2020.106962
      real cloudfrac                                                      ! cloud fraction in percentage
      integer xsrmi,xsrma,ysrmi,ysrma                                     ! limits of the loop valeur for the reflecting surfaces
      real rcloud                                                         ! cloud relfectance
      real azencl                                                         ! zenith angle from cloud to observer
      real icloud                                                         ! cloud reflected intensity
      real fcloud                                                         ! flux reaching the intrument from the cloud voxel
      real fccld                                                          ! correction for the FOV to the flux reaching the intrument from the cloud voxel
      real fctcld                                                         ! total flux from cloud at the sensor level
      real totlu(nzon)                                                    ! total flux of a source type
      real stoplim                                                        ! Stop computation when the new voxel contribution is less than 1/stoplim of the cumulated flux
      real ff,ff2,hh                                                          ! temporary obstacle filling factor and horizon blocking factor
      real cloudbase,cloudtop,cloudhei                                    ! cloud base and top altitude (m), cloud layer avg height (m)
      real distd                                                          ! distance to compute the scattering probability
      real volu                                                           ! volume of a voxel
      real scal                                                           ! stepping along the line of sight
      real scalo                                                          ! previous value of scal
      real siz                                                            ! resolution of the 2nd scat grid in meter
      real angvi1,angaz1,angze1                                           ! viewing angles in radian
      real ix,iy,iz                                                       ! base vector of the viewing (length=1)
      real dsc2,doc2                                                      ! square of the path lengths for the cloud contribution
      real azcl1,azcl2                                                    ! zenith angle from the (source, refl surface, or scattering voxel) to line of path and observer to line p.
      real dh,dho                                                         ! distance of the horizon limit
      integer n2nd                                                        ! desired number of voxel in the calculation of the 2nd scattering
      integer step                                                        ! skiping 2nd scat on 1 dim
      real omemax                                                         ! max solid angle allowed
      real tcloud                                                         ! low cloud transmission
      real rx_sp,ry_sp                                                    ! position of a low cloud pixel
      real flcld(width,width)                                             ! flux crossing a low cloud
      real ds1,ds2,ds3,dss                                                ! double scattering distances
      integer nss                                                         ! number of skipped 2nd scat elements
      integer ndi                                                         ! number of cell under ground
      integer nvol                                                        ! number of cell for second scat calc un full resolution
      real diamobj                                                        ! instrument objective diameter
      integer i,j,k,id,jd
      real tranam,tranaa                                                  ! atmospheric transmittancess of a path (molecular, aerosol)
      real zhoriz                                                         ! zenith angle of the horizon
      real direct                                                         ! direct radiance from sources on a surface normal to the line of sight (no scattering)
      real rdirect                                                        ! direct radiance from a reflecting surface on a surface normal to the line of sight (no scattering)
      real irdirect                                                       ! direct irradiance from sources on a surface normal to the line of sight (no scattering)
      real irrdirect                                                      ! direct irradiance from a reflecting surface on a surface normal to the line of sight (no scattering)
      real dang                                                           ! Angle between the line of sight and the direction of a source
      real dzen                                                           ! zenith angle of the source-observer line
      real ddir_obs                                                       ! distance between the source and the observer
      real rx,ry,rz                                                       ! driving vector for the calculation of the projection angle for direct radiance. It is 20km long
      real dfov                                                           ! field of view in degrees for the calculation of the direct radiance this number will be a kind of smoothing effect. The angular grid resolution to create a direct radiance panorama should be finer than that number
      real Fo                                                             ! flux correction factor for obstacles
      real thetali                                                        ! limit angle for the obstacles blocking of viirs
      integer viirs(width,width)                                          ! viirs flag 1=yes 0=no
      character*72 vifile                                                 ! name of the viirs flag file
      real dh0,dhmax                                                      ! horizontal distance along the line of sight and maximum distance before beeing blocked by topography
      character*72 layfile                                                ! filename of the optical properties of the particle layer
      real layaod                                                         ! 500 nm aod of the particle layer
      real layalp                                                         ! spectral exponent of the aod for the particle layer
      real hlay                                                           ! exponential vertical scale height of the particle layer
      real secdil                                                         ! scattering/extinction ratio for the particle layer
      real fdifl(181)                                                     ! scattering phase function of the particle layer
      real tranal                                                         ! top of atmos transmission of the particle layer
      real haer                                                           ! exponential vertical scale height of the background aerosol layer
      real distc,hcur                                                     ! distance to any cell and curvature  correction for the earth curvature
      real bandw                                                          ! bandwidth of the spectral bin
      real tabs                                                           ! TOA transmittance related to molecule absorption
      integer obsobs                                                      ! flag to activate the direct light obstacle blocking aroud the observer.
      verbose=1                                                           ! Very little printout=0, Many printout = 1, even more=2
      diamobj=1.                                                          ! A dummy value for the diameter of the objective of the instrument used by the observer.
      volu=0.
      zero=0.
      un=1.
      ff=0.
      ff2=0.
      step=1
      ncible=1024
      stepdi=1
      cloudslope=-0.013
      cloudfrac=100.
      if (verbose.ge.1) then
        print*,'Starting ILLUMINA computations...'
      endif
c reading of the fichier d'entree (illumina.in)
      print*,'Reading illumina.in input file'
      open(unit=1,file='illumina.in',status='old')
        read(1,*)
        read(1,*) basenm
        read(1,*) dx,dy
        read(1,*) diffil
        read(1,*) layfile, layaod, layalp, hlay
        read(1,*) ssswit
        read(1,*) fsswit
        read(1,*) lambda,bandw
        read(1,*) srefl
        read(1,*) pressi
        read(1,*) taua,alpha,haer
        read(1,*) ntype
        read(1,*) stoplim
        read(1,*)
        read(1,*) x_obs,y_obs,z_o
        read(1,*) obsobs
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
c conversion of the geographical viewing angles toward the cartesian
c angle we assume that the angle in the file illumina.in
c is consistent with the geographical definition
c geographical, azim=0 toward north, 90 toward east, 180 toward south
c etc
c cartesian, azim=0 toward east, 90 toward north, 180 toward west etc
      azim=90.-azim
      if (azim.lt.0.) azim=azim+360.
      if (azim.ge.360.) azim=azim-360.
      angvi1 = (pi*angvis)/180.
      angze1 = pi/2.-angvi1
      angaz1 = (pi*azim)/180.
      ix = ( sin((pi/2.)-angvi1) ) * (cos(angaz1))                        ! viewing vector components
      iy = ( sin((pi/2.)-angvi1) ) * (sin(angaz1))
      iz = (sin(angvi1))
      dfov=(dfov*pi/180.)/2.
      siz=2500.
      if (ssswit.eq.0) then
        effdif=0.
      else
        effdif=40000.
      endif
      scal=19.
      scalo=scal
      boxx=nint(reflsiz/dx)                                               ! Number of column to consider left/right of the source for the reflection.
      boxy=nint(reflsiz/dy)                                               ! Number of column to consider up/down of the source for the reflection.
c omemax: exclude calculations too close (<10m) this is a sustended angle of 1 deg.
c the calculated flux is highly sensitive to that number for a very high
c pixel resolution (a few 10th of meters). We assume anyway that somebody
c observing the sky will never lies closer than that distance to a
c light fixture. This number is however somehow subjective and that means
c that the value of sky brightness near sources will be affected by this
c choice
      omemax=1./((10.)**2.)
      if (verbose.gt.0) then
        print*,'2nd order scattering grid = ',siz,'m'
        print*,'2nd order scattering radius=',effdif,'m'
        print*,'Pixel size = ',dx,' x ',dy
        print*,'Maximum radius for reflection = ',reflsiz
      endif
c computing the actual AOD at the wavelength lambda
      if (verbose.ge.1) print*,'500nm AOD=',taua,'500nm angstrom coeff.=
     +',alpha
      taua=taua*(lambda/500.)**(-1.*alpha)
      layaod=layaod*(lambda/500.)**(-1.*layalp)
c  determine the Length of basenm
      lenbase=index(basenm,' ')-1
      mnaf=basenm(1:lenbase)//'_topogra.bin'                              ! determine the names of input and output files
      outfile=basenm(1:lenbase)//'.out'
      pclf=basenm(1:lenbase)//'_pcl.txt'
      pclimg=basenm(1:lenbase)//'_pcl.bin'
      pcwimg=basenm(1:lenbase)//'_pcw.bin'
      pclgp=basenm(1:lenbase)//'_pcl.gplot'
c opening output file
      open(unit=2,file=outfile,status='unknown')
        write(2,*) "ILLUMINA version __version__"
        write(2,*) 'FILE USED:'
        write(2,*) mnaf,diffil
        print*,'Wavelength (nm):',lambda,
     +       ' Aerosol optical depth:',taua
        write(2,*) 'Wavelength (nm):',lambda,
     +       ' Aerosol optical depth:',taua
        write(2,*) '2nd order scattering radius:',effdif,' m'
        print*,'2nd order scattering radius:',effdif,' m'
        write(2,*) 'Observer position (x,y,z)',x_obs,y_obs,z_o
        print*,'Observer position (x,y,z)',x_obs,y_obs,z_o
        write(2,*) 'Elevation angle:',angvis,' azim angle (counterclockwise
     +from east)',azim
        print*,'Elevation angle:',angvis,' azim angle (counterclockwise
     +from east)',azim
c Initialisation of the arrays and variables
        if (verbose.ge.1) print*,'Initializing variables...'
        if (cloudt.eq.0) then
          cloudbase=1000000000.
        endif
        prmaps=1
        iun=0
        ideux=1
        icloud=0.
        do i=1,width
          do j=1,width
            val2d(i,j)=0.
            altsol(i,j)=0.
            obsH(i,j)=0.
            viirs(i,j)=0
            ofill(i,j)=0.
            inclix(i,j)=0.
            incliy(i,j)=0.
            lpluto(i,j)=0.
            ITC(i,j)=0.
            FTC(i,j)=0.
            FCA(i,j)=0.
            flcld(i,j)=0.
            do k=1,nzon
              lamplu(i,j,k)=0.
              lampal(i,j)=0.
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
            zondif(i,j)=1.
          enddo
        enddo
        idif1=0.
        idif2=0.
        fdif2=0
        idif2p=0.
        fldir=0.
        flindi=0.
        fldiff=0.
        pdifdi=0.
        pdifin=0.
        pdifd1=0.
        pdifd2=0.
        intdir=0.
        intind=0.
        idiff2=0.
        angmin=0.
        isourc=0.
        itotty=0.
        itotci=0.
        itotrd=0.
        flcib=0.
        flrefl=0.
        irefl=0.
        irefl1=0.
        fldif1=0.
        fldif2=0.
        portio=0.
        fccld=0.
        fctcld=0.
        ometif=0.
        omefov=0.
        hh=1.
c determine the 2nd scattering zone
        if (ssswit.ne.0) then
          call zone_diffusion(effdif,
     +    zondif,ndiff,stepdi,siz)
          dss=1.*siz
          if (verbose.gt.0) then
            print*,'2nd order scattering grid points =',ndiff
            print*,'2nd order scattering smoothing radius =',dss,'m'
          endif
        endif
c determination of the vertical atmospheric transmittance
        call transtoa(lambda,bandw,taua,layaod,pressi,tranam,tranaa,      ! tranam and tranaa are the top of atmosphere transmittance (molecules and aerosols)
     +tranal,tabs)

c reading of the environment variables
c reading of the elevation file
        call twodin(nbx,nby,mnaf,altsol)
c computation of the tilt of the pixels along x and along y
        do i=1,nbx                                                        ! beginning of the loop over the column (longitude) of the domain.
          do j=1,nby                                                      ! beginning of the loop over the rows (latitu) of the domain.
            if (i.eq.1) then                                              ! specific case close to the border of the domain (vertical side left).
              inclix(i,j)=atan((altsol(i+1,j)-altsol(i,j))/real(dx))      ! computation of the tilt along x of the surface.
            elseif (i.eq.nbx) then                                        ! specific case close to the border of the domain (vertical side right).
              inclix(i,j)=atan((altsol(i-1,j)-altsol(i,j))/(real(dx)))    ! computation of the tilt along x of the surface.
            else
              inclix(i,j)=atan((altsol(i+1,j)-altsol(i-1,j))/(2.          ! computation of the tilt along x of the surface.
     1        *real(dx)))
            endif
            if (j.eq.1) then                                              ! specific case close to the border of the domain (horizontal side down).
              incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j))/(real(dy)))    ! computation of the tilt along y of the surface.
            elseif (j.eq.nby) then                                        ! specific case close to the border of the domain (horizontal side up).
              incliy(i,j)=atan((altsol(i,j-1)-altsol(i,j))/(real(dy)))    ! computation of the tilt along y of the surface.
            else
              incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j-1))/(2.          ! computation of the tilt along y of the surface
     1        *real(dy)))
            endif
          enddo                                                           ! end of the loop over the rows (latitu) of the domain
        enddo                                                             ! end of the loop over the column (longitude) of the domain
c reading of the values of P(theta), height, luminosities and positions
c of the sources, obstacle height and distance
        ohfile=basenm(1:lenbase)//'_obsth.bin'
        odfile=basenm(1:lenbase)//'_obstd.bin'
        alfile=basenm(1:lenbase)//'_altlp.bin'                            ! setting the file name of height of the sources lumineuse.
        offile=basenm(1:lenbase)//'_obstf.bin'
        vifile='origin.bin'
        dtheta=.017453293                                                 ! one degree
c reading lamp heights
        call twodin(nbx,nby,alfile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
          do j=1,nby                                                      ! beginning of the loop over all cells along y.
            lampal(i,j)=val2d(i,j)                                        ! filling of the array for the lamp stype
          enddo                                                           ! end of the loop over all cells along y.
        enddo                                                             ! end of the loop over all cells along x.
c reading subgrid obstacles average height
        call twodin(nbx,nby,ohfile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
          do j=1,nby                                                      ! beginning of the loop over all cells along y.
            obsH(i,j)=val2d(i,j)                                          ! filling of the array
          enddo                                                           ! end of the loop over all cells along y.
        enddo
c reading subgrid obstacles average distance
        call twodin(nbx,nby,odfile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
          do j=1,nby                                                      ! beginning of the loop over all cells along y.
            drefle(i,j)=val2d(i,j)/2.
            if (drefle(i,j).eq.0.) drefle(i,j)=dx                         ! when outside a zone, block to the size of the cell (typically 1km)
          enddo                                                           ! end of the loop over all cells along y.
        enddo
c reading subgrid obstacles filling factor
        call twodin(nbx,nby,offile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
          do j=1,nby                                                      ! beginning of the loop over all cells along y.
            ofill(i,j)=val2d(i,j)                                         ! Filling of the array 0-1
          enddo                                                           ! end of the loop over all cells along y.
        enddo
c reading viirs flag
        call twodin(nbx,nby,vifile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
          do j=1,nby                                                      ! beginning of the loop over all cells along y.
            viirs(i,j)=nint(val2d(i,j))                                   ! viirs flag array 0 or 1
          enddo                                                           ! end of the loop over all cells along y.
        enddo
c reading of the scattering parameters for background aerosols
        open(unit = 1, file = diffil,status= 'old')                       ! opening file containing the scattering parameters
          read(1,*)  secdif                                               ! the scattering / extinction ratio
          read(1,*)
          do i=1,181
            read(1,*) anglea(i), fdifa(i)                                 ! reading of the scattering functions
            fdifan(i)=fdifa(i)/pix4                                       ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle
          enddo
        close(1)
c reading scattering parameters of particle layer
         open(unit = 1, file = layfile,status= 'old')                     ! opening file containing the scattering parameters
          read(1,*)  secdil                                               ! the scattering / extinction ratio of particle layer
          read(1,*)
          do i=1,181
            read(1,*) anglea(i), fdifl(i)                                 ! reading of the scattering functions of the particle layer
            fdifl(i)=fdifl(i)/pix4                                        ! The integral of the imported phase fonction over sphere = 4 pi) We divide by 4 pi to get it per unit of solid angle
          enddo
        close(1)
c Some preliminary tasks
        do stype=1,ntype                                                  ! beginning of the loop 1 for the nzon types of sources.
          imin(stype)=nbx
          jmin(stype)=nby
          imax(stype)=1
          jmax(stype)=1
          pvalto=0.
          write(lampno, '(I3.3)' ) stype                                  ! support of nzon different sources (3 digits)
          pafile=basenm(1:lenbase)//'_fctem_'//lampno//'.dat'             ! setting the file name of angular photometry.
          lufile=basenm(1:lenbase)//'_lumlp_'//lampno//'.bin'             ! setting the file name of the luminosite of the cases.
c reading photometry files
          open(UNIT=1, FILE=pafile,status='OLD')                          ! opening file pa#.dat, angular photometry.
            do i=1,181                                                    ! beginning of the loop for the 181 data points
              read(1,*) pval(i,stype)                                     ! reading of the data in the array pval.
              pvalto=pvalto+pval(i,stype)*2.*pi*                          ! Sum of the values of the  photometric function
     a        sin(real(i-1)*dtheta)*dtheta                                ! (pvaleur x 2pi x sin theta x dtheta) (ou theta egale (i-1) x 1 degrees).
            enddo                                                         ! end of the loop over the 181 donnees of the fichier pa#.dat.
          close(1)                                                        ! closing file pa#.dat, angular photometry.
          do i=1,181
            if (pvalto.ne.0.) pvalno(i,stype)=pval(i,stype)/pvalto        ! Normalisation of the photometric function.
          enddo
c reading luminosity files
          call twodin(nbx,nby,lufile,val2d)
          do i=1,nbx                                                      ! beginning of the loop over all cells along x.
            do j=1,nby                                                    ! beginning of the loop over all cells along y.
              if (val2d(i,j).lt.0.) then                                  ! searching of negative fluxes
                print*,'***Negative lamp flux!, stopping execution'
                stop
              endif
            enddo                                                         ! end of the loop over all cells along y.
          enddo
          do i=1,nbx                                                      ! searching of the smallest rectangle containing the zone
            do j=1,nby                                                    ! of non-null luminosity to speedup the calculation
              if (val2d(i,j).ne.0.) then
                if (i-1.lt.imin(stype)) imin(stype)=i-2
                if (imin(stype).lt.1) imin(stype)=1
                goto 333
              endif
            enddo
          enddo
          imin(stype)=1
 333      do i=nbx,1,-1
            do j=1,nby
              if (val2d(i,j).ne.0.) then
                if (i+1.gt.imax(stype)) imax(stype)=i+2
                if (imax(stype).gt.nbx) imax(stype)=nbx
                goto 334
              endif
            enddo
          enddo
          imax(stype)=1
 334      do j=1,nby
            do i=1,nbx
              if (val2d(i,j).ne.0.) then
                if (j-1.lt.jmin(stype)) jmin(stype)=j-2
                if (jmin(stype).lt.1) jmin(stype)=1
                goto 335
              endif
            enddo
          enddo
          jmin(stype)=1
 335      do j=nby,1,-1
            do i=1,nbx
              if (val2d(i,j).ne.0.) then
                if (j+1.gt.jmax(stype)) jmax(stype)=j+2
                if (jmax(stype).gt.nby) jmax(stype)=nby
                goto 336
              endif
            enddo
          enddo
          jmax(stype)=1
 336      do i=1,nbx                                                      ! beginning of the loop over all cells along x.
            do j=1,nby                                                    ! beginning of the loop over all cells along y.
              lamplu(i,j,stype)=val2d(i,j)                                ! remplir the array of the lamp type: stype
c Atmospheric correction and obstacles masking corrections to the lamp
c flux arrays (lumlp)
              if (viirs(i,j).eq.1) then
                lamplu(i,j,stype)=lamplu(i,j,stype)/(tranam*tranaa*
     +          tranal)
                thetali=atan2(drefle(i,j),obsH(i,j))
                if (thetali .lt. 70.*pi/180.) then
                  Fo=(1.-cos(70.*pi/180.))/(1.-ofill(i,j)*cos(thetali)+
     +            (ofill(i,j)-1.)*cos(70.*pi/180.))
                  lamplu(i,j,stype)=lamplu(i,j,stype)*Fo
                else
                  Fo=1.
                endif
              endif
              totlu(stype)=totlu(stype)+lamplu(i,j,stype)                 ! the total lamp flux should be non-null to proceed to the calculations
            enddo                                                         ! end of the loop over all cells along y.
          enddo                                                           ! end of the loop over all cells along x.
        enddo                                                             ! end of the loop 1 over the nzon types of sources.
        dy=dx
        omefov=0.00000001                                                 ! solid angle of the spectrometer slit on the sky. Here we only need a small value
        z_obs=z_o+altsol(x_obs,y_obs)                                     ! z_obs = the local observer elevation plus the height of observation above ground (z_o)
        rx_obs=real(x_obs)*dx
        ry_obs=real(y_obs)*dy
        if (z_obs.eq.0.) z_obs=0.001
        largx=dx*real(nbx)                                                ! computation of the Width along x of the case.
        largy=dy*real(nby)                                                ! computation of the Width along y of the case.
        write(2,*) 'Width of the domain [NS](m):',largx,'#cases:',nbx
        write(2,*) 'Width of the domain [EO](m):',largy,'#cases:',nby
        write(2,*) 'Size of a cell (m):',dx,' X ',dy
        write(2,*) 'latitu center:',latitu






        direct=0.                                                         ! initialize the total direct radiance from sources to observer
        rdirect=0.                                                        ! initialize the total reflected radiance from surface to observer
        irdirect=0.                                                       ! initialize the total direct irradiance from sources to observer
        irrdirect=0.                                                      ! initialize the total reflected irradiance from surface to observer
c =================================
c Calculation of the direct radiances
c
        if (verbose.ge.1) print*,' Calculating obtrusive light...'
        do stype=1,ntype                                                  ! beginning of the loop over the source types.
          if (totlu(stype).ne.0.) then                                    ! check if there are any flux in that source type otherwise skip this lamp
            if (verbose.ge.1) print*,' Turning on lamps',stype
            if (verbose.ge.1) write(2,*) ' Turning on lamps',
     +      stype
            do x_s=imin(stype),imax(stype)                                ! beginning of the loop over the column (longitude the) of the domain.
            do y_s=jmin(stype),jmax(stype)                                ! beginning of the loop over the rows (latitud) of the domain.
              intdir=0.
              itotind=0.
              itodif=0.
              itotrd=0.
              isourc=0.
              rx_s=real(x_s)*dx
              ry_s=real(y_s)*dy
              if (lamplu(x_s,y_s,stype) .ne. 0.) then                     ! if the luminosite of the case is null, the program ignore this case.
                z_s=(altsol(x_s,y_s)+lampal(x_s,y_s))                     ! Definition of the position (metre) vertical of the source.
c
c *********************************************************************************************************
c calculation of the direct radiance of sources falling on a surface perpendicular
c to the viewing angle Units of W/nm/m2/sr
c *********************************************************************************************************
                rx=rx_obs+20000.*ix
                ry=ry_obs+20000.*iy
                rz=z_obs+20000.*iz
                dho=sqrt((rx_obs-rx_s)**2.
     +          +(ry_obs-ry_s)**2.)
                if ((dho.gt.0.).and.(z_s.ne.z_obs)) then
                  call anglezenithal(rx_obs,ry_obs,z_obs                  ! zenithal angle source-observer
     +            ,rx_s,ry_s,z_s,dzen)
                  call angleazimutal(rx_obs,ry_obs,rx_s,                  ! computation of the angle azimutal direct line of sight-source
     +            ry_s,angazi)
                  if (dzen.gt.pi/4.) then                                 ! 45deg. it is unlikely to have a 1km high mountain less than 1
                    call horizon(x_obs,y_obs,z_obs,dx,dy,
     +              altsol,angazi,zhoriz,dh)
                    if (dh.le.dho) then
                      if (dzen-zhoriz.lt.0.00001) then                    ! shadow the path line of sight-source is not below the horizon => we compute
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
               if (obsobs.eq.1) then
c sub-grid obstacles
                  if (dho.gt.drefle(x_obs,y_obs)+drefle(x_s,y_s)) then    ! light path to observer larger than the mean free path -> subgrid obstacles
                    angmin=pi/2.-atan2((altsol(x_obs,y_obs)+
     +              obsH(x_obs,y_obs)-z_obs),drefle(x_obs,
     +              y_obs))
                    if (dzen.lt.angmin) then                              ! condition sub-grid obstacles direct.
                      ff=0.
                    else
                      ff=ofill(x_obs,y_obs)
                    endif
                  endif                                                   ! end light path to the observer larger than mean free path
               endif
                  
 
 
 
                  
                  call anglezenithal(rx_s,ry_s,z_s                        ! zenithal angle source-observer
     +            ,rx_obs,ry_obs,z_obs,dzen)                  
                  ff2=0.
                  if (dho.gt.drefle(x_s,y_s)) then                        ! light path from source larger than the mean free path -> subgrid obstacles
                    angmin=pi/2.-atan2((altsol(x_s,y_s)+
     +              obsH(x_s,y_s)-z_s),drefle(x_s,
     +              y_s))
                    if (dzen.lt.angmin) then                              ! condition sub-grid obstacles direct.
                      ff2=0.
                    else
                      ff2=ofill(x_s,y_s)
                    endif
                  endif                                                   ! end light path to the observer larger than mean free path                  
                  call anglezenithal(rx_obs,ry_obs,z_obs                  ! zenithal angle source-observer
     +            ,rx_s,ry_s,z_s,dzen)                  
                  
                  
                  
                  
                  
c projection angle of line to the lamp and the viewing angle
                  call angle3points (rx_s,ry_s,z_s,rx_obs,                ! scattering angle.
     +            ry_obs,z_obs,rx,ry,rz,dang)
                  dang=pi-dang
c computation of the solid angle of the line of sight voxel seen from the source
                  anglez=nint(180.*(pi-dzen)/pi)+1
                  P_dir=pvalno(anglez,stype)
c computation of the flux direct reaching the line of sight voxel
                  if ((cos(dang).gt.0.).and.(dang.lt.pi/2.))
     +            then
                    ddir_obs=sqrt((rx_obs-rx_s)**2.+                      ! distance direct sight between source and observer
     +              (ry_obs-ry_s)**2.+(z_obs-z_s)**2.)
c computation of the solid angle 1m^2 at the observer as seen from the source
                    omega=1.*abs(cos(dang))/ddir_obs**2.
                    call transmitm(dzen,z_obs,z_s,ddir_obs,
     +              transm,tranam,tabs)
                    call transmita(dzen,z_obs,z_s,ddir_obs,
     +              haer,transa,tranaa)
                    call transmitl(dzen,z_obs,z_s,ddir_obs,
     +              hlay,transl,tranal)
                    if (dang.lt.dfov) then                                ! check if the reflecting surface enter the field of view of the observer
                      direct=direct+lamplu(x_s,y_s,stype)*
     +                transa*transm*transl*P_dir*omega*(1.-ff)*(1.-ff2)
     +                *hh/(pi*dfov**2.)                                      ! correction for obstacle filling factor
                    endif
                    irdirect=irdirect+lamplu(x_s,y_s,stype)*
     +              transa*transm*transl*P_dir*omega*(1.-ff)*(1.-ff2)*hh    ! correction for obstacle filling factor
                  endif
                endif
c
c **********************************************************************************
c * computation of the direct light toward the observer by the ground reflection   *
c **********************************************************************************
c
                xsrmi=x_s-boxx
                if (xsrmi.lt.1) xsrmi=1
                xsrma=x_s+boxx
                if (xsrma.gt.nbx) xsrma=nbx
                ysrmi=y_s-boxy
                if (ysrmi.lt.1) ysrmi=1
                ysrma=y_s+boxy
                if (ysrma.gt.nby) ysrma=nby
                do x_sr=xsrmi,xsrma                                       ! beginning of the loop over the column (longitude) reflecting.
                  rx_sr=real(x_sr)*dx
                  do y_sr=ysrmi,ysrma                                     ! beginning of the loop over the rows (latitu) reflecting.
                    ry_sr=real(y_sr)*dy
                    irefl=0.
                    z_sr=altsol(x_sr,y_sr)
                    if((x_sr.gt.nbx).or.(x_sr.lt.1).or.
     +              (y_sr.gt.nby).or.(y_sr.lt.1)) then
                        if (verbose.eq.2) then
                          print*,'Ground cell out of borders'
                        endif
                    else
                        if((x_s.eq.x_sr).and.(y_s.eq.y_sr)
     +                  .and.(z_s.eq.z_sr)) then
                          if (verbose.eq.2) then
                            print*,'Source pos = Ground cell'
                          endif
                        else
                          haut=-(rx_s-rx_sr)*tan(                         ! if haut is negative, the ground cell is lighted from below
     +                    inclix(x_sr,y_sr))-(ry_s-
     +                    ry_sr)*tan(incliy(x_sr,
     +                    y_sr))+z_s-z_sr
                          if (haut .gt. 0.) then                          ! Condition: the ground cell is lighted from above
c computation of the zenithal angle between the source and the surface reflectance
                            call anglezenithal(rx_s,ry_s,                 ! computation of the zenithal angle between the source and the line of sight voxel.
     +                      z_s,rx_sr,ry_sr,z_sr,                         ! end of the case "observer at the same latitu/longitude than the source".
     +                      angzen)
c computation of the transmittance between the source and the ground surface
                            distd=sqrt((rx_s-rx_sr)**2.
     +                      +(ry_s-ry_sr)**2.+
     +                      (z_s-z_sr)**2.)
                            call transmitm(angzen,z_s,
     +                      z_sr,distd,transm,tranam,tabs)
                            call transmita(angzen,z_s,
     +                      z_sr,distd,haer,transa,tranaa)
                            call transmitl(angzen,z_s,z_sr,distd,
     +                      hlay,transl,tranal)
c computation of the solid angle of the reflecting cell seen from the source
                            xc=dble(x_sr)*dble(dx)                        ! Position in meters of the observer voxel (longitude).
                            yc=dble(y_sr)*dble(dy)                        ! Position in meters of the observer voxel (latitu).
                            zc=dble(z_sr)                                 ! Position in meters of the observer voxel (altitude).
                            xn=dble(x_s)*dble(dx)                         ! Position in meters of the source (longitude).
                            yn=dble(y_s)*dble(dy)                         ! Position in meters of the source (latitu).
                            zn=dble(z_s)                                  ! Position in meters of the source (altitude).
                            epsilx=inclix(x_sr,y_sr)                      ! tilt along x of the ground reflectance
                            epsily=incliy(x_sr,y_sr)                      ! tilt along x of the ground reflectance
                            if (dx.gt.reflsiz) then                       ! use a sub-grid surface when the reflectance radius is smaller than the cell size
                              if ((x_sr.eq.x_s).and.(y_sr
     +                        .eq.y_s)) then
                                dxp=reflsiz
                              else
                                dxp=dx
                              endif
                            else
                              dxp=dx
                            endif
                            if (dy.gt.reflsiz) then
                              if ((x_sr.eq.x_s).and.(y_sr
     +                        .eq.y_s)) then
                                dyp=reflsiz
                              else
                                dyp=dy
                              endif
                            else
                              dyp=dy
                            endif
                            r1x=xc-dble(dxp)/2.-xn                        ! computation of the composante along x of the first vector.
                            r1y=yc+dble(dyp)/2.-yn                        ! computation of the composante along y of the first vector.
                            r1z=zc-tan(dble(epsilx))*
     +                      dble(dxp)/2.+tan(dble(epsily))                ! computation of the composante en z of the first vector.
     +                      *dble(dyp)/2.-zn
                            r2x=xc+dble(dxp)/2.-xn                        ! computation of the composante along x of the second vector.
                            r2y=yc+dble(dyp)/2.-yn                        ! computation of the composante along y of the second vector.
                            r2z=zc+tan(dble(epsilx))*
     +                      dble(dxp)/2.+tan(dble(epsily))                ! computation of the composante en z of the second vector.
     +                      *dble(dyp)/2.-zn
                            r3x=xc-dble(dxp)/2.-xn                        ! computation of the composante along x of the third vector.
                            r3y=yc-dble(dyp)/2.-yn                        ! computation of the composante along y of the third vector.
                            r3z=zc-tan(dble(epsilx))*
     +                      dble(dxp)/2.-tan(dble(epsily))                ! computation of the composante en z of the third vector.
     +                      *dble(dyp)/2.-zn
                            r4x=xc+dble(dxp)/2.-xn                        ! computation of the composante along x of the fourth vector.
                            r4y=yc-dble(dyp)/2.-yn                        ! computation of the composante along y of the fourth vector.
                            r4z=zc+tan(dble(epsilx))*
     +                      dble(dxp)/2.-tan(dble(epsily))                ! computation of the composante en z of the fourth vector.
     +                      *dble(dyp)/2.-zn
                            call anglesolide(omega,r1x,                   ! Call of the routine anglesolide to compute the angle solide.
     +                      r1y,r1z,r2x,r2y,r2z,r3x,r3y,
     +                      r3z,r4x,r4y,r4z)
         if (omega.lt.0.) then
           print*,'ERROR: Solid angle of the reflecting surface < 0.'
           stop
         endif
c estimation of the half of the underlying angle of the solid angle       ! this angle servira a obtenir un meilleur isime (moyenne) of
c                                                                         ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                            ouvang=sqrt(omega/pi)                         ! Angle in radian.
                            ouvang=ouvang*180./pi                         ! Angle in degrees.
c computation of the photometric function of the light fixture toward the reflection surface
c=======================================================================
c
                            anglez=nint(180.*angzen/pi)
                            if (anglez.lt.0)
     +                      anglez=-anglez
                            if (anglez.gt.180) anglez=360
     +                      -anglez
                            anglez=anglez+1                               ! Transform the angle in integer degree into the position in the array.
c average +- ouvang
                            naz=0
                            nbang=0.
                            P_indir=0.
                            do na=-nint(ouvang),nint(ouvang)
                              naz=anglez+na
                              if (naz.lt.0) naz=-naz
                              if (naz.gt.181) naz=362-naz                 ! symetric function
                              if (naz.eq.0) naz=1
                              P_indir=P_indir+pvalno(naz,
     +                        stype)*abs(sin(pi*real(naz)
     +                        /180.))/2.
                              nbang=nbang+1.*abs(sin(pi*
     +                        real(naz)/180.))/2.
                            enddo
                            P_indir=P_indir/nbang
c computation of the flux reaching the reflecting surface
                            flrefl=lamplu(x_s,y_s,stype)*
     +                      P_indir*omega*transm*transa*transl
c computation of the reflected intensity leaving the ground surface
                            irefl1=flrefl*srefl/pi                        ! The factor 1/pi comes from the normalisation of the fonction
c
c *********************************************************************************************************
c calculation of the direct radiance from reflection falling on a surface perpendicular
c to the viewing angle Units of W/nm/m2/sr
c *********************************************************************************************************
                            dho=sqrt((rx_obs-rx_sr)**2.
     +                      +(ry_obs-ry_sr)**2.)
                            if ((dho.gt.0.).and.(z_s.ne.z_obs)) then
                              call anglezenithal(rx_obs,ry_obs,z_obs      ! zenithal angle source-observer
     +                        ,rx_sr,ry_sr,z_sr,dzen)
                              call angleazimutal(rx_obs,ry_obs,rx_sr,     ! computation of the angle azimutal direct line of sight-source
     +                        ry_sr,angazi)
                              if (dzen.gt.pi/4.) then                     ! 45deg. it is unlikely to have a 1km high mountain less than 1
                                call horizon(x_obs,y_obs,z_obs,dx,dy,
     +                          altsol,angazi,zhoriz,dh)
                                if (dh.le.dho) then
                                  if (dzen-zhoriz.lt.0.00001) then        ! shadow the path line of sight-source is not below the horizon => we compute
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
c sub-grid obstacles
               ff=0.
               if (obsobs.eq.1) then
                              if (dho.gt.drefle(x_obs,y_obs)+
     +                        drefle(x_sr,y_sr)) then                     ! light path to observer larger than the mean free path -> subgrid obstacles
                                angmin=pi/2.-atan2((altsol(x_obs,y_obs)
     +                          +obsH(x_obs,y_obs)-z_obs),drefle(x_obs,
     +                          y_obs))
                                if (dzen.lt.angmin) then                  ! condition sub-grid obstacles direct.
                                  ff=0.
                                else
                                  ff=ofill(x_obs,y_obs)
                                endif
                              endif                                       ! end light path to the observer larger than mean free path
               endif
                              
                              
                              
                  call anglezenithal(rx_sr,ry_sr,z_sr                     ! zenithal angle surface-observer
     +            ,rx_obs,ry_obs,z_obs,dzen)                  
                  ff2=0.
                  if (dho.gt.drefle(x_sr,y_sr)) then                      ! light path from reflecting surface larger than the mean free path -> subgrid obstacles
                    angmin=pi/2.-atan2((altsol(x_sr,y_sr)+
     +              obsH(x_sr,y_sr)-z_sr),drefle(x_sr,
     +              y_sr))
                    if (dzen.lt.angmin) then                              ! condition sub-grid obstacles direct.
                      ff2=0.
                    else
                      ff2=ofill(x_sr,y_sr)
                    endif
                  endif                                                   ! end light path to the observer larger than mean free path                                                
                              call anglezenithal(rx_obs,ry_obs,z_obs      ! zenithal angle source-observer
     +                        ,rx_sr,ry_sr,z_sr,dzen)                              
                              
                              
                              
c projection angle of line to the lamp and the viewing angle
                              call angle3points (rx_sr,ry_sr,z_sr,        ! scattering angle.
     +                        rx_obs,ry_obs,z_obs,rx,ry,rz,dang)
                              dang=pi-dang

c computation of the flux direct reaching the line of sight voxel
                              if ((cos(dang).gt.0.).and.(dang.lt.pi/2.))
     +                        then
                                ddir_obs=sqrt((rx_obs-rx_sr)**2.+         ! distance direct sight between source and observer
     +                          (ry_obs-ry_sr)**2.+(z_obs-z_sr)**2.)
c computation of the solid angle of the line of sight voxel seen from the source
                                omega=1.*abs(cos(dang))/ddir_obs**2.
                                call transmitm(dzen,z_obs,z_sr,ddir_obs,
     +                          transm,tranam,tabs)
                                call transmita(dzen,z_obs,z_sr,ddir_obs,
     +                          haer,transa,tranaa)
                                call transmitl(dzen,z_obs,z_sr,ddir_obs,
     +                          hlay,transl,tranal)
                                if (dang.lt.dfov) then                    ! check if the reflecting surface enter the field of view of the observer
                                  rdirect=rdirect+irefl1*omega*transa*
     +                            transm*transl*hh*(1.-ff)*(1.-ff2)
     +                            /(pi*dfov**2.)
                                endif
                                irrdirect=irrdirect+irefl1*omega*transa*
     +                          transm*transl*hh*(1.-ff)*(1.-ff2)
                              endif

                            endif
                          endif
                        endif
                    endif
                  enddo
                enddo

              endif
            enddo
            enddo
          endif
        enddo
c
c End of calculation of the direct radiances
c =================================



        if (fsswit.ne.0) then
c =================================
c Calculation of the scattered radiances
c

c temporaire !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        cloudtop=100000.

        if ((z_obs.ge.cloudbase).and.(z_obs.le.cloudtop)) then
          print*,'The observer is inside the cloud! Abort computing.',
     +    z_obs,cloudbase
          stop
        endif
 1110   format(I4,1x,I4,1x,I4)
        fctcld=0.
        ftocap=0.                                                         ! Initialisation of the value of flux received by the sensor
        call horizon(x_obs,y_obs,z_obs,dx,dy,altsol,angaz1,zhoriz,        ! calculating the distance before the line of sight beeing blocked by topography
     +  dhmax)
        rx_c=real(x_obs)*dx-ix*scal/2.
        ry_c=real(y_obs)*dx-iy*scal/2.
        z_c=z_obs-iz*scal/2.
        do icible=1,ncible                                                ! beginning of the loop over the line of sight voxels
          rx_c=rx_c+ix*(scalo/2.+scal/2.)
          ry_c=ry_c+iy*(scalo/2.+scal/2.)
        dh0=sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2)
          if ((dh0.le.dhmax).or.((dh0.gt.dhmax).and.(angze1-zhoriz.lt.    ! the line of sight is not yet blocked by the topography
     +    0.00001))) then
          x_c=nint(rx_c/dx)
          if (x_c.lt.1) x_c=1
          if (x_c.gt.width) x_c=width
          y_c=nint(ry_c/dy)
          if (y_c.lt.1) y_c=1
          if (y_c.gt.width) y_c=width
          z_c=z_c+iz*(scalo/2.+scal/2.)
          if (z_c.gt.altsol(x_c,y_c)) then
          if ((fcapt.ge.ftocap/stoplim).and.(z_c.lt.cloudbase).and.       ! stop the calculation of the viewing line when the increment is lower than 1/stoplim
     +    (z_c.lt.35000.)) then                                           ! or when hitting a cloud or when z>40km (scattering probability =0 (given precision)
            fcapt=0.
            do i=1,nbx
              do j=1,nby
                FCA(i,j)=0.
              enddo
            enddo
c Calculate the solid angle of the line of sight voxel unit voxel
c (1 m^3) given the fixed FOV of the observer.
c For line of sight voxel near the observer
c we need to calculate the scattering on a part of the voxel. For far
c voxels we may be needed to increase the solid angle since the FOV can
c encompass more than the voxel size. This correction is done with the
c portio parameter calculated as the ratio of the solid angle of the
c observer FOV over the line of sight voxel solid angle as seen from the
c observer.
            distd=sqrt((rx_c-rx_obs)**2.+
     +      (ry_c-ry_obs)**2.+(z_c-z_obs)**2.)
c computation of the Solid angle of the line of sight voxel seen from the observer
            omega=1./distd**2.
            if (omega.gt.omemax) then
              omega=0.
              portio=0.
            else
              portio=(omefov/omega)
            endif
            itotci=0.                                                     ! Initialisation of the contribution of the line of sight at the sensor level
            do i=1,nbx
              do j=1,nby
                ITC(i,j)=0.
              enddo
            enddo
            if( (rx_c.gt.real(nbx*dx)).or.(rx_c.lt.dx).or.                ! Condition line of sight inside the modelling domain
     +      (ry_c.gt.(nby*dy)).or.(ry_c.lt.dy)) then
            else
              if (verbose.ge.1) print*,'================================
     +================'
      if (verbose.ge.1) print*,' Progression along the line of sight :'
     +,icible
      if (verbose.ge.1) print*,' Horizontal dist. line of sight =',
     +sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.),' m'
      if (verbose.ge.1) print*,' Vertical dist. line of sight =',
     +abs(z_c-z_obs),' m'
              if (verbose.ge.1) write(2,*) '========================
     +====================='
      if (verbose.ge.1) write(2,*) ' Progression along the line of sight
     + :',icible
      if (verbose.ge.1) write(2,*) ' Horizontal dist. line of sight =',
     +sqrt((rx_c-rx_obs)**2.+(ry_c-ry_obs)**2.),' m'
      if (verbose.ge.1) write(2,*) ' Vertical dist. line of sight =',
     +abs(z_c-z_obs),' m'
              dis_obs=sqrt((z_c-z_obs)**2.+(ry_c-ry_obs)**2.
     +          +(rx_c-rx_obs)**2.)
              if (dis_obs.eq.0.) then
                print*,'ERROR problem with dis_obs',dis_obs
                print*,rx_c,x_obs,y_c,y_obs,z_c,z_obs
                stop
              endif
              ometif=pi*(diamobj/2.)**2./dis_obs**2.
c beginning of the loop over the types of light sources
              do stype=1,ntype                                            ! beginning of the loop over the source types.
                if (totlu(stype).ne.0.) then                              ! check if there are any flux in that source type otherwise skip this lamp
                  if (verbose.ge.1) print*,' Turning on lamps',stype
                  if (verbose.ge.1) write(2,*) ' Turning on lamps',
     +            stype
                  itotty=0.                                               ! Initialisation of the contribution of a source types to
                  do x_s=1,nbx                                            ! the intensity toward the sensor by a line of sight voxel.
                    do y_s=1,nby
                      ITT(x_s,y_s,stype)=0.
                    enddo
                  enddo
                  do x_s=imin(stype),imax(stype)                          ! beginning of the loop over the column (longitude the) of the domain.
                    do y_s=jmin(stype),jmax(stype)                        ! beginning of the loop over the rows (latitud) of the domain.
                      intdir=0.
                      itotind=0.
                      itodif=0.
                      itotrd=0.
                      isourc=0.
                      rx_s=real(x_s)*dx
                      ry_s=real(y_s)*dy
                      if (lamplu(x_s,y_s,stype) .ne. 0.) then             ! if the luminosite of the case is null, the program ignore this case.
                        z_s=(altsol(x_s,y_s)+lampal(x_s,y_s))             ! Definition of the position (metre) vertical of the source.
c
c *********************************************************************************************************
c * computation of the scattered intensity toward the observer by a line of sight voxel from the source   *
c *********************************************************************************************************

                        dirck=0                                           ! Initialisation of the verification of the position of the source
                        if ((rx_s.eq.rx_c).and.(ry_s.eq.ry_c).and.            ! if the position of the source and the line of sight voxel are the
     +                  (z_s.eq.z_c))                                     ! same then...
     +                  then
                          dirck=1
                          if (verbose.ge.1) then
                            print*,'Source = line of sight'
                          endif
                        endif                                             ! end of the case positions x and y source and line of sight voxel identical.
                        if (dirck.ne.1) then                              ! the source is not at the line of sight voxel position
c computation of the zenithal angle between the source and the line of sight
c computation of the horizon for the resolved shadows direct              ! horizon resolution is 1 degree
                              distd=sqrt((rx_c-rx_s)**2.
     +                        +(ry_c-ry_s)**2.
     +                        +(z_c-z_s)**2.)
                              dho=sqrt((rx_c-rx_s)**2.
     +                        +(ry_c-ry_s)**2.)
                              call anglezenithal(rx_s,ry_s,z_s
     +                        ,rx_c,ry_c,z_c,angzen)                      ! computation of the zenithal angle between the source and the line of sight voxel.
                              call angleazimutal(rx_s,ry_s,rx_c,          ! computation of the angle azimutal direct line of sight-source
     +                        ry_c,angazi)
                              if (angzen.gt.pi/4.) then                   ! 45deg. it is unlikely to have a 1km high mountain less than 1
                                call horizon(x_s,y_s,z_s,dx,dy,altsol,
     +                          angazi,zhoriz,dh)
                                if (dh.le.dho) then
                                  if (angzen-zhoriz.lt.0.00001) then      ! shadow the path line of sight-source is not below the horizon => we compute
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
c sub-grid obstacles
                              ff=0.
                              if (dho.gt.drefle(x_s,y_s)) then            ! light path to observer larger than the mean free path -> subgrid obstacles
                                angmin=pi/2.-atan2((altsol(x_s,y_s)+
     +                          obsH(x_s,y_s)-z_s),drefle(x_s,y_s))
                                if (angzen.lt.angmin) then                ! condition sub-grid obstacles direct.
                                  ff=0.
                                else
                                  ff=ofill(x_s,y_s)
                                endif
                              endif
c computation of the transmittance between the source and the line of sight
                              call transmitm(angzen,z_s,z_c,distd,
     +                        transm,tranam,tabs)
                              call transmita(angzen,z_s,z_c,distd,
     +                        haer,transa,tranaa)
                              call transmitl(angzen,z_s,z_c,distd,
     +                        hlay,transl,tranal)
c computation of the solid angle of the line of sight voxel seen from the source
                              omega=1./distd**2.
                              if (omega.gt.omemax) omega=0.
                              anglez=nint(180.*angzen/pi)+1
                              P_dir=pvalno(anglez,stype)
c computation of the flux reaching the line of sight voxel
                              fldir=lamplu(x_s,y_s,stype)*P_dir*
     +                        omega*transm*transa*transl*(1.-ff)*hh       ! correction for obstacle filling factor
c computation of the scattering probability of the direct light
c distance pour traverser la cellule unitaire parfaitement orientée
                              if (omega.ne.0.) then
                                call angle3points (rx_s,ry_s,z_s,rx_c,    ! scattering angle.
     +                          ry_c,z_c,rx_obs,ry_obs,z_obs,
     +                          angdif)
                                call diffusion(angdif,                    ! scattering probability of the direct light. ############################################ secdif et un etaient inverses
     +                          tranam,tranaa,tranal,un,secdif,secdil,
     +                          fdifan,fdifl,haer,hlay,pdifdi,z_c)
                              else
                                pdifdi=0.
                              endif
c computation of the source contribution to the scattered intensity toward the sensor by a line of sight voxel
                              intdir=fldir*pdifdi
c contribution of the cloud reflection of the light coming directly from the source
                              if (cloudt.ne.0) then                       ! line of sight voxel = cloud
                                if (cloudbase-z_c.le.iz*scal) then
                                  call anglezenithal(rx_c,ry_c,z_c,
     +                            rx_obs,ry_obs,z_obs,azcl1)              ! zenith angle from cloud to observer
                                  call anglezenithal(rx_c,ry_c,z_c,
     +                            rx_s,ry_s,z_s,azcl2)                    ! zenith angle from source to cloud
                                  doc2=(rx_c-rx_obs)**2.+
     +                            (ry_c-ry_obs)**2.+(z_c-z_obs)**2.
                                  dsc2=(rx_s-rx_c)**2.+
     +                            (ry_s-ry_c)**2.+(z_s-z_c)**2.
                                  call cloudreflectance(angzen,           ! cloud intensity from direct illum
     +                            cloudt,rcloud)
                                  icloud=icloud+
     +                            fldir/omega*rcloud*doc2*omefov*
     +                            abs(cos(azcl2)/cos(azcl1))/dsc2/pi
                                endif
                              endif
                        else
                          intdir=0.
                        endif                                             ! end of the case Position Source is not equal to the line of sight voxel position
c end of the computation of the scattered intensity
c
c
c
c
c **********************************************************************************************************************
c * computation of the scattered light toward the observer by a line of sight voxel lighted by the ground reflection    *
c **********************************************************************************************************************
c etablissement of the conditions ands boucles
                            itotind=0.                                    ! Initialisation of the reflected intensity of the source
                            itotrd=0.
                            xsrmi=x_s-boxx
                            if (xsrmi.lt.1) xsrmi=1
                            xsrma=x_s+boxx
                            if (xsrma.gt.nbx) xsrma=nbx
                            ysrmi=y_s-boxy
                            if (ysrmi.lt.1) ysrmi=1
                            ysrma=y_s+boxy
                            if (ysrma.gt.nby) ysrma=nby
                            do x_sr=xsrmi,xsrma                           ! beginning of the loop over the column (longitude) reflecting.
                              rx_sr=real(x_sr)*dx
                              do y_sr=ysrmi,ysrma                         ! beginning of the loop over the rows (latitu) reflecting.
                                ry_sr=real(y_sr)*dy
                                irefl=0.
                                z_sr=altsol(x_sr,y_sr)
                                if((x_sr.gt.nbx).or.(x_sr.lt.1).or.
     +                          (y_sr.gt.nby).or.(y_sr.lt.1)) then
                                  if (verbose.eq.2) then
                                    print*,'Ground cell out of borders'
                                  endif
                                else
                                  if((x_s.eq.x_sr).and.(y_s.eq.y_sr)
     +                            .and.(z_s.eq.z_sr)) then
                                    if (verbose.eq.2) then
                                      print*,'Source pos = Ground cell'
                                    endif
                                  else
                                      haut=-(rx_s-rx_sr)*tan(             ! if haut is negative, the ground cell is lighted from below
     +                                inclix(x_sr,y_sr))-(ry_s-
     +                                ry_sr)*tan(incliy(x_sr,
     +                                y_sr))+z_s-z_sr
                                      if (haut .gt. 0.) then              ! Condition: the ground cell is lighted from above
c computation of the zenithal angle between the source and the surface reflectance
                                        call anglezenithal(rx_s,ry_s,     ! computation of the zenithal angle between the source and the line of sight voxel.
     +                                  z_s,rx_sr,ry_sr,z_sr,             ! end of the case "observer at the same latitu/longitude than the source".
     +                                  angzen)
c computation of the transmittance between the source and the ground surface
                                        distd=sqrt((rx_s-rx_sr)**2.
     +                                  +(ry_s-ry_sr)**2.+
     +                                  (z_s-z_sr)**2.)
                                        call transmitm(angzen,z_s,
     +                                  z_sr,distd,transm,tranam,tabs)
                                        call transmita(angzen,z_s,
     +                                  z_sr,distd,haer,transa,tranaa)
                                        call transmitl(angzen,z_s,z_sr,
     +                                  distd,hlay,transl,tranal)
c computation of the solid angle of the reflecting cell seen from the source
                                        xc=dble(x_sr)*dble(dx)            ! Position in meters of the observer voxel (longitude).
                                        yc=dble(y_sr)*dble(dy)            ! Position in meters of the observer voxel (latitu).
                                        zc=dble(z_sr)                     ! Position in meters of the observer voxel (altitude).
                                        xn=dble(x_s)*dble(dx)             ! Position in meters of the source (longitude).
                                        yn=dble(y_s)*dble(dy)             ! Position in meters of the source (latitu).
                                        zn=dble(z_s)                      ! Position in meters of the source (altitude).
                                        epsilx=inclix(x_sr,y_sr)          ! tilt along x of the ground reflectance
                                        epsily=incliy(x_sr,y_sr)          ! tilt along x of the ground reflectance
                                        if (dx.gt.reflsiz) then           ! use a sub-grid surface when the reflectance radius is smaller than the cell size
                                          if ((x_sr.eq.x_s).and.(y_sr
     +                                    .eq.y_s)) then
                                            dxp=reflsiz
                                          else
                                            dxp=dx
                                          endif
                                        else
                                          dxp=dx
                                        endif
                                        if (dy.gt.reflsiz) then
                                          if ((x_sr.eq.x_s).and.(y_sr
     +                                    .eq.y_s)) then
                                            dyp=reflsiz
                                          else
                                            dyp=dy
                                          endif
                                        else
                                          dyp=dy
                                        endif
                                        r1x=xc-dble(dxp)/2.-xn            ! computation of the composante along x of the first vector.
                                        r1y=yc+dble(dyp)/2.-yn            ! computation of the composante along y of the first vector.
                                        r1z=zc-tan(dble(epsilx))*
     +                                  dble(dxp)/2.+tan(dble(epsily))    ! computation of the composante en z of the first vector.
     +                                  *dble(dyp)/2.-zn
                                        r2x=xc+dble(dxp)/2.-xn            ! computation of the composante along x of the second vector.
                                        r2y=yc+dble(dyp)/2.-yn            ! computation of the composante along y of the second vector.
                                        r2z=zc+tan(dble(epsilx))*
     +                                  dble(dxp)/2.+tan(dble(epsily))    ! computation of the composante en z of the second vector.
     +                                  *dble(dyp)/2.-zn
                                        r3x=xc-dble(dxp)/2.-xn            ! computation of the composante along x of the third vector.
                                        r3y=yc-dble(dyp)/2.-yn            ! computation of the composante along y of the third vector.
                                        r3z=zc-tan(dble(epsilx))*
     +                                  dble(dxp)/2.-tan(dble(epsily))    ! computation of the composante en z of the third vector.
     +                                  *dble(dyp)/2.-zn
                                        r4x=xc+dble(dxp)/2.-xn            ! computation of the composante along x of the fourth vector.
                                        r4y=yc-dble(dyp)/2.-yn            ! computation of the composante along y of the fourth vector.
                                        r4z=zc+tan(dble(epsilx))*
     +                                  dble(dxp)/2.-tan(dble(epsily))    ! computation of the composante en z of the fourth vector.
     +                                  *dble(dyp)/2.-zn
                                        call anglesolide(omega,r1x,       ! Call of the routine anglesolide to compute the angle solide.
     +                                  r1y,r1z,r2x,r2y,r2z,r3x,r3y,
     +                                  r3z,r4x,r4y,r4z)
         if (omega.lt.0.) then
           print*,'ERROR: Solid angle of the reflecting surface < 0.'
           stop
         endif
c estimation of the half of the underlying angle of the solid angle       ! this angle servira a obtenir un meilleur isime (moyenne) of
c                                                                         ! P_dir for le cas of grans solid angles the , pvalno varie significativement sur +- ouvang.
                                        ouvang=sqrt(omega/pi)             ! Angle in radian.
                                        ouvang=ouvang*180./pi             ! Angle in degrees.
c computation of the photometric function of the light fixture toward the reflection surface
c=======================================================================
c
                                        anglez=nint(180.*angzen/pi)
                                        if (anglez.lt.0)
     +                                  anglez=-anglez
                                        if (anglez.gt.180) anglez=360
     +                                  -anglez
                                        anglez=anglez+1                   ! Transform the angle in integer degree into the position in the array.
c average +- ouvang
                                        naz=0
                                        nbang=0.
                                        P_indir=0.
                                        do na=-nint(ouvang),nint(ouvang)
                                          naz=anglez+na
                                          if (naz.lt.0) naz=-naz
                                          if (naz.gt.181) naz=362-naz     ! symetric function
                                          if (naz.eq.0) naz=1
                                          P_indir=P_indir+pvalno(naz,
     +                                    stype)*abs(sin(pi*real(naz)
     +                                    /180.))/2.
                                          nbang=nbang+1.*abs(sin(pi*
     +                                    real(naz)/180.))/2.
                                        enddo
                                        P_indir=P_indir/nbang
c computation of the flux reaching the reflecting surface
                                        flrefl=lamplu(x_s,y_s,stype)*
     +                                  P_indir*omega*transm*transa*
     +                                  transl
c computation of the reflected intensity leaving the ground surface
                                        irefl1=flrefl*srefl/pi            ! The factor 1/pi comes from the normalisation of the fonction
c
c
c
c
c
c **************************************************************************************
c * computation of the 2nd scattering contributions (2 order scattering and after reflection)
c **************************************************************************************
                                        if (effdif.gt.0.) then
      nss=0
      ndi=0
      do idi=1,ndiff                                                      ! beginning of the loop over the scattering voxels.
        rx_dif=zondif(idi,1)+(rx_s+rx_c)/2.
        x_dif=nint(rx_dif/dx)
        ry_dif=zondif(idi,2)+(ry_s+ry_c)/2.
        y_dif=nint(ry_dif/dy)
        z_dif=zondif(idi,3)+(z_s+z_c)/2.
        id=nint(rx_dif/dx)
        if (id.gt.width) id=width
        if (id.lt.1) id=1
        jd=nint(ry_dif/dy)
        if (jd.gt.width) jd=width
        if (jd.lt.1) jd=1
        if (z_dif-siz/2..le.altsol(id,jd).or.(z_dif.gt.35000.).or.
     +  (z_dif.gt.cloudbase)) then                                        ! beginning diffusing cell underground
          ndi=ndi+1
        else
          ds1=sqrt((rx_sr-rx_dif)**2.+(ry_sr-ry_dif)**2.+
     +    (z_sr-z_dif)**2.)
          ds2=sqrt((rx_c-rx_dif)**2.+(ry_c-ry_dif)**2.+
     +    (z_c-z_dif)**2.)
          ds3=sqrt((rx_s-rx_dif)**2.+(ry_s-ry_dif)**2.+
     +    (z_s-z_dif)**2.)
          if ((ds1.lt.dss).or.(ds2.lt.dss).or.(ds3.lt.dss)) then
            nss=nss+1
c       print*,ds1,ds2,ds3,'c',rx_c,ry_c,z_c,'d',rx_dif,ry_dif,z_dif
          else
            dho=sqrt((rx_dif-rx_sr)**2.+(ry_dif-ry_sr)**2.)
c computation of the zenithal angle between the reflection surface and the scattering voxel
c shadow reflection surface-scattering voxel
            call anglezenithal(rx_sr,ry_sr,
     +      z_sr,rx_dif,ry_dif,z_dif,angzen)                              ! computation of the zenithal angle reflection surface - scattering voxel.
            call angleazimutal(rx_sr,ry_sr,rx_dif,ry_dif,angazi)          ! computation of the angle azimutal line of sight-scattering voxel
c horizon blocking not a matte because dif are closeby and some downward
            hh=1.
c sub-grid obstacles
            ff=0.
            if (dho.gt.drefle(x_sr,y_sr)) then                            ! light path to observer larger than the mean free path -> subgrid obstacles
              angmin=pi/2.-atan2(obsH(x_sr,y_sr),drefle(x_sr,y_sr))
              if (angzen.lt.angmin) then                                  ! condition obstacle reflechi->scattered
                ff=0.
              else
                ff=ofill(x_sr,y_sr)
              endif
            endif
c computation of the transmittance between the reflection surface and the scattering voxel
            distd=sqrt((rx_dif-rx_sr)**2.+(ry_dif-ry_sr)**2.+
     +      (z_dif-z_sr)**2.)
            call transmitm(angzen,z_sr,z_dif,distd,transm,tranam,tabs)
            call transmita(angzen,z_sr,z_dif,distd,haer,transa,tranaa)
            call transmitl(angzen,z_sr,z_dif,
     +      distd,hlay,transl,tranal)
c computation of the solid angle of the scattering voxel seen from the reflecting surface
            omega=1./distd**2.
            if (omega.gt.omemax) omega=0.
c computing flux reaching the scattering voxel
            fldif2=irefl1*omega*transm*transa*transl*(1.-ff)*hh
c computing the scattering probability toward the line of sight voxel
c cell unitaire
            if (omega.ne.0.) then
              call angle3points (rx_sr,ry_sr,z_sr,rx_dif,ry_dif,z_dif,    ! scattering angle.
     +        rx_c,ry_c,z_c,angdif)
              call diffusion(angdif,tranam,tranaa,tranal,un,secdif,       ! scattering probability of the direct light.
     +        secdil,fdifan,fdifl,haer,hlay,pdifd1,z_dif)
            else
              pdifd1=0.
            endif
            volu=siz**3.
            if (volu.lt.0.) then
              print*,'ERROR, volume 2 is negative!'
              stop
            endif
c computing scattered intensity toward the line of sight voxel from the scattering voxel
            idif2=fldif2*pdifd1*volu
c computing zenith angle between the scattering voxel and the line of sight voxel
            call anglezenithal(rx_dif,ry_dif,z_dif,rx_c,ry_c,z_c,angzen)  ! computation of the zenithal angle between the scattering voxel and the line of sight voxel.
            call angleazimutal(rx_dif,ry_dif,rx_c,ry_c,angazi)            ! computation of the azimutal angle surf refl-scattering voxel
c subgrid obstacles
            if ((x_dif.lt.1).or.(x_dif.gt.nbx).or.(y_dif.lt.1).or.
     +      (y_dif.gt.nbx)) then
              ff=0.
            else
              dho=sqrt((rx_dif-rx_c)**2.+(ry_dif-ry_c)**2.)
              ff=0.
              if (dho.gt.drefle(x_dif,y_dif)) then                        ! light path to observer larger than the mean free path -> subgrid obstacles
                angmin=pi/2.-atan2((obsH(x_dif,y_dif)+
     +          altsol(x_dif,y_dif)-z_dif),drefle(x_dif,y_dif))
                if (angzen.lt.angmin) then                                ! condition subgrid obstacle scattering -> line of sight
                  ff=0.
                else
                  ff=ofill(x_dif,y_dif)
                endif
              endif
            endif
            hh=1.
c computing transmittance between the scattering voxel and the line of sight voxel
            distd=sqrt((rx_dif-rx_c)**2.+(ry_dif-ry_c)**2.+
     +      (z_dif-z_c)**2.)
            call transmitm(angzen,z_dif,z_c,distd,transm,tranam,tabs)
            call transmita(angzen,z_dif,z_c,distd,haer,transa,tranaa)
            call transmitl(angzen,z_dif,z_c,
     +      distd,hlay,transl,tranal)
c computing the solid angle of the line of sight voxel as seen from the scattering voxel
            omega=1./distd**2.
            if (omega.gt.omemax) omega=0.
c computation of the scattered flux reaching the line of sight voxel
            fdif2=idif2*omega*transm*transa*transl*(1.-ff)*hh
c cloud contribution for double scat from a reflecting pixel
            if (cloudt.ne.0) then                                         ! line of sight voxel = cloud
              if (cloudbase-z_c.le.iz*scal) then
                call anglezenithal(rx_c,ry_c,z_c,
     +          rx_obs,ry_obs,z_obs,azcl1)                                ! zenith angle from cloud to observer
                call anglezenithal(rx_c,ry_c,z_c,
     +          rx_dif,ry_dif,z_dif,azcl2)                                ! zenith angle from source to cloud
                doc2=(rx_c-rx_obs)**2.+
     +          (ry_c-ry_obs)**2.+(z_c-z_obs)**2.
                dsc2=(rx_dif-rx_c)**2.+
     +          (ry_dif-ry_c)**2.+(z_dif-z_c)**2.
                call cloudreflectance(angzen,                             ! cloud intensity from direct illum
     +          cloudt,rcloud)
                icloud=icloud+
     +          fldif2/omega*rcloud*doc2*omefov*
     +          abs(cos(azcl2)/cos(azcl1))/dsc2/pi
              endif
            endif
c computation of the scattering probability of the scattered light toward the observer voxel (exiting voxel_c)
            if (omega.ne.0.) then
              call angle3points(rx_dif,ry_dif,z_dif,rx_c,ry_c,z_c,        ! scattering angle.
     +        rx_obs,ry_obs,z_obs,angdif)
              call diffusion(angdif,tranam,tranaa,tranal,un,secdif,       ! scattering probability of the direct light.
     +        secdil,fdifan,fdifl,haer,hlay,pdifd2,z_c)
            else
              pdifd2=0.
            endif
c computing scattered intensity toward the observer from the line of sight voxel
            idif2p=fdif2*pdifd2
            idif2p=idif2p*real(stepdi)*real(ndiff-ndi)/
     +      real(ndiff-ndi-nss)                                           ! Correct the result for the skipping of 2nd scattering voxels to accelerate the calculation
            itotrd=itotrd+idif2p
c ********************************************************************************
c *  section for the calculation of the 2nd scat from the source without reflexion
c ********************************************************************************
            if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then                     ! beginning condition source = reflection for the computation of the source scat line of sight
c computation of the zenithal angle between the source and the scattering voxel
c shadow source-scattering voxel
              call anglezenithal(rx_s,ry_s,z_s,
     +        rx_dif,ry_dif,z_dif,angzen)                                 ! computation of the zenithal angle source-scattering voxel.
              call angleazimutal(rx_s,ry_s,rx_dif,                        ! computation of the angle azimutal line of sight-scattering voxel
     +        ry_dif,angazi)
c horizon blocking not a matter because some path are downward and most of them closeby
              hh=1.
              angmin=pi/2.-atan2((obsH(x_s,y_s)+
     +        altsol(x_s,y_s)-z_s),drefle(x_s,
     +        y_s))
              if (angzen.lt.angmin) then                                  ! condition obstacle source->scattering.
                ff=0.
              else
                ff=ofill(x_s,y_s)
              endif
c computation of the transmittance between the source and the scattering voxel
              distd=sqrt((rx_s-rx_dif)**2.
     +        +(ry_s-ry_dif)**2.
     +        +(z_s-z_dif)**2.)
              call transmitm(angzen,z_s,z_dif,
     +        distd,transm,tranam,tabs)
              call transmita(angzen,z_s,z_dif,
     +        distd,haer,transa,tranaa)
              call transmitl(angzen,z_s,z_dif,
     +        distd,hlay,transl,tranal)
c computation of the Solid angle of the scattering unit voxel seen from the source
              omega=1./distd**2.
              if (omega.gt.omemax) omega=0.
              anglez=nint(180.*angzen/pi)+1
              P_dif1=pvalno(anglez,stype)
c computing flux reaching the scattering voxel
              fldif1=lamplu(x_s,y_s,stype)*P_dif1*
     +        omega*transm*transa*transl*(1.-ff)*hh
c computing the scattering probability toward the line of sight voxel
              if (omega.ne.0.) then
                call angle3points (rx_s,ry_s,z_s,                         ! scattering angle.
     +          rx_dif,ry_dif,z_dif,rx_c,ry_c,z_c,
     +          angdif)
                call diffusion(angdif,                                    ! scattering probability of the direct light.
     +          tranam,tranaa,tranal,un,secdif,secdil,
     +          fdifan,fdifl,haer,hlay,pdifd1,z_dif)
              else
                pdifd1=0.
              endif
              volu=siz**3.
              if (volu.lt.0.) then
                print*,'ERROR, volume 1 is negative!'
                stop
              endif
c computing scattered intensity toward the line of sight voxel from the scattering voxel
              idif1=fldif1*pdifd1*volu
c computing zenith angle between the scattering voxel and the line of sight voxel
              call anglezenithal(rx_dif,ry_dif,
     +        z_dif,rx_c,ry_c,z_c,angzen)                                 ! computation of the zenithal angle between the scattering voxel and the line of sight voxel.
              call angleazimutal(rx_dif,ry_dif,                           ! computation of the azimutal angle surf refl-scattering voxel
     +        rx_c,ry_c,angazi)
c subgrid obstacles
            if ((x_dif.lt.1).or.(x_dif.gt.nbx).or.(y_dif.lt.1).or.
     +      (y_dif.gt.nbx)) then
              dho=sqrt((rx_dif-rx_c)**2.+(ry_dif-ry_c)**2.)
              ff=0.
               else
                 ff=0.
                 if (dho.gt.drefle(x_dif,y_dif)) then
                   angmin=pi/2.-atan2((obsH(x_dif,y_dif)
     +             +altsol(x_dif,y_dif)-z_dif),drefle(
     +             x_dif,y_dif))
                   if (angzen.lt.angmin) then                             ! condition obstacles scattering->line of sight
                     ff=0.
                   else
                     ff=ofill(x_dif,y_dif)
                   endif
                endif
              endif
              hh=1.
c Computing transmittance between the scattering voxel and the line of sight voxel
              distd=sqrt((rx_c-rx_dif)**2.
     +        +(ry_c-ry_dif)**2.
     +        +(z_c-z_dif)**2.)
              call transmitm(angzen,z_dif,z_c,
     +        distd,transm,tranam,tabs)
              call transmita(angzen,z_dif,z_c,
     +        distd,haer,transa,tranaa)
              call transmitl(angzen,z_dif,z_c,
     +        distd,hlay,transl,tranal)
c computing the solid angle of the line of sight voxel as seen from the scattering voxel
              omega=1./distd**2.
              if (omega.gt.omemax) omega=0.
c computation of the scattered flux reaching the line of sight voxel
              fldiff=idif1*omega*transm*transa*transl*(1.-ff)*hh
c cloud contribution to the double scattering from a source
              if (cloudt.ne.0) then                                       ! line of sight voxel = cloud
                if (cloudbase-z_c.le.iz*scal) then
                  call anglezenithal(rx_c,ry_c,z_c,
     +            rx_obs,ry_obs,z_obs,azcl1)                              ! zenith angle from cloud to observer
                  call anglezenithal(rx_c,ry_c,z_c,
     +            rx_dif,ry_dif,z_dif,azcl2)                              ! zenith angle from source to cloud
                  doc2=(rx_c-rx_obs)**2.+
     +            (ry_c-ry_obs)**2.+(z_c-z_obs)**2.
                  dsc2=(rx_dif-rx_c)**2.+
     +            (ry_dif-ry_c)**2.+(z_dif-z_c)**2.
                  call cloudreflectance(angzen,                           ! cloud intensity from direct illum
     +            cloudt,rcloud)
                  icloud=icloud+
     +            fldiff/omega*rcloud*doc2*omefov*
     +            abs(cos(azcl2)/cos(azcl1))/dsc2/pi
                endif
              endif
c computation of the scattering probability of the scattered light toward the observer voxel (exiting voxel_c)
              if (omega.ne.0.) then
                call angle3points(rx_dif,ry_dif,                          ! scattering angle.
     +          z_dif,rx_c,ry_c,z_c,rx_obs,ry_obs,
     +          z_obs,angdif)
                call diffusion(angdif,                                    ! scattering probability of the direct light.
     +          tranam,tranaa,tranal,un,secdif,secdil,
     +          fdifan,fdifl,haer,hlay,pdifd2,z_c)
              else
                pdifd2=0.
              endif
c computing scattered intensity toward the observer from the line of sight voxel
              idiff2=fldiff*pdifd2
              idiff2=idiff2*real(stepdi)*real(ndiff-ndi)/
     +        real(ndiff-ndi-nss)                                         ! Correct the result for the skipping of 2nd scattering voxels to accelerate the calculation
              itodif=itodif+idiff2                                        ! sum over the scattering voxels
            endif                                                         ! end condition source = reflection for the computation of the source scat line of sight
          endif                                                           ! end of the case scattering pos = Source pos or line of sight pos
        endif                                                             ! end diffusing celle underground
      enddo                                                               ! end of the loop over the scattering voxels.
                                        endif                             ! end of the condition ou effdif > 0
c End of 2nd scattered intensity calculations
c===================================================================
c
c
c
c **********************************************************************
c * section refected light with single scattering
c **********************************************************************
c verify if there is shadow between sr and line of sight voxel
                                        call anglezenithal(rx_sr,ry_sr,   ! zenithal angle between the reflecting surface and the line of sight voxel.
     +                                  z_sr,rx_c,ry_c,z_c,angzen)
                                        call angleazimutal(rx_sr,ry_sr,   ! computation of the azimutal angle reflect-line of sight
     +                                  rx_c,ry_c,angazi)
                                        distd=sqrt((rx_sr-rx_c)**2.
     +                                  +(ry_sr-ry_c)**2.
     +                                  +(z_sr-z_c)**2.)
                                        dho=sqrt((rx_sr-rx_c)**2.
     +                                  +(ry_sr-ry_c)**2.)
                                        if (angzen.gt.pi/4.) then         ! 45deg. it is unlikely to have a 1km high mountain less than 1
        call horizon(x_sr,y_sr,z_sr,dx,dy,altsol,angazi,zhoriz,dh)
                                          if (dh.le.dho) then
                                            if (angzen-zhoriz.lt.
     +                                      0.00001) then                 ! the path line of sight-reflec is not below the horizon => we compute
                                              hh=1.
                                            else
                                              hh=0.
                                            endif                         ! end condition reflecting surf. above horizon
                                          else
                                            hh=1.
                                          endif
                                        else
                                          hh=1.
                                        endif
                                        irefl=irefl1
c case: line of sight position = Position of reflecting cell
                                        if((rx_c.eq.rx_sr).and.(ry_c.eq.
     +                                  ry_sr).and.(z_c.eq.z_sr)) then
                                          intind=0.
                                        else
c obstacle
                                         dho=sqrt((rx_sr-rx_c)**2.
     +                                   +(ry_sr-ry_c)**2.)
                                         ff=0.
                                         if (dho.gt.drefle(x_sr,y_sr))
     +                                   then
                                           angmin=pi/2.-atan2(obsH(x_sr
     +                                     ,y_sr),drefle(x_sr,y_sr))
                                           if (angzen.lt.angmin) then     ! condition obstacle reflected.
                                             ff=0.
                                           else
                                             ff=ofill(x_sr,y_sr)
                                           endif
                                         endif
c computation of the transmittance between the ground surface and the line of sight voxel
                                          call transmitm(angzen,z_sr,
     +                                    z_c,distd,transm,tranam,tabs)
                                          call transmita(angzen,z_sr,
     +                                    z_c,distd,haer,transa,tranaa)
                                          call transmitl(angzen,z_sr,
     +                                    z_c,distd,hlay,transl,tranal)
c computation of the solid angle of the line of sight voxel seen from the reflecting cell
                                          omega=1./distd**2.
                                          if (omega.gt.omemax) omega=0.
c computation of the flux reflected reaching the line of sight voxel
                                          flindi=irefl*omega*transm*
     +                                    transa*transl*(1.-ff)*hh        ! obstacles correction
c cloud contribution to the reflected light from a ground pixel
                              if (cloudt.ne.0) then                       ! line of sight voxel = cloud
                                if (cloudbase-z_c.le.iz*scal) then
                                  call anglezenithal(rx_c,ry_c,z_c,
     +                            rx_obs,ry_obs,z_obs,azcl1)              ! zenith angle from cloud to observer
                                  call anglezenithal(rx_c,ry_c,z_c,
     +                            rx_sr,ry_sr,z_sr,azcl2)                 ! zenith angle from source to cloud
                                  doc2=(rx_c-rx_obs)**2.+
     +                            (ry_c-ry_obs)**2.+(z_c-z_obs)**2.
                                  dsc2=(rx_sr-rx_c)**2.+
     +                            (ry_sr-ry_c)**2.+(z_sr-z_c)**2.
                                  call cloudreflectance(angzen,           ! cloud intensity from direct illum
     +                            cloudt,rcloud)
                                  icloud=icloud+
     +                            flindi/omega*rcloud*doc2*omefov*
     +                            abs(cos(azcl2)/cos(azcl1))/dsc2/pi
                                endif
                              endif
c computation of the scattering probability of the reflected light
                                          if (omega.ne.0.) then
                                            call angle3points(rx_sr,      ! scattering angle.
     +                                      ry_sr,z_sr,rx_c,ry_c,z_c,
     +                                      rx_obs,ry_obs,z_obs,angdif)
                                            call diffusion(angdif,        ! scattering probability of the reflected light.
     +                                      tranam,tranaa,tranal,un,
     +                                      secdif,secdil,fdifan,fdifl,
     +                                      haer,hlay,pdifin,z_c)
                                          else
                                            pdifin=0.
                                          endif
c computation of the reflected intensity toward the sensor by a reflecting cell
                                          intind=flindi*pdifin*(1.-ff)
     +                                    *hh
                                        endif                             ! end of the case Posi reflecting cell =  line of sight voxel position
                                        itotind=itotind+intind            ! Sum of the intensities of each reflecting cell.
                                      endif                               ! end of the condition surface not lighted from the top.
                                  endif                                   ! end of the condition reflecting cell is not on the source.
                                endif                                     ! end of the condition surface of the domain.
                              enddo                                       ! end of the loop over the rows (latitu) reflecting.
                            enddo                                         ! end of the loop over the column (longitude) reflecting.
c   end of the computation of the reflected intensity
c
c**********************************************************************
c computation of the total intensity coming from a source to the line of sight voxel toward the sensor
c**********************************************************************
c In the order 1st scat; refl->1st scat; 1st scat->2nd scat,
c refl->1st scat->2nd scat
                            isourc=intdir+itotind+itodif+itotrd           ! Sum of the intensities of a given type of source reaching the line of sight voxel.
                            isourc=isourc*scal                            ! scaling the values according to the path length in the l. of sight voxel of 1m3
                            isourc=isourc*portio                          ! correct for the field of view of the observer
c include clouds in the total intensity
                            isourc=isourc+icloud


        if ((itodif.lt.0.).or.(itotrd.lt.0.)) then
          print*,intdir,itotind,itodif,itotrd
          stop
        endif



                            if (verbose.eq.2) then
       print*,' Total intensity per component for type ',ntype,':'
       print*,' source->scattering=',intdir
       print*,' source->reflexion->scattering=',itotind
       print*,' source->scattering->scattering=',itodif
       print*,' source->reflexion->scattering->scattering=',itotrd
       if (intdir*itotind*itodif*itotrd.lt.0.) then
         print*,'PROBLEM! Negative intensity.'
         stop
       endif
                            endif
c**********************************************************************
c computation of the total intensity coming from all the sources of a given type
c**********************************************************************
                            itotty=itotty+isourc                          ! Sum of the intensities all sources of the same typeand a given line of sight element
                            ITT(x_s,y_s,stype)=ITT(x_s,y_s,stype)+isourc  ! ITT stores itotty in a matrix
                        endif                                             ! end of the condition "the luminosity of the ground pixel x_s,y_s in not null".
                      enddo                                               ! end the loop over the lines (latitude) of the domain (y_s).
                    enddo                                                 ! end the loop over the column (longitude) of the domain (x_s).
c end of the computation of the intensity of one source type
                    itotci=itotci+itotty                                  ! Sum of the intensities all source all type to a line of sight element
                    do x_s=imin(stype),imax(stype)
                      do y_s=jmin(stype),jmax(stype)
                        ITC(x_s,y_s)=ITC(x_s,y_s)+ITT(x_s,y_s,stype)
                      enddo
                    enddo
c calculate total lamp flux matrix for all lamp types
                    do x_s=1,nbx
                      do y_s=1,nby
                        lpluto(x_s,y_s)=lpluto(x_s,y_s)+
     +                  lamplu(x_s,y_s,stype)
                      enddo
                    enddo
                  endif                                                   ! end of condition if there are any flux in that source type
              enddo                                                       ! end of the loop over the types of sources (stype).
c end of the computation of the intensity coming from a line of sight voxel toward the sensor
c
c
c***********************************************************************
c computation of the luminous flux reaching the observer
c***********************************************************************
c computation of the zenithal angle between the observer and the line of sight voxel
c=======================================================================
                call anglezenithal(rx_c,ry_c,z_c,rx_obs,ry_obs,z_obs,
     +          angzen)                                                   ! computation of the zenithal angle between the line of sight voxel and the observer.
c                                                                         ! end of the case "observer at the same latitu/longitude than the source".
c computation of the transmittance between the line of sight voxel and the observer
                                    distd=sqrt((rx_c-rx_obs)**2.
     +                              +(ry_c-ry_obs)**2.
     +                              +(z_c-z_obs)**2.)
                call transmitm(angzen,z_c,z_obs,distd,transm,tranam,
     +          tabs)
                call transmita(angzen,z_c,z_obs,distd,haer,transa,
     +          tranaa)
                call transmitl(angzen,z_c,z_obs,
     +          distd,hlay,transl,tranal)
c computation of the flux reaching the objective of the telescope from the line of sight voxel
                fcapt=itotci*ometif*transa*transm*transl                         ! computation of the flux reaching the intrument from the line of sight voxel
                do x_s=1,nbx
                  do y_s=1,nby
                    FCA(x_s,y_s)=ITC(x_s,y_s)*ometif*transa*transm*
     +              transl
                  enddo
                enddo
                if (cos(pi-angzen).eq.0.) then
                  print*,'ERROR perfectly horizontal sight is forbidden'
                  stop
                endif
c end of the computation of the flux reaching the observer voxel from the line of sight voxel
                ftocap=ftocap+fcapt                                       ! flux for all source all type all line of sight element
                do x_s=1,nbx
                  do y_s=1,nby
                    FTC(x_s,y_s)=FTC(x_s,y_s)+FCA(x_s,y_s)                ! FTC is the array of the flux total at the sensor to identify
                                                                          ! the contribution of each ground pixel to the total flux at the observer level
                                                                          ! The % is simply given by the ratio FTC/ftocap
                  enddo
                enddo
c correction for the FOV to the flux reaching the intrument from the cloud voxel
            if (cloudt.ne.0) then
c computation of the flux reaching the intrument from the cloud voxel
                fccld=icloud*ometif*transa*transm*transl
                fctcld=fctcld+fccld                                       ! cloud flux for all source all type all line of sight element
            endif
            if (verbose.ge.1) print*,'Added radiance =',
     +      fcapt/omefov/(pi*(diamobj/2.)**2.)
            if (verbose.ge.1) print*,'Radiance accumulated =',
     +      ftocap/omefov/(pi*(diamobj/2.)**2.)
            if (verbose.ge.1) write(2,*) 'Added radiance =',
     +      fcapt/omefov/(pi*(diamobj/2.)**2.)
            if (verbose.ge.1) write(2,*) 'Radiance accumulated =',
     +      ftocap/omefov/(pi*(diamobj/2.)**2.)
              endif                                                       ! end of the condition line of sight voxel inside the modelling domain
          endif                                                           ! end condition line of sight voxel 1/stoplim
c accelerate the computation as we get away from the sources
          scalo=scal
          if (scal.le.3000.)  scal=scal*1.12
          endif
        else
c           print*,'End of line of sight - touching the ground'
        endif                                                             ! line of sight not blocked by topography
        enddo                                                             ! end of the loop over the line of sight voxels.
        fctcld=fctcld*10**(0.4*(100.-cloudfrac)*cloudslope)               ! correction for the cloud fraction (defined from 0 to 100)
        if (prmaps.eq.1) then
c          open(unit=9,file=pclf,status='unknown')
            do x_s=1,nbx
              do y_s=1,nby
                FTC(x_s,y_s)=FTC(x_s,y_s)/ftocap                          ! Here FTC becomes the flux fraction of each pixel. The sum of FTC values over all pixels give the total flux
              enddo
            enddo
            if (verbose.eq.2) then
              print*,'Writing normalized contribution array'
              print*,'Warning Cloud contrib. excluded from that array.'
            endif
c            do x_s=1,nbx
c              do y_s=1,nby
c                write(9,*) x_s,y_s,FTC(x_s,y_s)                           ! emettrice au sol, c'est un % par unite of watt installes
c              enddo
c            enddo
            call twodout(nbx,nby,pclimg,FTC)
c          close(unit=9)
c creation of gnuplot file. To visualize, type gnuplot and then
c load 'BASENAME_pcl.gplot'
c          open(unit=9,file=pclgp,status='unknown')
c            write(9,*) 'sand dgrid3d',nbx,',',nby
c            write(9,*) 'sand hidden3d'
c            write(9,*) 'sand pm3d'
c            write(9,*) 'splot "'//basenm(1:lenbase)//'_pcl.txt"
c     +      with dots'
c          close(unit=9)
        endif                                                             ! end of condition for creating contrib and sensit maps
        endif                                                             ! end of scattered light
c
c End of calculation of the scattered radiances
c =================================

        if (verbose.ge.1) print*,'======================================
     +==============='
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
        print*,'            Diffuse radiance (W/str/m**2/nm)'
        write(*,2001) (ftocap+fctcld)/omefov/(pi*(diamobj/2.)**2.)
        if (verbose.ge.1) write(2,*) '==================================
     +================='
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
        write(2,*) '         Diffuse radiance (W/str/m**2/nm)          '
        write(2,2001) (ftocap+fctcld)/omefov/(pi*(diamobj/2.)**2.)
      close(2)
 2001 format('                   ',E10.3E2)
      stop
      end
c***********************************************************************************************************************
c*                                                                                                                     *
c*                                         end of the programme                                                        *
c*                                                                                                                     *
c***********************************************************************************************************************
