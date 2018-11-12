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
c ** Illumina en Fortran 77                                                                                           **
c ** Programmers in decreasing order of contribution  :                                                               **
c **                            Martin Aube, Loic Franchomme-Fosse,  Mathieu Provencher, Andre Morin                  **
c **                            Alex Neron, Etienne Rousseau                                                          ** 
c **                            William Desroches, Maxime Girardin, Tom Neron                                         **
c **                                                                                                                  **
c ** Illumina can be downloaded via:   hg clone  https://aubema@bitbucket.org/aubema/illumina                         **
c ** To compile:                                                                                                      **
c **    cd hg/illumina                                                                                                **
c **    mkdir bin                                                                                                     **
c **    bash makeILLUMINA                                                                                             **
c **                                                                                                                  **
c **  Current version features/limitations :                                                                          **
c **                                                                                                                  **
c **    - Calculation of flux entering a spectrometer in a given line of sight                                        **
c **    - Calculation of the sky spectral luminance in a given line of sight                                          **
c **    - Calculation of the atmospheric transmittance and 1st and 2nd order of scattering                            **
c **    - Lambertian reflexion on the ground                                                                          **
c **    - Terrain slope considered (apparent surface and shadows)                                                     **
c **    - Angular photometry of a lamp is considered uniform along the azimuth                                        **
c **    - Sub-grid obstacles considered (with the mean free path of light toward ground and mean obstacle height      **
c **    - Molecules and aerosol optics (phase function, scattering probability, aerosol absorption)                   **  
c **    - Exponential concentrations vertical profile (H aerosol= 2km, H molecules= 8km  )                            **
c **    - Exponential vertical resolution                                                                             **
c **    - Accounting for heterogeneity of ground reflectance, luminaires number, luminaires heights,                  **
c **      angular photometry                                                                                          **
c **    - Wavelength dependant                                                                                        **
c **    - Clouds models                                                                                               **
c **    - Ignore the flux scattered by the voxel occupied by the observer (voxelobs=voxelcible)                       **
c **    - Do not support direct observation of a source                                                               ** 
c **    - Direct observation of the ground not implemented                                                            **
c **    - Not accounting for molecular absorption                                                                     **
c **    - Do not consider earth curvature (i.e. local/regional model)                                                 **
c **    - No clouds                                                                                                   **
c **                                                                                                                  **
c ** Theoretical equations by Martin Aube, CEGEP of Sherbrooke (in french)                                            **
c **      http://cegepsherbrooke.qc.ca/~aubema/index.php/Prof/IllumEn?action=download&upname=intensity_lumineuse.pdf  **
c **                                                                                                                  **
c **********************************************************************************************************************
c   
c    Copyright (C) 2017 Martin Aube
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
c2345678901234567890123456789012345678901234567890123456789012345678901234
c
      program illumina                                                    ! Beginning
c
c=======================================================================
c     Variables declaration
c=======================================================================
c
      integer width,height                                                ! Matrix dimension in Length/width and height
      parameter (width=1024,height=100)
      integer iun,ideux
      real pi,pix4
      integer verbose                                                     ! verbose = 1 to have more print out, 0 for silent
      parameter (pi=3.1415926)
      parameter (pix4=4.*pi)
      real cthick(height)                                                 ! voxel thickness array (meter)
      real cellh(height)                                                  ! voxel height array (meter)
      real flcumu                                                         ! Accrued flux along the line of sight
      character*72 mnaf                                                   ! Terrain elevation file.
      character*72 reflf                                                  ! Reflectance file.
      character*72 diffil                                                 ! Aerosol file.
      character*72 outfile                                                ! Results file
      character*72 pclf,pcwf,pclgp,pcwgp                                  ! File containing contribution and sensitivity maps
      character*72 pclimg,pcwimg                                                                         
      character*72 basenm                                                 ! Base name of files.
      integer lenbase                                                     ! Length of the Base name of the experiment.
      real lambda,pressi,drefle(width,width)                              ! Wavelength (nanometer), atmospheric pressure (kPa), mean free path to the ground (meter).
      integer ntype                                                       ! Number of light source types or zones considered.
      real largx                                                          ! Width (x axis) of the modeling domain (meter).
      real largy                                                          ! Length (y axis) of the modeling domain (meter).
      integer nbx,nby                                                     ! Number of pixels in the modeling domain.  
      real val2d(width,width)                                             ! Temporary input array 2d
      real altsol(width,width)                                            ! Ground elevation (meter).
      real srefl(width,width)                                             ! Ground reflectance.
      real Hmin                                                           ! Minimum ground elevation of the modeling domain
      integer stype                                                       ! Source type or zone index
      character*72 pafile,lufile,alfile,ohfile,odfile,offile              ! Files related to light sources and obstacles (photometric function of the sources (sr-1), flux (W), height (m), obstacle height (m), obstacle distance (m), obstacle filling factor (0-1).    
      real lamplu(width,width,120)                                        ! Source fluxes
      real lampal(width,width)                                            ! Height of the light sources relative to the ground (meter).
      real pval(181,120),pvalto,pvalno(181,120)                           ! Values of the angular photometry functions (unnormalized, integral, normalized).
      real dtheta                                                         ! Angle increment of the photometric function of the sources 
      real dx,dy,dxp,dyp                                                  ! Width of the voxel (meter)
c      real pixsiz                                                         ! pixel size (meter)
      integer boxx,boxy                                                   ! reflection window size (pixels).
      real fdifa(181),fdifan(181)                                         ! Aerosol scattering functions (unnormalized and normalized).
      real extinc,scatte,anglea(181)                                      ! Aerosol cross sections (extinction and scattering), scattering angle (degree).
      real secdif                                                         ! Contribution of the scattering to the extinction
      real inclix(width,width)                                            ! tilt of the ground pixel along x (radian).
      real incliy(width,width)                                            ! tilt of the ground pixel along y (radian).   
      integer x_obs,y_obs,zcello                                          ! Position of the observer (INTEGER).
      real z_o                                                            ! observer height relative to the ground
      real z_obs                                                          ! Height of the observer (meter) to the vertical grid scale.
      integer lcible(width,3)                                             ! Array for the line of sight voxels along the line of sight.
      integer ncible,icible                                               ! Number of line of sight voxels, number loops over the voxels 
      integer x_c,y_c,zcellc                                              ! Position of the line of sight voxel (INTEGER).
      real z_c                                                            ! Height of the line of sight voxel (metre).
      real zcup,zcdown                                                    ! Lower and upper limits of the line of sight voxel.    
      integer dirck                                                       ! Test for the position of the source (case source=line of sight voxel).     
      integer x_s,y_s,x_sr,y_sr,x_dif,y_dif,zceldi                        ! Positions of the source, the reflecting surface, and the scattering voxels 
      real z_s,z_sr,z_dif                                                 ! Heights of the source, the reflecting surface, and the scattering voxel (metre).
      real angzen,ouvang                                                  ! Zenithal angle between two voxels (radians) and opening angle of the solid angle in degrees.
      integer anglez                                                      ! Emitting zenithal angle from the luminaire.      
      real P_dir,P_indir,P_dif1                                           ! photometric function of the light sources (direct,indirect,scattered) 
      real transa,transm                                                  ! Transmittance between two voxels (aerosols,molecules).
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
      real irefdi                                                         ! intensity from a voxel toward the sensor after reflexion and double scattering. 
      real flcib                                                          ! Flux reaching the observer voxel from a line of sight voxel.
      real fcapt                                                          ! Flux reaching the observer voxel from all FOV voxels in a given model level
      real ftocap                                                         ! Total flux reaching the observer voxel
      real haut                                                           ! Haut (negative indicate that the surface is lighted from inside the ground. I.e. not considered in the calculation
      real epsilx,epsily                                                  ! tilt of the ground pixel
      real flrefl                                                         ! flux reaching a reflecting surface (watts).
      real irefl,irefl1                                                   ! intensity leaving a reflecting surface toward the line of sight voxel.  
      real effdif                                                         ! Distance around the source voxel and line of sight voxel considered to compute the 2nd order of scattering.
      integer zondif(3000000,4)                                           ! Array for the scattering voxels, the 4th column represents the nearest integer value of the distance (en metre) to the line of single scattering.
      integer ndiff,idi                                                   ! Number of scattering voxels, counter of the loop over the scattering voxels
      integer stepdi                                                      ! scattering step to speedup the calculation e.g. if =2 one computation over two will be done
      integer nvis0                                                       ! starting value for the calculation along of the viewing line. 
c                                                                         ! by default the value is 1 but it can be larger  
c                                                                         ! when we resume a previous interrupted calculation.
      real fldif1                                                         ! flux reaching a scattering voxel.
      real idif1                                                          ! intensity toward a line of sight voxel from a scattering voxel.
      real portio                                                         ! ratio of voxel surface to the solid angle of the sensor field of view.
      real dis_obs                                                        ! Distance between the line of sight and the observer.
      real ometif                                                         ! Solid angle of the telescope objective as seen from the line of sight voxel
      real omefov                                                         ! Solid angle of the spectrometer slit.
      real lfente                                                         ! Width of the slit (or of the sensor area)
      real longfe                                                         ! Length of the slit (or of the sensor area)
      real focal                                                          ! Focal distance of the spectrometer objective  
      real angvis,azim                                                    ! viewing angles of the sensor.
      real projap                                                         ! Fraction of the reflecting surface relatively to the normal. 
c                                                                         ! Useful for the calculation of the lambertian reflectance.
      real nbang                                                          ! for the averaging of the photometric function
      real obsH(width,width),angmin                                       ! averaged height of the sub-grid obstacles, minimum angle under wich 
c                                                                         ! a light ray cannot propagate because it is blocked by a sub-grid obstable
      real ofill(width,width)                                             ! fill factor giving the probability to hit an obstacle when pointing in its direction integer 0-100

      integer naz,na 
      real ITT(width,width,120)                                           ! total intensity per type of lamp
      real ITC(width,width)                                               ! total intensity per line of sight voxel
      real FC(width,width)                                                ! line of sight flux
      real FTC(width,width)                                               ! fraction of the total flux at the sensor level 
      real FTCN(width,width)                                              ! fraction of the total flux at the sensor level normalized per unit of watt
      real FCA(width,width)                                               ! sensor flux array
      real lpluto(width,width)                                            ! total luminosity of the ground cell for all lamps
      real fctnto,ftcmax                                                  ! FTCN total for all the domain for all lamps
      character*3 lampno                                                  ! lamp number string
      integer imin(120),imax(120),jmin(120),jmax(120)                     ! x and y limits of the zone containing a type of lamp
c      real zhoriz                                                         ! horizon in rad over 360 deg, the first index of the array is for 0 deg while index 360 = 359 deg
      real angazi                                                         ! azimuth angle between two points in rad, max dist for the horizon determination
      real latitu                                                         ! approximate latitude of the domain center
      integer prmaps                                                      ! flag to enable the tracking of contribution and sensitivity maps
      integer cloudt                                                      ! cloud type 0=clear, 1=Thin Cirrus/Cirrostratus, 2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus, 4=Cumulus/Cumulonimbus, 5=Stratocumulus
      integer cloudh(5),cloudz                                            ! cloud base layer relative to the lower elevation 
      integer xsrmi,xsrma,ysrmi,ysrma                                     ! limits of the loop valeur for the reflecting surfaces
      real rcloud                                                         ! cloud relfectance 
      real azencl                                                         ! zenith angle from cloud to observer
      real icloud                                                         ! cloud reflected intensity
      real fcloud                                                         ! flux reaching the intrument from the cloud voxel
      real fccld                                                          ! correction for the FOV to the flux reaching the intrument from the cloud voxel
      real fctcld                                                         ! total flux from cloud at the sensor level
      real dsco                                                           ! distancesource-line of sight-observer
      real dminlp                                                         ! minimum distance between the observer and a lamp (m)
      real totlu(120)                                                     ! total flux of a source type
      real stoplim                                                        ! Stop computation when the new voxel contribution is less than 1/stoplim of the cumulated flux
      real zero
      real anaz
      real ff,hh                                                          ! temporary obstacle filling factor and horizon blocking factor
      data cloudh /73,73,65,54,53/                                        ! 9300.,9300.,4000.,1200.,1100.
      real dist,distm                                                     ! distance and minimal distance to find observer level
      real scalef                                                         ! scale factor to calculate the line of sight
      verbose=0
      zero=0.
      ff=0.
      print*,'Starting ILLUMINA computations...'
c
c determining the vertical scale
c
      call verticalscale(cthick,cellh)
c  
c=======================================================================
c        reading of the fichier d'entree (illumina.in)
c=======================================================================
      print*,'Reading illumina.in input file'
      open(unit=1,file='illumina.in',status='old')
       read(1,*)
       read(1,*) basenm
       read(1,*) nbx,nby
       read(1,*) dx,dy
       read(1,*) latitu
       read(1,*) diffil
       read(1,*) 
       read(1,*) effdif,stepdi
       if (verbose.eq.1) then
         print*,'2nd order scattering radius=',effdif,'m   1 voxel over
     a    ',stepdi
       endif
       read(1,*)
       read(1,*) lambda
       read(1,*) pressi
       read(1,*) taua,alpha
       read(1,*) ntype
       read(1,*) stoplim
       read(1,*)
       read(1,*) x_obs,y_obs,z_o,nvis0
       read(1,*)
       read(1,*) angvis,azim
       read(1,*) 
       read(1,*) lfente,longfe,focal,diamobj 
       read(1,*)
       read(1,*)
       read(1,*) cloudt  
       read(1,*) dminlp 
      close(1)
c
c computing the actual AOD at the wavelength lambda
c      
       print*,'500nm AOD=',taua,'500nm angstrom coeff.=',alpha
       taua=taua*(lambda/500.)**(-1.*alpha)
c
c  determine the Length of basenm
c 
      lenbase=index(basenm,' ')-1  
      mnaf=basenm(1:lenbase)//'_topogra.bin'                              ! determine the names of input and output files
      reflf=basenm(1:lenbase)//'_reflect.bin' 
      outfile=basenm(1:lenbase)//'.out'  
      pclf=basenm(1:lenbase)//'_pcl.txt'
      pcwf=basenm(1:lenbase)//'_pcw.txt'
      pclimg=basenm(1:lenbase)//'_pcl.bin'
      pcwimg=basenm(1:lenbase)//'_pcw.bin'
      pclgp=basenm(1:lenbase)//'_pcl.gplot'
      pcwgp=basenm(1:lenbase)//'_pcw.gplot'    
c  conversion of the geographical viewing angles toward the cartesian 
c  angle we assume that the angle in the file illumina.in
c  is consistent with the geographical definition 
c  geographical, azim=0 toward north, 90 toward east, 180 toward south 
c  etc
c  cartesian, azim=0 toward east, 90 toward north, 180 toward west etc
      azim=90.-azim
      if (azim.lt.0.) azim=azim+360.
      if (azim.ge.360.) azim=azim-360.
c  opening output file
      open(unit=2,file=outfile,status='unknown')      
       write(2,*) 'FILE USED:'
       write(2,*) mnaf,reflf,diffil
       print*,'Wavelength (nm):',lambda,
     +       ' Aerosol optical depth:',taua
       write(2,*) 'Wavelength (nm):',lambda,
     +       ' Aerosol optical depth:',taua
       write(2,*) '2nd order scattering radius:',effdif,' m'
       print*,'2nd order scattering radius:',effdif,' m'
       write(2,*) 'Scattering step:',stepdi
       print*,'Scattering step:',stepdi

       write(2,*) 'Observer position (x,y,z)',x_obs,y_obs,z_o
       print*,'Observer position (x,y,z)',x_obs,y_obs,z_o
       write(2,*) 'Elevation angle:',angvis,' azim angle (clockwise fro
     +m north)',azim     
       print*,'Elevation angle:',angvis,' azim angle (counterclockwise f
     +rom east)',azim 
c=======================================================================
c        Initialisation of the arrays and variables
c=======================================================================
       print*,'Initializing variables...'
       if (cloudt.eq.0) then
          cloudz=height
       else
          cloudz=cloudh(cloudt)
       endif
       prmaps=1
       iun=0
       ideux=1
       flcumu=0.
       icloud=0.
       do i=1,width
        do j=1,width
         val2d(i,j)=0.
         altsol(i,j)=0.
         srefl(i,j)=0.
         inclix(i,j)=0.
         incliy(i,j)=0.
         lpluto(i,j)=0.
         ITC(i,j)=0.
         FC(i,j)=0.
         FTC(i,j)=0.
         FTCN(i,j)=0.
         FCA(i,j)=0.
         do k=1,120
          lamplu(i,j,k)=0.
          lampal(i,j)=0.
          ITT(i,j,k)=0.
         enddo
        enddo
       enddo
       do i=1,181
        fdifa(i)=0.
        fdifan(i)=0.
        anglea(i)=0.
        do j=1,120
         pval(i,j)=0.
         pvalno(i,j)=0.
        enddo
       enddo  
       do i=1,width
        do j=1,3
         lcible(i,j)=1
        enddo
       enddo
       do i=1,3000000
        do j=1,4
         zondif(i,j)=1
        enddo
       enddo     
       irefdi=0.
       idif1=0.
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
       portio=0.
       fctnto=0.
       ftcmax=0.
       fccld=0.
       fctcld=0.
       ometif=0.
       omefov=0.
       hh=1.
c***********************************************************************
c        reading of the environment variables                          *
c***********************************************************************
c=======================================================================
c  reading of the elevation file
c=======================================================================
       call 2din(nbx,nby,mnaf,altsol)

       Hmin=3000000.
       do i=1,nbx                                                         ! beginning of the loop over all cells along x.
        do j=1,nby                                                        ! beginning of the loop over all cells along y.
c                                                                         ! searching of the Height minimale.
         if (Hmin.gt.altsol(i,j)) Hmin=altsol(i,j)
        enddo                                                             ! end of the loop over all cells along y.
       enddo 
       do i=1,nbx                                                         ! beginning of the loop over all cells along x.
        do j=1,nby                                                        ! beginning of the loop over all cells along y.
         altsol(i,j)=altsol(i,j)-Hmin                                     ! subtraction of the Minimum ground elevation
        enddo                                                             ! end of the loop over all cells along y.
       enddo
c=======================================================================
c reading reflectance file
c=======================================================================
       call 2din(nbx,nby,reflf,srefl)
       do i=1,nbx                                                         ! beginning of the loop over all cells along x.
        do j=1,nby                                                        ! beginning of the loop over all cells along y.
         if (srefl(i,j).lt.0.) then                                       ! searching of of the negative reflectances
           print*,'***,WARNING - Negative reflectance replacing by 0.!'
           srefl(i,j)=0.
         endif
        enddo                                                             ! end of the loop over all cells along y.
       enddo
c=======================================================================
c  reading of the values of P(theta), height, luminosities and positions 
c   of the sources, obstacle height and distance
c=======================================================================
c
       ohfile=basenm(1:lenbase)//'_obsth.bin'
       odfile=basenm(1:lenbase)//'_obstd.bin'
       alfile=basenm(1:lenbase)//'_altlp.bin'                             ! setting the file name of height of the sources lumineuse.
       offile=basenm(1:lenbase)//'_obstf.bin'
       dtheta=.017453293                                                  ! one degree
       do stype=1,ntype                                                   ! beginning of the loop 1 for the 120 types of sources.
        imin(stype)=nbx
        jmin(stype)=nby
        imax(stype)=1
        jmax(stype)=1       
        pvalto=0.
        write(lampno, '(I3.3)' ) stype                                    ! support of 120 different sources (3 digits)
        pafile=basenm(1:lenbase)//'_fctem_'//lampno//'.dat'               ! setting the file name of angular photometry.
        lufile=basenm(1:lenbase)//'_lumlp_'//lampno//'.bin'               ! setting the file name of the luminosite of the cases.
c    ===================================================================
c    reading photometry files
        open(UNIT=1, FILE=pafile,status='OLD')                            ! opening file pa#.dat, angular photometry.
        do i=1,181                                                        ! beginning of the loop for the 181 data points
         read(1,*) pval(i,stype)                                          ! reading of the data in the array pval.
         pvalto=pvalto+pval(i,stype)*2.*pi*                               ! Sum of the values of the  photometric function 
     a   sin(real(i-1)*dtheta)*dtheta                                     ! (pvaleur x 2pi x sin theta x dtheta) (ou theta egale 
c                                                                         ! (i-1) x 1 degrees).
        enddo                                                             ! end of the loop over the 181 donnees of the fichier pa#.dat.
        close(1)                                                          ! closing file pa#.dat, angular photometry.
        do i=1,181
         if (pvalto.ne.0.) pvalno(i,stype)=pval(i,stype)/pvalto           ! Normalisation of the photometric function.
        enddo   
c    ===================================================================
c    reading luminosity files
       call 2din(nbx,nby,lufile,val2d)
       do i=1,nbx                                                         ! beginning of the loop over all cells along x.
        do j=1,nby                                                        ! beginning of the loop over all cells along y.
         if (val2d(i,j).lt.0.) then                                       ! searching of negative fluxes
           print*,'***Negative lamp flux!, stopping execution'
           stop
         endif
        enddo                                                             ! end of the loop over all cells along y.
       enddo     
        do i=1,nbx                                                        ! searching of the smallest rectangle containing the zone
         do j=1,nby                                                       ! of non-null luminosity to speedup the calculation
          if (val2d(i,j).ne.0.) then
           if (i-1.lt.imin(stype)) imin(stype)=i-2
           if (imin(stype).lt.1) imin(stype)=1
           goto 333
          endif
         enddo 
        enddo
        imin(stype)=1   
 333    do i=nbx,1,-1
         do j=1,nby
          if (val2d(i,j).ne.0.) then
           if (i+1.gt.imax(stype)) imax(stype)=i+2    
           if (imax(stype).gt.nbx) imax(stype)=nbx
           goto 334
          endif
         enddo
        enddo
        imax(stype)=1
 334    do j=1,nby
         do i=1,nbx
          if (val2d(i,j).ne.0.) then
           if (j-1.lt.jmin(stype)) jmin(stype)=j-2 
           if (jmin(stype).lt.1) jmin(stype)=1
           goto 335
          endif
         enddo
        enddo 
        jmin(stype)=1
 335    do j=nby,1,-1
         do i=1,nbx
          if (val2d(i,j).ne.0.) then
           if (j+1.gt.jmax(stype)) jmax(stype)=j+2
           if (jmax(stype).gt.nby) jmax(stype)=nby
           goto 336
          endif
         enddo
        enddo  
        jmax(stype)=1
 336    do i=1,nbx                                                        ! beginning of the loop over all cells along x.
         do j=1,nby                                                       ! beginning of the loop over all cells along y.
          lamplu(i,j,stype)=val2d(i,j)                                    ! remplir the array of the lamp type: stype
          totlu(stype)=totlu(stype)+lamplu(i,j,stype)                     ! the total lamp flux should be non-null to proceed to the calculations
         enddo                                                            ! end of the loop over all cells along y.
        enddo                                                             ! end of the loop over all cells along x.
       enddo                                                              ! end of the loop 1 over the 120 types of sources. 
c    ==================================================================
c    reading lamp heights
         call 2din(nbx,nby,alfile,val2d)
         do i=1,nbx                                                       ! beginning of the loop over all cells along x.
           do j=1,nby                                                     ! beginning of the loop over all cells along y.
             lampal(i,j)=val2d(i,j)                                       ! filling of the array for the lamp stype
           enddo                                                          ! end of the loop over all cells along y.
         enddo                                                            ! end of the loop over all cells along x.
c    ==================================================================
c    reading subgrid obstacles average height
        call 2din(nbx,nby,ohfile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
         do j=1,nby                                                       ! beginning of the loop over all cells along y.
          obsH(i,j)=val2d(i,j)                                            ! filling of the array
         enddo                                                            ! end of the loop over all cells along y.
        enddo 
c    ==================================================================
c    reading subgrid obstacles average distance
        call 2din(nbx,nby,odfile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
         do j=1,nby                                                       ! beginning of the loop over all cells along y.
          drefle(i,j)=val2d(i,j)/2.                                       ! Filling of the array
          if (drefle(i,j).eq.0.) drefle(i,j)=dx                           ! when outside a zone, block to the size of the cell (typically 1km)
         enddo                                                            ! end of the loop over all cells along y.
        enddo 
c    ==================================================================
c    reading subgrid obstacles filling factor
        call 2din(nbx,nby,offile,val2d)
        do i=1,nbx                                                        ! beginning of the loop over all cells along x.
         do j=1,nby                                                       ! beginning of the loop over all cells along y.
          ofill(i,j)=val2d(i,j)                                           ! Filling of the array 0-1
         enddo                                                            ! end of the loop over all cells along y.
        enddo 
   
c=======================================================================
c        reading of the scattering parameters 
c=======================================================================
       open(unit = 1, file = diffil,status= 'old')                        ! opening file containing the parameters of scattering.
c                                                                         ! the scattering file is generated by the program imies 
c                                                                         ! of the progiciel AODSEM (Martin Aube).
        read(1,*)                                                           
        read(1,*)
        read(1,*)
        do i=1,181
         read(1,*) anglea(i), fdifa(i)                                    ! reading of the Scattering functions and the associate angle a 
c                                                                         ! this fonction of 0 a 180 degrees soit 181 lignes.
         fdifan(i)=fdifa(i)/pix4                                          ! Normalisation of the fonction a 4 pi (the integral of the 
c                                                                         ! fonction provided over all solid angles the doit etre egale a 4 pi).
c                                                                         ! in fact the file .mie.out is normalized ainsi (revefie par 
c                                                                         ! M. Aube en avril 2009)
        enddo
        do i = 1,7
         read(1,*)
        enddo
        read(1,*) extinc                                                  ! reading of the cross section extinction of the aerosols.
        read(1,*) scatte                                                  ! reading of the cross section of scattering of the aerosols.
       close(1)
       secdif=scatte/extinc                                               ! Rapport (sigmadif/sigmatotal).
c======================================================================
c        Quelques operations preparatoires
c======================================================================
       dy=dx                                                              ! we consider than the echelle is the same over the two axes

       z_obs=z_o+altsol(x_obs,y_obs)                                      ! z_obs = the local observer elevation plus the height of observation above ground (z_o)
c find nearest vertical grid
       distm=1000000000.
       do i=1,height
          dist=abs(cellh(i)-z_obs)
          if (dist.le.distm) then
             distm=dist
             zcello=i
          endif
       enddo
       z_obs=cellh(zcello)                                                ! Attribution of the value in meter to the position z of the observateur.
       largx=dx*real(nbx)                                                 ! computation of the Width along x of the case.
       largy=dy*real(nby)                                                 ! computation of the Width along y of the case.

       write(2,*) 'Width of the domain [NS](m):',largx,'#cases:',nbx
       write(2,*) 'Width of the domain [EO](m):',largy,'#cases:',nby
       write(2,*) 'Size of a cell (m):',dx,' X ',dy
       write(2,*) 'latitu center:',latitu
c=======================================================================
c        computation of the tilt of the cases along x and along y
c=======================================================================
       do i=1,nbx                                                         ! beginning of the loop over the column (longitude) of the domain.
        do j=1,nby                                                        ! beginning of the loop over the rows (latitu) of the domain.
         if (i.eq.1) then                                                 ! specific case close to the border of the domain (vertical side left).
          inclix(i,j)=atan((altsol(i+1,j)-altsol(i,j))/real(dx))          ! computation of the tilt along x of the surface.
         elseif (i.eq.nbx) then                                           ! specific case close to the border of the domain (vertical side right).
          inclix(i,j)=atan((altsol(i-1,j)-altsol(i,j))/(real(dx)))        ! computation of the tilt along x of the surface.
         else
          inclix(i,j)=atan((altsol(i+1,j)-altsol(i-1,j))/(2.              ! computation of the tilt along x of the surface.
     1    *real(dx)))
         endif
         if (j.eq.1) then                                                 ! specific case close to the border of the domain (horizontal side down).
          incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j))/(real(dy)))        ! computation of the tilt along y of the surface.
         elseif (j.eq.nby) then                                           ! specific case close to the border of the domain (horizontal side up).
          incliy(i,j)=atan((altsol(i,j-1)-altsol(i,j))/(real(dy)))        ! computation of the tilt along y of the surface.
         else
          incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j-1))/(2.              ! computation of the tilt along y of the surface
     1    *real(dy)))
         endif
        enddo                                                             ! end of the loop over the rows (latitu) of the domain
       enddo                                                              ! end of the loop over the column (longitude) of the domain
c=======================================================================
c        beginning of the loop over the line of sight voxels
c=======================================================================
       print*,'Determining voxels crossed by the line of sight...'
       scalef=1.00001
       call lignevisee(x_obs,y_obs,z_obs,dx,dy,angvis,                    ! Determination of the viewing line (line of sight voxels).
     + azim,nbx,nby,cloudz,lcible,ncible,scalef)
       print*,'Total final cells: ',ncible
       print*,'============'
       print*,'x   y   z'
       print*,'------------'
       do i=1,ncible
          write(*,1110) lcible(i,1),lcible(i,2),lcible(i,3)
       enddo
       if (zcello.ge.cloudz) then
          print*,'The observer is inside the cloud! Abort computing.'
          stop
       endif
 1110  format(I3,1x,I3,1x,I3)
       fctcld=0.
       ftocap=0.                                                          ! Initialisation of the value of flux received by the sensor
       fcapt=1.
       do icible=1,ncible                                                 ! beginning of the loop over the line of sight voxels
      if ((fcapt.ge.ftocap/stoplim).or.(cloudt.ne.0)) then                ! stop the calculation of the viewing line when the increment is lower than 1/stoplim
        if (fcapt.eq.1.) fcapt=0.
        if (icible.ge.nvis0) then                                         ! beginning condition for continuing of a computation stopped
         itotci=0.                                                        ! Initialisation of the contribution of the line of sight at the sensor level
         do i=1,nbx
          do j=1,nby
            ITC(i,j)=0.
          enddo
         enddo
         zcellc=lcible(icible,3)                                          ! Definition of the vertical position (voxel) of the line of sight
         z_c=cellh(zcellc)                                                ! Definition of the vertical position (meter) of the line of sight
         y_c=lcible(icible,2)                                             ! Definition of the position (voxel) of the line of sight
         x_c=lcible(icible,1)                                             ! Definition of the position (voxel) of the line of sight
         print*,'=================================================='
         print*,' Progression along the line of sight :',
     +   icible,'/',ncible,'(',x_c,',',y_c,')'
         print*,' Voxel height =',z_c,' m'
         print*,' Voxel thickness =',cthick(zcellc),' m'
         write(2,*) '=================================================='
         write(2,*) ' Progression along the line of sight :',
     +   icible,'/',ncible,'(',x_c,',',y_c,')'
         write(2,*) ' Voxel height =',z_c,' m'
      write(2,*) ' Voxel thickness =',cthick(zcellc),' m'
         if( (x_c.gt.nbx).or.(x_c.lt.1).or.(y_c.gt.nby).or.(y_c.lt.1)     ! Condition line of sight voxel inside the modelling domain
     +      .or.(zcellc.gt.height).or.(zcellc.lt.1) )then
         else
          if((x_c.eq.x_obs).and.(y_c.eq.y_obs).and.                       ! for the moment, if the line of sight voxel is the observer voxel, 
c                                                                         ! we do not compute the scattered flux
     +    (zcellc.eq.zcello))then
           if (verbose.eq.1) then
             print*,'Scat voxel = Observer voxel' 
           endif
          else
           dis_obs=sqrt((z_c-z_obs)**2.+((real(y_c-y_obs))*dy)**2.
     a     +((real(x_c-x_obs))*dx)**2.)
           if (dis_obs.eq.0.) then
            print*,'ERROR problem with dis_obs',dis_obs
            print*,x_c,x_obs,y_c,y_obs,z_c,z_obs
            stop
           endif
           ometif=pi*(diamobj/2.)**2./dis_obs**2.
           omefov=lfente*longfe/focal**2.                                 ! computation of the solid angle of the fente projete sur le ciel.

           zcdown=z_c-0.5*cthick(zcellc)                                  ! lower limit of the line of sight voxel.
           zcup=z_c+0.5*cthick(zcellc)                                    ! upper limit of the line of sight voxel.
c=======================================================================
c        beginning of the loop over the types of light sources
c=======================================================================
           do stype=1,ntype                                               ! beginning of the loop over the source types.
            if (totlu(stype).ne.0.) then                                  ! check if there are any flux in that source type
                                                                          ! otherwise skip this lamp
            print*,' Turning on zone',stype
            write(2,*) ' Turning on zone',stype
            itotty=0.                                                     ! Initialisation of the contribution of a source types to 
c                                                                         ! the intensity toward the sensor by a line of sight voxel.
            do x_s=1,nbx
             do y_s=1,nby
              ITT(x_s,y_s,stype)=0.
             enddo
            enddo     
            do x_s=imin(stype),imax(stype)                                ! beginning of the loop over the column (longitude the) of the domain.
             do y_s=jmin(stype),jmax(stype)                               ! beginning of the loop over the rows (latitud) of the domain.
              if (lamplu(x_s,y_s,stype) .ne. 0.) then                     ! if the luminosite of the case is null, the program ignore this case.
               z_s=(altsol(x_s,y_s)+lampal(x_s,y_s))                      ! Definition of the position (metre) vertical of the source.
c computation of the distance source-line of sight-observer if this distance is lower than  dx/2, pas of computation effectue
c the raison is que autrement on passe par of the voxels tres proches of the source and on is jamais dans of telles
c conditions lorsqu'on observe le ciel. C is un probleme cree par le fait than the sources and l observateur
c sont toujours considered au centre of the voxels.
               dsco=sqrt((real(x_s-x_c)*dx)**2.+(real(y_s-y_c)*dx)**2.+
     +         (z_s-z_c)**2.)+sqrt((real(x_obs-x_c)*dx)**2.+(real(y_obs
     +         -y_c)*dx)**2.+(z_obs-z_c)**2.)
          if (dsco.ge.dminlp) then                                        ! beginning condition distance source-line of sight-observer >= dx/2
c **********************************************************************************************************************
c *     computation of the direct intensity toward the sensor by a line of sight voxel en provenance of the source         *
c **********************************************************************************************************************         
               dirck=0                                                    ! Initialisation of the verification of the position of the source.
               if ( (x_s.eq.x_c).and.(y_s.eq.y_c).and.( abs(z_s-z_c)      ! if the positions x and y of the source and the line of sight voxel are the 
c                                                                         ! memes alors.
     +         .lt.(cthick(zcellc)/2.) ) )then
                dirck=1
                if (verbose.eq.1) then
                 print*,'Source insiof scat voxel' 
                endif
               endif                                                      ! end of the case positions x and y source and line of sight voxel identical.
               if (dirck.ne.1) then                                       ! the source is not at the line of sight voxel position
c=======================================================================
c        computation of the zenithal angle between the source and the line of sight
c=======================================================================
c
c computation of the horizon for the resolved shadows direct              ! horizon resolution is 1 degree
                call anglezenithal
     +          (x_s,y_s,z_s,x_c,y_c,z_c,dx,dy,angzen)                    ! computation of the zenithal angle between the source and the line of sight voxel.
                call angleazimutal(x_s,y_s,x_c,y_c,dx,dy,angazi)          ! computation of the angle azimutal direct line of sight-source
            if (angzen.gt.pi/4.) then                                     ! 45deg. it is unlikely to have a 1km high mountain less than 1
                call horizon(x_s,y_s,z_s,dx,dy,nbx,nby,altsol,
     +          latitu,angzen,angazi,zhoriz) 
                if (angzen.lt.zhoriz) then                                ! shadow the path line of sight-source is not below the horizon => we compute
                   hh=1.
                else
                   hh=0.
                endif
            else
               hh=1.
            endif

c                                                                         ! beginning condition above the horizon direct
c sub-grid obstacles             
                angmin=pi/2.-atan((altsol(x_s,y_s)+obsH(x_s,y_s)
     +          -z_s)/drefle(x_s,y_s))
                if (angzen.lt.angmin) then                                ! condition sub-grid obstacles direct.
                   ff=0.
                else 
                   ff=ofill(x_s,y_s)
                endif
c
c=======================================================================
c computation of the transmittance between the source and the line of sight
c=======================================================================
                  anaz=zero
                  call transmitm(angzen,anaz,x_s,y_s,z_s,x_c,y_c,z_c,
     +            lambda,dx,dy,pressi,transm)     
                  call transmita(angzen,anaz,x_s,y_s,z_s,x_c,y_c,z_c,
     +            dx,dy,taua,transa)
c=======================================================================
c computation of the Solid angle of the line of sight voxel seen from the source
c=======================================================================
                  xc=dble(x_c)*dble(dx)                                   ! Position in meters of the observer voxel (longitude).
                  yc=dble(y_c)*dble(dy)                                   ! Position in meters of the observer voxel (latitu).
                  zc=dble(z_c)                                            ! Position in meters of the observer voxel (altitude).
                  xn=dble(x_s)*dble(dx)                                   ! Position in meters of the source (longitude).
                  yn=dble(y_s)*dble(dy)                                   ! Position in meters of the source (latitu).
                  zn=dble(z_s)                                            ! Position in meters of the source (altitude).
c    ------------------------------------
c    solid angle for the central plane xy
c    ------------------------------------
                  if (z_c .ne. z_s) then
                   call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +             r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,
     +             r4z) 
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Call of the routine anglesolide to compute the solid angle 
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                   ! along the surface xy.
                   omega1 = omega
                  else
                   omega1=0.
                  endif
c     ------------------------------------
c     solid angle for the central plane zx
c     ------------------------------------
                  if (y_c .ne. y_s) then                                  ! if the latitu of the observer voxel is the same as
c                                                                         ! of the source voxel, we are not computing the angle solide
c                                                                         ! for the surface zx car il is egal a 0
                   call planzx(dx,xc,xn,yc,yn,zc,zn,cthick,
     +             zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y
     +             ,r4z)
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Call of the routine anglesolide to compute the solid angle 
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                   ! along the surface zx.
                  else
                   omega=0.
                  endif
                  if (omega.gt.0.) then
                   if (omega .gt. omega1) omega1 = omega                  ! On garof the solid angle le plus grand jusqu'a present.
                  endif
c     ------------------------------------
c     solid angle for the central plane yz
c     ------------------------------------
                  if (x_c .ne. x_s) then                                  ! if the longitude of the observer voxel is the same as
c                                                                         ! of the source voxel, we are not computing the angle solide
c                                                                         ! for the surface yz car il is egal a 0.
                   call planyz(dy,xc,xn,yc,yn,zc,zn,cthick,
     +             zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,
     +             r4z)
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Routine anglesolide to compute the solid angle along the surface yz.
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                  else 
                   omega=0.
                  endif
                  if (omega.gt.0.) then
                   if (omega .gt. omega1) omega1 = omega                  ! On garof the solid angle le plus grand
                  endif
                  omega=omega1
c=======================================================================
c    Estimation of the half of the underlying angle of the solid angle    ! this angle is used to get a better estimate (average) of 
c                                                                         ! P_dir for le cas of grans solid angles the ou pvalno 
c=======================================================================  ! varie significativement sur +- ouvang.
                  ouvang=sqrt(omega/pi)                                   ! Angle in radian.
                  ouvang=ouvang*180./pi                                   ! Angle in degrees.
c   
c=======================================================================
c computation of the photometric function of the source toward the line of sight voxel
c=======================================================================
c   
                  anglez=nint(180.*angzen/pi)
                  if (anglez.lt.0) anglez=-anglez
                  if (anglez.gt.180) anglez=360-anglez
                  anglez=anglez+1                                         ! Transform the angle in integer degree into the position in the array.
c
c  moyenner sur +- ouvang	
c
                  naz=0
                  nbang=0.
                  P_dir=0.
                  do na=-nint(ouvang),nint(ouvang)
                   naz=anglez+na+1
                   if (naz.lt.0) naz=-naz
                   if (naz.gt.181) naz=362-naz                            ! symetric function
                   if (naz.eq.0) naz=1
                   P_dir=P_dir+pvalno(naz,stype)
                   nbang=nbang+1.
                  enddo
                  P_dir=P_dir/nbang
c
c=======================================================================
c        computation of the flux direct reaching the line of sight voxel
c=======================================================================
                  fldir=lamplu(x_s,y_s,stype)*P_dir*omega*
     1            transm*transa*(1.-ff)*hh                                   ! correction for obstacle filling factor
c=======================================================================
c   computation of the scattering probability of the direct light
c=======================================================================
                  if (angzen.lt.(pi/2.)) then                             ! Attribution of the initial and final limit of the path of 
c                                                                         ! scattering dans the voxel.
                   zidif=zcdown
                   zfdif=zcup
                  else
                   zidif=zcup
                   zfdif=zcdown
                  endif
                  anaz=angazi
                  call transmitm (angzen,anaz,iun,iun,zidif,ideux,iun,    ! Transmittance molecular of the scattering voxel.
     +            zfdif,lambda,dx,dy,pressi,tran1m)
                  call transmita (angzen,anaz,iun,iun,zidif,ideux,iun,    ! Transmittance aerosols of the scattering voxel.
     +            zfdif,dx,dy,taua,tran1a)
                  call angle3points (x_s,y_s,z_s,x_c,y_c,z_c,x_obs,       ! scattering angle.
     +            y_obs,z_obs,dx,dy,angdif)
                  call diffusion(omega,angdif,tran1a,tran1m,              ! scattering probability of the direct light.     
     +            secdif,fdifan,pdifdi)
c=======================================================================
c   computation of the source contribution a the direct intensity toward the sensor by a line of sight voxel
c=======================================================================
                  intdir=fldir*pdifdi

                if (cloudt.ne.0) then                                     ! line of sight voxel = cloud
                  if (cloudh(cloudt).eq.zcellc) then
                     call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +               dx,dy,azencl)                                        ! zenith angle from cloud to observer                     
                     call cloudreflectance(angzen,cloudt,rcloud)          ! cloud intensity from direct illum
                     icloud=icloud+
     +               fldir*rcloud*abs(cos(azencl))/pi
                  endif
                endif
c                else
c                endif                                                     ! end condition below the horizon direct? 
               endif                                                      ! end of the case Position Source is not equal to the line of sight voxel position
c  end of the computation of the direct intensity
c **********************************************************************************************************************
c * computation of the reflected intensity toward the sensor by a line of sight voxel from the source           *
c **********************************************************************************************************************
c=======================================================================
c        etablissement of the conditions ands boucles
c=======================================================================
               itotind=0.                                                 ! Initialisation of the reflected intensity of the source
               itotrd=0.
       boxx=nint(drefle(x_s,y_s)/dx)                                      ! Number of column to consider left/right of the source 
c                                                                         ! for the reflexion.
       boxy=nint(drefle(x_s,y_s)/dy)                                      ! Number of column to consider up/down of the source for 
c                                                                         ! the reflexion.
       xsrmi=x_s-boxx
       if (xsrmi.lt.1) xsrmi=1
       xsrma=x_s+boxx
       if (xsrma.gt.nbx) xsrma=nbx
       ysrmi=y_s-boxy
       if (ysrmi.lt.1) ysrmi=1
       ysrma=y_s+boxy
       if (ysrma.gt.nby) ysrma=nby
               do x_sr=xsrmi,xsrma                                        ! beginning of the loop over the column (longitude) reflecting.
                do y_sr=ysrmi,ysrma                                       ! beginning of the loop over the rows (latitu) reflecting.
                 irefl=0.
                 z_sr=altsol(x_sr,y_sr)   
                  if( (x_sr.gt.nbx).or.(x_sr.lt.1).or.(y_sr.gt.nby)
     +            .or.(y_sr.lt.1) )then  
                   if (verbose.eq.1) then
                    print*,'Ground cell out of the borders'
                   endif
                  else  
                   if((x_s.eq.x_sr).and.(y_s.eq.y_sr).and.(z_s.eq.z_sr))
     +             then
                    if (verbose.eq.1) then
                     print*,'Source pos = Ground cell'
                    endif
                   else
                    if (srefl(x_sr,y_sr).ne.0.) then                      ! Condition: the surface reflectance is not null
                     haut=-real(x_s-x_sr)*dx*tan(inclix(x_sr,y_sr))       ! if haut is negative, the ground cell is lighted from below
     1               -real(y_s-y_sr)*dy                                 
     2               *tan(incliy(x_sr,y_sr))+z_s-z_sr
                     if (haut .gt. 0.) then                               ! Condition: the ground cell is lighted from above
c=======================================================================
c        computation of the zenithal angle between the source and the  surface reflectance
c=======================================================================
                      call anglezenithal(x_s,y_s,z_s,x_sr,y_sr,z_sr,dx,   ! computation of the zenithal angle between the source and the line of sight voxel.
     +                dy,angzen)                                          ! end of the case "observer at the same latitu/longitude than the source".

c=======================================================================
c        computation of the transmittance between the source and the ground surface
c=======================================================================
                      anaz=zero      
                      call transmitm(angzen,anaz,x_s,y_s,z_s,x_sr,y_sr,
     +                z_sr,lambda,dx,dy,pressi,transm)          
                      call transmita(angzen,anaz,x_s,y_s,z_s,x_sr,y_sr,
     +                z_sr,dx,dy,taua,transa)
c=======================================================================
c     computation of the Solid angle of the reflecting cell seen from the source
c=======================================================================
                      xc=dble(x_sr)*dble(dx)                              ! Position in meters of the observer voxel (longitude).
                      yc=dble(y_sr)*dble(dy)                              ! Position in meters of the observer voxel (latitu).
                      zc=dble(z_sr)                                       ! Position in meters of the observer voxel (altitude).
                      xn=dble(x_s)*dble(dx)                               ! Position in meters of the source (longitude).
                      yn=dble(y_s)*dble(dy)                               ! Position in meters of the source (latitu).
                      zn=dble(z_s)                                        ! Position in meters of the source (altitude).
                      epsilx=inclix(x_sr,y_sr)                            ! tilt along x of the ground reflectance
                      epsily=incliy(x_sr,y_sr)                            ! tilt along x of the ground reflectance
                      if (dx.gt.drefle(x_s,y_s)*2.) then                  ! use a sub-grid surface when the mean free path to the ground is smaller than the cell size
                       if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then
                        dxp=drefle(x_s,y_s)*2.
                       else
                        dxp=dx
                       endif
                      else
                       dxp=dx
                      endif
                      if (dy.gt.drefle(x_s,y_s)*2.) then
                       if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then         
                        dyp=drefle(x_s,y_s)*2.
                       else
                        dyp=dy
                       endif
                      else
                       dyp=dy
                      endif              
                      r1x=xc-dble(dxp)/2.-xn                              ! computation of the composante along x of the first vector.
                      r1y=yc+dble(dyp)/2.-yn                              ! computation of the composante along y of the first vector.
                      r1z=zc-tan(dble(epsilx))*dble(dxp)/2.+tan(dble(
     +                epsily))*dble(dyp)/2.-zn                            ! computation of the composante en z of the first vector.
                      r2x=xc+dble(dxp)/2.-xn                              ! computation of the composante along x of the second vector.
                      r2y=yc+dble(dyp)/2.-yn                              ! computation of the composante along y of the second vector.
                      r2z=zc+tan(dble(epsilx))*dble(dxp)/2.+tan(dble(
     +                epsily))*dble(dyp)/2.-zn                            ! computation of the composante en z of the second vector.
                      r3x=xc-dble(dxp)/2.-xn                              ! computation of the composante along x of the third vector.
                      r3y=yc-dble(dyp)/2.-yn                              ! computation of the composante along y of the third vector.
                      r3z=zc-tan(dble(epsilx))*dble(dxp)/2.-tan(
     +                dble(epsily))*dble(dyp)/2.-zn                       ! computation of the composante en z of the third vector.
                      r4x=xc+dble(dxp)/2.-xn                              ! computation of the composante along x of the fourth vector.
                      r4y=yc-dble(dyp)/2.-yn                              ! computation of the composante along y of the fourth vector.
                      r4z=zc+tan(dble(epsilx))*dble(dxp)/2.-tan(
     +                dble(epsily))*dble(dyp)/2.-zn                       ! computation of the composante en z of the fourth vector.
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Call of the routine anglesolide to compute the angle solide.
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z) 
c
c=======================================================================
c    isimation of the half of the underlying angle of the solid angle     ! this angle servira a obtenir un meilleur isime (moyenne) of 
c                                                                         ! P_dir for le cas of grans solid angles the ou pvalno
c=======================================================================  ! varie significativement sur +- ouvang.
                      ouvang=sqrt(omega/pi)                               ! Angle in radian.
                      ouvang=ouvang*180./pi                               ! Angle in degrees.
c   
c=======================================================================
c        computation of the photometric function of the lampadaire toward the  surface reflectance
c=======================================================================
c    
                      anglez=nint(180.*angzen/pi)
                      if (anglez.lt.0) anglez=-anglez
                      if (anglez.gt.180) anglez=360-anglez
                      anglez=anglez+1                                     ! Transform the angle in integer degree into the position in the array.
c
c  moyenner sur +- ouvang	
c
c
                      naz=0
                      nbang=0.
                      P_indir=0.
                      do na=-nint(ouvang),nint(ouvang)
                       naz=anglez+na+1
                       if (naz.lt.0) naz=-naz
                       if (naz.gt.181) naz=362-naz                        ! symetric function
                       if (naz.eq.0) naz=1
                       P_indir=P_indir+pvalno(naz,stype)
                       nbang=nbang+1.
                      enddo
                      P_indir=P_indir/nbang
c 
c=======================================================================
c        computation of the flux reaching the reflecting cell
c=======================================================================
                      flrefl=lamplu(x_s,y_s,stype)*P_indir*
     a                omega*transm*transa
c=======================================================================
c        computation of the intensity reflechie leaving the  surface reflectance
c=======================================================================
                      irefl1=flrefl*srefl(x_sr,y_sr)/pi                   ! The factor 1/pi comes from the normalisation of the fonction 
                      if (effdif.gt.(dx+dy)/2.) then 
                       call reflexdbledif (x_sr,y_sr,z_sr,x_c,y_c,
     +                 zcellc,dx,dy,effdif,nbx,nby,stepdi,
     +                 irefl1,lambda,pressi,taua,zcup,
     +                 zcdown,secdif,fdifan,x_obs,y_obs,z_obs,
     +                 epsilx,epsily,irefdi,drefle,obsH,ofill,
     +                 altsol,latitu,cloudt,cloudh,icloud)
                      endif
                      itotrd=itotrd+irefdi      
c
c  the projection apparente is calculee a partir of the produit scalaire of the vecteur normal a 
c  the reflecting surface and the line reflecting surface toward scattering voxel ou cible
c  it is the cosine correction for the lambertian reflectance for finite elements
c         
                      projap=(-tan(epsilx)*real(x_c-x_sr)*dx-
     +                tan(epsily)*real(y_c-y_sr)*dy+1.*(cellh(
     +                zcellc)-z_sr))/(sqrt(tan(epsilx)**2.+tan(epsily)
     +                **2.+1.)*sqrt((real(x_c-x_sr)*dx)**2.+(real(y_c-
     +                y_sr)*dy)**2.+(cellh(zcellc)-z_sr)**2.))
                      if (projap.lt.0.) projap=0.
c                                                                         ! no matter the direction we are taking the absolute value of cos theta 
c                                                                          
c verify if there is shadow between sr and line of sight voxel

                 call anglezenithal(x_sr,y_sr,z_sr,x_c,y_c,z_c,dx,        ! zenithal angle between the reflecting surface and the line of sight voxel.
     +           dy,angzen)     
                 call angleazimutal(x_sr,y_sr,x_c,y_c,dx,dy,angazi)       ! computation of the azimutal angle reflect-line of sight

            if (angzen.gt.pi/4.) then                                     ! 45deg. it is unlikely to have a 1km high mountain less than 1
                 call horizon(x_sr,y_sr,z_sr,dx,dy,nbx,nby,altsol,
     +           latitu,angzen,angazi,zhoriz) 

                 if (angzen.lt.zhoriz) then                               ! the path line of sight-reflec is not below the horizon => we compute
                    hh=1.
                 else
                    hh=0.
                 endif                                                    ! end condition reflecting surf. above horizon
            else
               hh=1.
            endif
c
c ????????????????????????????????
                 irefl=irefl1*projap
c ?????????????????????????????????
c
c=======================================================================
c        Case: line of sight position = Position of reflecting cell
c=======================================================================
                      if((x_c.eq.x_sr).and.(y_c.eq.y_sr).and.
     +                (z_c.eq.z_sr)) then
                       intind=irefl*(1.-ff)*hh
                      else
c
c            
c obstacle                 
                       angmin=pi/2.-atan(obsH(x_sr,y_sr)/
     +                 drefle(x_sr,y_sr))
                       if (angzen.lt.angmin) then                         ! condition obstacle reflected.
                          ff=0.
                       else 
                          ff=ofill(x_sr,y_sr)
                       endif
c
c=======================================================================
c        computation of the transmittance between the  ground surface and the line of sight voxel
c=======================================================================
                        anaz=zero
                        call transmitm(angzen,anaz,x_sr,y_sr,z_sr,x_c,
     +                  y_c,z_c,lambda,dx,dy,pressi,transm)        
                        call transmita(angzen,anaz,x_sr,y_sr,z_sr,x_c,
     +                  y_c,z_c,dx,dy,taua,transa)
c=======================================================================
c     computation of the Solid angle of the line of sight voxel seen from the reflecting cell
c=======================================================================
                        xc=dble(x_c)*dble(dx)                             ! Position in meters of the observer voxel (longitude).
                        yc=dble(y_c)*dble(dy)                             ! Position in meters of the observer voxel (latitu).
                        zc=dble(z_c)                                      ! Position in meters of the observer voxel (altitude).
                        xn=dble(x_sr)*dble(dx)                            ! Position in meters of the source (longitude).
                        yn=dble(y_sr)*dble(dy)                            ! Position in meters of the source (latitu).
                        zn=dble(z_sr)                                     ! Position in meters of the source (altitude).
c    ------------------------------------
c    solid angle for the central plane xy
c    ------------------------------------
                        if (z_c .ne. z_sr) then
                         call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                   r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Call of the routine anglesolide to compute the solid angle 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! along the surface xy.
                         omega1 = omega
                        else
                         omega1=0.
                        endif
c     ------------------------------------
c     solid angle for the central plane zx
c     ------------------------------------
                        if (y_c .ne. y_sr) then                           ! if the latitu of the observer voxel is the same as
c                                                                         ! of the source voxel, we are not computing the angle solide
c                                                                         ! for the surface zx car il is egal a 0.
                         call planzx(dx,xc,xn,yc,yn,zc,zn,
     +                   cthick,zcellc,r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Call of the routine anglesolide to compute the solid angle 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! along the surface zx.
                        else
                         omega=0.
                        endif
                        if (omega.gt.0.) then
                         if (omega.gt.omega1) omega1 = omega              ! On garof the solid angle le plus grand jusqu'a present.
                        endif
c     ------------------------------------
c     solid angle for the central plane yz
c     ------------------------------------
                        if (x_c.ne.x_sr) then                             ! if the longitude of the observer voxel is the same as
c                                                                         ! of the source voxel, we are not computing the angle solide
c                                                                         ! for the surface yz car il is egal a 0.
                         call planyz(dy,xc,xn,yc,yn,zc,zn,
     +                   cthick,zcellc,r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Call of the routine anglesolide to compute the solid angle 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! along the surface yz.
                        else 
                         omega=0.
                        endif
                        if (omega.gt.0.) then
                         if (omega.gt.omega1) omega1=omega                ! we keep the largest solid
                        endif
                        omega=omega1 

c=======================================================================
c        computation of the flux reflected reaching the line of sight voxel
c=======================================================================
                        flindi=irefl*omega*transm*
     +                  transa*(1.-ff)*hh                                 ! obstacles correction
                if (cloudt.ne.0) then                                     ! line of sight voxel = cloud
                  if (cloudh(cloudt).eq.zcellc) then
                     call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +               dx,dy,azencl)                                        ! zenith angle from cloud to observer                     
                     call cloudreflectance(angzen,cloudt,rcloud)          ! cloud intensity from reflected illum
                     icloud=icloud+
     +               flindi*rcloud*abs(cos(azencl))/pi
                  endif
                endif
c=======================================================================
c   computation of the scattering probability of the reflected light
c=======================================================================
                        if (angzen.lt.(pi/2.)) then                       ! Attribution of the initial and final limits of the path of 
c                                                                         ! scattering voxel.
                         zidif=zcdown
                         zfdif=zcup
                        else
                         zidif=zcup
                         zfdif=zcdown
                        endif 
                        anaz=angazi       
                        call transmitm(angzen,anaz,iun,iun,zidif,         ! Transmittance molecular of the scattering voxel.
     +                  ideux,iun,zfdif,lambda,dx,dy,pressi,tran1m)
                        call transmita(angzen,anaz,iun,iun,zidif,         ! Transmittance aerosols of the scattering voxel.
     +                  ideux,iun,zfdif,dx,dy,taua,tran1a)
                        call angle3points (x_sr,y_sr,z_sr,x_c,y_c,z_c,    ! scattering angle.
     +                  x_obs,y_obs,z_obs,dx,dy,angdif)
                        call diffusion(omega,angdif,tran1a,               ! scattering probability of the reflected light.
     +                  tran1m,secdif,fdifan,pdifin)
c=======================================================================
c   computation of the reflected intensity toward the sensor by a reflecting cell
c=======================================================================
                        intind=flindi*pdifin*(1.-ff)*hh 
                      endif                                               ! end of the case Posi reflecting cell =  line of sight voxel position                                 
                      itotind=itotind+intind                              ! Sum of the intensities of each reflecting cell.
                     endif                                                ! end of the condition surface not lighted from the top.
                    endif                                                 ! end of the condition non-zero reflectancee.
                   endif                                                  ! end of the condition reflecting cell is not on the source.
                  endif                                                   ! end of the condition surface of the domain.
                enddo                                                     ! end of the loop over the rows (latitu) reflecting.
               enddo                                                      ! end of the loop over the column (longitude) reflecting.
c   end of the computation of the reflected intensity        
c **********************************************************************************************************************
c * computation of the scattered intensity toward the sensor by a line of sight voxel coming from the source            *
c **********************************************************************************************************************
c
c=======================================================================
c    Determination of the scattering voxels en fonction of the source voxel and the line of sight voxel
c=======================================================================

               itodif=0.                                                  ! Initialisation of the scattered intensity by a source dans 
c                                                                         ! a line of sight voxel compute the double scattering only if 
               if (effdif.gt.(dx+dy)/2.) then                             ! radius of scattering is larger than the size of the voxels.

                call zone_diffusion(x_s,y_s,z_s,x_c,y_c,zcellc,dx,dy,
     +          effdif,nbx,nby,altsol,zondif,ndiff)
                do idi=1,ndiff,stepdi                                     ! beginning of the loop over the scattering voxels.
                 x_dif=zondif(idi,1)
                 y_dif=zondif(idi,2)
                 zceldi=zondif(idi,3)
                 z_dif=cellh(zceldi)             
                 if((x_dif.gt.nbx).or.(x_dif.lt.1).or.(y_dif.gt.nby).     ! Condition scattering voxel of the domain.
     +           or.(y_dif.lt.1)) then     
c
c !!!!!!!rien ici???????
c          
                 else
                  if ((x_s.eq.x_dif).and.(y_s.eq.y_dif).and.(z_s.eq. 
     +            z_dif)) then
                   if (verbose.eq.1) then
                     print*,'Scat voxel = Source position'
                   endif
                  elseif ((x_c.eq.x_dif).and.(y_c.eq.y_dif).and. 
     +                 (z_c.eq.z_dif)) then
                  else
c=======================================================================
c        computation of the zenithal angle between the source and the scattering voxel
c=======================================================================


c shadow source-scattering voxel

                   call anglezenithal(x_s,y_s,z_s,x_dif,y_dif,z_dif,dx,
     +             dy,angzen)                                             ! computation of the zenithal angle source-scattering voxel. 

                   call angleazimutal(x_s,y_s,x_dif,y_dif,dx,dy,          ! computation of the angle azimutal line of sight-scattering voxel
     +             angazi)

c                   call horizon(x_s,y_s,z_s,dx,dy,nbx,nby,altsol,
c     +             latitu,angzen,angazi,zhoriz) 
c                   if (angzen.lt.zhoriz) then                             ! beginning condition shadow source-diffusante
c                      hh=1.
c                   else
c                      hh=0.
c                   endif
                    hh=1.

c sub-grid obstacles               
                    angmin=pi/2.-atan((obsH(x_s,y_s)+
     +              altsol(x_s,y_s)-z_s)/drefle(x_s,y_s))
                    if (angzen.lt.angmin) then                            ! condition obstacle source->scattering.
                       ff=0.
                    else 
                       ff=ofill(x_s,y_s)
                    endif
c                                                                    
c=======================================================================
c        computation of the transmittance between the source and the scattering voxel
c=======================================================================
                     anaz=zero
                     call transmitm(angzen,anaz,x_s,y_s,z_s,x_dif,
     +               y_dif,z_dif,lambda,dx,dy,pressi,transm)
                     call transmita(angzen,anaz,x_s,y_s,z_s,x_dif,
     +               y_dif,z_dif,dx,dy,taua,transa) 
c=======================================================================
c     computation of the Solid angle of the par the scattering voxel seen from the source
c=======================================================================

                     xc=dble(x_dif)*dble(dx)                              ! Position in meters of the scattering voxel (longitude).
                     yc=dble(y_dif)*dble(dy)                              ! Position in meters of the scattering voxel (latitu).
                     zc=dble(z_dif)                                       ! Position in meters of the scattering voxel (altitude).
                     xn=dble(x_s)*dble(dx)                                ! Position in meters of the source (longitude).
                     yn=dble(y_s)*dble(dy)                                ! Position in meters of the source (latitu).
                     zn=dble(z_s)                                         ! Position in meters of the source (altitude).
c    ------------------------------------
c    solid angle for the central plane xy
c    ------------------------------------
                     if (z_dif .ne. z_s) then
                      call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                r1x,r1y,r1z,r2x,r2y,r2z,
     +                r3x,r3y,r3z,r4x,r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Calling the routine anglesolide for the computation of the solid angle 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! in the xy plane.
                      omega1 = omega
                     else
                      omega1=0.
                     endif
c     ------------------------------------
c     solid angle for the central plane zx
c     ------------------------------------
                     if (y_dif .ne. y_s) then                             ! if the latitude of the observer voxel is the same as the one
c                                                                         ! of the source voxel, we do not compute the solid angle
c                                                                         ! for le zx plane because is is null.
                      call planzx(dx,xc,xn,yc,yn,zc,zn,cthick,
     +                zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x
     +                ,r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Calling of the routine anglesolide for the computation of the solid angle 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! in the zx plane.
                     else
                      omega=0.
                     endif
                     if (omega.gt.0.) then
                      if (omega .gt. omega1) omega1 = omega               ! We keep the largest solid angle
                     endif
c     ------------------------------------
c     solid angle for the central plane yz
c     ------------------------------------
                     if (x_dif .ne. x_s) then                             ! if the longitude of the observer voxel is the same as the 
c                                                                         ! source voxel, we do not compute the solid angle
c                                                                         ! for the yz plane because it is null.
                      call planyz(dy,xc,xn,yc,yn,zc,zn,cthick,
     +                zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,
     +                r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Calling of the routine anglesolide for the computation of the solid angle 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! in the yz plane.       
                     else 
                      omega=0.
                     endif
                     if (omega.gt.0.) then
                      if (omega .gt. omega1) omega1 = omega               ! We keep the largest solid angle.
                     endif
                     omega=omega1
c=======================================================================
c estimation of the subtended angle of the solid angle                    ! this angle will allow a better estimate (average) of 
c                                                                         ! P_dir for the case of large solid angles when pvalno
c=======================================================================  ! vary significatively in +- ouvang.
                     ouvang=sqrt(omega/pi)                                ! Angle in radian.
                     ouvang=ouvang*180./pi                                ! Angle in degrees.
c 
c=======================================================================
c Computing emission function of the source toward the scattering voxel    
c=======================================================================
c    
                     anglez=nint(180.*angzen/pi)
                     if (anglez.lt.0) anglez=-anglez
                     if (anglez.gt.180) anglez=360-anglez
                     anglez=anglez+1                                      ! Transform the angle in degree integer into position inside the array.
c		
c  moyenner sur +- ouvang	
c
c
                     naz=0
                     nbang=0.
                     P_dif1=0.
                     do na=-nint(ouvang),nint(ouvang)
                      naz=anglez+na+1
                      if (naz.lt.0) naz=-naz
                      if (naz.gt.181) naz=362-naz                         ! symetric function
                      if (naz.eq.0) naz=1                     
                      P_dif1=P_dif1+pvalno(naz,stype)
                      nbang=nbang+1. 
                     enddo
                     P_dif1=P_dif1/nbang 
c
c=======================================================================
c Computing flux reaching the scattering voxel
c=======================================================================
                     fldif1=lamplu(x_s,y_s,stype)*P_dif1*
     +               omega*transm*transa*(1.-ff)*hh
c=======================================================================
c Computing the scattering probability toward the line of sight voxel
c=======================================================================
                     if (angzen.lt.(pi/2.)) then                          ! Attribution of the initial and final limits of the 
c                                                                         ! scattering path.
                      zidif=z_c-0.5*cthick(zceldi)
                      zfdif=z_c+0.5*cthick(zceldi)
                     else
                      zidif=z_c+0.5*cthick(zceldi)
                      zfdif=z_c-0.5*cthick(zceldi)
                     endif
                     anaz=angazi       
                     call transmitm(angzen,anaz,iun,iun,zidif,ideux,      ! Molecular transmittance of the scattering voxel.
     +               iun,zfdif,lambda,dx,dy,pressi,tran1m)
                     call transmita(angzen,anaz,iun,iun,zidif,ideux,      ! Aerosol transmittance of the scattering voxel.
     +               iun,zfdif,dx,dy,taua,tran1a)
                     call angle3points (x_s,y_s,z_s,x_dif,y_dif,z_dif,    ! scattering angle.
     +               x_c,y_c,z_c,dx,dy,angdif)
                     call diffusion(omega,angdif,tran1a,tran1m,           ! scattering probability of the direct light.
     +               secdif,fdifan,pdifd1)
c=======================================================================
c Computing scattered intensity toward the line of sight voxel from the scattering voxel  
c=======================================================================
                     idif1=fldif1*pdifd1
c=======================================================================
c Computing zenith angle between the scattering voxel and the line of sight voxel
c=======================================================================

                     call anglezenithal(x_dif,y_dif,z_dif,x_c,y_c,z_c,
     +               dx,dy,angzen)                                        ! computation of the zenithal angle between the scattering voxel and the 
c                                                                         ! line of sight voxel.
        call angleazimutal(x_dif,y_dif,x_c,y_c,dx,dy,angazi)              ! computation of the azimutal angle surf refl-scattering voxel

c        call horizon(x_dif,y_dif,z_dif,dx,dy,nbx,nby,altsol,
c     +  latitu,angzen,angazi,zhoriz) 
c        if (angzen.lt.zhoriz) then                                        ! beginning shadow condition diffuse-line of sight
c           hh=1.
c        else
c           hh=0.
c        endif
         hh=1.


c                                                                 
c subgrid obstacles                
                     angmin=pi/2.-atan((obsH(x_dif,y_dif)+
     +               altsol(x_dif,y_dif)-z_dif)/drefle(x_dif,y_dif))

                    if (angzen.lt.angmin) then                            ! condition obstacles scattering->line of sight
                       ff=0.
                    else 
                       ff=ofill(x_dif,y_dif)
                    endif
c                                                                   
c=======================================================================
c Computing transmittance between the scattering voxel and the line of sight voxel
c=======================================================================
                      anaz=zero
                      call transmitm(angzen,anaz,x_dif,y_dif,z_dif,x_c
     +                ,y_c,z_c,lambda,dx,dy,pressi,transm)
                      call transmita(angzen,anaz,x_dif,y_dif,z_dif,x_c
     +                ,y_c,z_c,dx,dy,taua,transa) 
c=======================================================================
c Computing the solid angle of the line of sight voxel as seen from the scattering voxel
c=======================================================================
                      xc=dble(x_c)*dble(dx)                               ! Position in meters of the line of sight voxel (longitude).
                      yc=dble(y_c)*dble(dy)                               ! Position in meters of the line of sight voxel (latitu).
                      zc=dble(z_c)                                        ! Position in meters of the line of sight voxel (altitude).
                      xn=dble(x_dif)*dble(dx)                             ! Position in meters of the scattering voxel (longitude).
                      yn=dble(y_dif)*dble(dy)                             ! Position in meters of the scattering voxel (latitu).
                      zn=dble(z_dif)                                      ! Position in meters of the scattering voxel (altitude).
c    ------------------------------------
c    solid angle for the central plane xy
c    ------------------------------------
                      if (z_c .ne. z_dif) then
                       call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                 r1x,r1y,r1z,r2x,r2y,r2z
     +                 ,r3x,r3y,r3z,r4x,r4y,r4z)
                       call anglesolide(omega,r1x,r1y,r1z,                ! Calling the routine anglesolide for the calculation of the solid angle 
     +                 r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)               ! in the xy plane.   
                       omega1 = omega
                      else
                       omega1=0.
                      endif
c     ------------------------------------
c     solid angle for the central plane zx
c     ------------------------------------
                      if (y_c .ne. y_dif) then                            ! if the latitude of the observer voxel is the same as the
c                                                                         ! source voxel, we do not calculate the solid angle
c                                                                         ! for the zx plane because it is null.
                       call planzx(dx,xc,xn,yc,yn,zc,zn,
     +                 cthick,zcellc,r1x,r1y,r1z,r2x,r2y,r2z,
     +                 r3x,r3y,r3z,r4x,r4y,r4z)
                       call anglesolide(omega,r1x,r1y,r1z,                ! Calling of the routine anglesolide for the calculation of the solid angle 
     +                 r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)               ! in the zx plane.
                      else
                       omega=0.
                      endif
                      if (omega.gt.0.) then
                       if (omega .gt. omega1) omega1 = omega              ! We keep the largest solid angle.
                      endif
c     ------------------------------------
c     solid angle for the central plane yz
c     ------------------------------------
                      if (x_c .ne. x_dif) then                            ! if the longitude of the observer voxel is the same as the
c                                                                         ! source voxel, we do not calculate the solid angle
c                                                                         ! for the yz plane because it is null.
                       call planyz(dy,xc,xn,yc,yn,zc,zn,
     +                 cthick,zcellc,r1x,r1y,r1z,r2x,r2y,r2z,
     +                 r3x,r3y,r3z,r4x,r4y,r4z)
                       call anglesolide(omega,r1x,r1y,r1z,                ! Calling the routine anglesolide for the calculation of the solid angle 
     +                 r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)               ! in the yz plane.
                      else 
                       omega=0.
                      endif
                      if (omega.gt.0.) then
                       if (omega .gt. omega1) omega1 = omega              ! We keep the largest solid angle.
                      endif
                      omega=omega1
c=======================================================================
c        computation of the scattered flux reaching the line of sight voxel
c=======================================================================
                      fldiff=idif1*omega*transm*
     +                transa*(1.-ff)*hh
                if (cloudt.ne.0) then                                     ! line of sight voxel = cloud
                  if (cloudh(cloudt).eq.zcellc) then
                     call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +               dx,dy,azencl)                                        ! zenith angle from cloud to observer                     
                     call cloudreflectance(angzen,cloudt,rcloud)          ! cloud intensity from direct illum
                     icloud=icloud+
     +               fldiff*rcloud*abs(cos(azencl))/pi
                  endif
                endif
c=======================================================================
c   computation of the scattering probability of the scattered light toward the observer voxel (exiting voxel_c)
c=======================================================================
                      if (angzen.lt.(pi/2.)) then                         ! Attribution of the initial and final limits of the path to the
c                                                                         ! scattering voxel.
                       zidif=zcdown
                       zfdif=zcup
                      else
                       zidif=zcup
                       zfdif=zcdown
                      endif
                      anaz=angazi
                      call transmitm(angzen,anaz,iun,iun,zidif,ideux,     ! Molecular transmittance of the scattering voxel.
     +                iun,zfdif,lambda,dx,dy,pressi,tran1m)
                      call transmita(angzen,anaz,iun,iun,zidif,ideux,     ! Aerosol transmittance of the scattering voxel.
     +                iun,zfdif,dx,dy,taua,tran1a)    
                      call angle3points (x_dif,y_dif,z_dif,x_c,y_c,       ! scattering angle.
     +                z_c,x_obs,y_obs,z_obs,dx,dy,angdif)
                      call diffusion(omega,angdif,tran1a,tran1m,          ! scattering probability of the direct light.
     +                secdif,fdifan,pdifd2)
c=======================================================================
c Computing scattered intensity toward the observer from the line of sight voxel
c=======================================================================
                      idiff2=fldiff*pdifd2
                      idiff2=
     +                idiff2*real(stepdi)                                 ! Correct the result for the skipping of 2nd scattering voxels to accelerate the calculation
                      itodif=                        
     +                itodif+idiff2
c                     endif                                               ! end condition obstacle scattering->line of sight
c        else
c        endif                                                             ! end condition shadow scattering-line of sight                     
c                    endif                                                ! end condition obstacle source->scattering.
c                   else
c                   endif                                                  ! end condition shadow source-scattering
                  endif                                                   ! end of the case scattering = Source or line of sight voxel
                 endif                                                    ! end of the condition "voxel of the domain".      
                enddo                                                     ! end of the loop over the scattering voxels.
               endif                                                      ! end of the condition ou effdif > dx.
c End of 2nd scattered intensity calculations    
c**********************************************************************
c        computation of the intensity coming from a source to the line of sight voxel toward the sensor
c**********************************************************************
               isourc=intdir+itotind+itodif+itotrd                        ! Sum of the intensities of each type of source  
c                                                                         ! reaching the line of sight voxel.
c                                                                         ! in the order 1st scat; refl->1st scat; 1st scat->2nd scat, refl->1st scat->2nd scat
               if (verbose.eq.1) then
                print*,' Total intensity components:'
                print*,' source->scattering=',intdir
                print*,' source->reflexion->scattering=',
     +          itotind
                print*,' source->scattering->scattering=',
     +          itodif
                print*,' source->reflexion->scattering->scattering=',
     a          itotrd  
               endif
c                   
c**********************************************************************
c        computation of the total intensity coming from all the sources of a given type
c**********************************************************************
               itotty=itotty
     +         +isourc                                                    ! Sum of the intensities of each source.
                                                                          ! ITT stores itotty in a matrix
               ITT(x_s,y_s,stype)=ITT(x_s,y_s,stype)+isourc

          endif                                                           ! end condition distance source-line of sight-observer <= dx/2
              endif                                                       ! end of the condition "the luminosity of the ground pixel x_s,y_s in not null".
             enddo                                                        ! end the loop over the lines (latitude) of the domain (y_s).
            enddo                                                         ! end the loop over the column (longitude) of the domain (x_s).
c
c   end of the computation of the intensity of one source type
            itotci=itotci+itotty                                          ! Sum of the intensities of each type to the line of sight voxel.
            do x_s=imin(stype),imax(stype)
             do y_s=jmin(stype),jmax(stype)
              ITC(x_s,y_s)=ITC(x_s,y_s)+ITT(x_s,y_s,stype)
             enddo   
            enddo  
c calculate lpluto 
            do x_s=1,nbx
             do y_s=1,nby
               lpluto(x_s,y_s)=lpluto(x_s,y_s)+   
     +         lamplu(x_s,y_s,stype)
             enddo
            enddo
           endif                                                          ! end of condition if there are any flux in that source type
           enddo                                                          ! end of the loop over the types of sources (stype).
c    end of the computation of the intensity coming from a line of sight voxel toward the sensor
c
c
c***********************************************************************
c        computation of the luminous flux reaching the observer voxel
c***********************************************************************
c
c=======================================================================
c        computation of the zenithal angle between the observer and the line of sight voxel
c=======================================================================
           call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,dx,dy,
     +     angzen)                                                        ! computation of the zenithal angle between the line of sight voxel and the observer.
c                                                                         ! end of the case "observer at the same latitu/longitude than the source".
c=======================================================================
c        computation of the transmittance between the line of sight voxel and the observer
c=======================================================================
           anaz=zero
           call transmitm(angzen,anaz,x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +     lambda,dx,dy,pressi,transm)
           call transmita(angzen,anaz,x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +     dx,dy,taua,transa)
c=======================================================================
c     computation of the Solid angle of the line of sight voxel seen fromthe observer
c=======================================================================
           xn=dble(x_obs)*dble(dx)                                        ! Position in meters of the observer voxel (longitude).
           yn=dble(y_obs)*dble(dy)                                        ! Position in meters of the observer voxel (latitu).
           zn=dble(z_obs)                                                 ! Position in meters of the observer voxel (altitude).
           xc=dble(x_c)*dble(dx)                                          ! Position in meters of the line of sight voxel (longitude).
           yc=dble(y_c)*dble(dy)                                          ! Position in meters of the line of sight voxel (latitu).
           zc=dble(z_c)                                                   ! Position in meters of the line of sight voxel (altitude).
c    ------------------------------------
c    solid angle for the central plane xy
c    ------------------------------------
           if (z_c .ne. z_obs) then
            call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)  
            call anglesolide(omega,r1x,r1y,r1z,                           ! Call of the routine anglesolide to compute the solid angle 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! along the surface xy.
            omega1 = omega
           else
            omega1=0.
           endif
c     ------------------------------------
c     solid angle for the central plane zx
c     ------------------------------------
           if (y_c .ne. y_obs) then                                       ! if the latitu of the observer voxel is the same as
c                                                                         ! of the source voxel, we are not computing the angle solide
c                                                                         ! for the surface zx car il is egal a 0.
            call planzx(dx,xc,xn,yc,yn,zc,zn,cthick,zcellc,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            call anglesolide(omega,r1x,r1y,r1z,                           ! Call of the routine anglesolide to compute the solid angle 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! along the surface zx.
           else
            omega=0.
           endif
           if (omega.gt.0.) then
            if (omega.gt.omega1) omega1 = omega                           ! We keep the largest solid angle.
           endif
c     ------------------------------------
c     solid angle for the central plane yz
c     ------------------------------------
           if (x_c.ne.x_obs) then                                         ! if the longitude of the observer voxel is is the same as
c                                                                         ! of the source voxel, we are not computing the angle solide
c                                                                         ! for the surface yz car il is egal a 0
            call planyz(dy,xc,xn,yc,yn,zc,zn,cthick,zcellc,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            call anglesolide(omega,r1x,r1y,r1z,                           ! Call of the routine anglesolide to compute the solid angle 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! along the surface yz.
           else 
            omega=0.
           endif
           if (omega.gt.0.) then
            if (omega.gt.omega1) omega1=omega                             ! On garof the solid angle le plus grand.
           endif
           omega=omega1
c=======================================================================
c        computation of the flux reaching the objective of the telescope from the line of sight voxel
c=======================================================================

           flcib=itotci*ometif*transa*transm                              ! computation of the flux reaching the intrument from the line of sight voxel

           do x_s=1,nbx
            do y_s=1,nby
             FC(x_s,y_s)=ITC(x_s,y_s)*ometif*transa*transm
            enddo
           enddo
           if (cos(pi-angzen).eq.0.) then 
            print*,'ERROR perfectly horizontal sight is forbidden!'
            stop
           else
             portio=omefov/omega                                          ! Fraction of the line of sight voxel covered by the fov (Fraction can be 
c                                                                         ! larger than 1). Le pi ici is du au fait
c                                                                         ! que angzen is calcule sur le path from the line of sight voxel toward the observer
          endif
           if (omega.eq.0.) then
            print*,'ERROR omega=0 (1)'
            stop
           endif
           fcapt=flcib*portio                                             ! correction for the FOV to the flux reaching the intrument from the line of sight voxel
           do x_s=1,nbx
            do y_s=1,nby
             FCA(x_s,y_s)=FC(x_s,y_s)*portio
            enddo
           enddo
c   end of the computation of the flux reaching the observer voxel from the line of sight voxel
           ftocap=ftocap+fcapt  
           do x_s=1,nbx
            do y_s=1,nby
             FTC(x_s,y_s)=FTC(x_s,y_s)+FCA(x_s,y_s)                       ! FTC is the array of the flux total at the sensor level permettant d'identifier
                                                                          ! the contribution of each voxel of the ground to the total flux at the observer level
                                                                          ! Le % is simplement donne par FTC/ftocap
             flcumu=flcumu+FCA(x_s,y_s)
            enddo
           enddo
          endif                                                           ! end of the condition line of sight voxel n'est pas observer voxel.
         endif                                                            ! end of the condition line of sight voxel inside the modelling domain
        endif                                                             ! end condition for continuing of a computation stopped.
c correction for the FOV to the flux reaching the intrument from the cloud voxel
           if (cloudt.ne.0) then
            if (cloudh(cloudt).eq.zcellc) then                            ! line of sight voxel = cloud
c=======================================================================
c  solid angle of the cloud pixel as seen from observer position
c=======================================================================
              xn=dble(x_obs)*dble(dx)                                     ! Position in meters of the observer voxel (longitude).
              yn=dble(y_obs)*dble(dy)                                     ! Position in meters of the observer voxel (latitu).
              zn=dble(z_obs)                                              ! Position in meters of the observer voxel (altitude).
              xc=dble(x_c)*dble(dx)                                       ! Position in meters of the line of sight voxel (longitude).
              yc=dble(y_c)*dble(dy)                                       ! Position in meters of the line of sight voxel (latitu).
              zc=dble(z_c)                                                ! Position in meters of the line of sight voxel (altitude).
c    ------------------------------------
c    solid angle for the central plane xy
c    ------------------------------------
              if (z_c .ne. z_obs) then
                 call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +           r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)  
                 call anglesolide(omega,r1x,r1y,r1z,                      ! Call of the routine anglesolide to compute the solid angle 
     +           r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                     ! along the surface xy.
              else
                 omega=0.
              endif
c computation of the flux reaching the intrument from the cloud voxel
              fcloud=icloud*ometif*transa*transm
              fccld=fcloud*omefov/omega
              fctcld=fctcld+fccld
            endif
           endif          
        print*,' Flux @ sensor (clear & cloudy)             =',fcapt,
     +  fccld
        print*,' Flux @ sensor accumulated (clear & cloudy) =',ftocap,
     +  fctcld
        write(2,*) ' Flux sensor (clear & cloudy)             =',fcapt
     +  ,fccld
        write(2,*) ' Flux sensor accumulated (clear & cloudy) =',ftocap
     +  ,fctcld   
      endif                                                               ! end condition line of sight voxel 1/stoplim
       enddo                                                              ! end of the loop over the line of sight voxels.
       if (prmaps.eq.1) then
          open(unit=9,file=pclf,status='unknown')
          open(unit=8,file=pcwf,status='unknown')
            fctnto=0.
            ftcmax=0.
            do x_s=1,nbx
               do y_s=1,nby
                  FTC(x_s,y_s)=FTC(x_s,y_s)/ftocap   
                  if (FTC(x_s,y_s).gt.ftcmax) ftcmax=FTC(x_s,y_s)
                  if (lpluto(x_s,y_s).ne.0.) then
                  fctnto=fctnto+FTC(x_s,y_s)/lpluto(x_s,y_s)
                  endif
               enddo
            enddo
            if (verbose.eq.1) then
               print*,'Writing normalized contribution matrix'
            endif
            do x_s=1,nbx
               do y_s=1,nby
                  if (lpluto(x_s,y_s).ne.0.) then
                     FTCN(x_s,y_s)=(FTC(x_s,y_s)/lpluto(x_s,y_s))
     +               /fctnto  
                  else 
                     FTCN(x_s,y_s)=0.
                  endif                                                   ! FTCN is le % par unite of luminosite of the voxel
                  write(9,*) x_s,y_s,FTC(x_s,y_s)                         ! emettrice au sol, c'est un % par unite of watt installes
                  write(8,*) x_s,y_s,FTCN(x_s,y_s)
               enddo
            enddo
            call 2dout(nbx,nby,pclimg,FTC)
            call 2dout(nbx,nby,pcwimg,FTCN)     
          close(unit=8)
          close(unit=9)
c creation of files gnuplot for the visualiser il faut betweenr gnuplot
c puis load 'fichier.gplot'
          open(unit=9,file=pclgp,status='unknown')
          open(unit=8,file=pcwgp,status='unknown')
            write(9,*) 'sand dgrid3d',nbx,',',nby
            write(9,*) 'sand hidden3d'
            write(9,*) 'sand pm3d'
            write(9,*) 'splot "'//basenm(1:lenbase)//'_pcl.txt"
     +      with dots'
            write(8,*) 'sand dgrid3d',nbx,',',nby
            write(8,*) 'sand hidden3d'
            write(8,*) 'sand pm3d'
            write(8,*) 'splot "'//basenm(1:lenbase)//'_pcw.txt"
     +      with dots'    
          close(unit=8)
          close(unit=9) 
       endif                                                              ! end of condition for creating contrib and sensit maps
          print*,'====================================================='
          print*,'          Total flux entering instrument (W)'
          write(*,2001) ftocap+fctcld  
          print*,'              Sky radiance (W/str/m**2)'       
          write(*,2001) (ftocap+fctcld)/(lfente*
     +          longfe/focal**2.)/(pi*(diamobj/2.)**2.)
       print*,'  '
       print*,' Interpolation flux error= ',
     +          ftocap-flcumu
       print*,'======================================================='
       write(2,*) '==================================================='
       write(2,*) '          Total flux entering instrument (W)'
       write(2,2001) ftocap+fctcld
        write(2,*) '            Sky radiance (W/str/m**2)          '      
       write(2,2001) (ftocap+fctcld)/(lfente*
     +          longfe/focal**2.)/(pi*(diamobj/2.)**2.)
       write(2,*) '  '                                                
       write(2,*) 'Interpolation flux errror= ',
     +          ftocap-flcumu
       write(2,*) '==================================================='
      close(2)
 2001 format('                   ',E10.3E2)
      stop
      end
c***********************************************************************************************************************
c*                                                                                                                     *
c*                                         end of the programme                                                            *
c*                                                                                                                     *
c***********************************************************************************************************************
