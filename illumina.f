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
c **    - Calculation of the atmospheric transmittance and 1st and 2nd order of scatte                            **
c **    - Lambertian reflexion on the ground                                                                          **
c **    - Terrain slope considered (apparent surface and shadows)                                                     **
c **    - Angular photometry of a lamp is considered uniform along the azimuth                                        **
c **    - Sub-grid obstacles considered (with the mean free path of light toward ground and mean obstacle height      **
c **      but these parameters are fixed on the modelling domain                                                      **
c **    - Molecules and aerosol optics (phase function, scatte probability, aerosol absorption)                   **  
c **    - Exponential concentrations vertical profile (H aerosol= 2km, H molecules= 8km                               **
c **    - Exponential vertical resolution (max height= 30 km)                                                         **
c **    - Accounting for hetecloudh(cloudt)rogeneity of ground reflectance, luminaires number, luminaire height, angular photometry **
c **    - Wavelength dependant                                                                                        **
c **    - Ignore the flux scattered by the voxel occupied by the observer (cellobs=cellcible)                         **
c **    - Do not support direct observation of a source                                                               ** 
c **    - Direct observation of the ground not implemented                                                            **
c **    - Not accounting for molecular absorption                                                                     **
c **    - Do not consider earth curvature (i.e. local/regional model)                                                 **
c **    - No clouds                                                                                                   **
c **                                                                                                                  **
c ** Theoritical equations by Martin Aube, CEGEP de Sherbrooke (in french)                                            **
c **      http://cegepsherbrooke.qc.ca/~aubema/index.php/Prof/IllumEn?action=download&upname=intensite_lumineuse.pdf  **
c **                                                                                                                  **
c **********************************************************************************************************************
c   
c    Copyright (C) 2012 Martin Aube
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
      parameter (width=1024,height=50)
      integer iun,ideux
      real pi,pix4
      integer verbose                                                     ! verbose = 1 to have more print out, 0 for silent
      parameter (pi=3.1415926)
      parameter (pix4=4.*pi)
      real cthick(50)                                                     ! Cell thickness array (meter)

      real cellh(50)                                                      ! Cell height array (meter)

      real flcumu                                                         ! Flux cumule en fonction du deplacement le long de la ligne de vise
      character*72 mnaf                                                   ! Fichier relief.
      character*72 reflf                                                  ! Fichier reflexion.
      character*72 diffil                                                 ! Fichier diffusion.
      character*72 outfile                                                ! Fichier de resultat.
      character*72 pclf,pcwf,pclgp,pcwgp                                  ! fichiers de matrice de % du flux total normalise par unite de flux 
      character*72 pclimg,pcwimg
                                                                          ! et par flux et wattage
      character*72 basenm                                                 ! Nom de base des fichiers.
      character*12 nom                                                    ! Nom de la variable 2d lue ou ecrite
      integer lenbase                                                     ! Longueur du nom de base de l'experimentation.
      real lambda,pressi,drefle                                           ! Longueur d'onde (nanometre), pressi atmospherique (kPa), zone 
c                                                                         ! de portee de la reflexion (metre).
      integer ntype                                                       ! Nombre de type de sources lumineuses consideres.
      real largx                                                          ! Largeur (valeur x) du domaine de l'environnement (metre).
      real largy                                                          ! Longueur (valeur y) du domaine de l'environnement (metre).
      integer nbx,nby                                                     ! Nombre de cellule dans le domaine.  
      real val2d(width,width)                                             ! Matrice temporaire pour les intrants 2d
      real altsol(width,width)                                            ! Altitude du sol (metre).
      real srefl(width,width)                                             ! Reflectance des surfaces.
      real Hmin                                                           ! Hauteur minimale du domaine.
      real xcell0                                                         ! Longitude sud-ouest du domaine.
      real ycell0                                                         ! latitu sud-ouest du domaine.
      real gain                                                           ! Coefficient multiplicatif des valeurs des fichiers source.
      real offset                                                         ! Constante d'addition des valeurs des fichiers source.
      integer valmax                                                      ! Maximum value of the output pgms
      integer stype                                                       ! Identification du type de source.
      character*72 pafile,lufile,alfile                                   ! Fichiers relatifs aux sources lumineuses (fonction d'emission 
c                                                                         ! des sources (lampadaires), luminosite (DIMENSIONS??), altitude (metre).    
      real lamplu(width,width,9)                                          ! Luminosite des sources dans chaque case (DIMENSION???).
      real lampal(width,width,9)                                         ! Altitude des sources lumineuses par rapport au sol (metre).
      real pval(181,9),pvalto,pvalno(181,9)                               ! Valeurs des fonctions angulaires d'emission (arbitraires, 
c                                                                         ! totale non-matrice, normalisees).
      real dtheta                                                         ! Increment d'angle de la fonction d'emission des sources 
      real dx,dy,dxp,dyp,pixsiz                                           ! Largeur d'une cellule (metre), nombre de cellule dans la portee 
c                                                                         ! de reflexion (cellule), largeur d'une cellule.
c                                                                         ! pour verif fichiers compatibles, taille des sous cellules pour 
c                                                                         ! l'angle solide de la portio de surface qui reflechit.
      integer boxx,boxy                                                   ! Portee de reflexion (cellule).
      real fdifa(181),fdifan(181)                                         ! Fonction de diffusion (arbitraire et normalisee) des aerosols.
      real extinc,scatte,anglea(181)                                      ! Sections efficaces d'extinc et de diffusion, angle de 
c                                                                         ! diffusion (degre).
      real secdif                                                         ! Contribution de la diffusion a l'extinc lumineuse.
      real inclix(width,width)                                            ! Inclinaison de la case en x (DIMENSIONS).
      real incliy(width,width)                                            ! Inclinaison de la case en y (DIMENSIONS).   
      integer x_obs,y_obs,zcello                                          ! Position de l'observateur (cellule).
      real z_obs                                                          ! Hauteur de l'observateur (metre).
      integer lcible(width,3)                                             ! Matrice des cellules cibles.
      integer ncible,icible                                               ! Nombre de cellules cibles, compteur de la boucle sur les cellules 
c                                                                         ! cibles.     
      integer x_c,y_c,zcellc                                              ! Position de la cellule cible (cellule).
      real z_c                                                            ! Hauteur de la cellule cible (metre).
      real zcup,zcdown                                                    ! Limites inferieure et superieure de la cellule cible.    
      integer dirck                                                       ! Verificateur de position de la source (cas source=cell cible).     
      integer x_s,y_s,x_sr,y_sr,x_dif,y_dif,zceldi                        ! Positions de la source, surface reflectrice, cellules diffusantes 
c                                                                         ! (cellule).
      real z_s,z_sr,z_dif                                                 ! Hauteur de la source, surface reflectrice, cellule diffusante (metre).
      real angzen,ouvang                                                  ! Angle zenithal entre deux cellules (radians) et angle d'ouverture 
c                                                                         ! du cone d'angle solide en degres.
      integer anglez                                                        
c                                                                         ! d'emission des lampadaires.      
      real P_dir,P_indir,P_dif1                                           ! Fonction d'emission (directe,indirecte,diffusee) des sources 
c                                                                         ! lumineuses.
      real transa,transm                                                  ! Transmittance entre deux cellules (aerosols,molecules).
      real tran1a,tran1m                                                  ! Transmittance a l'interieur d'une cellule (aerosols,molecules).
      real taua                                                           ! epaisseur optique des aerosols.
      real*8 xc,yc,zc,xn,yn,zn                                            ! Position (metre) des elements (arrivee, depart) pour le calcul 
c                                                                         ! de l'angle solide.  
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! Composantes des vecteurs utilises dans la routine angle solide.
      real omega,omega1,omega2                                            ! Angle solide couvert par une cellule vue d'une autre, angle de 
c                                                                         ! comparaison.
      real fldir                                                          ! Flux provenant d'une cellule source dans une cellule cible (watt).
      real flindi                                                         ! Flux provenant d'une cellule reflectrice dans une cellule cible (watt).
      real fldiff                                                         ! Flux provenant d'une cellule diffusante dans une cellule cible (watt).
      real zidif,zfdif                                                    ! Limites initiale et finale du parcours de diffusion dans une cellule.
      real angdif                                                         ! Angle de diffusion.
      real pdifdi,pdifin,pdifd1,pdifd2                                    ! Probabilite de diffusion (directe,indirecte,premiere et deuxieme 
c                                                                         ! diffusions.
      real intdir                                                         ! Contribution d'une source a l'intensite directe dirigee vers 
c                                                                         ! le capteur par une cellule cible.
      real intind                                                         ! Contribution d'une cellule reflectrice a l'intensite indirecte 
c                                                                         ! dirigee vers le capteur.
      real itotind                                                        ! Contribution totale d'une source a l'intensite indirecte dirigee 
c                                                                         ! vers le capteur.
      real idiff2                                                         ! Contribution d'une cellule diffusante a l'intensite diffusee 
c                                                                         ! dirigee vers le capteur.
      real itodif                                                         ! Contribution totale d'une source a l'intensite diffusee dirigee 
c                                                                         ! vers le capteur.
      real isourc                                                         ! Contribution totale d'une source a l'intensite dirigee par une 
c                                                                         ! cellule cible vers le capteur.
      real itotty                                                         ! Contribution totale d'un type de source a l'intensite dirigee par 
c                                                                         ! une cellule cible vers le capteur.
      real itotci                                                         ! Intensite totale dirigee par une cellule cible vers le capteur.
      real itotrd                                                         ! Intensite totale dirigee par une cellule vers le capteur apres 
c                                                                         ! reflexion et double diffusion.
      real irefdi                                                         ! Intensite dirigee par une cellule vers le capteur apres reflexion 
c                                                                         ! et double diffusion. 
      real flcib                                                          ! Flux parvenant dans la cellule observatrice par une cellule cible.
      real fcapt                                                          ! Flux parvenant dans la cellule observatrice par toutes les cellules 
c                                                                         ! cibles d'un niveau contenues dans le champ de vision du capteur.
      real ftocap                                                         ! Flux total parvenant dans la cellule receptrice.
      real haut                                                           ! Haut negatif indique que la surface est eclairee par en dessous
c                                                                         ! (on ne calcule alors pas).
      real epsilx,epsily                                                  ! Inclinaison de la surface reflectrice.
      real flrefl                                                         ! Flux atteignant une surface reflectrice (watts).
      real irefl,irefl1                                                   ! Intensite quittant une surface reflectrice en direction d'une 
c                                                                         ! cellule cible.  
      real effdif                                                         ! Distance autour des cellule source et cible qui seront considerees  
c                                                                         ! pour calculer la double diffusion.
      integer zondif(3000000,4)                                           ! Matrice des cellules diffusantes la 4e colonne represente la valeur 
c                                                                         ! entiere la plus proche de la distance (en metre) a la ligne de 
c                                                                         ! diffusion simple.
      integer ndiff,idi                                                   ! Nombre de cellules diffusantes, compteur de la boucle sur les cellules
c                                                                         ! diffusantes. 
      integer stepdi                                                      ! Saut de diffusion pour accelerer le calcul e.g. si =2 il fera un  
      integer redudi                                                      ! calcul sur deux flag indiquant si le rayon de diffusion a ete reduit
      integer nvis0                                                       ! valeur de depart pour le calcul le long de la ligne de visee. 
c                                                                         ! Par defaut cette valeur est de 1 mais elle peut etre plus grande  
c                                                                         ! lorsque l'on reprend un calcul precedent qui a ete arrete anormalement.
      real fldif1                                                         ! Flux atteignant une cellule diffusante.
      real idif1                                                          ! Intensite dirigee vers une cellule cible par une cellule diffusante.
      real portio                                                         ! portio de cellule observee dans le champ de vision du capteur.
      real dis_obs                                                        ! Distance d'observation entre la cible et l'observateur.
      real ometif                                                         ! Angle solide couvert par l'objectif du telescope vu de la cellule 
c                                                                         ! cible.
      real omefov                                                         ! Angle solide couvert sur le ciel vu par la fente.
      real lfente                                                         ! Largeur de la fente (ou de la surface du capteur)
      real longfe                                                         ! Longueur de la fente (ou de la surface du capteur)
      real focal                                                          ! Distance focale effective objectif (selon le Throughput de l'appareil)  
      real angvis,azim                                                    ! Angles d'observation du recepteur.
      real projap                                                         ! Fraction de la surface reflectrice vue par rapport a la direction 
c                                                                         ! normale. Utile pour le calcul de la reflectance lambertienne.
      real nbang                                                          ! pour le moyennage de la fonction d'emission
      real obsH,angmin                                                    ! Hauteur moyenne des obstacle sous maille, angle minimum en dessous 
c                                                                         ! duquel un faisceau ne peut etre propage en raison de l'obstacle 
c                                                                         ! sous-maille.
      integer naz,na 
      real ITT(width,width,9)                                             ! Intensite totale type en forme matricielle
      real ITC(width,width)                                               ! Intensite totale cible en forme matricielle
      real FC(width,width)                                                ! Flux cible en forme matricielle
      real FTC(width,width)                                               ! Fraction du Flux total au capteur en forme matricielle
      real FTCN(width,width)                                              ! Fraction Flux total au capteur normalise par unite de wattage au sol
      real FCA(width,width)                                               ! Flux capteur sous forme matricielle
      real lpluto(width,width)                                            ! Luminosite totale de la cellule de surface toutes lampes confondues
      real fctnto,ftcmax                                                  ! FTCN total pour tout le domaine pour toutes lampes
      character*3 lampno                                                  ! mot de trois lettres contenant le numero de lampe
      integer imin(9),imax(9),jmin(9),jmax(9),step(9)                     ! bornes x et y de la zone contenant un type de lampe
      real defval                                                         ! valeur ignoree dans l'interpolation
      real dat(1024,1024)                                                 ! matrice a interpoler
      integer autom,intype,ii,jj                                          ! switch manuel automatique pour l'interpolation; interpolation type
      real window                                                         ! interpolation diameter
      real zhoriz(360)                                                    ! horizon en rad sur 360 deg, la posi 1 = 0 deg et 360 = 359deg
      real angazi,d2                                                      ! angle azimutal entre deux points en rad, dist max pour lhorizon
      integer az                                                          ! azimut de l'horizon
      real latitu                                                         ! latitu approx du centre du domaine pour l'instant on utilises ycell0
      integer vistep                                                      ! line of sight step for low elevation angles vistep=ncells_along_sight/50
      integer prmaps                                                      ! frag to enable saving contribution and sensitivity maps
      integer cloudt                                                      ! cloud type 0=clear, 1=Thin Cirrus/Cirrostratus, 2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus, 4=Cumulus/Cumulonimbus, 5=Stratocumulus
      integer cloudh(5),cloudz                                            ! cloud base layer relative to the lower elevation 
      real rcloud                                                         ! cloud relfectance 
      real azencl                                                         ! zenith angle from cloud to observer
      real icloud                                                         ! cloud reflected intensity
      real fcloud                                                         ! flux reaching the intrument from the cloud cell
      real fccld                                                          ! correction for the FOV to the flux reaching the intrument from the cloud cell
      real fctcld                                                         ! total flux from cloud at the sensor level
      real dsco                                                           ! distance source-cible-observateur
      real dminlp                                                         ! minimum distance between the observer and a lamp (m)
      data cthick /0.5,0.6,0.72,0.86,1.04,1.26,1.52,1.84,2.22,            ! Epaisseur des niveaux.
     a 2.68,3.24,3.92,4.74,5.72,6.9,8.34,10.08,12.18,14.72,17.78,21.48,
     b 25.94,31.34,37.86,45.74,55.26,66.76,80.64,97.42,117.68,142.16,
     c 171.72,207.44,250.58,302.7,365.66,441.72,533.6,644.58,778.66,
     d 940.62,1136.26,1372.6,1658.1,2002.98,2419.6,2922.88,3530.84,
     e 4265.26,5152.44/
      data cellh /0.25,0.8,1.46,2.25,3.2,4.35,5.74,7.42,9.45,             ! Hauteur du centre de chaque niveau.
     a 11.9,14.86,18.44,22.77,28.,34.31,41.93,51.14,62.27,75.72,91.97,
     b 111.6,135.31,163.95,198.55,240.35,290.85,351.86,425.56,514.59,
     c 622.14,752.06,909.,1098.58,1327.59,1604.23,1938.41,2342.1,
     d 2829.76,3418.85,4130.47,4990.11,6028.55,7282.98,8798.33,
     e 10628.87,12840.16,15511.4,18738.26,22636.31,27345.16/
      data cloudh /44,44,40,33,33/                                        ! 9300.,9300.,4000.,1200.,1100.
      verbose=0
c  
c=======================================================================
c        Lecture du fichier d'entree (illumina.in)
c=======================================================================
      open(unit=1,file='illumina.in',status='old')
       read(1,*)
       read(1,*) basenm
       read(1,*) 
       read(1,*) dx,dy
       read(1,*)
       read(1,*) diffil
       read(1,*) 
       read(1,*) effdif,stepdi
       if (verbose.eq.1) then
         print*,'Rayon de diffusion=',effdif,'m   1 calcul sur ',
     a   stepdi
       endif
       read(1,*)
       read(1,*) lambda
       read(1,*) pressi
       read(1,*) taua
       read(1,*) ntype
       read(1,*) 
       read(1,*) drefle,obsH
       read(1,*)
       read(1,*) x_obs,y_obs,zcello,nvis0
       read(1,*)
       read(1,*) angvis,azim
       read(1,*)  
       read(1,*) lfente,longfe,focal,diamobj 
       read(1,*)
       read(1,*)
       read(1,*) cloudt  
       read(1,*) dminlp 
       if (verbose.eq.1) then
         print*,'Minimum distance to the source=',dminlp    
       endif
      close(1)
c
c  determiner la longueur du nom
c 
      lenbase=index(basenm,' ')-1  
      mnaf=basenm(1:lenbase)//'_topogra.pgm'                              ! Determiner les noms de fichiers d'entree et de sortie
      reflf=basenm(1:lenbase)//'_reflect.pgm' 
      outfile=basenm(1:lenbase)//'.out'  
      pclf=basenm(1:lenbase)//'_pcl.txt'
      pcwf=basenm(1:lenbase)//'_pcw.txt'
      pclimg=basenm(1:lenbase)//'_pcl.pgm'
      pcwimg=basenm(1:lenbase)//'_pcw.pgm'
      pclgp=basenm(1:lenbase)//'_pcl.gplot'
      pcwgp=basenm(1:lenbase)//'_pcw.gplot'    
c  conversion des angles de visee geographique vers angle mathematique
c  on presume que maintenant l'angle mis dans le fichier illumina.in
c  est en reference a la definition geographique
c  geographie, azim=0 au nord, 90 a l'est, 180 au sud etc
c  mathematiques, azim=0 a l'est, 90 au nord, 180 a l'ouest etc
      azim=90.-azim
      if (azim.lt.0.) azim=azim+360.
      if (azim.ge.360.) azim=azim-360.
c  ouverture en ecriture du fichier de sortie
      open(unit=2,file=outfile,status='unknown')      
       write(2,*) 'FICHIERS UTILISES'
       write(2,*) mnaf,reflf,diffil
       print*,'Longueur d''onde (nm):',lambda,
     +       ' epaisseur optique:',taua
       write(2,*) 'Longueur d''onde (nm):',lambda,
     +       ' epaisseur optique:',taua
       write(2,*) 'Zone de diffusion:',effdif,' m'
       print*,'Zone de diffusion:',effdif,' m'
       write(2,*) 'Saut de diffusion:',stepdi
       print*,'Saut de diffusion:',stepdi
       write(2,*) 'Portee de la reflexion:',drefle,' m'
       write(2,*) 'Position observateur(x,y,z)',x_obs,y_obs,zcello
       print*,'Position observateur(x,y,z)',x_obs,y_obs,zcello
       write(2,*) 'Angle d elevation:',angvis,' angle azim (horaire
     + p/r au nord)',azim     
       print*,'Angle d elevation:',angvis,' angle azim (anti-horair
     +e p/r a l est)',azim 
c=======================================================================
c        Initialisation des matrices et variables
c=======================================================================
       if (cloudt.eq.0) then
          cloudz=50
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
         do k=1,9
          lamplu(i,j,k)=0.
          lampal(i,j,k)=0.
          ITT(i,j,k)=0.
         enddo
        enddo
       enddo
       do i=1,181
        do j=1,9
         pval(i,j)=0.
         pvalno(i,j)=0.
        enddo
       enddo
       do i=1,181
        fdifa(i)=0.
        fdifan(i)=0.
        anglea(i)=0.
       enddo     
       do i=1,1024
        do j=1,3
         lcible(i,j)=1
        enddo
       enddo
       do i=1,3000000
        do j=1,4
         zondif(i,j)=1
        enddo
       enddo     
       redudi=1
       irefdi=0.
       angmin=0.
       vistep=1.
c***********************************************************************
c        Lecture des variables propres a l'environnement               *
c***********************************************************************
c=======================================================================
c  Lecture du fichier de relief
c=======================================================================
       nom='relief'
       call intrants2d(mnaf,altsol,nom,xcell0,ycell0,pixsiz,
     + nbx,nby)

       latitu=ycell0

       Hmin=3000000.
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
c                                                                         ! Recherche de la hauteur minimale.
         if (Hmin.gt.altsol(i,j)) Hmin=altsol(i,j)
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
       enddo 
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
         altsol(i,j)=altsol(i,j)-Hmin                                     ! Soustraction de la hauteur minimale du domaine.
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
       enddo
c=======================================================================
c Lecture du fichier de reflexion
c=======================================================================
       nom='reflexion'
       call intrants2d(reflf,srefl,nom,xcell0,ycell0,
     + pixsiz,nbx,nby)
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
         if (srefl(i,j).lt.0.) then                                       ! Recherche de reflectances negatives
           print*,'***Negative reflectance!, stopping execution'
           stop
         endif
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
       enddo
c========================================================================
c  Lecture des valeurs de P(theta), luminosites et positions des sources
c========================================================================
c
       dtheta=.017453293                                                  ! Un degre.
       do stype=1,ntype                                                   ! Debut de la boucle pour les 99 types de sources.
        imin(stype)=nbx
        jmin(stype)=nby
        imax(stype)=1
        jmax(stype)=1       
        pvalto=0.
c        if (stype.lt.10) then
c          lampno=char(48)//char(48+stype)
c        endif
        write(lampno, '(I3.3)' ) stype                                    ! support de 999 sources differentes (3 digits)
        pafile=basenm(1:lenbase)//'_fctem_'//lampno//'.dat'               ! Attribution du nom du fichier  fonction angulaire d'emission.
        lufile=basenm(1:lenbase)//'_lumlp_'//lampno//'.pgm'               ! Attribution du nom du fichier de la luminosite des cases.
        alfile=basenm(1:lenbase)//'_altlp_'//lampno//'.pgm'               ! Attribution du nom du fichier  hauteur des sources lumineuse.
        open(UNIT=1, FILE=pafile,status='OLD')                            ! Ouverture du fichier pa#.dat, fonction angulaire d'emission.
        do i=1,181                                                        ! Debut de la boucle pour les 181 donnees.
         read(1,*) pval(i,stype)                                          ! Lecture des donnees qui sont inscrites dans la matrice pval.
         pvalto=pvalto+pval(i,stype)*2.*pi*                               ! Sommation de la valeur tot  fonction d'emission prenormalisee 
     a   sin(real(i-1)*dtheta)*dtheta                                     ! (pvaleur x 2pi x sin theta x dtheta) (ou theta egale 
c                                                                         ! (i-1) x 1 degres).
        enddo                                                             ! Fin de la boucle sur les 181 donnees du fichier pa#.dat.
        close(1)                                                          ! Fermeture du fichier pa#.dat, fonction angulaire d'emission.
        do i=1,181
         pvalno(i,stype)=pval(i,stype)/pvalto                             ! Normalisation de la fonction d'emission.
        enddo                                                       
c    ===================================================================
        nom='luminosite'
        call intrants2d  (lufile,val2d,nom,xcell0,ycell0,pixsiz,
     +  nbx,nby)
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
         if (val2d(i,j).lt.0.) then                                       ! Recherche de luminosites negatives
           print*,'***Negative lamp luminosity!, stopping execution'
           stop
         endif
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
       enddo     
        do i=1,nbx                                                        ! recherche du plus petit rectangle encadrant la zone
         do j=1,nby                                                       ! de luminosite non nulle pour accelerer le calcul
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
 336    do i=1,nbx                                                        ! Debut de la boucle sur toutes les cases en x.
         do j=1,nby                                                       ! Debut de la boucle sur toutes les cases en y.
          lamplu(i,j,stype)=val2d(i,j)                                    ! remplir la matrice du type de lampe stype
         enddo                                                            ! Fin de la boucle sur toutes les cases en y.
        enddo                                                             ! Fin de la boucle sur toutes les cases en x.
c
c determination de la distance du centroide d un type de lampes (utiles lorsque une lampe est definie par 
c zone geographique regroupee comme une ile ou un ville
c cette fonction est maintenant commentee car nous utilisons une algorithme externe
c pour faire une selection statistique des points a calculer
c ainsi nous definirons la saut spatial a 1 pour toutes les zones 
c les lignes ci-dessous seront commentees en consequence
c        dzone=nint(sqrt(real(((imax(stype)+imin(stype))/2-x_obs)**2+
c     +  ((jmax(stype)+jmin(stype))/2-y_obs)**2)))
        step(stype)=1
c        if (dzone.gt.50) step(stype)=2                                   ! saut pour le calcul des cellules au sol ce saut augmente
c        if (dzone.gt.100) step(stype)=3                                  ! avec la distance
c        if (dzone.gt.150) step(stype)=4 
c        if (dzone.gt.200) step(stype)=5
c    ==================================================================
        nom='hauteurs'
        call intrants2d(alfile,val2d,nom,xcell0,ycell0,pixsiz,nbx,
     +  nby)
        do i=1,nbx                                                        ! Debut de la boucle sur toutes les cases en x.
         do j=1,nby                                                       ! Debut de la boucle sur toutes les cases en y.
          lampal(i,j,stype)=val2d(i,j)                                    ! Remplissage de la matrice pour la lampe stype
         enddo                                                            ! Fin de la boucle sur toutes les cases en y.
        enddo                                                             ! Fin de la boucle sur toutes les cases en x.
       enddo                                                              ! Fin de la boucle sur les 99 types de sources.    
c=======================================================================
c        Lecture des parametres de diffusion
c=======================================================================
       open(unit = 1, file = diffil,status= 'old')                        ! Ouverture du fichier contenant les parametres de diffusion.
c                                                                         ! Le fichier de diffusion est genere par le programme imies 
c                                                                         ! du progiciel AODSEM (Martin Aube).
        read(1,*)                                                           
        read(1,*)
        read(1,*)
        do i=1,181
         read(1,*) anglea(i), fdifa(i)                                    ! Lecture de la fonction de diffusion et de l'angle associe a 
c                                                                         ! cette fonction de 0 a 180 degres soit 181 lignes.
         fdifan(i)=fdifa(i)/pix4                                          ! Normalisation de la fonction a 4 pi (l'integrale de la 
c                                                                         ! fonction fournie sur tous les angles solides doit etre egale a 4 pi).
c                                                                         ! En effet le fichier .mie.out est normalisé ainsi (revefie par 
c                                                                         ! M. Aube en avril 2009)
        enddo
        do i = 1,7
         read(1,*)
        enddo
        read(1,*) extinc                                                  ! Lecture de la section efficace d'extinc des aerosols.
        read(1,*) scatte                                                  ! Lecture de la section efficace de diffusion des aerosols.
       close(1)
       secdif=scatte/extinc                                               ! Rapport (sigmadif/sigmatotal).
c======================================================================
c        Quelques operations preparatoires
c======================================================================
       dy=dx                                                              ! Pour le moment on considere que l'echelle est la meme sur les deux axes
       z_obs=cellh(zcello)                                                ! Attribution d'une valeur en metre a la position z de l'observateur.
       largx=dx*real(nbx)                                                 ! Calcul de la largeur en x d'une case.
       largy=dy*real(nby)                                                 ! Calcul de la largeur en y d'une case.
       boxx=nint(drefle/dx)                                               ! Nombre de colonnes a considerer a gauche/droite de la source 
c                                                                         ! pour la reflexion.
       boxy=nint(drefle/dy)                                               ! Nombre de colonnes a considerer en haut/bas de la source pour 
c                                                                         ! la reflexion.
       write(2,*) 'Largeur du domaine [NS](m):',largx,'#cases:',nbx
       write(2,*) 'Largeur du domaine [EO](m):',largy,'#cases:',nby
       write(2,*) 'Taille d''une cellule (m):',dx,' X ',dy
       write(2,*) 'latitu sud-ouest:',ycell0,' Longitude sud-ouest:',
     + xcell0
c=======================================================================
c        Calcul de l'inclinaison des cases en x et en y
c=======================================================================
       do i=1,nbx                                                         ! Debut de la boucle sur les colonnes (longitude) du domaine.
        do j=1,nby                                                        ! Debut de la boucle sur les ranges (latitu) du domaine.
         if (i.eq.1) then                                                 ! Cas particulier au bord du domaine (cote vertical gauche).
          inclix(i,j)=atan((altsol(i+1,j)-altsol(i,j))/real(dx))          ! Calcul de l'inclinaison en x de la surface.
         elseif (i.eq.nbx) then                                           ! Cas particulier au bord du domaine (cote vertical droit).
          inclix(i,j)=atan((altsol(i-1,j)-altsol(i,j))/(real(dx)))        ! Calcul de l'inclinaison en x de la surface.
         else
          inclix(i,j)=atan((altsol(i+1,j)-altsol(i-1,j))/(2.              ! Calcul de l'inclinaison en x de la surface.
     1    *real(dx)))
         endif
         if (j.eq.1) then                                                 ! Cas particulier au bord du domaine (cote horizontal en bas).
          incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j))/(real(dy)))        ! Calcul de l'inclinaison en y de la surface.
         elseif (j.eq.nby) then                                           ! Cas particulier au bord du domaine (cote horizontal en haut).
          incliy(i,j)=atan((altsol(i,j-1)-altsol(i,j))/(real(dy)))        ! Calcul de l'inclinaison en y de la surface.
         else
          incliy(i,j)=atan((altsol(i,j+1)-altsol(i,j-1))/(2.              ! Calcul de l'inclinaison en y de la surface.
     1    *real(dy)))
         endif
        enddo                                                             ! Fin de la boucle sur les ranges (latitu) du domaine.
       enddo                                                              ! Fin de la boucle sur les colonnes (longitude) du domaine.
c=======================================================================
c        Debut de la boucle sur les cellules cibles
c=======================================================================
       call lignevisee(x_obs,y_obs,z_obs,dx,dy,angvis,                    ! Determination de la ligne de visee (cellules cibles).
     + azim,nbx,nby,vistep,cloudz,lcible,ncible)
       fctcld=0.
       ftocap=0.                                                          ! Initialisation de la valeur du flux recu par le capteur.
       fcapt=1.
       do icible=1,ncible                                                 ! Debut de la boucle sur les cellules cibles.
      if ((fcapt.ge.ftocap/5000.).or.(cloudt.ne.0)) then                  ! arreter de calculer la ligne de visee lorsque l'increment est inferieur a 1/5000
        if (icible.ge.nvis0) then                                         ! Debut condition pour la reprise d'un calcul arrete   .
         itotci=0.                                                        ! Initialisation de la contribution d'une cible au capteur.
         do i=1,nbx
          do j=1,nby
            ITC(i,j)=0.
          enddo
         enddo
         zcellc=lcible(icible,3)                                          ! Definition de la position (cellule) verticale de la cible.
         z_c=cellh(zcellc)                                                ! Definition de la position (metre) verticale de la cible.
         y_c=lcible(icible,2)                                             ! Definition de la position (cellule) de la cible.
         x_c=lcible(icible,1)                                             ! Definition de la position (cellule) de la cible.
         print*,'=================================================='
         print*,' Avancement le long de la ligne de visee :',
     +   icible,'/',ncible,'(',x_c,',',y_c,')'
         print*,' Hauteur de la cible   =',z_c,' m'
         print*,' Epaisseur de la cible =',cthick(zcellc),' m'
         write(2,*) '=================================================='
         write(2,*) ' Avancement le long de la ligne de visee :',
     +   icible,'/',ncible,'(',x_c,',',y_c,')'
         write(2,*) ' Hauteur de la cible   =',z_c,' m'
      write(2,*) ' Epaisseur de la cible =',cthick(zcellc),' m'
         if( (x_c.gt.nbx).or.(x_c.lt.1).or.(y_c.gt.nby).or.(y_c.lt.1)     ! Condition Cellule cible dans le domaine de calcul.
     +      .or.(zcellc.gt.50).or.(zcellc.lt.1) )then
         else
          if((x_c.eq.x_obs).and.(y_c.eq.y_obs).and.                       ! Pour le moment, si la cellule cible est la cellule observatrice, 
c                                                                         ! on ne calcule pas le flux diffuse.
     +    (zcellc.eq.zcello))then
           if (verbose.eq.1) then
             print*,'Cellule Cible = Cellule Observateur'
           endif
          else
           zcdown=z_c-0.5*cthick(zcellc)                                  ! Limite inferieure de la cellule cible.
           zcup=z_c+0.5*cthick(zcellc)                                    ! Limite superieure de la cellule cible.
c=======================================================================
c        Debut de la boucle sur les types de sources lumineuses
c=======================================================================
           do stype=1,ntype                                               ! Debut de la boucle sur les type de source.
            print*,' Turning on lamp',stype
            write(2,*) ' Turning on lamp',stype
            itotty=0.                                                     ! Initialisation de la contribution d'un type de source a 
c                                                                         ! l'intensite dirigee vers le capteur par une cellule cible.
            do x_s=1,nbx
             do y_s=1,nby
              ITT(x_s,y_s,stype)=0.
             enddo
            enddo     
            do x_s=imin(stype),imax(stype),step(stype)                    ! Debut de la boucle sur les colonnes (longitudes) du domaine.
             do y_s=jmin(stype),jmax(stype),step(stype)                   ! Debut de la boucle sur les rangees (latitud) du domaine.
              if (lamplu(x_s,y_s,stype) .ne. 0.) then                     ! Si la luminosite d'une case est nulle, le programme ignore cette case.
               z_s=(altsol(x_s,y_s)+lampal(x_s,y_s,stype))                ! Definition de la position (metre) verticale de la source.








c calcul de la distance source-cible-observateur si cette distance est inferieure a dx/2, pas de calcul effectue
c la raison est que autrement on passe par des cellules tres proches de la source et on est jamais dans de telles
c conditions lorsqu'on observe le ciel. C est un probleme cree par le fait que les sources et l observateur
c sont toujours consideres au centre des cellules.
               dsco=sqrt((real(x_s-x_c)*dx)**2.+(real(y_s-y_c)*dx)**2.+
     +         (z_s-z_c)**2.)+sqrt((real(x_obs-x_c)*dx)**2.+(real(y_obs
     +         -y_c)*dx)**2.+(z_obs-z_c)**2.)
          if (dsco.ge.dminlp) then                                    ! debut condition distance source-cible-observateur >= dx/2









c **********************************************************************************************************************
c *     Calcul de l'intensite directe dirigee vers le capteur par une cellule cible en provenance d'une source         *
c **********************************************************************************************************************         
               dirck=0                                                    ! Initialisation de la verification de la position de la source.
               if ( (x_s.eq.x_c).and.(y_s.eq.y_c).and.( abs(z_s-z_c)      ! Si les positions x et y de la source et de la cible sont les 
c                                                                         ! memes alors.
     +         .lt.(cthick(zcellc)/2.) ) )then
                dirck=1
                if (verbose.eq.1) then
                 print*,'Source dans Cellule cible'
                endif
               endif                                                      ! Fin du cas positions x et y source et cible identiques.
               if (dirck.ne.1) then                                       ! Cas ou la source n'est pas dans la cellule cible.
c=======================================================================
c        Calcul de l'angle zenithal entre la source et la cible
c=======================================================================
c
c calcul de la ligne d'horizon pour les ombrages resolus direct           ! Il y a une ligne d'horizon par cellule cible a une resolution de 1 deg
         
                d2=sqrt((real(x_s-x_c)*dx)**2.+(real(y_s-y_c)*dy)**2.)    ! dist max pour l'horizon (i.e. l horizon passe la source ne compte pas
                call horizon(x_s,y_s,z_s,d2,altsol,nbx,nby,dx,dy,
     +          zhoriz,latitu)
                call anglezenithal
     +          (x_s,y_s,z_s,x_c,y_c,z_c,dx,dy,angzen)                    ! Calcul de l'angle zenithal entre la source et la cellule cible.
                call angleazimutal(x_s,y_s,x_c,y_c,dx,dy,angazi)          ! calcul de l'angle azimutal direct cible-source
                az=nint(angazi*180./pi)+1
                if ((angzen).lt.zhoriz(az)) then                          ! la ligne cible-source n'est pas sous l'horizon => on calcule
c                                                                         ! debut condition sous l'horizon direct
c obstacle sous maille             
                 angmin=pi/2.-atan((altsol(x_s,y_s)+obsH-z_s)/
     +           drefle)
                 if (angzen.lt.angmin) then                               ! Debut condition obstacle sous maille direct.
c
c=======================================================================
c        Calcul de la transmittance entre la source et la cible
c=======================================================================
                  call transmitm (angzen,x_s,y_s,z_s,x_c,y_c,z_c,
     +            lambda,dx,dy,pressi,transm)     
                  call transmita (angzen,x_s,y_s,z_s,x_c,y_c,z_c,
     +            dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cible vue de la source
c=======================================================================

c                omega2=0.

                  xc=dble(x_c)*dble(dx)                                   ! Position en metres de la cellule observatrice (longitude).
                  yc=dble(y_c)*dble(dy)                                   ! Position en metres de la cellule observatrice (latitu).
                  zc=dble(z_c)                                            ! Position en metres de la cellule observatrice (altitude).
                  xn=dble(x_s)*dble(dx)                                   ! Position en metres de la source (longitude).
                  yn=dble(y_s)*dble(dy)                                   ! Position en metres de la source (latitu).
                  zn=dble(z_s)                                            ! Position en metres de la source (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                  if (z_c .ne. z_s) then
                   call planxy(dx,dy,xc,xn,yc,yn,zc,zn,cthick,
     +             zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,
     +             r4z) 
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Appel de la routine anglesolide qui calcule l'angle solide 
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                   ! selon le plan xy.
                   omega1 = omega
                  else
                   omega1=0.
                  endif


c                  omega2=omega2+omega1


c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
                  if (y_c .ne. y_s) then                                  ! Si la latitu de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0
                   call planzx(dx,xc,xn,yc,yn,zc,zn,cthick,
     +             zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y
     +             ,r4z)
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Appel de la routine anglesolide qui calcule l'angle solide 
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                   ! selon le plan zx.
                  else
                   omega=0.
                  endif
                  if (omega.gt.0.) then
                   if (omega .gt. omega1) omega1 = omega                  ! On garde l'angle solide le plus grand jusqu'a present.
                  endif

c                  omega2=omega2+omega1


c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
                  if (x_c .ne. x_s) then                                  ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0.
                   call planyz(dy,xc,xn,yc,yn,zc,zn,cthick,
     +             zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,
     +             r4z)
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Routine anglesolide qui calcule l'angle solide selon le plan yz.
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                  else 
                   omega=0.
                  endif
                  if (omega.gt.0.) then
                   if (omega .gt. omega1) omega1 = omega                  ! On garde l'angle solide le plus grand
                  endif
                  omega=omega1


c                  omega2=omega2+omega1
c                  omega=omega2


c=======================================================================
c    Estimation du demi angle directeur de l'angle solide                 ! Cet angle servira a obtenir un meilleur estime (moyenne) de 
c                                                                         ! P_dir pour le cas de grans angles solides ou pvalno 
c=======================================================================  ! varie significativement sur +- ouvang.
                  ouvang=sqrt(omega/pi)                                   ! Angle en radian.
                  ouvang=ouvang*180./pi                                   ! Angle en degres.
c   
c=======================================================================
c        Calcul de la fonction d'emission de la source vers la cible
c=======================================================================
c   
                  anglez=nint(angzen/pi*180.)
                  if (anglez.lt.0) anglez=-anglez
                  if (anglez.gt.180) anglez=360-anglez
                  anglez=anglez+1                                         ! Transformer l'angle en deg entier en position dans la matrice.
c
c  moyenner sur +- ouvang	
c
                  naz=0

                  nbang=0.
                  P_dir=0.
                  do na=-nint(ouvang),nint(ouvang)
                   naz=anglez+na
                   if (naz.lt.0) naz=-naz
                   if (naz.gt.181) naz=362-naz                            ! La fonction est symetrique.
                   P_dir=P_dir+pvalno(naz,stype)
                   nbang=nbang+1.
                  enddo
                  P_dir=P_dir/nbang
c
c=======================================================================
c        Calcul du flux direct atteignant la cellule cible
c=======================================================================
                  fldir=lamplu(x_s,y_s,stype)*P_dir*omega*
     1            transm*transa
c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere directe
c=======================================================================
                  if (angzen.lt.(pi/2.)) then                             ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                   zidif=zcdown
                   zfdif=zcup
                  else
                   zidif=zcup
                   zfdif=zcdown
                  endif
                  call transmitm (angzen,iun,iun,zidif,ideux,ideux,       ! Transmittance moleculaire a l'interieur de la cellule diffusante.
     +            zfdif,lambda,dx,dy,pressi,tran1m)
                  call transmita (angzen,iun,iun,zidif,ideux,ideux,       ! Transmittance aerosols a l'interieur de la cellule diffusante.
     +            zfdif,dx,dy,taua,tran1a)
                  call angle3points (x_s,y_s,z_s,x_c,y_c,z_c,x_obs,       ! Angle de diffusion.
     +            y_obs,z_obs,dx,dy,angdif)
                  call diffusion(omega,angdif,tran1a,tran1m,              ! Probabilite de diffusion de la lumiere directe.     
     +            secdif,fdifan,pdifdi)
c=======================================================================
c   Calcul de la contribution d'une source a l'intensite directe dirigee vers le capteur par une cellule cible
c=======================================================================
                  intdir=fldir*pdifdi








c                print*,'balise1',cloudh(cloudt),zcellc

                if ((cloudt.ne.0).and.(cloudh(cloudt).eq.zcellc)) then    ! target cell = cloud
                     call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +               dx,dy,azencl)                                        ! zenith angle from cloud to observer                     
                     call cloudreflectance(angzen,cloudt,rcloud)          ! cloud intensity from direct illum
                     icloud=icloud+
     +               fldir*rcloud*abs(cos(azencl))/pi
c        print*,'balise1b',icloud
                endif













                 else 
                  intdir=0.                                      
                 endif                                                    ! Fin condition obstacle sous maille direct.
                else
                endif                                                     ! Fin condition sous l'horizon direct? 
               endif                                                      ! Fin du cas Position Source n'egale pas Position Cible.


c        print*,'deb3',intdir,fldir,pdifdi,
c     +  fldir,lamplu(x_s,y_s,stype),P_dir,omega,
c     +  transm,transa
c        print*,omega,angdif,tran1a,tran1m,secdif,pdifdi    
c        print*,x_s,y_s,z_s,x_c,y_c,z_c,x_obs,y_obs,z_obs 
c ok=lamplu,fldir
c        if (icible.eq.23) stop



c  fin du calcul de l'intensite directe
c **********************************************************************************************************************
c * Calcul de l'intensite indirecte dirigee vers le capteur par une cellule cible en provenance d'une source           *
c **********************************************************************************************************************
c=======================================================================
c        etablissement des conditions et des boucles
c=======================================================================
               itotind=0.                                                 ! Initialisation de l'intensite indirecte d'une source a une cible.    
               itotrd=0.
               do x_sr=x_s-boxx,x_s+boxx                                  ! Debut de la boucle sur les colonnes (longitude) reflectrices.
                do y_sr=y_s-boxy,y_s+boxy                                 ! Debut de la boucle sur les ranges (latitu) relfectrices.
                 irefl=0.
                 z_sr=altsol(x_sr,y_sr)   
                  if( (x_sr.gt.nbx).or.(x_sr.lt.1).or.(y_sr.gt.nby)
     +            .or.(y_sr.lt.1) )then  
                   if (verbose.eq.1) then
                    print*,'Cell Reflectrice a l''exterieur du domaine'
                   endif
                  else  
                   if((x_s.eq.x_sr).and.(y_s.eq.y_sr).and.(z_s.eq.z_sr))
     +             then
                    if (verbose.eq.1) then
                     print*,'Pos Source = Pos Cell Reflectrice'
                    endif
                   else
                    if (srefl(x_sr,y_sr).ne.0.) then                      ! Condition: la surface reflectrice n'a pas une reflectance nulle.
                     haut=-real(x_s-x_sr)*dx*tan(inclix(x_sr,y_sr))       ! Si la variable haut a une valeur negative, cela signifie que la surface
     1               -real(y_s-y_sr)*dy                                   ! reflectrice est eclairee par en-dessous.
     2               *tan(incliy(x_sr,y_sr))+z_s-z_sr
                     if (haut .gt. 0.) then                               ! Condition: la surface est eclairee par au-dessus.
c=======================================================================
c        Calcul de l'angle zenithal entre la source et la surface reflectrice
c=======================================================================
                      call anglezenithal(x_s,y_s,z_s,x_sr,y_sr,z_sr,dx,   ! Calcul de l'angle zenithal entre la source et la cellule cible.
     +                dy,angzen)                                          ! Fin du cas "observateur a la meme latitu/longitude que la source".

c=======================================================================
c        Calcul de la transmittance entre la source et la surface reflectrice
c=======================================================================
                      call transmitm (angzen,x_s,y_s,z_s,x_sr,y_sr,
     +                z_sr,lambda,dx,dy,pressi,transm)          
                      call transmita (angzen,x_s,y_s,z_s,x_sr,y_sr,
     +                z_sr,dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cellule reflectrice vue de la source
c=======================================================================
                      xc=dble(x_sr)*dble(dx)                              ! Position en metres de la cellule observatrice (longitude).
                      yc=dble(y_sr)*dble(dy)                              ! Position en metres de la cellule observatrice (latitu).
                      zc=dble(z_sr)                                       ! Position en metres de la cellule observatrice (altitude).
                      xn=dble(x_s)*dble(dx)                               ! Position en metres de la source (longitude).
                      yn=dble(y_s)*dble(dy)                               ! Position en metres de la source (latitu).
                      zn=dble(z_s)                                        ! Position en metres de la source (altitude).
                      epsilx=inclix(x_sr,y_sr)                            ! Inclinaison en x de la surface reflectrice.
                      epsily=incliy(x_sr,y_sr)                            ! Inclinaison en x de la surface reflectrice.
                      if (dx.gt.drefle*2.) then                           ! Utiliser une surface sous-grille lors que la portee de la 
c                                                                         ! reflexion est inferieure a la taille de cellule.
                       if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then
                        dxp=drefle*2.
                       endif
                      else
                       dxp=dx
                      endif
                      if (dy.gt.drefle*2.) then
                       if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then         
                        dyp=drefle*2.
                       endif
                      else
                       dyp=dy
                      endif              
                      r1x=xc-dble(dxp)/2.-xn                              ! Calcul de la composante en x du premier vecteur.
                      r1y=yc+dble(dyp)/2.-yn                              ! Calcul de la composante en y du premier vecteur.
                      r1z=zc-tan(dble(epsilx))*dble(dxp)/2.+tan(dble(
     +                epsily))*dble(dyp)/2.-zn                            ! Calcul de la composante en z du premier vecteur.
                      r2x=xc+dble(dxp)/2.-xn                              ! Calcul de la composante en x du deuxieme vecteur.
                      r2y=yc+dble(dyp)/2.-yn                              ! Calcul de la composante en y du deuxieme vecteur.
                      r2z=zc+tan(dble(epsilx))*dble(dxp)/2.+tan(dble(
     +                epsily))*dble(dyp)/2.-zn                            ! Calcul de la composante en z du deuxieme vecteur.
                      r3x=xc-dble(dxp)/2.-xn                              ! Calcul de la composante en x du troisieme vecteur.
                      r3y=yc-dble(dyp)/2.-yn                              ! Calcul de la composante en y du troisieme vecteur.
                      r3z=zc-tan(dble(epsilx))*dble(dxp)/2.-tan(
     +                dble(epsily))*dble(dyp)/2.-zn                       ! Calcul de la composante en z du troisieme vecteur.
                      r4x=xc+dble(dxp)/2.-xn                              ! Calcul de la composante en x du quatrieme vecteur.
                      r4y=yc-dble(dyp)/2.-yn                              ! Calcul de la composante en y du quatrieme vecteur.
                      r4z=zc+tan(dble(epsilx))*dble(dxp)/2.-tan(
     +                dble(epsily))*dble(dyp)/2.-zn                       ! Calcul de la composante en z du quatrieme vecteur.
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Appel de la routine anglesolide qui calcule l'angle solide.
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z) 
c
c=======================================================================
c    Estimation du demi angle directeur de l'angle solide                 ! Cet angle servira a obtenir un meilleur estime (moyenne) de 
c                                                                         ! P_dir pour le cas de grans angles solides ou pvalno
c=======================================================================  ! varie significativement sur +- ouvang.
                      ouvang=sqrt(omega/pi)                               ! Angle en radian.
                      ouvang=ouvang*180./pi                               ! Angle en degres.
c   
c=======================================================================
c        Calcul de la fonction d'emission du lampadaire vers la surface reflectrice
c=======================================================================
c    
                      anglez=nint(angzen/pi*180.)
                      if (anglez.lt.0) anglez=-anglez
                      if (anglez.gt.180) anglez=360-anglez
                      anglez=anglez+1                                     ! Transformer l'angle en deg entier en position dans la matrice.
c
c  moyenner sur +- ouvang	
c
c
                      nbang=0.
                      P_indir=0.
                      do na=-nint(ouvang),nint(ouvang)
                       naz=anglez+na
                       if (naz.lt.0) naz=-naz
                       if (naz.gt.181) naz=362-naz                        ! La fonction est symetrique.
                       P_indir=P_indir+pvalno(naz,stype)
                       nbang=nbang+1.
                      enddo
                      P_indir=P_indir/nbang
c 
c=======================================================================
c        Calcul du flux atteignant la cellule reflectrice
c=======================================================================
                      flrefl=lamplu(x_s,y_s,stype)*P_indir*
     a                omega*transm*transa
c=======================================================================
c        Calcul de l'intensite reflechie quittant la surface reflectrice
c=======================================================================
                      irefl1=flrefl*srefl(x_sr,y_sr)/pi                   ! Le facteur 1/pi vient de la normalisation de la fonction 
                      if (effdif.gt.(dx+dy)/2.) then 

c                  print*,'balise3',cloudt,cloudh,icloud

                       call reflexdbledif (x_sr,y_sr,z_sr,x_c,y_c,
     +                 zcellc,dx,dy,effdif,nbx,nby,stepdi,
     +                 irefl1,lambda,pressi,taua,zcup,
     +                 zcdown,secdif,fdifan,x_obs,y_obs,z_obs,
     +                 epsilx,epsily,irefdi,
     +                 drefle,obsH,altsol,latitu,cloudt,cloudh,icloud)


c                 print*,'balise3b',icloud

                      endif
                      itotrd=itotrd+irefdi      
c
c  la projection apparente est calculee a partir du produit scalaire du vecteur normal a 
c  la surface reflechissante et la ligne surface reflechissante vers cellule diffusante ou cible
c  it is the cosine correction for the lambertian reflectance for finite elements
c         
                      projap=(-tan(epsilx)*real(x_c-x_sr)*dx-
     +                tan(epsily)*real(y_c-y_sr)*dy+1.*(cellh(
     +                zcellc)-z_sr))/(sqrt(tan(epsilx)**2.+tan(epsily)
     +                **2.+1.)*sqrt((real(x_c-x_sr)*dx)**2.+(real(y_c-
     +                y_sr)*dy)**2.+(cellh(zcellc)-z_sr)**2.))
c                                                                         ! Peu importe la direction c'est la valeur absolue du cos theta 
c                                                                         ! qui compte.   
c verifier s'il y a ombrage entre sr et cible ici 
                 d2=sqrt((real(x_sr-x_c)*dx)**2.+(real(y_sr-y_c)*         ! dist max pour l'horizon (i.e. l horizon passe la source ne compte pas
     +           dy)**2.)                 
                 call horizon(x_sr,y_sr,z_sr,d2,altsol,nbx,nby,dx,dy,
     +           zhoriz,latitu) 
                 call anglezenithal(x_sr,y_sr,z_sr,x_c,y_c,z_c,dx,        ! Angle zenithal entre la surface reflechissante et la cellule cible.
     +           dy,angzen)     
                 call angleazimutal(x_sr,y_sr,x_c,y_c,dx,dy,angazi)       ! calcul de l'angle azimutal reflect-cible
                 az=nint(angazi*180./pi)+1          
                 if ((angzen).lt.zhoriz(az)) then                         ! la ligne cible-reflec n'est pas sous l'horizon => on calcule
                 
                 
                      if (projap.lt.0.) projap=0.
                      irefl=irefl1*
     +                projap
                 else
c                  print*,'ombrage reflexion',x_sr,y_sr,z_sr,x_c,y_c,z_c
                 endif                                                    ! fin condition surf. reflectrice au-dessus horizon
c=======================================================================
c        Cas Position Cible = Position Cellule reflectrice
c=======================================================================
                      if((x_c.eq.x_sr).and.(y_c.eq.y_sr).and.
     +                (z_c.eq.z_sr)) then
                       intind=irefl
                      else
c
c            
c obstacle                 
                       angmin=pi/2.-atan(obsH/drefle)
                       if (angzen.lt.angmin) then                         ! Debut condition obstacle indirect.
c
c=======================================================================
c        Calcul de la transmittance entre la surface reflectrice et la cellule cible
c=======================================================================
                        call transmitm (angzen,x_sr,y_sr,z_sr,x_c,
     +                  y_c,z_c,lambda,dx,dy,pressi,transm)        
                        call transmita (angzen,x_sr,y_sr,z_sr,x_c,
     +                  y_c,z_c,dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cible vue de la cellule reflectrice
c=======================================================================


c                omega2=0.

                        xc=dble(x_c)*dble(dx)                             ! Position en metres de la cellule observatrice (longitude).
                        yc=dble(y_c)*dble(dy)                             ! Position en metres de la cellule observatrice (latitu).
                        zc=dble(z_c)                                      ! Position en metres de la cellule observatrice (altitude).
                        xn=dble(x_sr)*dble(dx)                            ! Position en metres de la source (longitude).
                        yn=dble(y_sr)*dble(dy)                            ! Position en metres de la source (latitu).
                        zn=dble(z_sr)                                     ! Position en metres de la source (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                        if (z_c .ne. z_sr) then
                         call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                   cthick,zcellc,r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! selon le plan xy.
                         omega1 = omega
                        else
                         omega1=0.
                        endif


c           omega2=omega2+omega1



c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
                        if (y_c .ne. y_sr) then                           ! Si la latitu de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
                         call planzx(dx,xc,xn,yc,yn,zc,zn,
     +                   cthick,zcellc,r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! selon le plan zx.
                        else
                         omega=0.
                        endif
                        if (omega.gt.0.) then
                         if (omega.gt.omega1) omega1 = omega              ! On garde l'angle solide le plus grand jusqu'a present.
                        endif


c           omega2=omega2+omega1



c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
                        if (x_c.ne.x_sr) then                             ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0.
                         call planyz(dy,xc,xn,yc,yn,zc,zn,
     +                   cthick,zcellc,r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! selon le plan yz.
                        else 
                         omega=0.
                        endif
                        if (omega.gt.0.) then
                         if (omega.gt.omega1) omega1=omega                ! On garde l'angle solide le plus grand.
                        endif
                        omega=omega1 


c           omega2=omega2+omega1
c           omega=omega2

    
c=======================================================================
c        Calcul du flux indirect atteignant la cellule cible
c=======================================================================
                        flindi=irefl*omega*transm*
     +                  transa     
















c                print*,'balise2',cloudh(cloudt),zcellc

                if ((cloudt.ne.0).and.(cloudh(cloudt).eq.zcellc)) then    ! target cell = cloud
                     call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +               dx,dy,azencl)                                        ! zenith angle from cloud to observer                     
                     call cloudreflectance(angzen,cloudt,rcloud)          ! cloud intensity from indirect illum
                     icloud=icloud+
     +               flindi*rcloud*abs(cos(azencl))/pi
c                  print*,'balise2',icloud
                endif























 
c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere indirecte
c=======================================================================
                        if (angzen.lt.(pi/2.)) then                       ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                         zidif=zcdown
                         zfdif=zcup
                        else
                         zidif=zcup
                         zfdif=zcdown
                        endif        
                        call transmitm (angzen,iun,iun,zidif,ideux,       ! Transmittance moleculaire a l'interieur de la cell diffusante.
     +                  ideux,zfdif,lambda,dx,dy,pressi,tran1m)
                        call transmita (angzen,iun,iun,zidif,ideux,       ! Transmittance aerosols a l'interieur de la cell diffusante.
     +                  ideux,zfdif,dx,dy,taua,tran1a)
                        call angle3points (x_sr,y_sr,z_sr,x_c,y_c,z_c,    ! Angle de diffusion.
     +                  x_obs,y_obs,z_obs,dx,dy,angdif)
                        call diffusion(omega,angdif,tran1a,               ! Probabilite de diffusion de la lumiere indirecte.
     +                  tran1m,secdif,fdifan,pdifin)
c=======================================================================
c   Calcul de l'intensite indirecte dirigee vers le capteur par une cellule reflectrice
c=======================================================================
                        intind=flindi*pdifin 
                       else
                        intind=0.
                       endif                                              ! Fin condition obstacle indirect.
                      endif                                               ! Fin du cas Posi Cell Reflectrice = Position Cible.                                 
                      itotind=
     a                itotind+intind                                      ! Somme des intensites de chaque cell reflectrices propres 
c                                                                         ! a une source.
                     endif                                                ! Fin de la condition surface non-eclairee par le haut.
                    endif                                                 ! Fin de la condition reflectance non-nulle.
                   endif                                                  ! Fin de la condition cellule reflectrice n'est pas une source.
                  endif                                                   ! Fin de la condition surface a l'interieur du domaine.
                enddo                                                     ! Fin de la boucle sur les ranges (latitu) relfectrices.
               enddo                                                      ! Fin de la boucle sur les colonnes (longitude) reflectrices.
c   fin du calcul de l'intensite indirecte        
c **********************************************************************************************************************
c * Calcul de l'intensite diffusee dirigee vers le capteur par une cellule cible en provenance d'une source            *
c **********************************************************************************************************************
c
c=======================================================================
c    Determination des cellules diffusantes en fonction de la cellule source et de la cellule cible
c=======================================================================

               itodif=0.                                                  ! Initialisation de l'intensite diffusee par une source dans 
c                                                                         ! une cellule cible calculer le double diffusion seulement si 
               if (effdif.gt.(dx+dy)/2.) then                             ! le rayon de diffusion est superieur a la taille des cellules.
                call zone_diffusion(x_s,y_s,z_s,x_c,y_c,zcellc,dx,dy,
     +          effdif,nbx,nby,altsol,zondif,ndiff)
                do idi=1,ndiff,stepdi                                     ! Debut de la boucle sur les cellules diffusantes.
                 x_dif=zondif(idi,1)
                 y_dif=zondif(idi,2)
                 zceldi=zondif(idi,3)
                 z_dif=cellh(zceldi)             
                 if((x_dif.gt.nbx).or.(x_dif.lt.1).or.(y_dif.gt.nby).     ! Condition cellule diffusante a l'interieur du domaine.
     +           or.(y_dif.lt.1)) then     
c
c !!!!!!!rien ici???????
c          
                 else
                  if ((x_s.eq.x_dif).and.(y_s.eq.y_dif).and.(z_s.eq. 
     +            z_dif)) then
                   if (verbose.eq.1) then
                     print*,'Position Cell Diffusante = Position Source'
                   endif
                  elseif ((x_c.eq.x_dif).and.(y_c.eq.y_dif).and. 
     +                 (z_c.eq.z_dif)) then
                  else
c=======================================================================
c        Calcul de l'angle zenithal entre la source et la cellule diffusante
c=======================================================================


c ombrage source-cellule diffusante
                   d2=sqrt((real(x_dif-x_s)*dx)**2.+(real(y_dif-y_s)      ! dist max pour l'horizon (i.e. l horiz passe la cell-diff ne compte pas)
     +             *dy)**2.)
                   call horizon(x_s,y_s,z_s,d2,altsol,nbx,nby,
     +             dx,dy,zhoriz,latitu)
                   call anglezenithal(x_s,y_s,z_s,x_dif,y_dif,z_dif,dx,
     +             dy,angzen)                                             ! Calcul de l'angle zenithal source-cell diffusante. 
                   call angleazimutal(x_s,y_s,x_dif,y_dif,dx,dy,          ! calcul de l'angle azimutal cible-cell diffusante
     +             angazi)
                   az=nint(angazi*180./pi)+1
                   if ((angzen).lt.zhoriz(az)) then                       ! debut condition ombrage source-diffusante
c                                                                   
c obstacle sous maille               
                    angmin=pi/2.-atan((obsH+altsol(x_s,y_s)-z_s
     +              )/drefle)
                    if (angzen.lt.angmin) then                            ! Debut condition obstacle source->diffuse.
c                                                                    
c=======================================================================
c        Calcul de la transmittance entre la source et la cellule diffusante
c=======================================================================
                     call transmitm (angzen,x_s,y_s,z_s,x_dif,y_dif,
     +               z_dif,lambda,dx,dy,pressi,transm)
                     call transmita (angzen,x_s,y_s,z_s,x_dif,y_dif,
     +               z_dif,dx,dy,taua,transa) 
c=======================================================================
c     Calcul de l'angle solide couvert par la cellule diffusante vue de la source
c=======================================================================


c                omega2=0.


                     xc=dble(x_dif)*dble(dx)                              ! Position en metres de la cellule diffusante (longitude).
                     yc=dble(y_dif)*dble(dy)                              ! Position en metres de la cellule diffusante (latitu).
                     zc=dble(z_dif)                                       ! Position en metres de la cellule diffusante (altitude).
                     xn=dble(x_s)*dble(dx)                                ! Position en metres de la source (longitude).
                     yn=dble(y_s)*dble(dy)                                ! Position en metres de la source (latitu).
                     zn=dble(z_s)                                         ! Position en metres de la source (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                     if (z_dif .ne. z_s) then
                      call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                cthick,zcellc,r1x,r1y,r1z,r2x,r2y,r2z,
     +                r3x,r3y,r3z,r4x,r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! selon le plan xy.
                      omega1 = omega
                     else
                      omega1=0.
                     endif


c           omega2=omega2+omega1



c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
                     if (y_dif .ne. y_s) then                             ! Si la latitu de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
                      call planzx(dx,xc,xn,yc,yn,zc,zn,cthick,
     +                zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x
     +                ,r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! selon le plan zx.
                     else
                      omega=0.
                     endif
                     if (omega.gt.0.) then
                      if (omega .gt. omega1) omega1 = omega               ! On garde l'angle solide le plus grand jusqu'a present.
                     endif


c           omega2=omega2+omega1



c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
                     if (x_dif .ne. x_s) then                             ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0.
                      call planyz(dy,xc,xn,yc,yn,zc,zn,cthick,
     +                zcellc,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,
     +                r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! selon le plan yz.       
                     else 
                      omega=0.
                     endif
                     if (omega.gt.0.) then
                      if (omega .gt. omega1) omega1 = omega               ! On garde l'angle solide le plus grand.
                     endif
                     omega=omega1


c           omega2=omega2+omega1
c           omega=omega2


c=======================================================================
c Estimation of the subtended angle of the solid angle                    ! Cet angle servira a obtenir un meilleur estime (moyenne) de 
c                                                                         ! P_dir pour le cas de grands angles solides ou pvalno
c=======================================================================  ! varie significativement sur +- ouvang.
                     ouvang=sqrt(omega/pi)                                ! Angle en radian.
                     ouvang=ouvang*180./pi                                ! Angle en degres.
c 
c=======================================================================
c Computing emission function of the source toward the scattering cell    
c=======================================================================
c    
                     anglez=nint(angzen/pi*180.)
                     if (anglez.lt.0) anglez=-anglez
                     if (anglez.gt.180) anglez=360-anglez
                     anglez=anglez+1                                      ! Transformer l'angle en deg entier en position dans la matrice.
c		
c  moyenner sur +- ouvang	
c
c
                     nbang=0.
                     P_dif1=0.
                     do na=-nint(ouvang),nint(ouvang)
                      naz=anglez+na
                      if (naz.lt.0) naz=-naz
                      if (naz.gt.181) naz=362-naz                         ! La fonction est symetrique.
                      P_dif1=P_dif1+pvalno(naz,stype)
                      nbang=nbang+1. 
                     enddo
                     P_dif1=P_dif1/nbang 
c
c=======================================================================
c Computing flux reaching the scattering cell
c=======================================================================
                     fldif1=lamplu(x_s,y_s,stype)*P_dif1*
     +               omega*transm*transa
c=======================================================================
c Computing the scattering probability toward the line of sight cell
c=======================================================================
                     if (angzen.lt.(pi/2.)) then                          ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                      zidif=z_c-0.5*cthick(zceldi)
                      zfdif=z_c+0.5*cthick(zceldi)
                     else
                      zidif=z_c+0.5*cthick(zceldi)
                      zfdif=z_c-0.5*cthick(zceldi)
                     endif       
                     call transmitm (angzen,iun,iun,zidif,ideux,          ! Transmittance moleculaire a l'interieur de la cellule diffusante.
     +               ideux,zfdif,lambda,dx,dy,pressi,tran1m)
                     call transmita (angzen,iun,iun,zidif,ideux,          ! Transmittance aerosols a l'interieur de la cellule diffusante.
     +               ideux,zfdif, dx,dy,taua,tran1a)
                     call angle3points (x_s,y_s,z_s,x_dif,y_dif,z_dif,    ! Angle de diffusion.
     +               x_c,y_c,z_c,dx,dy,angdif)
                     call diffusion(omega,angdif,tran1a,tran1m,           ! Probabilite de diffusion de la lumiere directe.
     +               secdif,fdifan,pdifd1)
c=======================================================================
c Computing scattered intensity toward the line of sight cell from the scattering cell  
c=======================================================================
                     idif1=fldif1*pdifd1
c=======================================================================
c Computing zenith angle between the scattering cell and the line of sight cell
c=======================================================================
        d2=sqrt((real(x_dif-x_c)*dx)**2.+(real(y_dif-y_c)*dy)**2.)        ! dist max pour l'horiz (i.e. l horizon passe la cell-diff ne compte pas)
        call horizon(x_dif,y_dif,z_dif,d2,altsol,nbx,nby,dx,dy,
     +  zhoriz,latitu)
                     call anglezenithal(x_dif,y_dif,z_dif,x_c,y_c,z_c,
     +               dx,dy,angzen)                                        ! Calcul de l'angle zenithal entre la cellule diffusante et la 
c                                                                         ! cellule cible.
        call angleazimutal(x_dif,y_dif,x_c,y_c,dx,dy,angazi)              ! calcul de l'angle azimutal surf refl-cell diffusante
        az=nint(angazi*180./pi)+1
        if ((angzen).lt.zhoriz(az)) then                                  ! debut condition ombrage diffuse-cible  
c                                                                 
c subgrid obstacles                
                     angmin=pi/2.-atan((obsH+altsol(x_dif,
     +               y_dif)-z_dif)/drefle)
                     if (angzen.lt.angmin) then                           ! debut condition obstacle sous maille diffuse->cible.
c                                                                   
c=======================================================================
c Computing transmittance between the scattering cell and the line of sight cell
c=======================================================================
                      call transmitm (angzen,x_dif,y_dif,z_dif,x_c,
     +                y_c,z_c,lambda,dx,dy,pressi,transm)
                      call transmita (angzen,x_dif,y_dif,z_dif,x_c,
     +                y_c,z_c,dx,dy,taua,transa) 
c=======================================================================
c Computing the solid angle of the line of sight cell as seen from the scattering cell
c=======================================================================
                      xc=dble(x_c)*dble(dx)                               ! Position en metres de la cellule cible (longitude).
                      yc=dble(y_c)*dble(dy)                               ! Position en metres de la cellule cible (latitu).
                      zc=dble(z_c)                                        ! Position en metres de la cellule cible (altitude).
                      xn=dble(x_dif)*dble(dx)                             ! Position en metres de la diffusante (longitude).
                      yn=dble(y_dif)*dble(dy)                             ! Position en metres de la diffusante (latitu).
                      zn=dble(z_dif)                                      ! Position en metres de la diffusante (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                      if (z_c .ne. z_dif) then
                       call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                 cthick,zcellc,r1x,r1y,r1z,r2x,r2y,r2z
     +                 ,r3x,r3y,r3z,r4x,r4y,r4z)
                       call anglesolide(omega,r1x,r1y,r1z,                ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                 r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)               ! selon le plan xy.   
                       omega1 = omega
                      else
                       omega1=0.
                      endif
c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
                      if (y_c .ne. y_dif) then                            ! Si la latitu de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
                       call planzx(dx,xc,xn,yc,yn,zc,zn,
     +                 cthick,zcellc,r1x,r1y,r1z,r2x,r2y,r2z,
     +                 r3x,r3y,r3z,r4x,r4y,r4z)
                       call anglesolide(omega,r1x,r1y,r1z,                ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                 r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)               ! selon le plan zx.
                      else
                       omega=0.
                      endif
                      if (omega.gt.0.) then
                       if (omega .gt. omega1) omega1 = omega              ! On garde l'angle solide le plus grand jusqu'a present.
                      endif
c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
                      if (x_c .ne. x_dif) then                            ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0.
                       call planyz(dy,xc,xn,yc,yn,zc,zn,
     +                 cthick,zcellc,r1x,r1y,r1z,r2x,r2y,r2z,
     +                 r3x,r3y,r3z,r4x,r4y,r4z)
                       call anglesolide(omega,r1x,r1y,r1z,                ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                 r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)               ! selon le plan yz.
                      else 
                       omega=0.
                      endif
                      if (omega.gt.0.) then
                       if (omega .gt. omega1) omega1 = omega              ! On garde l'angle solide le plus grand.
                      endif
                      omega=omega1
c=======================================================================
c        Calcul du flux diffuse atteignant la cellule cible
c=======================================================================
                      fldiff=idif1*omega*transm*
     +                transa















c verifie mais je crosi que le calcul de azencl est facultatif ici
                if ((cloudt.ne.0).and.(cloudh(cloudt).eq.zcellc)) then           ! target cell = cloud
                     call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +               dx,dy,azencl)                                        ! zenith angle from cloud to observer                     
                     call cloudreflectance(angzen,cloudt,rcloud)                 ! cloud intensity from direct illum
                     icloud=icloud+
     +               fldiff*rcloud*abs(cos(azencl))/pi
c                    print*,'balise4',icloud
                endif






















c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere diffuse vers la cellule observatrice(SORTANT de cell_c)
c=======================================================================
                      if (angzen.lt.(pi/2.)) then                         ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                       zidif=zcdown
                       zfdif=zcup
                      else
                       zidif=zcup
                       zfdif=zcdown
                      endif
                      call transmitm (angzen,iun,iun,zidif,ideux,         ! Transmittance moleculaire a l'interieur de la cellule diffusante.
     +                ideux,zfdif,lambda,dx,dy,pressi,tran1m)
                      call transmita (angzen,iun,iun,zidif,ideux,         ! Transmittance aerosols a l'interieur de la cellule diffusante.
     +                ideux,zfdif,dx,dy,taua,tran1a)    
                      call angle3points (x_dif,y_dif,z_dif,x_c,y_c,       ! Angle de diffusion.
     +                z_c,x_obs,y_obs,z_obs,dx,dy,angdif)
                      call diffusion(omega,angdif,tran1a,tran1m,          ! Probabilite de diffusion de la lumiere directe.
     +                secdif,fdifan,pdifd2)
c=======================================================================
c Computing scattered intensity toward the observer from the line of sight cell
c=======================================================================
                      idiff2=fldiff*pdifd2
                      idiff2=
     +                idiff2*real(stepdi)                                 ! Corriger le resultat pour le fait d'avoir passe des cellules
                      itodif=                                             ! afin d'accelerer le calcul.
     +                itodif+idiff2
                     endif                                                ! Fin condition obstacle diffuse->cible.
       else
c          print*,'ombrage diff-cible1',x_dif,y_dif,z_dif,x_c,y_c,z_c
        endif                                                             ! fin condition ombrage diffuse-cible                     
                    endif                                                 ! Fin condition obstacle source->diffuse.
                   else
c        print*,'ombrage source-diff',x_s,y_s,z_s,x_dif,y_dif,z_dif
                   endif                                                  ! Fin condition ombrage source-diffusante
                  endif                                                   ! Fin du cas Diffusante = Source ou Cible.
                 endif                                                    ! Fin de la condition "cellule a l'interieur du domaine".      
                enddo                                                     ! Fin de la boucle sur les cellules diffusante.
               endif                                                      ! fin de la condition ou effdif > dx.
c End of 2nd scattered intensity calculations    
c**********************************************************************
c        Calcul de l'intensite provenant d'une source dans la cible vers le capteur
c**********************************************************************
               isourc=intdir
     a         +itotind+itodif 
     +         +itotrd                                                    ! Somme des intensites de chaque type propre a une source  
c                                                                         ! atteignant une cellule cible.
c                                                                         ds l ordre 1st scat; refl->1st scat; 1st scat->2nd scat, refl->1st scat->2nd scat
               if (verbose.eq.1) then
                print*,' Composantes de l''intensite totale:'
                print*,' source->diffuse=',intdir
                print*,' source->reflexion->diffuse=',
     +          itotind
                print*,' source->diffuse->diffuse=',
     +          itodif
                print*,' source->reflexion->diffuse->diffuse=',
     a          itotrd  
               endif
c                   
c**********************************************************************
c        Calcul de l'intensite totale provenant de toutes les sources d'un type dans la cible vers le capteur 
c**********************************************************************
               itotty=itotty
     +         +isourc*real(step(stype)*step(stype))                      ! Somme des intensites de chaque source a une cellule cible.
                                                                          ! Vers le calcul du poids de chaque cellule source i.e. ITT (equivaut 
               ITT(x_s,y_s,stype)=ITT(x_s,y_s,stype)+isourc               ! a itotty mais sous forme matricielle 







          endif                                                           ! fin condition distance source-cible-observateur <= dx/2








              endif                                                       ! Fin de la condition "La luminosite de la case x_s,y_s n'est pas nulle".
             enddo                                                        ! Fin la boucle sur les rangees (latitus) du domaine (y_s).
            enddo                                                         ! Fin la boucle sur les colonnes (longitudes) du domaine (x_s).
c
c   fin du calcul de l'intensite par un type de source
            itotci=itotci
     1      + itotty                                                      ! Somme des intensites de chaque type a une cellule cible.
c interpoler ITT pour combler le step(stype)
            if (step(stype).gt.1) then
             defval=0.
             autom=0
             intype=0
             window=real(step(stype))
             do ii=1,nbx
              do jj=1,nby
               dat(ii,jj)=0.
              enddo
             enddo
             do ii=imin(stype),imax(stype)
              do jj=jmin(stype),jmax(stype)
               dat(ii,jj)=ITT(ii,jj,stype)
              enddo
             enddo
             call interpmatrix(dat,imin(stype),imax(stype),jmin(stype),
     +       jmax(stype),intype,window,autom,defval)
             do ii=imin(stype),imax(stype)
              do jj=jmin(stype),jmax(stype)
               if (lamplu(ii,jj,stype).ne.0.) then
                ITT(ii,jj,stype)=dat(ii,jj)
               endif
              enddo
             enddo            
            endif
            do x_s=imin(stype),imax(stype)
             do y_s=jmin(stype),jmax(stype)
              ITC(x_s,y_s)=ITC(x_s,y_s)+ITT(x_s,y_s,stype)
             enddo   
            enddo  
c calculer lpluto 
            do x_s=1,nbx
             do y_s=1,nby
               lpluto(x_s,y_s)=lpluto(x_s,y_s)+   
     +         lamplu(x_s,y_s,stype)
             enddo
            enddo
           enddo                                                          ! Fin de la boucle sur les types de sources (stype).
c    fin du calcul de l'intensite provenant d'une cellule cible se dirigeant vers le capteur
c
c
c***********************************************************************
c        Calcul du flux lumineux atteignant la cellule ou se trouve le capteur
c***********************************************************************
c
c=======================================================================
c        Calcul de l'angle zenithal entre l'observateur et la cible
c=======================================================================
           call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,dx,dy,
     +     angzen)                                                        ! Calcul de l'angle zenithal entre la cellule cible et l'observateur.
c                                                                         ! Fin du cas "observateur a la meme latitu/longitude que la source".
c=======================================================================
c        Calcul de la transmittance entre la cible et l'observateur
c=======================================================================
           call transmitm (angzen,x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +     lambda,dx,dy,pressi,transm)
           call transmita (angzen,x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +     dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cible vu par l'observateur
c=======================================================================


c             omega2=0.


           xn=dble(x_obs)*dble(dx)                                        ! Position en metres de la cellule observatrice (longitude).
           yn=dble(y_obs)*dble(dy)                                        ! Position en metres de la cellule observatrice (latitu).
           zn=dble(z_obs)                                                 ! Position en metres de la cellule observatrice (altitude).
           xc=dble(x_c)*dble(dx)                                          ! Position en metres de la cible (longitude).
           yc=dble(y_c)*dble(dy)                                          ! Position en metres de la cible (latitu).
           zc=dble(z_c)                                                   ! Position en metres de la cible (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
           if (z_c .ne. z_obs) then
            call planxy(dx,dy,xc,xn,yc,yn,zc,zn,cthick,zcellc,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)  
            call anglesolide(omega,r1x,r1y,r1z,                           ! Appel de la routine anglesolide qui calcule l'angle solide 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! selon le plan xy.
            omega1 = omega
           else
            omega1=0.
           endif


c           omega2=omega2+omega1


c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
           if (y_c .ne. y_obs) then                                       ! Si la latitu de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
            call planzx(dx,xc,xn,yc,yn,zc,zn,cthick,zcellc,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            call anglesolide(omega,r1x,r1y,r1z,                           ! Appel de la routine anglesolide qui calcule l'angle solide 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! selon le plan zx.
           else
            omega=0.
           endif
           if (omega.gt.0.) then
            if (omega.gt.omega1) omega1 = omega                           ! On garde l'angle solide le plus grand jusqu'a present.
           endif


c           omega2=omega2+omega1



c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
           if (x_c.ne.x_obs) then                                         ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0
            call planyz(dy,xc,xn,yc,yn,zc,zn,cthick,zcellc,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            call anglesolide(omega,r1x,r1y,r1z,                           ! Appel de la routine anglesolide qui calcule l'angle solide 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! selon le plan yz.
           else 
            omega=0.
           endif
           if (omega.gt.0.) then
            if (omega.gt.omega1) omega1=omega                             ! On garde l'angle solide le plus grand.
           endif
           omega=omega1


c           omega2=omega2+omega1
c           omega=omega2


c=======================================================================
c        Calcul du flux atteignant l'objectif du telescope en provenance de la cellule cible
c=======================================================================
           dis_obs=sqrt((z_c-z_obs)**2.+((real(y_c-y_obs))*dy)**2.
     a     +((real(x_c-x_obs))*dx)**2.)
           ometif=pi*(diamobj/2.)**2./dis_obs**2.
           if (dis_obs.eq.0.) then
            print*,'ERREUR probleme avec dis_obs',dis_obs
            stop
           endif
           flcib=itotci*ometif*transa*transm                              ! computation of the flux reaching the intrument from the line of sight cell
           do x_s=1,nbx
            do y_s=1,nby
             FC(x_s,y_s)=ITC(x_s,y_s)*ometif*transa*transm
            enddo
           enddo
           omefov=lfente*longfe/focal**2.                                 ! Calcul de l'angle solide de la fente projete sur le ciel.
           if (cos(pi-angzen).eq.0.) then 
            print*,'ERREUR la visee est pafaitement horizontale!'
            stop
           else


             portio=omefov/omega
c            portio=(omefov*dis_obs*dis_obs)/(cos(pi-angzen)*dx*dy)       ! Fraction de la cellule cible vue par le fov (Fraction peut 
c                                                                         ! etre superieure a 1). Le pi ici est du au fait
c                                                                         ! que angzen est calcule sur le trajet cible vers l'observateur
           endif
           if (omega.eq.0.) then
            print*,'ERREUR omega=0 (1)'
            stop
           endif
           fcapt=flcib*portio                                             ! correction for the FOV to the flux reaching the intrument from the line of sight cell












           do x_s=1,nbx
            do y_s=1,nby
             FCA(x_s,y_s)=FC(x_s,y_s)*portio
            enddo
           enddo
c   fin du calcul du flux atteignant la cellule observatrice en provenance d'une cellule cible
           ftocap=ftocap+fcapt  

           do x_s=1,nbx
            do y_s=1,nby
             FTC(x_s,y_s)=FTC(x_s,y_s)+FCA(x_s,y_s)                       ! FTC est la matrice du flux total au capteur permettant d'identifier
                                                                          ! la contribution de chaque cellule du sol au flux total au capteur
                                                                          ! Le % est simplement donne par FTC/ftocap
             flcumu=flcumu+FCA(x_s,y_s)
            enddo
           enddo
          endif                                                           ! Fin de la condition cellule cible n'est pas cellule observatrice.
         endif                                                            ! Fin de la condition cellule cible dans le domaine de calcul.
        endif                                                             ! Fin condition pour la reprise d'un calcul arrete.









c correction for the FOV to the flux reaching the intrument from the cloud cell
           if ((cloudt.ne.0).and.(cloudh(cloudt).eq.zcellc)) then         ! target cell = cloud
c=======================================================================
c  solid angle of the cloud pixel as seen from observer position
c=======================================================================
              xn=dble(x_obs)*dble(dx)                                     ! Position en metres de la cellule observatrice (longitude).
              yn=dble(y_obs)*dble(dy)                                     ! Position en metres de la cellule observatrice (latitu).
              zn=dble(z_obs)                                              ! Position en metres de la cellule observatrice (altitude).
              xc=dble(x_c)*dble(dx)                                       ! Position en metres de la cible (longitude).
              yc=dble(y_c)*dble(dy)                                       ! Position en metres de la cible (latitu).
              zc=dble(z_c)                                                ! Position en metres de la cible (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
              if (z_c .ne. z_obs) then
                 call planxy(dx,dy,xc,xn,yc,yn,zc,zn,cthick,zcellc,
     +           r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)  
                 call anglesolide(omega,r1x,r1y,r1z,                      ! Appel de la routine anglesolide qui calcule l'angle solide 
     +           r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                     ! selon le plan xy.
              else
                 omega=0.
              endif
c computation of the flux reaching the intrument from the cloud cell
              fcloud=icloud*ometif*transa*transm
              fccld=fcloud*omefov/omega
              fctcld=fctcld+fccld
              print*,'tot',icloud,fcloud,fccld,fctcld
            endif          







        print*,' Flux capteur (air/cloud)         =',fcapt,fccld
        print*,' Flux capteur cumule (air/cloud)  =',ftocap,fctcld
        write(2,*) ' Flux capteur (air/cloud)         =',fcapt,fccld
        write(2,*) ' Flux capteur cumule (air/cloud)  =',ftocap,fctcld   
      endif                                                               ! fin condition cellule cibles 1/5000
       enddo                                                              ! Fin de la boucle sur les cellules cibles.
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
                  endif                                                   ! FTCN est le % par unite de luminosite de la cellule
                  write(9,*) x_s,y_s,FTC(x_s,y_s)                         ! emettrice au sol, c'est un % par unite de watt installes
                  write(8,*) x_s,y_s,FTCN(x_s,y_s)
               enddo
            enddo
            nom='Grid weight '
            valmax=65535       
            gain=ftcmax/real(valmax)
            offset=0.
            call extrants2d (pclimg,FTC,nom,xcell0,ycell0,pixsiz,
     +      gain,offset,nbx,nby,valmax)
            nom='NormGrid wgt'
            call extrants2d (pcwimg,FTCN,nom,xcell0,ycell0,pixsiz,
     +      gain,offset,nbx,nby,valmax)     
          close(unit=8)
          close(unit=9)
c creation des fichiers gnuplot pour les visualiser il faut entrer gnuplot
c puis load 'fichier.gplot'
          open(unit=9,file=pclgp,status='unknown')
          open(unit=8,file=pcwgp,status='unknown')
            write(9,*) 'set dgrid3d',nbx,',',nby
            write(9,*) 'set hidden3d'
            write(9,*) 'set pm3d'
            write(9,*) 'splot "'//basenm(1:lenbase)//'_pcl.txt"
     +      with dots'
            write(8,*) 'set dgrid3d',nbx,',',nby
            write(8,*) 'set hidden3d'
            write(8,*) 'set pm3d'
            write(8,*) 'splot "'//basenm(1:lenbase)//'_pcw.txt"
     +      with dots'    
          close(unit=8)
          close(unit=9) 
       endif                                                              ! end of condition for creating contrib and sensit maps
          print*,'====================================================='
          print*,'          Total flux entering instrument (W)'
          write(*,2001) ftocap*real(vistep)+fctcld  
          print*,'              Sky luminance (W/str/m**2)'       
          write(*,2001) (ftocap+fctcld)/(lfente*
     +          longfe/focal**2.)/(pi*(diamobj/2.)**2.)*
     +          real(vistep)
       print*,'  '
       print*,' Interpolation flux error= ',
     +          ftocap-flcumu
       print*,'======================================================='
       write(2,*) '==================================================='
       write(2,*) '          Total flux entering instrument (W)'
       write(2,2001) ftocap*real(vistep)+fctcld
        write(2,*) '            Sky luminance (W/str/m**2)          '      
       write(2,2001) (ftocap+fctcld)/(lfente*
     +          longfe/focal**2.)/(pi*(diamobj/2.)**2.)*
     +          real(vistep)
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
c*                                         fin du programme                                                            *
c*                                                                                                                     *
c***********************************************************************************************************************
