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
c ** Illumina can be downloaded via:   svn checkout http://illumina.googlecode.com/svn/trunk/ illumina                **
c ** To compile, run bash makeILLUMINA         makeILLUMINA                                                           **
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
c **    - Exponential concentrations vertical profile                                                                 **
c **    - Exponential vertical resolution                                                                             **
c **    - Accounting for heterogeneity of ground reflectance, luminaires number, luminaire height, angular photometry **
c **      de fut.                                                                                                     **
c **    - Wavelength dependant                                                                                        **
c **    - Ignore the flux scattered by the voxel occupied by the observer (cellobs=cellcible)                         **
c **    - Do not support direct observation of a source                                                               ** 
c **    - Direct observation of the ground not implemented                                                            **
c **    - Not accounting for molecular absorption                                                                     **
c **    - Do not consider earth curvature (i.e. local/regional model)                                                 **
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
      real cell_thickness(50)                                             ! Cell thickness array (meter)

      real cell_height(50)                                                ! Cell height array (meter)

      real fluxcumule                                                     ! Flux cumule en fonction du deplacement le long de la ligne de vise
      character*72 relief_file                                            ! Fichier relief.
      character*72 reflex_file                                            ! Fichier reflexion.
      character*72 dif_file                                               ! Fichier diffusion.
      character*72 outfile                                                ! Fichier de resultat.
      character*72 pclfile,pcwfile,pclgplot,pcwgplot                      ! fichiers de matrice de % du flux total normalise par unite de flux 
      character*72 pclimg,pcwimg
                                                                          ! et par flux et wattage
      character*72 basename                                               ! Nom de base des fichiers.
      character*12 nom                                                    ! Nom de la variable 2d lue ou ecrite
      integer lenbase                                                     ! Longueur du nom de base de l'experimentation.
      real lambda,pression,portee_reflex                                  ! Longueur d'onde (nanometre), pression atmospherique (kPa), zone 
c                                                                         ! de portee de la reflexion (metre).
      integer ntype                                                       ! Nombre de type de sources lumineuses consideres.
      real largeur_x                                                      ! Largeur (valeur x) du domaine de l'environnement (metre).
      real largeur_y                                                      ! Longueur (valeur y) du domaine de l'environnement (metre).
      integer nbx,nby                                                     ! Nombre de cellule dans le domaine.  
      real valeur2d(width,width)                                          ! Matrice temporaire pour les intrants 2d
      real alt_sol(width,width)                                           ! Altitude du sol (metre).
      real surface_reflect(width,width)                                   ! Reflectance des surfaces.
      real Hmin                                                           ! Hauteur minimale du domaine.
      real xcell0                                                         ! Longitude sud-ouest du domaine.
      real ycell0                                                         ! Latitude sud-ouest du domaine.
      real gain                                                           ! Coefficient multiplicatif des valeurs des fichiers source.
      real offset                                                         ! Constante d'addition des valeurs des fichiers source.
      integer valmax                                                      ! Maximum value of the output pgms
      integer stype                                                       ! Identification du type de source.
      character*72 pafile,lufile,alfile                                   ! Fichiers relatifs aux sources lumineuses (fonction d'emission 
c                                                                         ! des sources (lampadaires), luminosite (DIMENSIONS??), altitude (metre).    
      real lamp_lumin(width,width,9)                                      ! Luminosite des sources dans chaque case (DIMENSION???).
      real lamp_alt(width,width,9)                                        ! Altitude des sources lumineuses par rapport au sol (metre).
      real p_valeur(181,9),p_valtot,p_valnorm(181,9)                      ! Valeurs des fonctions angulaires d'emission (arbitraires, 
c                                                                         ! totale non-matrice, normalisees).
      real dtheta                                                         ! Increment d'angle de la fonction d'emission des sources 
      real dx,dy,dxp,dyp,pixsiz                                           ! Largeur d'une cellule (metre), nombre de cellule dans la portee 
c                                                                         ! de reflexion (cellule), largeur d'une cellule.
c                                                                         ! pour verif fichiers compatibles, taille des sous cellules pour 
c                                                                         ! l'angle solide de la portion de surface qui reflechit.
      integer boxx,boxy                                                   ! Portee de reflexion (cellule).
      real foncdif_a(181),foncdif_anorm(181)                              ! Fonction de diffusion (arbitraire et normalisee) des aerosols.
      real extinction,scattering,angle_a(181)                             ! Sections efficaces d'extinction et de diffusion, angle de 
c                                                                         ! diffusion (degre).
      real secdif                                                         ! Contribution de la diffusion a l'extinction lumineuse.
      real inclin_x(width,width)                                          ! Inclinaison de la case en x (DIMENSIONS).
      real inclin_y(width,width)                                          ! Inclinaison de la case en y (DIMENSIONS).   
      integer x_obs,y_obs,zcell_obs                                       ! Position de l'observateur (cellule).
      real z_obs                                                          ! Hauteur de l'observateur (metre).
      integer lignecible(width,3)                                         ! Matrice des cellules cibles.
      integer ncible,icible                                               ! Nombre de cellules cibles, compteur de la boucle sur les cellules 
c                                                                         ! cibles.     
      integer x_c,y_c,zcell_c                                             ! Position de la cellule cible (cellule).
      real z_c                                                            ! Hauteur de la cellule cible (metre).
      real zcup,zcdown                                                    ! Limites inferieure et superieure de la cellule cible.    
      integer dircheck                                                    ! Verificateur de position de la source (cas source=cell cible).     
      integer x_s,y_s,x_sr,y_sr,x_dif,y_dif,zcell_dif                     ! Positions de la source, surface reflectrice, cellules diffusantes 
c                                                                         ! (cellule).
      real z_s,z_sr,z_dif                                                 ! Hauteur de la source, surface reflectrice, cellule diffusante (metre).
      real anglezen,ouvang                                                ! Angle zenithal entre deux cellules (radians) et angle d'ouverture 
c                                                                         ! du cone d'angle solide en degres.
      integer anglez                                                        
c                                                                         ! d'emission des lampadaires.      
      real P_dir,P_indir,P_dif1                                           ! Fonction d'emission (directe,indirecte,diffusee) des sources 
c                                                                         ! lumineuses.
      real transa,transm                                                  ! Transmittance entre deux cellules (aerosols,molecules).
      real trans_1a,trans_1m                                              ! Transmittance a l'interieur d'une cellule (aerosols,molecules).
      real taua                                                           ! epaisseur optique des aerosols.
      real*8 xc,yc,zc,xn,yn,zn                                            ! Position (metre) des elements (arrivee, depart) pour le calcul 
c                                                                         ! de l'angle solide.  
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! Composantes des vecteurs utilises dans la routine angle solide.
      real omega,omega1                                                   ! Angle solide couvert par une cellule vue d'une autre, angle de 
c                                                                         ! comparaison.
      real flux_direct                                                    ! Flux provenant d'une cellule source dans une cellule cible (watt).
      real flux_indirect                                                  ! Flux provenant d'une cellule reflectrice dans une cellule cible (watt).
      real flux_diffuse                                                   ! Flux provenant d'une cellule diffusante dans une cellule cible (watt).
      real zidif,zfdif                                                    ! Limites initiale et finale du parcours de diffusion dans une cellule.
      real angle_dif                                                      ! Angle de diffusion.
      real probdif_dir,probdif_indir,probdif_dif1,probdif_dif2            ! Probabilite de diffusion (directe,indirecte,premiere et deuxieme 
c                                                                         ! diffusions.
      real intensite_directe                                              ! Contribution d'une source a l'intensite directe dirigee vers 
c                                                                         ! le capteur par une cellule cible.
      real intensite_indirecte                                            ! Contribution d'une cellule reflectrice a l'intensite indirecte 
c                                                                         ! dirigee vers le capteur.
      real intensite_totale_indirecte                                     ! Contribution totale d'une source a l'intensite indirecte dirigee 
c                                                                         ! vers le capteur.
      real intensite_diffusee2                                            ! Contribution d'une cellule diffusante a l'intensite diffusee 
c                                                                         ! dirigee vers le capteur.
      real intensite_totale_diffusee                                      ! Contribution totale d'une source a l'intensite diffusee dirigee 
c                                                                         ! vers le capteur.
      real intensite_source                                               ! Contribution totale d'une source a l'intensite dirigee par une 
c                                                                         ! cellule cible vers le capteur.
      real intensite_totale_type                                          ! Contribution totale d'un type de source a l'intensite dirigee par 
c                                                                         ! une cellule cible vers le capteur.
      real intensite_totale_cible                                         ! Intensite totale dirigee par une cellule cible vers le capteur.
      real intensite_totale_reflechie_diffusee                            ! Intensite totale dirigee par une cellule vers le capteur apres 
c                                                                         ! reflexion et double diffusion.
      real intensite_reflechie_diffusee                                   ! Intensite dirigee par une cellule vers le capteur apres reflexion 
c                                                                         ! et double diffusion. 
      real flux_cible                                                     ! Flux parvenant dans la cellule observatrice par une cellule cible.
      real flux_capteur                                                   ! Flux parvenant dans la cellule observatrice par toutes les cellules 
c                                                                         ! cibles d'un niveau contenues dans le champ de vision du capteur.
      real flux_total_capteur                                             ! Flux total parvenant dans la cellule receptrice.
      real haut                                                           ! Haut negatif indique que la surface est eclairee par en dessous
c                                                                         ! (on ne calcule alors pas).
      real epsilx,epsily                                                  ! Inclinaison de la surface reflectrice.
      real flux_reflechi                                                  ! Flux atteignant une surface reflectrice (watts).
      real intensite_reflechie,intensite_reflechie1                       ! Intensite quittant une surface reflectrice en direction d'une 
c                                                                         ! cellule cible.  
      real effetdif                                                       ! Distance autour des cellule source et cible qui seront considerees  
c                                                                         ! pour calculer la double diffusion.
      integer zonedif(3000000,4)                                          ! Matrice des cellules diffusantes la 4e colonne represente la valeur 
c                                                                         ! entiere la plus proche de la distance (en metre) a la ligne de 
c                                                                         ! diffusion simple.
      integer ndiff,idi                                                   ! Nombre de cellules diffusantes, compteur de la boucle sur les cellules
c                                                                         ! diffusantes. 
      integer stepdif                                                     ! Saut de diffusion pour accelerer le calcul e.g. si =2 il fera un  
      integer reducdif                                                    ! calcul sur deux flag indiquant si le rayon de diffusion a ete reduit
      integer nvis0                                                       ! valeur de depart pour le calcul le long de la ligne de visee. 
c                                                                         ! Par defaut cette valeur est de 1 mais elle peut etre plus grande  
c                                                                         ! lorsque l'on reprend un calcul precedent qui a ete arrete anormalement.
      real flux_dif1                                                      ! Flux atteignant une cellule diffusante.
      real intensite_diffusee1                                            ! Intensite dirigee vers une cellule cible par une cellule diffusante.
      real portion                                                        ! Portion de cellule observee dans le champ de vision du capteur.
      real dis_obs                                                        ! Distance d'observation entre la cible et l'observateur.
      real omegtif                                                        ! Angle solide couvert par l'objectif du telescope vu de la cellule 
c                                                                         ! cible.
      real omegafov                                                       ! Angle solide couvert sur le ciel vu par la fente.
      real largefente                                                     ! Largeur de la fente (ou de la surface du capteur)
      real longfente                                                      ! Longueur de la fente (ou de la surface du capteur)
      real distfocal                                                      ! Distance focale effective objectif (selon le Throughput de l'appareil)  
      real angle_visee,angle_azimut                                       ! Angles d'observation du recepteur.
      real projapparente                                                  ! Fraction de la surface reflectrice vue par rapport a la direction 
c                                                                         ! normale. Utile pour le calcul de la reflectance lambertienne.
      real nbang                                                          ! pour le moyennage de la fonction d'emission
      real obstacleH,anglemin                                             ! Hauteur moyenne des obstacle sous maille, angle minimum en dessous 
c                                                                         ! duquel un faisceau ne peut etre propage en raison de l'obstacle 
c                                                                         ! sous-maille.
      integer naz,na 
      real ITT(width,width,9)                                             ! Intensite totale type en forme matricielle
      real ITC(width,width)                                               ! Intensite totale cible en forme matricielle
      real FC(width,width)                                                ! Flux cible en forme matricielle
      real FTC(width,width)                                               ! Fraction du Flux total au capteur en forme matricielle
      real FTCN(width,width)                                              ! Fraction Flux total au capteur normalise par unite de wattage au sol
      real FCA(width,width)                                               ! Flux capteur sous forme matricielle
      real lamp_lum_tot(width,width)                                      ! Luminosite totale de la cellule de surface toutes lampes confondues
      real FTCNTOT,ftcmax                                                 ! FTCN total pour tout le domaine pour toutes lampes
      character*3 lampno                                                  ! mot de trois lettres contenant le numero de lampe
      integer imin(9),imax(9),jmin(9),jmax(9),step(9)                     ! bornes x et y de la zone contenant un type de lampe
      real defval                                                         ! valeur ignoree dans l'interpolation
      real dat(1024,1024)                                                 ! matrice a interpoler
      integer autom,intype,ii,jj                                          ! switch manuel automatique pour l'interpolation; interpolation type
      real window                                                         ! interpolation diameter
      real zen_horiz(360)                                                 ! horizon en rad sur 360 deg, la posi 1 = 0 deg et 360 = 359deg
      real angleazi,d2                                                    ! angle azimutal entre deux points en rad, dist max pour lhorizon
      integer az                                                          ! azimut de l'horizon
      real latitude                                                       ! latitude approx du centre du domaine pour l'instant on utilises ycell0
      integer vistep                                                      ! line of sight step for low elevation angles vistep=ncells_along_sight/50
      integer prmaps                                                      ! frag to enable saving contribution and sensitivity maps
      data cell_thickness /0.5,0.6,0.72,0.86,1.04,1.26,1.52,1.84,2.22,    ! Epaisseur des niveaux.
     a 2.68,3.24,3.92,4.74,5.72,6.9,8.34,10.08,12.18,14.72,17.78,21.48,
     b 25.94,31.34,37.86,45.74,55.26,66.76,80.64,97.42,117.68,142.16,
     c 171.72,207.44,250.58,302.7,365.66,441.72,533.6,644.58,778.66,
     d 940.62,1136.26,1372.6,1658.1,2002.98,2419.6,2922.88,3530.84,
     e 4265.26,5152.44/
      data cell_height /0.25,0.8,1.46,2.25,3.2,4.35,5.74,7.42,9.45,       ! Hauteur du centre de chaque niveau.
     a 11.9,14.86,18.44,22.77,28.,34.31,41.93,51.14,62.27,75.72,91.97,
     b 111.6,135.31,163.95,198.55,240.35,290.85,351.86,425.56,514.59,
     c 622.14,752.06,909.,1098.58,1327.59,1604.23,1938.41,2342.1,
     d 2829.76,3418.85,4130.47,4990.11,6028.55,7282.98,8798.33,
     e 10628.87,12840.16,15511.4,18738.26,22636.31,27345.16/
c  
c=======================================================================
c        Lecture du fichier d'entree (illumina.in)
c=======================================================================
      open(unit=1,file='illumina.in',status='old')
       read(1,*)
       read(1,*) basename
       read(1,*) 
       read(1,*) dx,dy
       read(1,*)
       read(1,*) dif_file
       read(1,*) 
       read(1,*) effetdif,stepdif
       if (verbose.eq.1) then
         print*,'Rayon de diffusion=',effetdif,'m   1 calcul sur ',
     a   stepdif
       endif
       read(1,*)
       read(1,*) lambda
       read(1,*) pression
       read(1,*) taua
       read(1,*) ntype
       read(1,*) 
       read(1,*) portee_reflex,obstacleH
       read(1,*)
       read(1,*) x_obs,y_obs,zcell_obs,nvis0
       read(1,*)
       read(1,*) angle_visee,angle_azimut
       read(1,*)  
       read(1,*) largefente,longfente,distfocal,diamobj     
      close(1)
c
c  determiner la longueur du nom
c 
      lenbase=index(basename,' ')-1  
      relief_file=basename(1:lenbase)//'_topogra.pgm'                     ! Determiner les noms de fichiers d'entree et de sortie
      reflex_file=basename(1:lenbase)//'_reflect.pgm' 
      outfile=basename(1:lenbase)//'.out'  
      pclfile=basename(1:lenbase)//'_pcl.txt'
      pcwfile=basename(1:lenbase)//'_pcw.txt'
      pclimg=basename(1:lenbase)//'_pcl.pgm'
      pcwimg=basename(1:lenbase)//'_pcw.pgm'
      pclgplot=basename(1:lenbase)//'_pcl.gplot'
      pcwgplot=basename(1:lenbase)//'_pcw.gplot'    
c  conversion des angles de visee geographique vers angle mathematique
c  on presume que maintenant l'angle mis dans le fichier illumina.in
c  est en reference a la definition geographique
c  geographie, azim=0 au nord, 90 a l'est, 180 au sud etc
c  mathematiques, azim=0 a l'est, 90 au nord, 180 a l'ouest etc
      angle_azimut=90.-angle_azimut
      if (angle_azimut.lt.0.) angle_azimut=angle_azimut+360.
      if (angle_azimut.ge.360.) angle_azimut=angle_azimut-360.
c  ouverture en ecriture du fichier de sortie
      open(unit=2,file=outfile,status='unknown')      
       write(2,*) 'FICHIERS UTILISES'
       write(2,*) relief_file,reflex_file,dif_file
       print*,'Longueur d''onde (nm):',lambda,
     +       ' epaisseur optique:',taua
       write(2,*) 'Longueur d''onde (nm):',lambda,
     +       ' epaisseur optique:',taua
       write(2,*) 'Zone de diffusion:',effetdif,' m'
       print*,'Zone de diffusion:',effetdif,' m'
       write(2,*) 'Saut de diffusion:',stepdif
       print*,'Saut de diffusion:',stepdif
       write(2,*) 'Portee de la reflexion:',portee_reflex,' m'
       write(2,*) 'Position observateur(x,y,z)',x_obs,y_obs,zcell_obs
       print*,'Position observateur(x,y,z)',x_obs,y_obs,zcell_obs
       write(2,*) 'Angle d elevation:',angle_visee,' angle azim (horaire
     + p/r au nord)',angle_azimut     
       print*,'Angle d elevation:',angle_visee,' angle azim (anti-horair
     +e p/r a l est)',angle_azimut 
c=======================================================================
c        Initialisation des matrices et variables
c=======================================================================
       verbose=0
       prmaps=1
       iun=0
       ideux=1
       fluxcumule=0.
       do i=1,width
        do j=1,width
         valeur2d(i,j)=0.
         alt_sol(i,j)=0.
         surface_reflect(i,j)=0.
         inclin_x(i,j)=0.
         inclin_y(i,j)=0.
         lamp_lum_tot(i,j)=0.
         ITC(i,j)=0.
         FC(i,j)=0.
         FTC(i,j)=0.
         FTCN(i,j)=0.
         FCA(i,j)=0.
         do k=1,9
          lamp_lumin(i,j,k)=0.
          lamp_alt(i,j,k)=0.
          ITT(i,j,k)=0.
         enddo
        enddo
       enddo
       do i=1,181
        do j=1,9
         p_valeur(i,j)=0.
         p_valnorm(i,j)=0.
        enddo
       enddo
       do i=1,181
        foncdif_a(i)=0.
        foncdif_anorm(i)=0.
        angle_a(i)=0.
       enddo     
       do i=1,1024
        do j=1,3
         lignecible(i,j)=1
        enddo
       enddo
       do i=1,3000000
        do j=1,4
         zonedif(i,j)=1
        enddo
       enddo     
       reducdif=1
       intensite_reflechie_diffusee=0.
       anglemin=0.
       vistep=1.
c***********************************************************************
c        Lecture des variables propres a l'environnement               *
c***********************************************************************
c=======================================================================
c  Lecture du fichier de relief
c=======================================================================
       nom='relief'
       call intrants2d(relief_file,alt_sol,nom,xcell0,ycell0,pixsiz,
     + nbx,nby)

       latitude=ycell0

       Hmin=3000000.
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
c                                                                         ! Recherche de la hauteur minimale.
         if (Hmin.gt.alt_sol(i,j)) Hmin=alt_sol(i,j)
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
       enddo 
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
         alt_sol(i,j)=alt_sol(i,j)-Hmin                                   ! Soustraction de la hauteur minimale du domaine.
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
       enddo
c=======================================================================
c Lecture du fichier de reflexion
c=======================================================================
       nom='reflexion'
       call intrants2d(reflex_file,surface_reflect,nom,xcell0,ycell0,
     + pixsiz,nbx,nby)
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
         if (surface_reflect(i,j).lt.0.) then                             ! Recherche de reflectances negatives
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
        p_valtot=0.
c        if (stype.lt.10) then
c          lampno=char(48)//char(48+stype)
c        endif
        write(lampno, '(I3.3)' ) stype                                    ! support de 999 sources differentes (3 digits)
        pafile=basename(1:lenbase)//'_fctem_'//lampno//'.dat'             ! Attribution du nom du fichier  fonction angulaire d'emission.
        lufile=basename(1:lenbase)//'_lumlp_'//lampno//'.pgm'             ! Attribution du nom du fichier de la luminosite des cases.
        alfile=basename(1:lenbase)//'_altlp_'//lampno//'.pgm'             ! Attribution du nom du fichier  hauteur des sources lumineuse.
        open(UNIT=1, FILE=pafile,status='OLD')                            ! Ouverture du fichier pa#.dat, fonction angulaire d'emission.
        do i=1,181                                                        ! Debut de la boucle pour les 181 donnees.
         read(1,*) p_valeur(i,stype)                                      ! Lecture des donnees qui sont inscrites dans la matrice p_valeur.
         p_valtot=p_valtot+p_valeur(i,stype)*2.*pi*                       ! Sommation de la valeur tot  fonction d'emission prenormalisee 
     a   sin(real(i-1)*dtheta)*dtheta                                     ! (pvaleur x 2pi x sin theta x dtheta) (ou theta egale 
c                                                                         ! (i-1) x 1 degres).
        enddo                                                             ! Fin de la boucle sur les 181 donnees du fichier pa#.dat.
        close(1)                                                          ! Fermeture du fichier pa#.dat, fonction angulaire d'emission.
        do i=1,181
         p_valnorm(i,stype)=p_valeur(i,stype)/p_valtot                    ! Normalisation de la fonction d'emission.
        enddo                                                       
c    ===================================================================
        nom='luminosite'
        call intrants2d  (lufile,valeur2d,nom,xcell0,ycell0,pixsiz,
     +  nbx,nby)
       do i=1,nbx                                                         ! Debut de la boucle sur toutes les cases en x.
        do j=1,nby                                                        ! Debut de la boucle sur toutes les cases en y.
         if (valeur2d(i,j).lt.0.) then                                    ! Recherche de luminosites negatives
           print*,'***Negative lamp luminosity!, stopping execution'
           stop
         endif
        enddo                                                             ! Fin de la boucle sur toutes les cases en y.
       enddo     
        do i=1,nbx                                                        ! recherche du plus petit rectangle encadrant la zone
         do j=1,nby                                                       ! de luminosite non nulle pour accelerer le calcul
          if (valeur2d(i,j).ne.0.) then
           if (i-1.lt.imin(stype)) imin(stype)=i-2
           if (imin(stype).lt.1) imin(stype)=1
           goto 333
          endif
         enddo 
        enddo
        imin(stype)=1   
 333    do i=nbx,1,-1
         do j=1,nby
          if (valeur2d(i,j).ne.0.) then
           if (i+1.gt.imax(stype)) imax(stype)=i+2    
           if (imax(stype).gt.nbx) imax(stype)=nbx
           goto 334
          endif
         enddo
        enddo
        imax(stype)=1
 334    do j=1,nby
         do i=1,nbx
          if (valeur2d(i,j).ne.0.) then
           if (j-1.lt.jmin(stype)) jmin(stype)=j-2 
           if (jmin(stype).lt.1) jmin(stype)=1
           goto 335
          endif
         enddo
        enddo 
        jmin(stype)=1
 335    do j=nby,1,-1
         do i=1,nbx
          if (valeur2d(i,j).ne.0.) then
           if (j+1.gt.jmax(stype)) jmax(stype)=j+2
           if (jmax(stype).gt.nby) jmax(stype)=nby
           goto 336
          endif
         enddo
        enddo  
        jmax(stype)=1
 336    do i=1,nbx                                                        ! Debut de la boucle sur toutes les cases en x.
         do j=1,nby                                                       ! Debut de la boucle sur toutes les cases en y.
          lamp_lumin(i,j,stype)=valeur2d(i,j)                             ! remplir la matrice du type de lampe stype
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
c        if (dzone.gt.50) step(stype)=2                                    ! saut pour le calcul des cellules au sol ce saut augmente
c        if (dzone.gt.100) step(stype)=3                                   ! avec la distance
c        if (dzone.gt.150) step(stype)=4 
c        if (dzone.gt.200) step(stype)=5
c    ==================================================================
        nom='hauteurs'
        call intrants2d(alfile,valeur2d,nom,xcell0,ycell0,pixsiz,nbx,
     +  nby)
        do i=1,nbx                                                        ! Debut de la boucle sur toutes les cases en x.
         do j=1,nby                                                       ! Debut de la boucle sur toutes les cases en y.
          lamp_alt(i,j,stype)=valeur2d(i,j)                               ! Remplissage de la matrice pour la lampe stype
         enddo                                                            ! Fin de la boucle sur toutes les cases en y.
        enddo                                                             ! Fin de la boucle sur toutes les cases en x.
       enddo                                                              ! Fin de la boucle sur les 99 types de sources.    
c=======================================================================
c        Lecture des parametres de diffusion
c=======================================================================
       open(unit = 1, file = dif_file,status= 'old')                      ! Ouverture du fichier contenant les parametres de diffusion.
c                                                                         ! Le fichier de diffusion est genere par le programme imies 
c                                                                         ! du progiciel AODSEM (Martin Aube).
        read(1,*)                                                           
        read(1,*)
        read(1,*)
        do i=1,181
         read(1,*) angle_a(i), foncdif_a(i)                               ! Lecture de la fonction de diffusion et de l'angle associe a 
c                                                                         ! cette fonction de 0 a 180 degres soit 181 lignes.
         foncdif_anorm(i)=foncdif_a(i)/pix4                               ! Normalisation de la fonction a 4 pi (l'integrale de la 
c                                                                         ! fonction fournie sur tous les angles solides doit etre egale a 4 pi).
c                                                                         ! En effet le fichier .mie.out est normalisé ainsi (revefie par 
c                                                                         ! M. Aube en avril 2009)
        enddo
        do i = 1,7
         read(1,*)
        enddo
        read(1,*) extinction                                              ! Lecture de la section efficace d'extinction des aerosols.
        read(1,*) scattering                                              ! Lecture de la section efficace de diffusion des aerosols.
       close(1)
       secdif=scattering/extinction                                       ! Rapport (sigmadif/sigmatotal).
c======================================================================
c        Quelques operations preparatoires
c======================================================================
       dy=dx                                                              ! Pour le moment on considere que l'echelle est la meme sur les deux axes
       z_obs=cell_height(zcell_obs)                                       ! Attribution d'une valeur en metre a la position z de l'observateur.
       largeur_x=dx*real(nbx)                                             ! Calcul de la largeur en x d'une case.
       largeur_y=dy*real(nby)                                             ! Calcul de la largeur en y d'une case.
       boxx=nint(portee_reflex/dx)                                        ! Nombre de colonnes a considerer a gauche/droite de la source 
c                                                                         ! pour la reflexion.
       boxy=nint(portee_reflex/dy)                                        ! Nombre de colonnes a considerer en haut/bas de la source pour 
c                                                                         ! la reflexion.
       write(2,*) 'Largeur du domaine [NS](m):',largeur_x,'#cases:',nbx
       write(2,*) 'Largeur du domaine [EO](m):',largeur_y,'#cases:',nby
       write(2,*) 'Taille d''une cellule (m):',dx,' X ',dy
       write(2,*) 'Latitude sud-ouest:',ycell0,' Longitude sud-ouest:',
     + xcell0
c=======================================================================
c        Calcul de l'inclinaison des cases en x et en y
c=======================================================================
       do i=1,nbx                                                         ! Debut de la boucle sur les colonnes (longitude) du domaine.
        do j=1,nby                                                        ! Debut de la boucle sur les ranges (latitude) du domaine.
         if (i.eq.1) then                                                 ! Cas particulier au bord du domaine (cote vertical gauche).
          inclin_x(i,j)=atan((alt_sol(i+1,j)-alt_sol(i,j))/real(dx))      ! Calcul de l'inclinaison en x de la surface.
         elseif (i.eq.nbx) then                                           ! Cas particulier au bord du domaine (cote vertical droit).
          inclin_x(i,j)=atan((alt_sol(i-1,j)-alt_sol(i,j))/(real(dx)))    ! Calcul de l'inclinaison en x de la surface.
         else
          inclin_x(i,j)=atan((alt_sol(i+1,j)-alt_sol(i-1,j))/(2.          ! Calcul de l'inclinaison en x de la surface.
     1    *real(dx)))
         endif
         if (j.eq.1) then                                                 ! Cas particulier au bord du domaine (cote horizontal en bas).
          inclin_y(i,j)=atan((alt_sol(i,j+1)-alt_sol(i,j))/(real(dy)))    ! Calcul de l'inclinaison en y de la surface.
         elseif (j.eq.nby) then                                           ! Cas particulier au bord du domaine (cote horizontal en haut).
          inclin_y(i,j)=atan((alt_sol(i,j-1)-alt_sol(i,j))/(real(dy)))    ! Calcul de l'inclinaison en y de la surface.
         else
          inclin_y(i,j)=atan((alt_sol(i,j+1)-alt_sol(i,j-1))/(2.          ! Calcul de l'inclinaison en y de la surface.
     1    *real(dy)))
         endif
        enddo                                                             ! Fin de la boucle sur les ranges (latitude) du domaine.
       enddo                                                              ! Fin de la boucle sur les colonnes (longitude) du domaine.
c=======================================================================
c        Debut de la boucle sur les cellules cibles
c=======================================================================
       call lignevisee(x_obs,y_obs,z_obs,dx,dy,angle_visee,               ! Determination de la ligne de visee (cellules cibles).
     + angle_azimut,nbx,nby,lignecible,ncible,vistep)
       flux_total_capteur=0.                                              ! Initialisation de la valeur du flux recu par le capteur.
       do icible=1,ncible                                                 ! Debut de la boucle sur les cellules cibles.
        if (icible.ge.nvis0) then                                         ! Debut condition pour la reprise d'un calcul arrete   .
         intensite_totale_cible=0.                                        ! Initialisation de la contribution d'une cible au capteur.
         do i=1,nbx
          do j=1,nby
            ITC(i,j)=0.
          enddo
         enddo
         zcell_c=lignecible(icible,3)                                     ! Definition de la position (cellule) verticale de la cible.
         z_c=cell_height(zcell_c)                                         ! Definition de la position (metre) verticale de la cible.
         y_c=lignecible(icible,2)                                         ! Definition de la position (cellule) de la cible.
         x_c=lignecible(icible,1)                                         ! Definition de la position (cellule) de la cible.
         print*,'=================================================='
         print*,' Avancement le long de la ligne de visee :',
     +   icible,'/',ncible
         print*,' Hauteur de la cible   =',z_c,' m'
         print*,' Epaisseur de la cible =',cell_thickness(zcell_c),' m'
         write(2,*) '=================================================='
         write(2,*) ' Avancement le long de la ligne de visee :',
     +   icible,'/',ncible
         write(2,*) ' Hauteur de la cible   =',z_c,' m'
      write(2,*) ' Epaisseur de la cible =',cell_thickness(zcell_c),' m'
         if( (x_c.gt.nbx).or.(x_c.lt.1).or.(y_c.gt.nby).or.(y_c.lt.1)     ! Condition Cellule cible dans le domaine de calcul.
     +      .or.(zcell_c.gt.50).or.(zcell_c.lt.1) )then
         else
          if((x_c.eq.x_obs).and.(y_c.eq.y_obs).and.                       ! Pour le moment, si la cellule cible est la cellule observatrice, 
c                                                                         ! on ne calcule pas le flux diffuse.
     +    (zcell_c.eq.zcell_obs))then
           if (verbose.eq.1) then
             print*,'Cellule Cible = Cellule Observateur'
           endif
          else
           zcdown=z_c-0.5*cell_thickness(zcell_c)                         ! Limite inferieure de la cellule cible.
           zcup=z_c+0.5*cell_thickness(zcell_c)                           ! Limite superieure de la cellule cible.
c=======================================================================
c        Debut de la boucle sur les types de sources lumineuses
c=======================================================================
           do stype=1,ntype                                               ! Debut de la boucle sur les type de source.
            print*,' Turning on lamp',stype
            write(2,*) ' Turning on lamp',stype
            intensite_totale_type=0.                                      ! Initialisation de la contribution d'un type de source a 
c                                                                         ! l'intensite dirigee vers le capteur par une cellule cible.
            do x_s=1,nbx
             do y_s=1,nby
              ITT(x_s,y_s,stype)=0.
             enddo
            enddo     
            do x_s=imin(stype),imax(stype),step(stype)                    ! Debut de la boucle sur les colonnes (longitudes) du domaine.
             do y_s=jmin(stype),jmax(stype),step(stype)                   ! Debut de la boucle sur les rangees (latitudes) du domaine.
              if (lamp_lumin(x_s,y_s,stype) .ne. 0.) then                 ! Si la luminosite d'une case est nulle, le programme ignore cette case.
               z_s=(alt_sol(x_s,y_s)+lamp_alt(x_s,y_s,stype))             ! Definition de la position (metre) verticale de la source.

c **********************************************************************************************************************
c *     Calcul de l'intensite directe dirigee vers le capteur par une cellule cible en provenance d'une source         *
c **********************************************************************************************************************         
               dircheck=0                                                 ! Initialisation de la verification de la position de la source.
               if ( (x_s.eq.x_c).and.(y_s.eq.y_c).and.( abs(z_s-z_c)      ! Si les positions x et y de la source et de la cible sont les 
c                                                                         ! memes alors.
     +         .lt.(cell_thickness(zcell_c)/2.) ) )then
                dircheck=1
                if (verbose.eq.1) then
                 print*,'Source dans Cellule cible'
                endif
               endif                                                      ! Fin du cas positions x et y source et cible identiques.
               if (dircheck.ne.1) then                                    ! Cas ou la source n'est pas dans la cellule cible.
c=======================================================================
c        Calcul de l'angle zenithal entre la source et la cible
c=======================================================================
c
c calcul de la ligne d'horizon pour les ombrages resolus direct           ! Il y a une ligne d'horizon par cellule cible a une resolution de 1 deg
         
                d2=sqrt((real(x_s-x_c)*dx)**2.+(real(y_s-y_c)*dy)**2.)    ! dist max pour l'horizon (i.e. l horizon passe la source ne compte pas
                call horizon(x_s,y_s,z_s,d2,alt_sol,nbx,nby,dx,dy,
     +          zen_horiz,latitude)
                call anglezenithal
     +          (x_s,y_s,z_s,x_c,y_c,z_c,dx,dy,anglezen)                  ! Calcul de l'angle zenithal entre la source et la cellule cible.
                call angleazimutal(x_s,y_s,x_c,y_c,dx,dy,angleazi)        ! calcul de l'angle azimutal direct cible-source
                az=nint(angleazi*180./pi)+1
                if ((anglezen).lt.zen_horiz(az)) then                     ! la ligne cible-source n'est pas sous l'horizon => on calcule
c                                                                         ! debut condition sous l'horizon direct
c obstacle sous maille             
                 anglemin=pi/2.-atan((alt_sol(x_s,y_s)+obstacleH-z_s)/
     +           portee_reflex)
                 if (anglezen.lt.anglemin) then                           ! Debut condition obstacle sous maille direct.
c
c=======================================================================
c        Calcul de la transmittance entre la source et la cible
c=======================================================================
                  call transmitm (anglezen,x_s,y_s,z_s,x_c,y_c,z_c,
     +            lambda,dx,dy,pression,transm)     
                  call transmita (anglezen,x_s,y_s,z_s,x_c,y_c,z_c,
     +            dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cible vue de la source
c=======================================================================
                  xc=dble(x_c)*dble(dx)                                   ! Position en metres de la cellule observatrice (longitude).
                  yc=dble(y_c)*dble(dy)                                   ! Position en metres de la cellule observatrice (latitude).
                  zc=dble(z_c)                                            ! Position en metres de la cellule observatrice (altitude).
                  xn=dble(x_s)*dble(dx)                                   ! Position en metres de la source (longitude).
                  yn=dble(y_s)*dble(dy)                                   ! Position en metres de la source (latitude).
                  zn=dble(z_s)                                            ! Position en metres de la source (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                  if (z_c .ne. z_s) then
                   call planxy(dx,dy,xc,xn,yc,yn,zc,zn,cell_thickness,
     +             zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,
     +             r4z) 
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Appel de la routine anglesolide qui calcule l'angle solide 
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                   ! selon le plan xy.
                   omega1 = omega
                  else
                   omega1=0.
                  endif
c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
                  if (y_c .ne. y_s) then                                  ! Si la latitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0
                   call planzx(dx,xc,xn,yc,yn,zc,zn,cell_thickness,
     +             zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y
     +             ,r4z)
                   call anglesolide(omega,r1x,r1y,r1z,                    ! Appel de la routine anglesolide qui calcule l'angle solide 
     +             r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                   ! selon le plan zx.
                  else
                   omega=0.
                  endif
                  if (omega.gt.0.) then
                   if (omega .gt. omega1) omega1 = omega                  ! On garde l'angle solide le plus grand jusqu'a present.
                  endif
c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
                  if (x_c .ne. x_s) then                                  ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0.
                   call planyz(dy,xc,xn,yc,yn,zc,zn,cell_thickness,
     +             zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,
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
c=======================================================================
c    Estimation du demi angle directeur de l'angle solide                 ! Cet angle servira a obtenir un meilleur estime (moyenne) de 
c                                                                         ! P_dir pour le cas de grans angles solides ou p_valnorm 
c=======================================================================  ! varie significativement sur +- ouvang.
                  ouvang=sqrt(omega/pi)                                   ! Angle en radian.
                  ouvang=ouvang*180./pi                                   ! Angle en degres.
c   
c=======================================================================
c        Calcul de la fonction d'emission de la source vers la cible
c=======================================================================
c   
                  anglez=nint(anglezen/pi*180.)
                  if (anglez.lt.0) anglez=-anglez
                  if (anglez.gt.180) anglez=360-anglez
                  anglez=anglez+1                                         ! Transformer l'angle en deg entier en position dans la matrice.
c
c  moyenner sur +- ouvang	
c
                  nbang=0.
                  P_dir=0.
                  do na=-nint(ouvang),nint(ouvang)
                   naz=anglez+na
                   if (naz.lt.0) naz=-naz
                   if (naz.gt.181) naz=362-naz                            ! La fonction est symetrique.
                   P_dir=P_dir+p_valnorm(naz,stype)
                   nbang=nbang+1.
                  enddo
                  P_dir=P_dir/nbang
c
c=======================================================================
c        Calcul du flux direct atteignant la cellule cible
c=======================================================================
                  flux_direct=lamp_lumin(x_s,y_s,stype)*P_dir*omega*
     1            transm*transa
c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere directe
c=======================================================================
                  if (anglezen.lt.(pi/2.)) then                           ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                   zidif=zcdown
                   zfdif=zcup
                  else
                   zidif=zcup
                   zfdif=zcdown
                  endif
                  call transmitm (anglezen,iun,iun,zidif,ideux,ideux,     ! Transmittance moleculaire a l'interieur de la cellule diffusante.
     +            zfdif,lambda,dx,dy,pression,trans_1m)
                  call transmita (anglezen,iun,iun,zidif,ideux,ideux,     ! Transmittance aerosols a l'interieur de la cellule diffusante.
     +            zfdif,dx,dy,taua,trans_1a)
                  call angle3points (x_s,y_s,z_s,x_c,y_c,z_c,x_obs,       ! Angle de diffusion.
     +            y_obs,z_obs,dx,dy,angle_dif)
                  call diffusion(omega,angle_dif,trans_1a,trans_1m,       ! Probabilite de diffusion de la lumiere directe.     
     +            secdif,foncdif_anorm,probdif_dir)
c=======================================================================
c   Calcul de la contribution d'une source a l'intensite directe dirigee vers le capteur par une cellule cible
c=======================================================================
                  intensite_directe=flux_direct*probdif_dir
                 else 
                  intensite_directe=0.                                      
                 endif                                                    ! Fin condition obstacle sous maille direct.
                else
c                 print*,'ombrage direct'
                endif                                                     ! Fin condition sous l'horizon direct? 
               endif                                                      ! Fin du cas Position Source n'egale pas Position Cible.
c  fin du calcul de l'intensite directe
c **********************************************************************************************************************
c * Calcul de l'intensite indirecte dirigee vers le capteur par une cellule cible en provenance d'une source           *
c **********************************************************************************************************************
c=======================================================================
c        etablissement des conditions et des boucles
c=======================================================================
               intensite_totale_indirecte=0.                              ! Initialisation de l'intensite indirecte d'une source a une cible.    
               intensite_totale_reflechie_diffusee=0.
               do x_sr=x_s-boxx,x_s+boxx                                  ! Debut de la boucle sur les colonnes (longitude) reflectrices.
                do y_sr=y_s-boxy,y_s+boxy                                 ! Debut de la boucle sur les ranges (latitude) relfectrices.
                 intensite_reflechie=0.
                 z_sr=alt_sol(x_sr,y_sr)   
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
                    if (surface_reflect(x_sr,y_sr).ne.0.) then            ! Condition: la surface reflectrice n'a pas une reflectance nulle.
                     haut=-real(x_s-x_sr)*dx*tan(inclin_x(x_sr,y_sr))     ! Si la variable haut a une valeur negative, cela signifie que la surface
     1               -real(y_s-y_sr)*dy                                   ! reflectrice est eclairee par en-dessous.
     2               *tan(inclin_y(x_sr,y_sr))+z_s-z_sr
                     if (haut .gt. 0.) then                               ! Condition: la surface est eclairee par au-dessus.
c=======================================================================
c        Calcul de l'angle zenithal entre la source et la surface reflectrice
c=======================================================================
                      call anglezenithal(x_s,y_s,z_s,x_sr,y_sr,z_sr,dx,   ! Calcul de l'angle zenithal entre la source et la cellule cible.
     +                dy,anglezen)                                        ! Fin du cas "observateur a la meme latitude/longitude que la source".

c=======================================================================
c        Calcul de la transmittance entre la source et la surface reflectrice
c=======================================================================
                      call transmitm (anglezen,x_s,y_s,z_s,x_sr,y_sr,
     +                z_sr,lambda,dx,dy,pression,transm)          
                      call transmita (anglezen,x_s,y_s,z_s,x_sr,y_sr,
     +                z_sr,dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cellule reflectrice vue de la source
c=======================================================================
                      xc=dble(x_sr)*dble(dx)                              ! Position en metres de la cellule observatrice (longitude).
                      yc=dble(y_sr)*dble(dy)                              ! Position en metres de la cellule observatrice (latitude).
                      zc=dble(z_sr)                                       ! Position en metres de la cellule observatrice (altitude).
                      xn=dble(x_s)*dble(dx)                               ! Position en metres de la source (longitude).
                      yn=dble(y_s)*dble(dy)                               ! Position en metres de la source (latitude).
                      zn=dble(z_s)                                        ! Position en metres de la source (altitude).
                      epsilx=inclin_x(x_sr,y_sr)                          ! Inclinaison en x de la surface reflectrice.
                      epsily=inclin_y(x_sr,y_sr)                          ! Inclinaison en x de la surface reflectrice.
                      if (dx.gt.portee_reflex*2.) then                    ! Utiliser une surface sous-grille lors que la portee de la 
c                                                                         ! reflexion est inferieure a la taille de cellule.
                       if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then
                        dxp=portee_reflex*2.
                       endif
                      else
                       dxp=dx
                      endif
                      if (dy.gt.portee_reflex*2.) then
                       if ((x_sr.eq.x_s).and.(y_sr.eq.y_s)) then         
                        dyp=portee_reflex*2.
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
c                                                                         ! P_dir pour le cas de grans angles solides ou p_valnorm
c=======================================================================  ! varie significativement sur +- ouvang.
                      ouvang=sqrt(omega/pi)                               ! Angle en radian.
                      ouvang=ouvang*180./pi                               ! Angle en degres.
c   
c=======================================================================
c        Calcul de la fonction d'emission du lampadaire vers la surface reflectrice
c=======================================================================
c    
                      anglez=nint(anglezen/pi*180.)
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
                       P_indir=P_indir+p_valnorm(naz,stype)
                       nbang=nbang+1.
                      enddo
                      P_indir=P_indir/nbang
c 
c=======================================================================
c        Calcul du flux atteignant la cellule reflectrice
c=======================================================================
                      flux_reflechi=lamp_lumin(x_s,y_s,stype)*P_indir*
     a                omega*transm*transa
c=======================================================================
c        Calcul de l'intensite reflechie quittant la surface reflectrice
c=======================================================================
                      intensite_reflechie1=flux_reflechi*                  ! Le facteur 1/pi vient de la normalisation de la fonction 
     +                surface_reflect(x_sr,y_sr)/pi
                      if (effetdif.gt.(dx+dy)/2.) then 
                       call reflexdbledif (x_sr,y_sr,z_sr,x_c,y_c,
     +                 zcell_c,dx,dy,effetdif,nbx,nby,stepdif,
     +                 intensite_reflechie1,lambda,pression,taua,zcup,
     +                 zcdown,secdif,foncdif_anorm,x_obs,y_obs,z_obs,
     +                 epsilx,epsily,intensite_reflechie_diffusee,
     +                 portee_reflex,obstacleH,alt_sol,latitude)
                      endif
                      intensite_totale_reflechie_diffusee=
     +                intensite_totale_reflechie_diffusee+
     +                intensite_reflechie_diffusee      
c
c  la projection apparente est calculee a partir du produit scalaire du vecteur normal a 
c  la surface reflechissante et la ligne surface reflechissante vers cellule diffusante ou cible
c         
                      projapparente=(-tan(epsilx)*real(x_c-x_sr)*dx-
     +                tan(epsily)*real(y_c-y_sr)*dy+1.*(cell_height(
     +                zcell_c)-z_sr))/(sqrt(tan(epsilx)**2.+tan(epsily)
     +                **2.+1.)*sqrt((real(x_c-x_sr)*dx)**2.+(real(y_c-
     +                y_sr)*dy)**2.+(cell_height(zcell_c)-z_sr)**2.))
c                                                                         ! Peu importe la direction c'est la valeur absolue du cos theta 
c                                                                         ! qui compte.   
c verifier s'il y a ombrage entre sr et cible ici 
                 d2=sqrt((real(x_sr-x_c)*dx)**2.+(real(y_sr-y_c)*         ! dist max pour l'horizon (i.e. l horizon passe la source ne compte pas
     +           dy)**2.)                 
                 call horizon(x_sr,y_sr,z_sr,d2,alt_sol,nbx,nby,dx,dy,
     +           zen_horiz,latitude) 
                 call anglezenithal(x_sr,y_sr,z_sr,x_c,y_c,z_c,dx,        ! Angle zenithal entre la surface reflechissante et la cellule cible.
     +           dy,anglezen)     
                 call angleazimutal(x_sr,y_sr,x_c,y_c,dx,dy,angleazi)     ! calcul de l'angle azimutal reflect-cible
                 az=nint(angleazi*180./pi)+1          
                 if ((anglezen).lt.zen_horiz(az)) then                    ! la ligne cible-reflec n'est pas sous l'horizon => on calcule
                 
                 
                      if (projapparente.lt.0.) projapparente=0.
                      intensite_reflechie=intensite_reflechie1*
     +                projapparente
                 else
c                  print*,'ombrage reflexion',x_sr,y_sr,z_sr,x_c,y_c,z_c
                 endif                                                    ! fin condition surf. reflectrice au-dessus horizon
c=======================================================================
c        Cas Position Cible = Position Cellule reflectrice
c=======================================================================
                      if((x_c.eq.x_sr).and.(y_c.eq.y_sr).and.
     +                (z_c.eq.z_sr)) then
                       intensite_indirecte=intensite_reflechie
                      else
c
c            
c obstacle                 
                       anglemin=pi/2.-atan(obstacleH/portee_reflex)
                       if (anglezen.lt.anglemin) then                     ! Debut condition obstacle indirect.
c
c=======================================================================
c        Calcul de la transmittance entre la surface reflectrice et la cellule cible
c=======================================================================
                        call transmitm (anglezen,x_sr,y_sr,z_sr,x_c,
     +                  y_c,z_c,lambda,dx,dy,pression,transm)        
                        call transmita (anglezen,x_sr,y_sr,z_sr,x_c,
     +                  y_c,z_c,dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cible vue de la cellule reflectrice
c=======================================================================
                        xc=dble(x_c)*dble(dx)                             ! Position en metres de la cellule observatrice (longitude).
                        yc=dble(y_c)*dble(dy)                             ! Position en metres de la cellule observatrice (latitude).
                        zc=dble(z_c)                                      ! Position en metres de la cellule observatrice (altitude).
                        xn=dble(x_sr)*dble(dx)                            ! Position en metres de la source (longitude).
                        yn=dble(y_sr)*dble(dy)                            ! Position en metres de la source (latitude).
                        zn=dble(z_sr)                                     ! Position en metres de la source (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                        if (z_c .ne. z_sr) then
                         call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                   cell_thickness,zcell_c,r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! selon le plan xy.
                         omega1 = omega
                        else
                         omega1=0.
                        endif
c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
                        if (y_c .ne. y_sr) then                           ! Si la latitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
                         call planzx(dx,xc,xn,yc,yn,zc,zn,
     +                   cell_thickness,zcell_c,r1x,r1y,r1z,r2x,r2y,
     +                   r2z,r3x,r3y,r3z,r4x,r4y,r4z)
                         call anglesolide(omega,r1x,r1y,r1z,              ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                   r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)             ! selon le plan zx.
                        else
                         omega=0.
                        endif
                        if (omega.gt.0.) then
                         if (omega.gt.omega1) omega1 = omega              ! On garde l'angle solide le plus grand jusqu'a present.
                        endif
c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
                        if (x_c.ne.x_sr) then                             ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0.
                         call planyz(dy,xc,xn,yc,yn,zc,zn,
     +                   cell_thickness,zcell_c,r1x,r1y,r1z,r2x,r2y,
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
c=======================================================================
c        Calcul du flux indirect atteignant la cellule cible
c=======================================================================
                        flux_indirect=intensite_reflechie*omega*transm*
     +                  transa      
c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere indirecte
c=======================================================================
                        if (anglezen.lt.(pi/2.)) then                     ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                         zidif=zcdown
                         zfdif=zcup
                        else
                         zidif=zcup
                         zfdif=zcdown
                        endif        
                        call transmitm (anglezen,iun,iun,zidif,ideux,     ! Transmittance moleculaire a l'interieur de la cell diffusante.
     +                  ideux,zfdif,lambda,dx,dy,pression,trans_1m)
                        call transmita (anglezen,iun,iun,zidif,ideux,     ! Transmittance aerosols a l'interieur de la cell diffusante.
     +                  ideux,zfdif,dx,dy,taua,trans_1a)
                        call angle3points (x_sr,y_sr,z_sr,x_c,y_c,z_c,    ! Angle de diffusion.
     +                  x_obs,y_obs,z_obs,dx,dy,angle_dif)
                        call diffusion(omega,angle_dif,trans_1a,          ! Probabilite de diffusion de la lumiere indirecte.
     +                  trans_1m,secdif,foncdif_anorm,probdif_indir)
c=======================================================================
c   Calcul de l'intensite indirecte dirigee vers le capteur par une cellule reflectrice
c=======================================================================
                        intensite_indirecte=flux_indirect*probdif_indir 
                       else
                        intensite_indirecte=0.
                       endif                                              ! Fin condition obstacle indirect.
                      endif                                               ! Fin du cas Posi Cell Reflectrice = Position Cible.                                 
                      intensite_totale_indirecte=
     a                intensite_totale_indirecte+intensite_indirecte      ! Somme des intensites de chaque cell reflectrices propres 
c                                                                         ! a une source.
                     endif                                                ! Fin de la condition surface non-eclairee par le haut.
                    endif                                                 ! Fin de la condition reflectance non-nulle.
                   endif                                                  ! Fin de la condition cellule reflectrice n'est pas une source.
                  endif                                                   ! Fin de la condition surface a l'interieur du domaine.




                enddo                                                     ! Fin de la boucle sur les ranges (latitude) relfectrices.
               enddo                                                      ! Fin de la boucle sur les colonnes (longitude) reflectrices.
c   fin du calcul de l'intensite indirecte        
c **********************************************************************************************************************
c * Calcul de l'intensite diffusee dirigee vers le capteur par une cellule cible en provenance d'une source            *
c **********************************************************************************************************************
c
c=======================================================================
c    Determination des cellules diffusantes en fonction de la cellule source et de la cellule cible
c=======================================================================

               intensite_totale_diffusee=0.                               ! Initialisation de l'intensite diffusee par une source dans 
c                                                                         ! une cellule cible calculer le double diffusion seulement si 
               if (effetdif.gt.(dx+dy)/2.) then                           ! le rayon de diffusion est superieur a la taille des cellules.
                call zone_diffusion(x_s,y_s,z_s,x_c,y_c,zcell_c,dx,dy,
     +          effetdif,nbx,nby,alt_sol,zonedif,ndiff)
                do idi=1,ndiff,stepdif                                 ! Debut de la boucle sur les cellules diffusantes.
                 x_dif=zonedif(idi,1)
                 y_dif=zonedif(idi,2)
                 zcell_dif=zonedif(idi,3)
                 z_dif=cell_height(zcell_dif)             
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
                   call horizon(x_s,y_s,z_s,d2,alt_sol,nbx,nby,
     +             dx,dy,zen_horiz,latitude)
                   call anglezenithal(x_s,y_s,z_s,x_dif,y_dif,z_dif,dx,
     +             dy,anglezen)                                           ! Calcul de l'angle zenithal source-cell diffusante. 
                   call angleazimutal(x_s,y_s,x_dif,y_dif,dx,dy,          ! calcul de l'angle azimutal cible-cell diffusante
     +             angleazi)
                   az=nint(angleazi*180./pi)+1
                   if ((anglezen).lt.zen_horiz(az)) then                  ! debut condition ombrage source-diffusante
c                                                                   
c obstacle sous maille               
                    anglemin=pi/2.-atan((obstacleH+alt_sol(x_s,y_s)-z_s
     +              )/portee_reflex)
                    if (anglezen.lt.anglemin) then                        ! Debut condition obstacle source->diffuse.
c                                                                    
c=======================================================================
c        Calcul de la transmittance entre la source et la cellule diffusante
c=======================================================================
                     call transmitm (anglezen,x_s,y_s,z_s,x_dif,y_dif,
     +               z_dif,lambda,dx,dy,pression,transm)
                     call transmita (anglezen,x_s,y_s,z_s,x_dif,y_dif,
     +               z_dif,dx,dy,taua,transa) 
c=======================================================================
c     Calcul de l'angle solide couvert par la cellule diffusante vue de la source
c=======================================================================
                     xc=dble(x_dif)*dble(dx)                              ! Position en metres de la cellule diffusante (longitude).
                     yc=dble(y_dif)*dble(dy)                              ! Position en metres de la cellule diffusante (latitude).
                     zc=dble(z_dif)                                       ! Position en metres de la cellule diffusante (altitude).
                     xn=dble(x_s)*dble(dx)                                ! Position en metres de la source (longitude).
                     yn=dble(y_s)*dble(dy)                                ! Position en metres de la source (latitude).
                     zn=dble(z_s)                                         ! Position en metres de la source (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                     if (z_dif .ne. z_s) then
                      call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                cell_thickness,zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,
     +                r3x,r3y,r3z,r4x,r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! selon le plan xy.
                      omega1 = omega
                     else
                      omega1=0.
                     endif
c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
                     if (y_dif .ne. y_s) then                             ! Si la latitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
                      call planzx(dx,xc,xn,yc,yn,zc,zn,cell_thickness,
     +                zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x
     +                ,r4y,r4z)
                      call anglesolide(omega,r1x,r1y,r1z,                 ! Appel de la routine anglesolide qui calcule l'angle solide 
     +                r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                ! selon le plan zx.
                     else
                      omega=0.
                     endif
                     if (omega.gt.0.) then
                      if (omega .gt. omega1) omega1 = omega               ! On garde l'angle solide le plus grand jusqu'a present.
                     endif
c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
                     if (x_dif .ne. x_s) then                             ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0.
                      call planyz(dy,xc,xn,yc,yn,zc,zn,cell_thickness,
     +                zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,
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
c=======================================================================
c    Estimation du demi angle directeur de l'angle solide                 ! Cet angle servira a obtenir un meilleur estime (moyenne) de 
c                                                                         ! P_dir pour le cas de grands angles solides ou p_valnorm
c=======================================================================  ! varie significativement sur +- ouvang.
                     ouvang=sqrt(omega/pi)                                ! Angle en radian.
                     ouvang=ouvang*180./pi                                ! Angle en degres.
c 
c=======================================================================
c        Calcul de la fonction d'emission de la source vers la cellule diffusante
c=======================================================================
c    
                     anglez=nint(anglezen/pi*180.)
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
                      P_dif1=P_dif1+p_valnorm(naz,stype)
                      nbang=nbang+1. 
                     enddo
                     P_dif1=P_dif1/nbang 
c
c=======================================================================
c        Calcul du flux atteignant la cellule diffusante
c=======================================================================
                     flux_dif1=lamp_lumin(x_s,y_s,stype)*P_dif1*
     +               omega*transm*transa
c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere diffuse vers la cellule cible
c=======================================================================
                     if (anglezen.lt.(pi/2.)) then                        ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                      zidif=z_c-0.5*cell_thickness(zcell_dif)
                      zfdif=z_c+0.5*cell_thickness(zcell_dif)
                     else
                      zidif=z_c+0.5*cell_thickness(zcell_dif)
                      zfdif=z_c-0.5*cell_thickness(zcell_dif)
                     endif       
                     call transmitm (anglezen,iun,iun,zidif,ideux,        ! Transmittance moleculaire a l'interieur de la cellule diffusante.
     +               ideux,zfdif,lambda,dx,dy,pression,trans_1m)
                     call transmita (anglezen,iun,iun,zidif,ideux,        ! Transmittance aerosols a l'interieur de la cellule diffusante.
     +               ideux,zfdif, dx,dy,taua,trans_1a)
                     call angle3points (x_s,y_s,z_s,x_dif,y_dif,z_dif,    ! Angle de diffusion.
     +               x_c,y_c,z_c,dx,dy,angle_dif)
                     call diffusion(omega,angle_dif,trans_1a,trans_1m,    ! Probabilite de diffusion de la lumiere directe.
     +               secdif,foncdif_anorm,probdif_dif1)
c=======================================================================
c   Calcul de l'intensite diffusee dirigee vers la cellule cible en provenance de la cellule diffusante
c=======================================================================
                     intensite_diffusee1=flux_dif1*probdif_dif1
c=======================================================================
c        Calcul de l'angle zenithal entre la cellule diffusante et la cellule cible
c=======================================================================


        d2=sqrt((real(x_dif-x_c)*dx)**2.+(real(y_dif-y_c)*dy)**2.)        ! dist max pour l'horiz (i.e. l horizon passe la cell-diff ne compte pas)
        call horizon(x_dif,y_dif,z_dif,d2,alt_sol,nbx,nby,dx,dy,
     +  zen_horiz,latitude)


                     call anglezenithal(x_dif,y_dif,z_dif,x_c,y_c,z_c,
     +               dx,dy,anglezen)                                      ! Calcul de l'angle zenithal entre la cellule diffusante et la 
c                                                                         ! cellule cible.


        call angleazimutal(x_dif,y_dif,x_c,y_c,dx,dy,angleazi)            ! calcul de l'angle azimutal surf refl-cell diffusante
        az=nint(angleazi*180./pi)+1
        if ((anglezen).lt.zen_horiz(az)) then                             ! debut condition ombrage diffuse-cible  


c                                                                 
c obstacle sous maille                
                     anglemin=pi/2.-atan((obstacleH+alt_sol(x_dif,
     +               y_dif)-z_dif)/portee_reflex)
                     if (anglezen.lt.anglemin) then                       ! debut condition obstacle sous maille diffuse->cible.
c                                                                   
c=======================================================================
c        Calcul de la transmittance entre la cellule diffusante et la cellule cible
c=======================================================================
                      call transmitm (anglezen,x_dif,y_dif,z_dif,x_c,
     +                y_c,z_c,lambda,dx,dy,pression,transm)
                      call transmita (anglezen,x_dif,y_dif,z_dif,x_c,
     +                y_c,z_c,dx,dy,taua,transa) 
c=======================================================================
c     Calcul de l'angle solide couvert par la cellule cible vue de la cellule diffusante
c=======================================================================
                      xc=dble(x_c)*dble(dx)                               ! Position en metres de la cellule cible (longitude).
                      yc=dble(y_c)*dble(dy)                               ! Position en metres de la cellule cible (latitude).
                      zc=dble(z_c)                                        ! Position en metres de la cellule cible (altitude).
                      xn=dble(x_dif)*dble(dx)                             ! Position en metres de la diffusante (longitude).
                      yn=dble(y_dif)*dble(dy)                             ! Position en metres de la diffusante (latitude).
                      zn=dble(z_dif)                                      ! Position en metres de la diffusante (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
                      if (z_c .ne. z_dif) then
                       call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +                 cell_thickness,zcell_c,r1x,r1y,r1z,r2x,r2y,r2z
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
                      if (y_c .ne. y_dif) then                            ! Si la latitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
                       call planzx(dx,xc,xn,yc,yn,zc,zn,
     +                 cell_thickness,zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,
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
     +                 cell_thickness,zcell_c,r1x,r1y,r1z,r2x,r2y,r2z,
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
                      flux_diffuse=intensite_diffusee1*omega*transm*
     +                transa
c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere diffuse vers la cellule observatrice(SORTANT de cell_c)
c=======================================================================
                      if (anglezen.lt.(pi/2.)) then                       ! Attribution des limites initiale et finale du parcours de 
c                                                                         ! diffusion dans la cellule.
                       zidif=zcdown
                       zfdif=zcup
                      else
                       zidif=zcup
                       zfdif=zcdown
                      endif
                      call transmitm (anglezen,iun,iun,zidif,ideux,       ! Transmittance moleculaire a l'interieur de la cellule diffusante.
     +                ideux,zfdif,lambda,dx,dy,pression,trans_1m)
                      call transmita (anglezen,iun,iun,zidif,ideux,       ! Transmittance aerosols a l'interieur de la cellule diffusante.
     +                ideux,zfdif,dx,dy,taua,trans_1a)    
                      call angle3points (x_dif,y_dif,z_dif,x_c,y_c,       ! Angle de diffusion.
     +                z_c,x_obs,y_obs,z_obs,dx,dy,angle_dif)
                      call diffusion(omega,angle_dif,trans_1a,trans_1m,   ! Probabilite de diffusion de la lumiere directe.
     +                secdif,foncdif_anorm,probdif_dif2)
c=======================================================================
c   Calcul de l'intensite diffusee dirigee vers l'observateur en provenance de la cellule cible
c=======================================================================
                      intensite_diffusee2=flux_diffuse*probdif_dif2
                      intensite_diffusee2=
     +                intensite_diffusee2*real(stepdif)                   ! Corriger le resultat pour le fait d'avoir passe des cellules
                      intensite_totale_diffusee=                          ! afin d'accelerer le calcul.
     +                intensite_totale_diffusee+intensite_diffusee2
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
               endif                                                      ! fin de la condition ou effetdif > dx.
c    fin du calcul de l'intensite diffusee    
c**********************************************************************
c        Calcul de l'intensite provenant d'une source dans la cible vers le capteur
c**********************************************************************
               intensite_source=intensite_directe
     a         +intensite_totale_indirecte+intensite_totale_diffusee 
     +         +intensite_totale_reflechie_diffusee                       ! Somme des intensites de chaque type propre a une source 
c                                                                         ! atteignant une cellule cible.
               if (verbose.eq.1) then
                print*,' Composantes de l''intensite totale:'
                print*,' source->diffuse=',intensite_directe
                print*,' source->reflexion->diffuse=',
     +          intensite_totale_indirecte
                print*,' source->diffuse->diffuse=',
     +          intensite_totale_diffusee
                print*,' source->reflexion->diffuse->diffuse=',
     a          intensite_totale_reflechie_diffusee  
               endif
c                   
c**********************************************************************
c        Calcul de l'intensite totale provenant de toutes les sources d'un type dans la cible vers le capteur 
c**********************************************************************
               intensite_totale_type=intensite_totale_type
     +         +intensite_source*real(step(stype)*step(stype))            ! Somme des intensites de chaque source a une cellule cible.
                                                                          ! Vers le calcul du poids de chaque cellule source i.e. ITT (equivaut 
               ITT(x_s,y_s,stype)=ITT(x_s,y_s,stype)+intensite_source     ! a intensite_totale_type mais sous forme matricielle 
              endif                                                       ! Fin de la condition "La luminosite de la case x_s,y_s n'est pas nulle".
             enddo                                                        ! Fin la boucle sur les rangees (latitudes) du domaine (y_s).
            enddo                                                         ! Fin la boucle sur les colonnes (longitudes) du domaine (x_s).
c
c   fin du calcul de l'intensite par un type de source
            intensite_totale_cible=intensite_totale_cible
     1      + intensite_totale_type                                       ! Somme des intensites de chaque type a une cellule cible.
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
               if (lamp_lumin(ii,jj,stype).ne.0.) then
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
c calculer lamp_lum_tot 
            do x_s=1,nbx
             do y_s=1,nby
               lamp_lum_tot(x_s,y_s)=lamp_lum_tot(x_s,y_s)+   
     +         lamp_lumin(x_s,y_s,stype)
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
     +     anglezen)                                                      ! Calcul de l'angle zenithal entre la cellule cible et l'observateur.
c                                                                         ! Fin du cas "observateur a la meme latitude/longitude que la source".
c=======================================================================
c        Calcul de la transmittance entre la cible et l'observateur
c=======================================================================
           call transmitm (anglezen,x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +     lambda,dx,dy,pression,transm)
           call transmita (anglezen,x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +     dx,dy,taua,transa)
c=======================================================================
c     Calcul de l'angle solide couvert par la cible vu par l'observateur
c=======================================================================
           xn=dble(x_obs)*dble(dx)                                        ! Position en metres de la cellule observatrice (longitude).
           yn=dble(y_obs)*dble(dy)                                        ! Position en metres de la cellule observatrice (latitude).
           zn=dble(z_obs)                                                 ! Position en metres de la cellule observatrice (altitude).
           xc=dble(x_c)*dble(dx)                                          ! Position en metres de la cible (longitude).
           yc=dble(y_c)*dble(dy)                                          ! Position en metres de la cible (latitude).
           zc=dble(z_c)                                                   ! Position en metres de la cible (altitude).
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
           if (z_c .ne. z_obs) then
            call planxy(dx,dy,xc,xn,yc,yn,zc,zn,cell_thickness,zcell_c,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)  
            call anglesolide(omega,r1x,r1y,r1z,                           ! Appel de la routine anglesolide qui calcule l'angle solide 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! selon le plan xy.
            omega1 = omega
           else
            omega1=0.
           endif
c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
           if (y_c .ne. y_obs) then                                       ! Si la latitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan zx car il est egal a 0.
            call planzx(dx,xc,xn,yc,yn,zc,zn,cell_thickness,zcell_c,
     +      r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            call anglesolide(omega,r1x,r1y,r1z,                           ! Appel de la routine anglesolide qui calcule l'angle solide 
     +      r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                          ! selon le plan zx.
           else
            omega=0.
           endif
           if (omega.gt.0.) then
            if (omega.gt.omega1) omega1 = omega                           ! On garde l'angle solide le plus grand jusqu'a present.
           endif
c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
           if (x_c.ne.x_obs) then                                         ! Si la longitude de la cellule observatrice est la meme que celle
c                                                                         ! de la cellule source, on ne calcule pas l'angle solide
c                                                                         ! pour le plan yz car il est egal a 0
            call planyz(dy,xc,xn,yc,yn,zc,zn,cell_thickness,zcell_c,
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
c=======================================================================
c        Calcul du flux atteignant l'objectif du telescope en provenance de la cellule cible
c=======================================================================
           dis_obs=sqrt((z_c-z_obs)**2.+((real(y_c-y_obs))*dy)**2.
     a     +((real(x_c-x_obs))*dx)**2.)
           omegtif=pi*(diamobj/2.)**2./dis_obs**2.
           if (dis_obs.eq.0.) then
            print*,'ERREUR probleme avec dis_obs',dis_obs
            stop
           endif
           flux_cible=intensite_totale_cible*omegtif*transa*transm
           do x_s=1,nbx
            do y_s=1,nby
             FC(x_s,y_s)=ITC(x_s,y_s)*omegtif*transa*transm
            enddo
           enddo
           omegafov=largefente*longfente/distfocal**2.                   ! Calcul de l'angle solide de la fente projete sur le ciel.
c      portion=omegafov/omega                                             ! Methode obsolete pour calculer la fraction de la cellule 
c                                                                         ! cible vue par le fov (Fraction peut etre superieure a 1 si 
c                                                                         ! omegafov.gt.omega).
           if (cos(pi-anglezen).eq.0.) then 
            print*,'ERREUR la visee est pafaitement horizontale!'
            stop
           else
            portion=(omegafov*dis_obs*dis_obs)/(cos(pi-anglezen)*dx*dy)   ! Fraction de la cellule cible vue par le fov (Fraction peut 
c                                                                         ! etre superieure a 1). Le pi ici est du au fait
c                                                                         ! que anglezen est calcule sur le trajet cible vers l'observateur
           endif
           if (omega.eq.0.) then
            print*,'ERREUR omega=0 (1)'
            stop
           endif
           flux_capteur=flux_cible*portion
           do x_s=1,nbx
            do y_s=1,nby
             FCA(x_s,y_s)=FC(x_s,y_s)*portion
            enddo
           enddo
c   fin du calcul du flux atteignant la cellule observatrice en provenance d'une cellule cible
           flux_total_capteur=flux_total_capteur+flux_capteur         
           do x_s=1,nbx
            do y_s=1,nby
             FTC(x_s,y_s)=FTC(x_s,y_s)+FCA(x_s,y_s)                       ! FTC est la matrice du flux total au capteur permettant d'identifier
                                                                          ! la contribution de chaque cellule du sol au flux total au capteur
                                                                          ! Le % est simplement donne par FTC/flux_total_capteur
             fluxcumule=fluxcumule+FCA(x_s,y_s)
            enddo
           enddo
          endif                                                           ! Fin de la condition cellule cible n'est pas cellule observatrice.
         endif                                                            ! Fin de la condition cellule cible dans le domaine de calcul.
        endif                                                             ! Fin condition pour la reprise d'un calcul arrete.
        print*,' Flux capteur          =',flux_capteur
        print*,' Flux capteur cumule   =',flux_total_capteur
        write(2,*) ' Flux capteur          =',flux_capteur
        write(2,*) ' Flux capteur cumule   =',flux_total_capteur    
       enddo                                                              ! Fin de la boucle sur les cellules cibles.
       if (prmaps.eq.1) then
          open(unit=9,file=pclfile,status='unknown')
          open(unit=8,file=pcwfile,status='unknown')
            FTCNTOT=0.
            ftcmax=0.
            do x_s=1,nbx
               do y_s=1,nby
                  FTC(x_s,y_s)=FTC(x_s,y_s)/flux_total_capteur   
                  if (FTC(x_s,y_s).gt.ftcmax) ftcmax=FTC(x_s,y_s)
                  if (lamp_lum_tot(x_s,y_s).ne.0.) then
                  FTCNTOT=FTCNTOT+FTC(x_s,y_s)/lamp_lum_tot(x_s,y_s)
                  endif
               enddo
            enddo
            if (verbose.eq.1) then
               print*,'Writing normalized contribution matrix'
            endif
            do x_s=1,nbx
               do y_s=1,nby
                  if (lamp_lum_tot(x_s,y_s).ne.0.) then
                     FTCN(x_s,y_s)=(FTC(x_s,y_s)/lamp_lum_tot(x_s,y_s))
     +               /FTCNTOT  
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
          open(unit=9,file=pclgplot,status='unknown')
          open(unit=8,file=pcwgplot,status='unknown')
            write(9,*) 'set dgrid3d',nbx,',',nby
            write(9,*) 'set hidden3d'
            write(9,*) 'set pm3d'
            write(9,*) 'splot "'//basename(1:lenbase)//'_pcl.txt"
     +      with dots'
            write(8,*) 'set dgrid3d',nbx,',',nby
            write(8,*) 'set hidden3d'
            write(8,*) 'set pm3d'
            write(8,*) 'splot "'//basename(1:lenbase)//'_pcw.txt"
     +      with dots'    
          close(unit=8)
          close(unit=9) 
       endif                                                              ! end of condition for creating contrib and sensit maps
          print*,'====================================================='
          print*,'          Total flux entering instrument (W)'
          print*,'                   ',flux_total_capteur*real(vistep)  
          print*,'              Sky luminance (W/str/m**2)'       
          print*,'                   ',flux_total_capteur/(largefente*
     +          longfente/distfocal**2.)/(pi*(diamobj/2.)**2.)*
     +          real(vistep)
       print*,'  '
       print*,' Interpolation flux error= ',
     +          flux_total_capteur-fluxcumule
       print*,'======================================================='
       write(2,*) '==================================================='
       write(2,*) '          Total flux entering instrument (W)'
       write(2,*) '                   ',flux_total_capteur*real(vistep)
        write(2,*) '            Sky luminance (W/str/m**2)          '      
       write(2,*) '                   ',flux_total_capteur/(largefente*
     +          longfente/distfocal**2.)/(pi*(diamobj/2.)**2.)*
     +          real(vistep)
       write(2,*) '  '                                                
       write(2,*) 'Interpolation flux errror= ',
     +          flux_total_capteur-fluxcumule
       write(2,*) '==================================================='
      close(2)
      stop
      end
c***********************************************************************************************************************
c*                                                                                                                     *
c*                                         fin du programme                                                            *
c*                                                                                                                     *
c***********************************************************************************************************************
