c **********************************************************************************************************************
c * Calcul de l'intensite diffusee dirigee vers le capteur par une cellule cible en provenance d'une cellule reflectrice *
c **********************************************************************************************************************
c
c=======================================================================
c    Determination des cellules diffusantes en fonction de la cellule reflechissante et de la cellule cible
c=======================================================================
c   
c    Copyright (C) 2010  Martin Aube
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
             subroutine reflexdbledif (x_sr,y_sr,z_sr,x_c,y_c,zcell_c,
     +       dx,dy,effdif,nbx,nby,stpdif,irefl,lambda
     +       ,pressi,taua,zcup,zcdown,secdif,fdifan,
     +       x_obs,y_obs,z_obs,epsilx,epsily,
     +       irefdi,drefl,hobst,ofill,
     +       alt_sol,latitu,cloudt,cloudh,icloud,stype)
c
c   declarations de variables
c  
      integer width
      parameter (width=1024)  
      integer x_sr,y_sr,x_dif,y_dif,zcell_dif                             ! Positions source, surface reflectrice, celldiffusantes (cellule)
      integer x_c,y_c,zcell_c,nbx,nby,stype
      real z_sr,z_dif,dx,dy                             
      real cell_t(50),alt_sol(width,width)                                ! Matrice de l'epaisseur des niveaux (metre)

      real cell_h(50)                                                     ! Matrice de la hauteur du centre de chaque niveau (metre)

      real effdif                                                         ! Distance autour des cellule source et cible qui seront considerees pour calculer la double diffusion
      integer zondif(3000000,4)                                           ! Matrice des cellules diffusantes
      integer ndiff,idi                                                   ! Nombre de cell diffusantes, compteur de boucle sur cell-diffus
      integer stpdif                                                      ! saut de diffus pour accelerer calcul e.g. si =2 fera un calcul/2
      integer iun,ideux
      real flux_dif1                                                      ! Flux atteignant une cellule diffusante
      real irefl
      real pi,angzen
      real lambda,pressi
      real transm,transa,taua
      real*8 xc,yc,zc,xn,yn,zn
      real*8 r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z              ! Composantes des vecteurs utilises dans la routine angle solide
      real irefdi
      real omega,omega1                                                   ! Angle solide couvert par une cellule vue d'une autre, angle de comparaison
      real zidif,zfdif                                                    ! Limites initiale et finale du parcours de diffusion dans une cellule
      real tran1a,tran1m                                                  ! Transmittance a l'interieur d'une cellule (aerosols,molecules)
      real angdif, secdif
      real pdif_dif1,pdifd2                                               ! Probabilite de diffusion (directe,indirecte,premiere et deuxieme diffusions
      real fdifan(181)                                                    ! Fonction de diffusion (arbitraire et normalisee) des aerosols
      real projap                                                         ! taille apparente de la surface
      real idiff1                                                         ! Intensite dirigee vers une cellule cible par une cellule diffusante
      real zcup,zcdown                                                    ! Limites inferieure et superieure de la cellule cible        
      real idiff2                                                         ! Contribution d'une cellule diffusante a l'intensite diffusee dirigee vers le capteur
      real fdiff                                                          ! Flux provenant d'une cellule diffusante dans une cellule cible (watt)
      integer x_obs,y_obs                                                 ! Position de l'observateur (cellule)
      real z_obs
      real epsilx,epsily                                                  ! inclinaison de la surface reflechissante
      real angmin,hobst(width,width),drefl(width,width)
      real ofill(width,width)
      parameter (pi=3.1415926)
      real zenhor(360),d2,angazi                                          ! angle zenithal de l'horizon,distance horizon, angle azimut
      integer az
      real latitu
      integer cloudt                                                      ! cloud type 0=clear, 1=Thin Cirrus/Cirrostratus, 2=Thick Cirrus/Cirrostratus, 3=Altostratus/Altocumulus, 4=Cumulus/Cumulonimbus, 5=Stratocumulus
      integer cloudh(5)                                                   ! cloud base layer relative to the lower elevation 
      real rcloud                                                         ! cloud relfecicloudtance 
      real azencl                                                         ! zenith angle from cloud to observer
      real icloud                                                         ! cloud reflected intensity
      real zero,anaz
      real monte                                                          ! random number between 0 and 1

      data cell_t /0.5,0.6,0.72,0.86,1.04,1.26,1.52,1.84,2.22,            ! Epaisseur des niveaux
     a 2.68,3.24,3.92,4.74,5.72,6.9,8.34,10.08,12.18,14.72,17.78,21.48,
     b 25.94,31.34,37.86,45.74,55.26,66.76,80.64,97.42,117.68,142.16,
     c 171.72,207.44,250.58,302.7,365.66,441.72,533.6,644.58,778.66,
     d 940.62,1136.26,1372.6,1658.1,2002.98,2419.6,2922.88,3530.84,
     e 4265.26,5152.44/
      data cell_h /0.25,0.8,1.46,2.25,3.2,4.35,5.74,7.42,9.45,            ! Hauteur du centre de chaque niveau
     a 11.9,14.86,18.44,22.77,28.,34.31,41.93,51.14,62.27,75.72,91.97,
     b 111.6,135.31,163.95,198.55,240.35,290.85,351.86,425.56,514.59,
     c 622.14,752.06,909.,1098.58,1327.59,1604.23,1938.41,2342.1,
     d 2829.76,3418.85,4130.47,4990.11,6028.55,7282.98,8798.33,
     e 10628.87,12840.16,15511.4,18738.26,22636.31,27345.16/
      iun=1
      ideux=2
      zero=0.
      call zone_diffusion(x_sr,y_sr,z_sr,x_c,y_c,zcell_c,                ! Determiner la zone de diffusion
     +dx,dy,effdif,nbx,nby,alt_sol,zondif,ndiff)
      z_c=cell_h(zcell_c)
      irefdi=0.                                                          ! Initialisation de l'intensite diffus par une source ds 1 cell cible         
      do idi=1,ndiff,stpdif                                              ! Debut de la boucle sur les cellules diffusantes
       x_dif=zondif(idi,1)
       y_dif=zondif(idi,2)
       zcell_dif=zondif(idi,3)
       z_dif=cell_h(zcell_dif)
c
c  la projection apparente est calculee a partir du produit scalaire du vecteur normal a 
c  la surface reflechissante et la ligne surface reflechissante vers cellule diffusante ou cible
c         
       projap=(-tan(epsilx)*real(x_dif-x_sr)*dx-tan(epsily)*
     + real(y_dif-y_sr)*dy+1.*(z_dif-z_sr))/(sqrt(
     + tan(epsilx)**2.+tan(epsily)**2.+1.)*sqrt((real(x_dif-x_sr)
     + *dx)**2.+(real(y_dif-y_sr)*dy)**2.+(z_dif-z_sr)**2.))
       if (projap.lt.0.) projap=0.    
        if((x_dif.gt.nbx).or.(x_dif.lt.1).or.(y_dif.gt.nby).or.           ! Condition cellule diffusante a l'interieur du domaine
     +  (y_dif.lt.1)) then     
        else
         if(((x_sr.eq.x_dif).and.(y_sr.eq.y_dif).and.
     +   (z_sr.eq.z_dif)) .or.
     +   ((x_c.eq.x_dif).and.(y_c.eq.y_dif).and. 
     +   (z_c.eq.z_dif))) then
         else
c ombrage s_reflechissante-diffusante
          call anglezenithal(x_sr,y_sr,z_sr,x_dif,y_dif,z_dif,dx,dy,      ! Calcul de l'angle zenithal entre la surf reflechissante et la cell diff
     +    angzen)                                                       
          call angleazimutal(x_sr,y_sr,x_dif,y_dif,dx,dy,angazi)          ! calcul de l'angle azimutal surf refl-cell diffusante
          az=nint(angazi*180./pi)+1
          d2=sqrt((real(x_dif-x_sr)*dx)**2.+(real(y_dif-y_sr)*dy)**2.)    ! dist max pour l'horiz (i.e. l horizon passe la cell-diff ne compte pas)
          call horizon(x_sr,y_sr,z_sr,d2,alt_sol,nbx,nby,dx,dy,
     +    zenhor,latitu,angazi) 
          if ((angzen).lt.zenhor(az)) then                                ! debut condition ombrage surface refl - diffuse        
c MA j'ai verifie que angzen ne depasse  jamais pi ou jamais moins que 0
                                                                          ! Fin du cas "observateur a la meme latitu/longitude que la source"
c obstacle sous maille
           angmin=pi/2.-atan(hobst(x_sr,y_sr)/
     +     drefl(x_sr,y_sr))
           monte=rand()
           if ((monte.gt.ofill(x_sr,y_sr)).or.
     +     (angzen.lt.angmin)) then                                       ! debut condition obstacle reflechi->diffuse               
c=======================================================================
c        Calcul de la transmittance entre la surface reflechissane et la cellule diffusante
c=======================================================================
            anaz=zero
            call transmitm(angzen,anaz,x_sr,y_sr,z_sr,x_dif,y_dif,
     +      z_dif,lambda,dx,dy,pressi,transm)
c MA j'ai verifie que transm est > 0 et <=1
            call transmita(angzen,anaz,x_sr,y_sr,z_sr,x_dif,y_dif,
     +      z_dif,dx,dy,taua,transa) 
c MA j'ai verifie que transa est > 0 et <=1
c=======================================================================
c     Calcul de l'angle solide couvert par la cellule diffusante vue de la surface reflechissante
c=======================================================================
            xc=dble(x_dif)*dble(dx)                                       ! Position en metres de la cellule diffusante (longitude)
            yc=dble(y_dif)*dble(dy)                                       ! Position en metres de la cellule diffusante (latitu)
            zc=dble(z_dif)                                                ! Position en metres de la cellule diffusante (altitude)
            xn=dble(x_sr)*dble(dx)                                        ! Position en metres de la surface (longitude)
            yn=dble(y_sr)*dble(dy)                                        ! Position en metres de la surface (latitu)
            zn=dble(z_sr)                                                 ! Position en metres de la surface (altitude)
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
            if (z_dif .ne. z_sr) then
             call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +       r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                  
             call anglesolide(omega,r1x,r1y,r1z,                          ! Appel de la rout. anglesolide qui calcule l'ang. solide selon le plan xy
     +       r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
             omega1 = omega
            else
             omega1=0.
            endif
c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
            if (y_dif .ne. y_sr) then                                     ! Si la latitu de la cellule observatrice est la meme que celle
                                                                          ! de la cellule source, on ne calcule pas l'angle solide
             call planzx(dx,xc,xn,yc,yn,zc,zn,cell_t,zcell_c,
     +       r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
             call anglesolide(omega,r1x,r1y,r1z,                          ! Appel de la rout. anglesolide qui calcule l'ang solide selon le plan zx
     +       r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            else
             omega=0.
            endif
            if (omega.gt.0.) then
             if (omega .gt. omega1) omega1 = omega                        ! On garde l'angle solide le plus grand jusqu'a present
            endif
c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
            if (x_dif .ne. x_sr) then                                     ! Si la longitude de la cellule observatrice est la meme que celle
                                                                          ! de la cellule source, on ne calcule pas l'angle solide
                                                                          ! pour le plan yz car il est egal a 0
             call planyz(dy,xc,xn,yc,yn,zc,zn,cell_t,zcell_c,
     +       r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
             call anglesolide(omega,r1x,r1y,r1z,                          ! Appel de la rout anglesolide qui calcule l'angle solide selon le plan yz
     +       r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            else 
             omega=0.
            endif
            if (omega.gt.0.) then
             if (omega .gt. omega1) omega1 = omega                        ! On garde l'angle solide le plus grand
            endif
            omega=omega1
c oups omega depasse pi et va meme jusqu a 6.26 ->ok c'est normal puisque on observe a peu pres la demi sphere
            if (omega.gt.2.*pi) then
             print*,'omega=',omega
             stop 
            elseif (omega.lt.0.) then
             print*,'omega=',omega
             stop 
            endif
c=======================================================================
c        Calcul du flux atteignant la cellule diffusante
c=======================================================================
            flux_dif1=irefl*projap*omega*transm*
     +      transa
c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere diffuse vers la cellule cible
c=======================================================================
            if (angzen.lt.(pi/2.)) then                                   ! Attribution des lim init et finale du parcours de diffu ds la cellule
             zidif=z_c-0.5*cell_t(zcell_dif)
             zfdif=z_c+0.5*cell_t(zcell_dif)
            else
             zidif=z_c+0.5*cell_t(zcell_dif)
             zfdif=z_c-0.5*cell_t(zcell_dif)
            endif 
            anaz=angazi
            call transmitm(angzen,anaz,iun,iun,zidif,ideux,iun,zfdif,         ! Transmittance moleculaire a l'interieur de la cellule diffusante
     +      lambda,dx,dy,pressi,tran1m)
            call transmita(angzen,anaz,iun,iun,zidif,ideux,iun,zfdif,         ! Transmittance aerosols a l'interieur de la cellule diffusante
     +      dx,dy,taua,tran1a)
            call angle3points(x_sr,y_sr,z_sr,x_dif,y_dif,z_dif,x_c,       ! Angle de diffusion
     +      y_c,z_c,dx,dy,angdif)
            call diffusion(omega,angdif,tran1a,tran1m,secdif,             ! Probabilite de diffusion de la lumiere directe      
     +      fdifan,pdif_dif1)
            if (flux_dif1.lt.0.) print*,'FLUX_DIF1=',flux_dif1,
     +      irefl,projap,omega,transm,transa
c=======================================================================
c   Calcul de l'intensite diffusee dirigee vers la cellule cible en provenance de la cellule diffusante
c=======================================================================     
            idiff1=flux_dif1*pdif_dif1
c=======================================================================
c        Calcul de l'angle zenithal entre la cellule diffusante et la cellule cible
c=======================================================================
c ombrage s_reflechissante-diffusante

     
     

            call anglezenithal(x_dif,y_dif,z_dif,x_c,y_c,z_c,dx,dy,
     +      angzen)                                                       ! Calcul de l'ang zenit  cellule diffus. - la cellule cible   
     
     
        call angleazimutal(x_dif,y_dif,x_c,y_c,dx,dy,angazi)              ! calcul de l'angle azimutal surf refl-cell diffusante
        az=nint(angazi*180./pi)+1
        d2=sqrt((real(x_dif-x_c)*dx)**2.+(real(y_dif-y_c)*dy)**2.)        ! dist max pour l'horiz (i.e. l horizon passe la cell-diff ne compte pas)
        call horizon(x_dif,y_dif,z_dif,d2,alt_sol,nbx,nby,dx,dy,
     +  zenhor,latitu,angazi)  

        if ((angzen).lt.zenhor(az)) then                                  ! debut condition ombrage diffuse-cible           
     
     
c obstacle sous maille
            angmin=pi/2.-atan((hobst(x_dif,y_dif)+
     +      alt_sol(x_dif,y_dif)-z_dif)/drefl(x_dif,y_dif))
            monte=rand()
            if ((monte.gt.ofill(x_dif,y_dif)).or.
     +      (angzen.lt.angmin)) then                                      ! debut condition obstacle sous maille diffuse->cible
                                                                          ! Fin du cas "observateur a la meme latitu/longitude que la source"
c=======================================================================
c        Calcul de la transmittance entre la cellule diffusante et la cellule cible
c=======================================================================
            anaz=zero
            call transmitm(angzen,anaz,x_dif,y_dif,z_dif,x_c,y_c,z_c,
     +      lambda,dx,dy,pressi,transm)
            call transmita(angzen,anaz,x_dif,y_dif,z_dif,x_c,y_c,z_c,
     +      dx,dy,taua,transa) 
c=======================================================================
c     Calcul de l'angle solide couvert par la cellule cible vue de la cellule diffusante
c=======================================================================
            xc=dble(x_c)*dble(dx)                                         ! Position en metres de la cellule cible (longitude)
            yc=dble(y_c)*dble(dy)                                         ! Position en metres de la cellule cible (latitu)
            zc=dble(z_c)                                                  ! Position en metres de la cellule cible (altitude)
            xn=dble(x_dif)*dble(dx)                                       ! Position en metres de la diffusante (longitude)
            yn=dble(y_dif)*dble(dy)                                       ! Position en metres de la diffusante (latitu)
            zn=dble(z_dif)                                                ! Position en metres de la diffusante (altitude)
c    ------------------------------------
c    Angle solide pour le plan central xy
c    ------------------------------------
            if (z_c .ne. z_dif) then
             call planxy(dx,dy,xc,xn,yc,yn,zc,zn,
     +       r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z) 
             call anglesolide(omega,r1x,r1y,r1z,                          ! Appel de la rout anglesolide qui calcule l'angle solide selon le plan xy
     +       r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
             omega1 = omega
            else
             omega1=0.
            endif
c     ------------------------------------
c     Angle solide pour le plan central zx
c     ------------------------------------
            if (y_c .ne. y_dif) then                                      ! Si la latitu de la cellule observatrice est la meme que celle
                                                                          ! de la cellule source, on ne calcule pas l'angle solide
                                                                          ! pour le plan zx car il est egal a 0
             call planzx(dx,xc,xn,yc,yn,zc,zn,cell_t,zcell_c,
     +       r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)                                                                         
             call anglesolide(omega,r1x,r1y,r1z,                          ! Appel de la rout anglesolide qui calcule l'angle solide selon le plan zx
     +       r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            else
             omega=0.
            endif
            if (omega.gt.0.) then
             if (omega .gt. omega1) omega1 = omega                        ! On garde l'angle solide le plus grand jusqu'a present
            endif
c     ------------------------------------
c     Angle solide pour le plan central yz
c     ------------------------------------
            if (x_c .ne. x_dif) then                                      ! Si la longitude de la cellule observatrice est la meme que celle
                                                                          ! de la cellule source, on ne calcule pas l'angle solide
                                                                          ! pour le plan yz car il est egal a 0
             call planyz(dy,xc,xn,yc,yn,zc,zn,cell_t,zcell_c,
     +       r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
             call anglesolide(omega,r1x,r1y,r1z,                          ! Appel de la rout anglesolide qui calcule l'angle solide selon le plan yz
     +       r2x,r2y,r2z,r3x,r3y,r3z,r4x,r4y,r4z)
            else 
             omega=0.
            endif
            if (omega.gt.0.) then
             if (omega .gt. omega1) omega1 = omega                        ! On garde l'angle solide le plus grand
            endif
            omega=omega1
c=======================================================================
c        Calcul du flux diffuse atteignant la cellule cible
c=======================================================================
            fdiff=idiff1*omega*transm*transa

c verifie mais je crosi que le calcul de azencl est facultatif ici
                if (cloudt.ne.0) then                                     ! target cell = cloud
                  if (cloudh(cloudt).eq.zcellc) then
                     call anglezenithal(x_c,y_c,z_c,x_obs,y_obs,z_obs,
     +               dx,dy,azencl)                                        ! zenith angle from cloud to observer                     
                     call cloudreflectance(angzen,cloudt,rcloud)          ! cloud intensity from direct illum
                     icloud=icloud+
     +               fdiff*rcloud*abs(cos(azencl))/pi
                  endif
                endif



c=======================================================================
c   Calcul de la probabilite de diffusion de la lumiere diffuse vers la cellule observatrice(SORTANT de cell_c)
c=======================================================================
            if (angzen.lt.(pi/2.)) then                                   ! Attribution des limites init et finale du parcours de diffus ds la cell
             zidif=zcdown
             zfdif=zcup
            else
             zidif=zcup
             zfdif=zcdown
            endif 
            anaz=angazi
            call transmitm(angzen,anaz,iun,iun,zidif,ideux,iun,zfdif,         ! Transmittance moleculaire a l'interieur de la cellule diffusante
     +      lambda,dx,dy,pressi,tran1m)
            call transmita(angzen,anaz,iun,iun,zidif,ideux,iun,zfdif,         ! Transmittance aerosols a l'interieur de la cellule diffusante
     +      dx,dy,taua,tran1a)    
            call angle3points (x_dif,y_dif,z_dif,x_c,y_c,z_c,x_obs,       ! Angle de diffusion
     +      y_obs,z_obs,dx,dy,angdif)
            call diffusion(omega,angdif,tran1a,tran1m,secdif,             ! Probabilite de diffusion de la lumiere directe
     +      fdifan,pdifd2)
c=======================================================================
c   Calcul de l'intensite diffusee dirigee vers l'observateur en provenance de la cellule cible
c=======================================================================
            idiff2=fdiff*pdifd2*real(stpdif)                              ! corriger le result pr avoir passe des cell afin d'accel le calcul
            irefdi=irefdi+idiff2      
           endif                                                          ! fin condition obstacle sous maille diffuse->cible 
        endif                                                             ! fin condition ombrage diffuse-cible
          endif                                                           ! fin condition obstacle reflechie->diffuse    
         endif                                                            ! fin  condition ombrage surface refl - diffuse 
        endif                                                             ! Fin du cas Diffusante = Source ou Cible        
       endif                                                              ! Fin de la condition "cellule a l'interieur du domaine"             
      enddo                                                               ! Fin de la boucle sur les cellules diffusante
c     fin du calcul de l'intensite diffusee
      return
      end      
