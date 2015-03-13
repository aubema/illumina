c------------------------------------------------------------------------
c
c=======================================================================
c  Routine diffusion (Alex Neron 2004) (Modifie par Martin Aube, Valerie Houle, Philippe Robert-Staehler)
c
c  Determine la probabilite de diffusion de la lumiere par unite d'angle solide
c  dans la direction angle_dif. Les parametres de diffusion sont donnes
c  par secdif (rapport de la section efficace de diffusion sur la section
c  efficace d'extinction totale), fonc_anorm (fonction de diffusion
c  normalisees des aerosols)
c  Retourne la probabilite de diffusion probdif
c
c  pour utilisation avec Illumina
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2009  Martin Aube
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
      subroutine diffusion (omega,angle_dif,trans_a,trans_m,secdif,
     +   fonc_a,probdif)
      real angle_dif,probdif,prob_a,prob_m,trans_a,trans_m,secdif 
      real fonc_mol,pi,fonc_a(181),fonc_ae
      real angle_deg
      real omega,ouvang,nbang
      integer rang,na,naz
      parameter (pi=3.1415926)
c--------------------------------------------------------
c   Calcul et normalisation des fonctions de diffusion
c--------------------------------------------------------      
      if (angle_dif.lt.0.) angle_dif=-angle_dif      
      if (angle_dif-pi.gt.0.00001) angle_dif=pi
      angle_deg=((angle_dif*180.)/pi)
c
c=======================================================================
c    Estimation du demi angle directeur de l'angle solide                 ! Cet angle servira a obtenir un meilleur estime (moyenne) de 
c                                                                         ! P_dir pour le cas de grands angles solides ou p_valnorm
c=======================================================================  ! varie significativement sur +- ouvang.
                ouvang=sqrt(omega/pi)                                     ! Angle en radian.
                ouvang=ouvang*180./pi                                     ! Angle en degres.
c   
c=======================================================================
c        Calcul de la fonction d'emission de la source vers la cible
c=======================================================================
c                                                                         ! Transformer l'angle en deg entier en position dans la matrice.
c  moyenner sur +- ouvang       
c
                rang=int(angle_deg)+1
                nbang=0.
                fonc_ae=0.
                do na=-nint(ouvang),nint(ouvang)
                  naz=rang+na
                  if (naz.le.0) naz=-naz+2
                  if (naz.gt.181) naz=362-naz                             ! La fonction est symetrique.
                  fonc_ae=fonc_ae+fonc_a(naz)
                  nbang=nbang+1.
                enddo
                fonc_ae=fonc_ae/nbang 
                nbang=0.
                fonc_mol=0.
                do na=-nint(ouvang),nint(ouvang)
       fonc_mol=fonc_mol+0.75*(1.+((cos(angle_dif+real(na)*pi/180.))
     +                  **2.))/(4.*pi)                                    ! L'integrale de la fonction de diffusion des molecules sur tous 
c                                                                         ! les angles solides est de 4pi  fonc_ae=fonc_ae+fonc_a(naz).
                  nbang=nbang+1.
                enddo
                fonc_mol=fonc_mol/nbang                      
c      
c----------------------------------------
c  Calcul des probabilites de diffusion par unite d'angle solide 
c----------------------------------------    
      prob_a=(1.-trans_a)*secdif*fonc_ae                                  ! Les fonctions utilisees ici sont deja normalisees
      prob_m=(1.-trans_m)*fonc_mol                                        ! Fonc_ae normalisee dans le MAIN, fonc_mol dans la routine (voir 
c                                                                         ! la division par 4 pi).
      probdif = prob_a+prob_m                                             ! Ce calcul est approximatif et bon seulement si 1-trans_a et
c                                                                         ! 1-trans_m sont tres petits.
      return
      end
