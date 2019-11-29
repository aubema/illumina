c------------------------------------------------------------------------
c
c=======================================================================
c  Routine diffusion (Alex Neron 2004) (Modifie par Martin Aube, Valerie Houle, Philippe Robert-Staehler)
c
c  Determine la probabilite de diffusion de la lumiere par unite d'angle solide
c  dans la direction angdif. Les parametres de diffusion sont donnes
c  par secdif (rapport de la section efficace de diffusion sur la section
c  efficace d'extinction totale), fonc_anorm (fonction de diffusion
c  normalisees des aerosols)
c  Retourne la probabilite de diffusion pdif
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
      subroutine diffusion (angdif,tranam,tranaa,un,secdif,
     +   fonc_a,pdif,altit)
      real angdif,pdif,prob_a,prob_m,secdif 
      real fctmol,pi,fonc_a(181),fonc_ae
      real angdeg,tranam,tranaa
      real altit,un
      integer rang,na,naz
      parameter (pi=3.1415926)
c--------------------------------------------------------
c   Calcul et normalisation des fonctions de diffusion
c--------------------------------------------------------      
      if (angdif.lt.0.) angdif=-angdif      
      if (angdif-pi.gt.0.00001) angdif=pi
      angdeg=((angdif*180.)/pi)
      rang=int(angdeg)+1
c=======================================================================
c        Calcul de la fonction d'emission de la source vers la cible
c=======================================================================
      fonc_ae=fonc_a(rang)
      fctmol=0.75*(1.+((cos(angdif))**2.))/(4.*pi)
c----------------------------------------
c  Calcul des probabilites de diffusion par unite d'angle solide 
c----------------------------------------    
c      prob_a=(1.-tranaa)*exp(-1.*altit/2000.)*distd*secdif*fonc_ae        ! Les fonctions utilisees ici sont deja normalisees
c     +/2000.
c      prob_m=(1.-tranam)*exp(-1.*altit/8000.)*distd*fctmol/8000.               ! Fonc_ae normalisee dans le MAIN, fctmol dans la routine (voir 
c                                                                         ! la division par 4 pi).

      prob_a=(1.-exp(log(tranaa)*exp(-1.*altit/2000.)*un/2000.))*        ! Les fonctions utilisees ici sont deja normalisees
     +secdif*fonc_ae
      prob_m=(1.-exp(log(tranam)*exp(-1.*altit/8000.)*un/8000.))*               ! Fonc_ae normalisee dans le MAIN, fctmol dans la routine (voir 
     +fctmol                                                                    ! la division par 4 pi).

      pdif = prob_a+prob_m                                                ! Ce calcul est approximatif et bon seulement si 1-transa et
                                                                          ! 1-transm sont tres petits.
      if (prob_a.gt.1.) then
         print*,'prob_a>1.'
         stop
      endif
      if (prob_a.lt.0.) then
         print*,'prob_a<0..'
         stop
      endif
      if (prob_m.gt.1.) then
         print*,'prob_m>1.'
         stop
      endif
      if (prob_m.lt.0.) then
         print*,'prob_m`¸^<0..'
         stop
      endif      
      if (pdif.gt.1.) then
         print*,'prob>1.',pdif,prob_a,prob_m,tranaa,tranam,altit,distd,
     +omega,omega*prob_a
         stop
      endif
      if (pdif.lt.0.) then
         print*,'prob<0.',pdif,prob_a,prob_m
         stop
      endif
       
c      if (pdif.gt.1.) pdif=1.

      return
      end
