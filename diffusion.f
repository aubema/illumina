c------------------------------------------------------------------------
c
c=======================================================================
c  Routine diffusion
c
c Determine the probability of light scattering per unit of solid angle
c in the angdif direction. The scattering parameters are given
c by secdif (ratio of the effective section of scattering on the section
c of extinction), fonc_anorm (scattering phase function
c
c Returns the scattering probability pdif
c
c an Illumina routine
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2021  Martin Aube
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
      subroutine diffusion(angdif,tranam,tranaa,tranal,un,secdif,secdil,
     +   fonc_a,fonc_l,haer,hlay,pdif,altit)
      real angdif,pdif,prob_a,prob_m,prob_l,secdif,secdil
      real fctmol,pi,fonc_a(181),fonc_l(181),fonc_ae,fonc_le
      real angdeg,tranam,tranaa,tranal
      real altit,un,hlay,haer
      integer rang,na,naz
      parameter (pi=3.1415926)
c--------------------------------------------------------      
      if (angdif.lt.0.) angdif=-angdif      
      if (angdif-pi.gt.0.00001) angdif=pi
      angdeg=((angdif*180.)/pi)
      rang=int(angdeg)+1   
c----------------------------------------
c  Calculate scattering probability per unit of steradian                 ! The probability is for a voxel of 1x1x1m refer to equation 1 in Aubé et al.
c
c  Aubé, M., Simoneau, A., Muñoz-Tuñón, C., Díaz-Castro, J., & 
c  Serra-Ricart, M. (2020). Restoring the night sky darkness at 
c  Observatorio del Teide: First application of the model Illumina 
c  version 2. Monthly Notices of the Royal Astronomical Society, 
c  497(3), 2501-2516.
c---------------------------------------- 
      if ((tranaa.le.1.).and.(tranaa.gt.0.)) then
         fonc_ae=fonc_a(rang)                                             ! value of the aerosol phase function
         prob_a=(1.-exp(log(tranaa)*exp(-1.*altit/haer)*un/haer))*        ! Functions are normalized in the main code. See their division by 4pi
     +   secdif*fonc_ae
      else
         prob_a=0.
      endif      
      if ((tranal.le.1.).and.(tranal.gt.0.)) then
         fonc_le=fonc_l(rang)                                             ! value of the layer phase function
         prob_l=(1.-exp(log(tranal)*exp(-1.*altit/hlay)*un/hlay))*           
     +   secdil*fonc_le
      else
         prob_l=0.
      endif
      fctmol=0.75*(1.+((cos(angdif))**2.))/(4.*pi)                        ! value of the molecule phase function
      prob_m=(1.-exp(log(tranam)*exp(-1.*altit/8000.)*un/8000.))*
     +fctmol

      pdif = prob_a+prob_m+prob_l                                         ! This is an approximation valide if 1-transa,
                                                                          ! 1-transm et 1-transl are small
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
      if (prob_l.gt.1.) then
         print*,'prob_l>1.'
         stop
      endif
      if (prob_l.lt.0.) then
         print*,'prob_l`¸^<0..'
         stop
      endif
      if (pdif.gt.1.) then
         pdif=1.
      endif
      if (pdif.lt.0.) then
         pdif=0.
      endif
      return
      end
