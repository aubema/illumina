c-----------------------------------------------------------------------
c
c=======================================================================
c Routine transmitl 
c
c Determine la transmittance des aerosols sur un parcours entre les 
c cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
c Est fonction de l'epaisseur optique des aerosols
c Recoit des longueurs d'ondes en nanometre et les transforme en microns.
c La pression doit etre en KPa
c Retourne la transmittance transl
c
c
c pour utilisation avec Illumina 
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
      subroutine transmitl(angz,z_i,z_f,distd,hlay,transl,tranal)
      real transl                                                         ! Declaration des variables.
      real tranal                                                         ! vertical transmittance of the complete atmosphere (aerosols)
      real angz,hlay
      real z_i,z_f,z1,z2
      real distd


      if (z_i.gt.z_f) then
        z2=z_i
        z1=z_f
      else
        z1=z_i
        z2=z_f
      endif
      if (z1.ne.z2) then    
        transl=exp((log(tranal)/abs(cos(angz)))*(exp(-1.*z1/hlay)-
     +  exp(-1.*z2/hlay)))
      else
        transl=exp((log(tranal))*exp(-1.*z1/hlay)*distd/hlay)  
      endif
         if (transl.eq.0.) then
            print*,'ERREUR transl - no transmission',z_i,z_f,angz
         endif
         if (transl.gt.1.) then 
            print*,'ERREUR avec transa',transl,z_i,z_f,angz
            stop
         endif
      return
      end
