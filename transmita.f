c-----------------------------------------------------------------------
c
c=======================================================================
c Routine transmita (Andre Morin, Alex Neron, Etienne Rousseau 2004)
c
c Determine la transmittance des aerosols sur un parcours entre les 
c cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
c Est fonction de l'epaisseur optique des aerosols
c Recoit des longueurs d'ondes en nanometre et les transforme en microns.
c La pression doit etre en KPa
c Retourne la transmittance transa
c
c  *** J'ai valide le calcul zenith tout atm avec modtran et 
c      cela concorde M. Aubé mars 2010
c
c pour utilisation avec Illumina 
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
      subroutine transmita(angz,z_i,z_f,distd,transa,tranaa)
      real angle,transa,transa1,transa2                                   ! Declaration des variables.
      real tranaa                                                         ! vertical transmittance of the complete atmosphere (aerosols)
      real angz
      real z_i,z_f,z1,z2
      real airm
      real distd
      airm=1.
         if (cos(angz).eq.0.) then
           print*,'prob w airm 2',angz
           stop
         endif
      if (abs(angz-1.5707963).gt.0.00001) airm=1./abs(cos(angz))
         if (z_i.gt.z_f) then
            z2=z_i
            z1=z_f
         else
            z1=z_i
            z2=z_f
         endif           
         transa1=1.-((tranaa-1.)*(exp(-1.*z2/2000.)-exp(-1.*z1
     +   /2000.)))
         transa1=transa1**airm
         transa2=1.-((1.-tranaa)/2000.*distd*exp(-(z1+z2)/(2.*2000.)))
      if (abs(distd-airm*abs(z1-z2)).le.0.1) then
        transa=transa1
c        print*,'verti',transa,distd,airm*abs(z1-z2),z1,z2
      else
        transa=transa2
c        print*,'horiz',transa,distd,airm*abs(z1-z2),z1,z2
      endif
         if (transa.eq.0.) then
            print*,'ERREUR transa - no transmission',z_i,z_f,airm
     +      ,angz
            stop
         endif
         if (transa.gt.1.) then 
            print*,'ERREUR avec transa',transa,z_i,z_f,airm,angz
            stop
         endif
      return
      end
