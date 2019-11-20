c-----------------------------------------------------------------------
c
c=======================================================================
c Routine transmitm (Andre Morin, Alex Neron, Etienne Rousseau 2004)
c debuggee par Martin Aube 2004
c Determine la transmittance des molecules atmospheriques sur un parcours
c entre les cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
c Recoit des longueurs d'ondes en nanometre et les transforme en microns.
c La pressi doit etre en KPa
c Retourne la transmittance transm
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
      subroutine transmitm(angz,z_i,z_f,distd,transm,tranam)
      real transm,transm1,transm2                                                         ! Declaration des variables.
      real tranam                                                         ! vertical transmittance of the complete atmosphere (molecules)
      real angz
      real distd
      real airm
      real z_i,z_f,z1,z2
      airm=1.
      if (cos(angz).eq.0.) then
        print*,'prob w airm 1',angz
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
      transm1=1.-((tranam-1.)*(exp(-1.*z2/8000.)-exp(-1.*z1/8000.)))
      transm1=transm1**airm
      transm2=1.-((1.-tranam)/8000.*distd*exp(-(z1+z2)/(2.*8000.)))
      if (abs(distd-airm*abs(z1-z2)).le.0.1) then                            ! si la distance fournie et celle calculée avec la différence de z fois la masse d'air est faible. Donc longue distance
        transm=transm1
      else
        transm=transm2
      endif
      if ((transm.lt.0.).or.(transm.gt.1.)) then
        print*,'ERREUR avec transm',tranam,z_f,z_i,distd,angle,
     +  lamdm,lambda
        stop
      endif
      return
      end
