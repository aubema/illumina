c-----------------------------------------------------------------------
c
c=======================================================================
c Routine distance (Martin Aube 2019)
c
c Determine the traveled distance 
c
c pour utilisation avec Illumina 
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2019  Martin Aube
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
      subroutine distance(angle,anaz,x_i,y_i,z_i,x_f,y_f,z_f,
     + dx,dy,distd,axe1)
      real angle,pi                                                       ! Declaration des variables.
      real z_i,z_f,dx,dy,distd,z1,z2                                      ! angle is the zenith angle
      integer x_i,y_i,x_f,y_f
      real anaz    
      integer axe1
      pi=3.1415926                                                        ! Attribution d'une valeur a la constante pi.
      if (axe1.eq.0) then
        distd=abs(z_f-z_i)/abs(cos(angle))
        if (cos(angle).eq.0.) then
          print*,'ERREUR cos(angle)=0 (1b), anaz=',angle
          print*,x_i,y_i,z_i,zinf,x_f,y_f,z_f,zsup
          stop
        endif   
      else
        distd=sqrt((real(x_f-x_i)*dx)**2.+(real(y_f-y_i)*dy)**2.)
        if (distd.lt.dx) then
          anaz=abs(anaz-real(nint(anaz/(pi/2.)))*(pi/2.))
           
          distd=dx/abs(sin(angle))/abs(cos(anaz))
        endif
        if (sin(angle).eq.0.) then
          print*,'ERREUR sin(angle)=0 (1b), angle=',angle
          print*,x_i,y_i,z_i,zinf,x_f,y_f,z_f,zsup
          stop
        endif 
        if (cos(anaz).eq.0.) then
          print*,'ERREUR cos(anaz)=0 (1b), anaz=',anaz
          print*,x_i,y_i,z_i,zinf,x_f,y_f,z_f,zsup
          stop
        endif   
      endif
      if (distd.lt.0.) then
        print*,'Distd lt 0.', distd
        stop
      endif
      return
      end
