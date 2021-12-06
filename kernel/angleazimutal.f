c=======================================================================
c  Routine angleazimutal (Martin Aube 2010)
c
c
c  Determine l'angle azimutal entre les points (x1,y1,z1) et (x2,y2,z2)
c  Retourne l'angle anglezen en radians
c
c  pour utilisation avec Illumina 
c-----------------------------------------------------------------------
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
      subroutine angleazimutal(x1,y1,x2,y2,angazi) 
      real x1,y1,x2,y2
      real pi,angazi                                                   
      parameter (pi=3.141592654)   
      if (x2-x1.ne.0.) angazi=abs(atan((y2-y1)/(x2-x1)))
      if (((x2-x1).eq.0.).and.((y2-y1).eq.0.)) then
         angazi=0.
      else
        if (x2-x1.gt.0.) then
         if (y2-y1.lt.0.) then
           angazi=2.*pi-angazi 
         endif 
        elseif (x2-x1.lt.0.) then
         if (y2-y1.lt.0.) then
           angazi=angazi+pi 
         else
           angazi=pi-angazi
         endif
        else                                                              ! i.e. x2-x1=0.
         if (y2.gt.y1) angazi=pi/2.
         if (y2.lt.y1) angazi=3.*pi/2.
        endif
        if ((angazi.lt.0.).or.(angazi.gt.2.*pi)) then
          print*,'ERREUR angazi=',angazi,x1,y2,x2,y2
          stop
        endif
c        if ((x1.eq.x2).and.(y1.eq.y2)) then
c          print*,'ERREUR cant compute angle between identical points!'
c          stop
c        endif
      endif
      return
        end 
