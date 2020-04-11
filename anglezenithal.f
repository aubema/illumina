c=======================================================================
c  Routine anglezenithal (Andre Morin 2004)
c
c  debugge par Martin Aube 2004 (cette routine ne calculait pas du tout
c  l'angle zenithal
c
c  Determine l'angle zenithal entre les points (x1,y1,z1) et (x2,y2,z2)
c  Retourne l'angle angzen en radians
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
      subroutine anglezenithal(x1,y1,z1,x2,y2,z2,angzen) 
      real x1,y1,x2,y2
      real z1,z2,pi,angzen                                                    
      parameter (pi=3.141592654)   
      hdist=sqrt((x2-x1)**2.
     a    +(y2-y1)**2.)
      if (z2-z1.ne.0.) then
         angzen=atan(hdist/abs(z2-z1))
         if (z2-z1.lt.0.) angzen=pi-angzen
      else
         angzen=pi/2.
      endif
      if ((angzen.lt.0.).or.(angzen.gt.pi)) then
        print*,'ERREUR angzen2=',angzen
        stop
      endif
      return
      end 
