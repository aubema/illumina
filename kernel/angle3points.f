c-----------------------------------------------------------------------
c
c=======================================================================
c  Routine angle3points (Andre Morin 2004)
c  debuggee et modifiee par Martin Aube 2004
c
c  Determine l'angle entre 3 points (x1,y1,z1), (x2,y2,z2) et (x3,y3,z3)
c  dont le sommet est au point 2
c  Retourne l'angle angle en radians
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
      subroutine angle3points(x1,y1,z1,x2,y2,z2,x3,y3,z3,an3pts)
      real x1,y1,x2,y2,x3,y3
      real z1,z2,z3,an3pts,argume
      real xu,yu,zu,xv,yv,zv,dx,dy,pi,comp
      parameter (pi = 3.141592654)
      xu=x2-x1                                                            ! Voici les composantes du vecteur u.
      yu=y2-y1
      zu=z2-z1
      xv=x3-x2                                                            ! Voici les composantes du vecteur v.
      yv=y3-y2
      zv=z3-z2
      if ((xv.eq.0.).and.(yv.eq.0.).and.(zv.eq.0.)) then
         print*,'ERREUR vecteur sortie nul'
         print*,x1,y1,z1,x2,y2,z2,x3,y3,z3
         stop
      endif
      if ((xu.eq.0.).and.(yu.eq.0.).and.(zu.eq.0.)) then
         print*,'ERREUR vecteur d entree nul'
         print*,x1,y1,z1,x2,y2,z2,x3,y3,z3
         stop
      endif
      argume=(xu*xv+yu*yv+zu*zv)/(sqrt(xu**2.+yu**2.+zu**2.)*
     a        sqrt(xv**2.+yv**2.+zv**2.))
      if (argume.ge.1.) then
         an3pts=0.
      elseif (argume.le.-1.) then
         an3pts=pi
      else
         an3pts=acos(argume)
      endif 
      if (an3pts.lt.0.) then
         print*,'ERREUR an3pts < 0'
         stop
      endif
      return
      end
